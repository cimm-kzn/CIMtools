# -*- coding: utf-8 -*-
#
#  Copyright 2016, 2017 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of MODtools.
#
#  MODtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
from gzip import open as gzip_open
from tempfile import mkdtemp
from threading import Thread, active_count
from time import sleep
from subprocess import call
from os import W_OK, access, mkdir
from os.path import join, exists, isdir, dirname
from dill import load, dump
from copy import deepcopy
from collections import OrderedDict
from functools import partial
from itertools import product, cycle
from sortedcontainers import SortedListWithKey
from CGRtools.files.RDFrw import RDFread
from CGRtools.files.SDFrw import SDFread
from .config import GACONF
from .descriptors.cxcalc import Pkab
from .descriptors.descriptoragregator import DescriptorsDict, DescriptorsChain
from .descriptors.eed import Eed
from .descriptors.fragmentor import Fragmentor
from .estimators.svmodel import SVModel
from .parsers import MBparser


def desc_starter(func, in_file, out_file, fformat, header, is_reaction):
    with open(in_file) as f:
        inp = (RDFread(f) if is_reaction else SDFread(f)).read()

    dsc = func(structures=inp, parsesdf=True)
    if dsc and not (dsc['X'].isnull().values.any() or dsc['Y'].isnull().any()):
        fformat(out_file, dsc['X'], dsc['Y'], header=header)
        return True

    print('BAD Descriptor generator params')
    return False


class ModelBuilder(MBparser):
    def __init__(self, **kwargs):
        if not kwargs['output'] and (not exists(kwargs['description']) or isdir(kwargs['description'])):
            raise Exception('path to model description file invalid')

        description = self.parse_model_description(kwargs['description'])
        description['type'] = 2 if kwargs['isreaction'] else 1
        if any(x not in description for x in ('nlim', 'name', 'example', 'description')):
            raise Exception('Description Invalid')

        self.__options = kwargs
        self.__description = description
        self.__generators = self.__descriptors_config(kwargs)
        self.__estimators = []

    def run(self):
        if not self.__options['output']:
            if isdir(self.__options['model']) or \
               (exists(self.__options['model']) and not access(self.__options['model'], W_OK)) or \
               isdir(self.__options['model'] + '.save') or \
               (exists(self.__options['model'] + '.save') and
                    not (access(self.__options['model'] + '.save', W_OK) or
                         self.__options['model'] + '.save' == self.__options['reload'])) or \
               not access(dirname(self.__options['model']), W_OK):
                print('path for model saving not writable')
                return
            self.prepare_estimators()
            self.fit()
        else:
            self.__gen_desc(self.__options['output'], fformat=self.__options['format'], header=True)

    def prepare_estimators(self):
        svm = {'svr', 'svc'}.intersection(self.__options['estimator']).pop()
        rf = False  # todo: implement RF
        if self.__options['reload']:
            tmp = load(gzip_open(self.__options['reload'], 'rb'))
            self.__clean_desc_gens(tmp['descgens'])
            if svm:
                est_svm_params = tmp['svm']
        else:
            for_save = dict(descgens=list(range(len(self.__generators))))
            if svm:
                if not self.__options['svm']:
                    svm_cfg, cleared = self.__dragos_svm_fit(svm)
                    self.__clean_desc_gens(cleared)
                    for_save['descgens'] = cleared
                else:
                    svm_cfg = self.get_svm_param(self.__options['svm'])

                est_svm_params = self.__chk_est(svm_cfg)
                for_save['svm'] = est_svm_params

            ''' save configuration
            '''
            dump(for_save, gzip_open(self.__options['model'] + '.save', 'wb'))
            ''' end save
            '''

        if svm:
            self.__estimators.append((partial(SVModel, estimator=svm, probability=self.__options['probability'],
                                      max_iter=self.__options['max_iter']), est_svm_params))

    @staticmethod
    def __descriptors_config(options):
        descgenerator = OrderedDict()
        if options['fragments']:
            descgenerator['F'] = [partial(Fragmentor, is_reaction=options['isreaction'], **x)
                                  for x in MBparser.parse_fragmentor_opts(options['fragments'])]

        if options['extension']:
            descgenerator['E'] = [partial(DescriptorsDict, **MBparser.parse_ext(options['extension']))]

        if options['eed']:
            descgenerator['D'] = [partial(Eed, is_reaction=options['isreaction'], **x)
                                  for x in MBparser.parse_fragmentor_opts(options['eed'])]

        if options['pka']:
            descgenerator['P'] = [partial(Pkab, is_reaction=options['isreaction'], **x)
                                  for x in MBparser.parse_fragmentor_opts(options['pka'])]

        if options['chains']:
            if options['ad'] and len(options['ad']) != len(options['chains']):
                raise Exception('number of generators chains should be equal to number of ad modifiers')

            descgens = []
            for ch, ad in zip(options['chains'], options['ad'] or cycle([None])):
                gen_chain = [x for x in ch.split(':') if x in descgenerator]
                if ad:
                    ad_marks = [x in ('y', 'Y', '1', 'True', 'true') for x in ad.split(':')]
                    if len(ad_marks) != len(gen_chain):
                        raise Exception('length of generators chain should be equal to length of ad modifier')
                else:
                    ad_marks = cycle([True])

                ad_chain = OrderedDict()
                for k, v in zip(gen_chain, ad_marks):
                    ad_chain.setdefault(k, []).append(v)

                combo = []
                for k, v in ad_chain.items():
                    try:
                        if len(v) > 1:
                            if len(descgenerator[k]) != len(v):
                                raise Exception('length of same generators chain should be equal '
                                                'to number of same generators')
                            combo.append([list(zip(descgenerator[k], v))])
                        else:
                            combo.append(list(zip(descgenerator[k], cycle(v))))
                    except:
                        raise Exception('Invalid chain. check configured descriptors generators')
                descgens.extend([DescriptorsChain(*[(g(), a) for gs in c
                                                    for g, a in (gs if isinstance(gs, list) else [gs])])
                                 for c in product(*combo)])
        else:
            descgens = [g() for x in descgenerator.values() for g in x]

        if not descgens:
            raise Exception('Descriptor generators not configured')

        return descgens

    def __order(self, model):
        s = (1 if self.__options['fit'] == 'rmse' else -1) * model.get_model_stats()[self.__options['fit']]
        v = self.__options['dispcoef'] * model.get_model_stats()['%s_var' % self.__options['fit']]
        return s + v

    def fit(self):
        models = SortedListWithKey(key=self.__order)
        with open(self.__options['input']) as f:
            inp = (RDFread(f) if self.__options['isreaction'] else SDFread(f)).read()

        for g, e in self.__estimators:
            for x, y in zip(self.__generators, e):
                models.add(g(x, list(y.values()), inp, parsesdf=True,
                           dispcoef=self.__options['dispcoef'], fit=self.__options['fit'],
                           scorers=self.__options['scorers'],
                           n_jobs=self.__options['n_jobs'], nfold=self.__options['nfold'],
                           rep_boost=self.__options['rep_boost'], repetitions=self.__options['repetition'],
                           normalize='scale' in y or self.__options['normalize']))

                if len(models) > self.__options['consensus']:
                    models.pop()

        if 'tol' not in self.__description:
            self.__description['tol'] = models[0].get_model_stats()['dragostolerance']

        print('name', self.__description['name'])
        print('description', self.__description['description'])
        print('tol', self.__description['tol'])
        print('nlim', self.__description.get('nlim'))
        dump(dict(models=models, config=self.__description), gzip_open(self.__options['model'], 'wb'))

    def __chk_est(self, estimatorparams):
        if not estimatorparams or 1 < len(estimatorparams) < len(self.__generators) or \
                        len(estimatorparams) > len(self.__generators):
            print('NUMBER of estimator params files SHOULD BE EQUAL to '
                  'number of descriptor generator params files or to 1')
            raise Exception('SVM estimators not configured')

        if len(estimatorparams) == 1:
            tmp = []
            for i in range(len(self.__generators)):
                tmp.append(deepcopy(estimatorparams[0]))
            estimatorparams = tmp
        return estimatorparams

    def __gen_desc(self, output, fformat='svm', header=False):
        queue = enumerate(self.__generators, start=1)
        workpath = mkdtemp(prefix='svm_', dir=self.__options['workpath'])
        while True:
            if active_count() < self.__options['n_jobs']:
                tmp = next(queue, None)
                if tmp:
                    n, dgen = tmp
                    subworkpath = join(workpath, str(n))
                    mkdir(subworkpath)
                    dgen.set_work_path(subworkpath)
                    t = Thread(target=desc_starter,
                               args=[dgen.get, self.__options['input'], '%s.%d' % (output, n),
                                     (self.save_svm if fformat == 'svm' else self.save_csv),
                                     header, self.__options['isreaction']])
                    t.start()
                else:
                    while active_count() > 1:
                        sleep(2)
                    break
            sleep(2)

        return True

    def __dragos_svm_fit(self, _type):
        """ files - basename for descriptors.
        """
        workpath = mkdtemp(prefix='gac_', dir=self.__options['workpath'])
        files = join(workpath, 'drag')
        dragos_work = join(workpath, 'work')

        execparams = [GACONF, workpath, _type, str(self.__options['ga_maxconfigs']),
                      str(self.__options['repetition']), str(self.__options['nfold'])]
        if self.__gen_desc(files):
            if call(execparams) == 0:
                best = {}
                with open(join(dragos_work, 'best_pop')) as f:
                    for n, line in enumerate(f):
                        if len(best) == self.__options['best_pop']:
                            break  # get only first n's params.
                        dset, normal, *_, attempt, _, _ = line.split()
                        best.setdefault(int(dset[5:]), (normal, attempt))

                cleared, svmpar, scale = [], [], []
                for k, (nv, av) in best.items():
                    cleared.append(k - 1)
                    svmpar.append(join(dragos_work, av, 'svm.pars'))
                    scale.append(nv)

                svm = []
                svmpar = self.get_svm_param(svmpar)
                if len(svmpar) == len(scale):
                    for x, y in zip(svmpar, scale):
                        svm.append({'scale' if y == 'scaled' else 'orig': list(x.values())[0]})
                    return svm, cleared
        raise Exception('GAConf failed')

    def __clean_desc_gens(self, select):
        self.__generators = [self.__generators[x] for x in select]
