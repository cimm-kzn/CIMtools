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
from io import StringIO
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
from .parsers import MBparser, argparser


def desc_starter(func, file_dump, out_file, fformat, header, is_reaction):
    with StringIO(file_dump) as f:
        inp = (RDFread(f) if is_reaction else SDFread(f)).read()

    dsc = func(structures=inp, parsesdf=True)
    if dsc and not (dsc['X'].isnull().values.any() or dsc['Y'].isnull().any()):
        fformat(out_file, dsc['X'], dsc['Y'], header=header)
        return True

    print('BAD Descriptor generator params')
    return False


class ModelBuilder(MBparser):
    def __init__(self, description, workpath='.', is_reaction=True, model=None, reload=None, output=None,
                 out_format='svm', fragments=None, extension=None, eed=None, pka=None, chains=None, ad=None,
                 estimator='svr', svm=None, max_iter=100000, probability=False, fit='rmse', dispcoef=0, best_pop=20,
                 n_jobs=2, nfold=5, repetition=1, rep_boost=25, consensus=10, normalize=False, ga_maxconfigs=3000,
                 scorers=('rmse', 'r2')):

        if not (exists(workpath) and access(workpath, W_OK)):
            raise Exception('work path not writable')

        if not output:
            if not exists(description) or isdir(description):
                raise Exception('path to model description file invalid')

            if not model:
                raise Exception('path to model save file invalid')

            model_save = '%s.save' % model
            if isdir(model) or not access(dirname(model) or '.', W_OK) or (exists(model) and not access(model, W_OK)) \
                    or not reload and (isdir(model_save) or (exists(model_save) and not access(model_save, W_OK))):
                raise Exception('path for model saving not writable or incorrect')

            if reload and (not exists(reload) or isdir(reload)):
                raise Exception('reload file path incorrect')

            if reload:
                self.__reload = load(gzip_open(reload, 'rb'))
            else:
                self.__reload = None
                self.__svm = self.get_svm_param(svm) if svm else None

            self.__output = False
            self.__model = model
            self.__model_save = model_save
        else:
            if not access(dirname(output) or '.', W_OK) or \
                    (exists('%s.1.svm' % output) and not access('%s.1.svm' % output, W_OK)) or \
                    (exists('%s.1.csv' % output) and not access('%s.1.csv' % output, W_OK)):
                raise Exception('path for descriptors saving not writable or incorrect')

            self.__output = output
            self.__format = out_format

        description = self.parse_model_description(description)
        description['type'] = 2 if is_reaction else 1
        if any(x not in description for x in ('nlim', 'name', 'example', 'description')):
            raise Exception('Description Invalid')

        self.__workpath = workpath
        self.__estimator = estimator
        self.__description = description
        self.__max_iter = max_iter
        self.__probability = probability
        self.__fit = fit
        self.__disp_coef = dispcoef
        self.__is_reaction = is_reaction
        self.__best_pop = best_pop
        self.__n_jobs = n_jobs
        self.__nfold = nfold
        self.__repetition = repetition
        self.__rep_boost = rep_boost
        self.__consensus = consensus
        self.__normalize = normalize
        self.__ga_maxconfigs = ga_maxconfigs
        self.__scorers = scorers
        self.__generators = self.__descriptors_config(is_reaction, fragments=fragments, extension=extension, eed=eed,
                                                      pka=pka, chains=chains, ad=ad)

    __estimators = []

    def run(self, input_file):
        if not self.__output:
            self.prepare_estimators(input_file)
            self.fit(input_file)
        else:
            self.__gen_desc(input_file, self.__output, fformat=self.__format, header=True)

    def prepare_estimators(self, input_file):
        svm = {'svr', 'svc'}.intersection(self.__estimator).pop()
        rf = False  # todo: implement RF
        if self.__reload:
            self.__clean_desc_gens(self.__reload['descgens'])
            if svm:
                self.__svm = self.__reload['svm']
        else:
            for_save = dict(descgens=list(range(len(self.__generators))))
            if svm:
                if not self.__svm:
                    self.__svm, cleared = self.__dragos_svm_fit(input_file, svm)
                    self.__clean_desc_gens(cleared)
                    for_save['descgens'] = cleared

                self.__chk_est()
                for_save['svm'] = self.__svm

            ''' save configuration
            '''
            dump(for_save, gzip_open(self.__model_save, 'wb'))
            ''' end save
            '''

        if svm:
            self.__estimators.append((partial(SVModel, estimator=svm, probability=self.__probability,
                                      max_iter=self.__max_iter), self.__svm))

    @staticmethod
    def __descriptors_config(is_reaction, fragments=None, extension=None, eed=None, pka=None, chains=None, ad=None):
        descgenerator = OrderedDict()
        if fragments:
            descgenerator['F'] = [partial(Fragmentor, is_reaction=is_reaction, **x)
                                  for x in MBparser.parse_fragmentor_opts(fragments)]

        if extension:
            descgenerator['E'] = [partial(DescriptorsDict, **MBparser.parse_ext(extension))]

        if eed:
            descgenerator['D'] = [partial(Eed, is_reaction=is_reaction, **x)
                                  for x in MBparser.parse_fragmentor_opts(eed)]

        if pka:
            descgenerator['P'] = [partial(Pkab, is_reaction=is_reaction, **x)
                                  for x in MBparser.parse_fragmentor_opts(pka)]

        if chains:
            if ad and len(ad) != len(chains):
                raise Exception('number of generators chains should be equal to number of ad modifiers')

            descgens = []
            for ch, ad in zip(chains, ad or cycle([None])):
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
        s = (1 if self.__fit == 'rmse' else -1) * model.get_model_stats()[self.__fit]
        v = self.__disp_coef * model.get_model_stats()['%s_var' % self.__fit]
        return s + v

    def fit(self, input_file):
        models = SortedListWithKey(key=self.__order)
        with open(input_file) as f:
            data = f.read()

        for g, e in self.__estimators:
            for x, y in zip(self.__generators, e):
                inp = (RDFread(StringIO(data)) if self.__is_reaction else SDFread(StringIO(data))).read()
                models.add(g(x, list(y.values()), inp, parsesdf=True, dispcoef=self.__disp_coef, fit=self.__fit,
                             scorers=self.__scorers, n_jobs=self.__n_jobs, nfold=self.__nfold,
                             rep_boost=self.__rep_boost, repetitions=self.__repetition,
                             normalize='scale' in y or self.__normalize))

                if len(models) > self.__consensus:
                    models.pop()

        if 'tol' not in self.__description:
            self.__description['tol'] = models[0].get_model_stats()['dragostolerance']

        print('name', self.__description['name'])
        print('description', self.__description['description'])
        print('tol', self.__description['tol'])
        print('nlim', self.__description.get('nlim'))
        dump(dict(models=models, config=self.__description), gzip_open(self.__model, 'wb'))

    def __chk_est(self):
        est_params = self.__svm
        if not est_params or 1 < len(est_params) < len(self.__generators) or len(est_params) > len(self.__generators):
            print('NUMBER of estimator params files SHOULD BE EQUAL to '
                  'number of descriptor generator params files or to 1')
            raise Exception('SVM estimators not configured')

        if len(est_params) == 1:
            tmp = []
            for i in range(len(self.__generators)):
                tmp.append(deepcopy(est_params[0]))
            self.__svm = tmp

    def __clean_desc_gens(self, select):
        self.__generators = [self.__generators[x] for x in select]

    def __gen_desc(self, input_file, output, fformat='svm', header=False):
        with open(input_file) as f:
            data = f.read()

        queue = enumerate(self.__generators, start=1)
        workpath = mkdtemp(prefix='svm_', dir=self.__workpath)
        while True:
            if active_count() < self.__n_jobs:
                tmp = next(queue, None)
                if tmp:
                    n, dgen = tmp
                    subworkpath = join(workpath, str(n))
                    mkdir(subworkpath)
                    dgen.set_work_path(subworkpath)
                    t = Thread(target=desc_starter,
                               args=[dgen.get, data, '%s.%d' % (output, n),
                                     (self.save_svm if fformat == 'svm' else self.save_csv),
                                     header, self.__is_reaction])
                    t.start()
                else:
                    while active_count() > 1:
                        sleep(2)
                    break
            sleep(2)

        return True

    def __dragos_svm_fit(self, input_file, _type):
        """ files - basename for descriptors.
        """
        workpath = mkdtemp(prefix='gac_', dir=self.__workpath)
        files = join(workpath, 'drag')
        dragos_work = join(workpath, 'work')

        execparams = [GACONF, workpath, _type, str(self.__ga_maxconfigs),
                      str(self.__repetition), str(self.__nfold)]

        if self.__gen_desc(input_file, files):
            if call(execparams) == 0:
                best = {}
                with open(join(dragos_work, 'best_pop')) as f:
                    for n, line in enumerate(f):
                        if len(best) == self.__best_pop:
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


def launcher():
    args = argparser()
    input_file = args.pop('input')
    main = ModelBuilder(args)
    main.run(input_file)
