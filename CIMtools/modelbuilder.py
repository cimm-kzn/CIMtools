# -*- coding: utf-8 -*-
#
#  Copyright 2016, 2017 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of CIMtools.
#
#  CIMtools is free software; you can redistribute it and/or modify
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
from CGRtools.files.RDFrw import RDFread
from CGRtools.files.SDFrw import SDFread
from collections import OrderedDict
from copy import deepcopy
from dill import load, dump
from functools import partial
from gzip import open as gzip_open
from io import StringIO
from itertools import product, cycle
from multiprocess import Queue, Process
from os import W_OK, access, mkdir
from os.path import join, exists, isdir, dirname
from pandas import concat
from shutil import rmtree
from sortedcontainers import SortedListWithKey
from subprocess import call
from tempfile import mkdtemp
from .config import GACONF
from .descriptors import DescriptorsChain
from .descriptors.cxcalc import Pkab
from .descriptors.descriptorsdict import DescriptorsDict
from .descriptors.eed import Eed
from .descriptors.fragmentor import Fragmentor
from .estimators.svmodel import SVModel
from .mbparser import MBparser


def save_svm(outputfile, x, y, header=True):
    with open(outputfile + '.svm', 'w', encoding='utf-8') as f:
        if header:
            f.write(' '.join(['Property'] + ['%s:%s' % i for i in enumerate(x.columns, start=1)]) + '\n')

        for i, j in zip(x.values, y):
            f.write(' '.join(['%s ' % j] + ['%s:%s' % x for x in enumerate(i, start=1) if x[1] != 0]) + '\n')


def save_csv(outputfile, x, y, header=True):
    concat([y, x], axis=1).to_csv(outputfile + '.csv', index=False, header=header)


def desc_starter(gen, file_dump, out_file, fformat, header, is_reaction):
    with StringIO(file_dump) as f:
        inp = (RDFread(f) if is_reaction else SDFread(f)).read()

    dsc = gen.get(structures=inp, in_structures=True)

    if dsc and not (dsc['X'].isnull().values.any() or dsc['Y'].isnull().any()):
        (save_svm if fformat == 'svm' else save_csv)(out_file, dsc['X'], dsc['Y'], header=header)
        return True

    print('BAD Descriptor generator params')
    return False


def worker(input_queue, output_queue):
    for args in iter(input_queue.get, 'STOP'):
        result = desc_starter(*args[1:])
        output_queue.put((args[0], result))


class ModelBuilder(MBparser):
    def __init__(self, description, workpath='.', is_reaction=True, model=None, reload=None, resume=None, output=None,
                 out_format='svm', fragments=None, extension=None, eed=None, pka=None, chains=None, ad=None,
                 estimator='svr', svm=None, max_iter=100000, probability=False, fit='rmse', dispcoef=0,
                 n_jobs=2, nfold=5, repetition=1, rep_boost=25, consensus=10, normalize=False, ga_maxconfigs=3000,
                 scorers=('rmse', 'r2')):
        clean_descgens = False
        if not (exists(workpath) and access(workpath, W_OK)):
            raise Exception('work path not writable')

        if not output:
            if not exists(description) or isdir(description):
                raise Exception('path to model description file invalid')

            description = self.parse_model_description(description)
            description['type'] = 2 if is_reaction else 1
            if any(x not in description for x in ('nlim', 'name', 'example', 'description')):
                raise Exception('Description Invalid')

            print('model description loaded\nname: {name}\ndescription: {description}\n'
                  'nlim: {nlim}'.format(**description))

            if not model:
                raise Exception('path to model save file invalid')

            model_save = '%s.save' % model
            if isdir(model) or not access(dirname(model) or '.', W_OK) or (exists(model) and not access(model, W_OK)) \
                    or not reload and (isdir(model_save) or (exists(model_save) and not access(model_save, W_OK))):
                raise Exception('path for model saving not writable or incorrect')

            if reload and (not exists(reload) or isdir(reload)):
                raise Exception('reload file path incorrect')

            if not reload and resume and not isdir(join(resume, 'work')):
                raise Exception('resume dir path incorrect')

            if reload:
                tmp = load(gzip_open(reload, 'rb'))
                clean_descgens = tmp.pop('descgens')
                print('reloaded save')
                if 'svm' in tmp:  # for svm. todo: for rf etc.
                    self.__svm = tmp['svm']
                    print('found SVM params in save')

            self.__model = model
            self.__model_save = model_save
        else:
            if not access(dirname(output) or '.', W_OK) or \
                    (exists('%s.1.svm' % output) and not access('%s.1.svm' % output, W_OK)) or \
                    (exists('%s.1.csv' % output) and not access('%s.1.csv' % output, W_OK)):
                raise Exception('path for descriptors saving not writable or incorrect')

            self.__output = output
            self.__format = out_format

        self.__generators = self.__descriptors_config(is_reaction, workpath, fragments=fragments, extension=extension,
                                                      eed=eed, pka=pka, chains=chains, ad=ad)

        if clean_descgens:
            self.__clean_desc_gens(clean_descgens)
        else:
            if svm:  # for svm. todo: for rf etc.
                self.__svm = self.__chk_est(self.get_svm_param(svm))
                print('SVM params loaded')

        self.__workpath = workpath
        self.__estimator = estimator
        self.__description = description
        self.__max_iter = max_iter
        self.__probability = probability
        self.__fit = fit
        self.__disp_coef = dispcoef
        self.__is_reaction = is_reaction
        self.__n_jobs = n_jobs
        self.__nfold = nfold
        self.__repetition = repetition
        self.__rep_boost = rep_boost
        self.__consensus = consensus
        self.__normalize = normalize
        self.__ga_maxconfigs = ga_maxconfigs
        self.__scorers = scorers
        self.__estimators = []
        self.__resume = join(resume, 'work') if resume else None

    __output = False
    __svm = None

    def run(self, input_file):
        if not self.__output:
            self.prepare_estimators(input_file)
            self.fit(input_file)
        else:
            self.__gen_desc(input_file, self.__output, fformat=self.__format, header=True)

            for dgen in self.__generators:
                if hasattr(dgen, 'delete_work_path'):
                    dgen.delete_work_path()

    def prepare_estimators(self, input_file):
        svm = {'svr', 'svc'}.intersection(self.__estimator).pop()

        if svm:
            if not self.__svm:
                svm_params, cleared = self.__parse_dragos_results(self.__resume) if self.__resume else \
                    self.__dragos_svm_fit(input_file, svm)
                self.__clean_desc_gens(cleared)
                self.__svm = self.__chk_est(svm_params)
                print('SVM params loaded')

                for_save = dict(descgens=cleared, svm=self.__svm)
                dump(for_save, gzip_open(self.__model_save, 'wb'))
                print('configuration saved')

            self.__estimators.append((partial(SVModel, estimator=svm, probability=self.__probability,
                                      max_iter=self.__max_iter), self.__svm))

    @staticmethod
    def __descriptors_config(is_reaction, workpath, fragments=None, extension=None, eed=None, pka=None,
                             chains=None, ad=None):
        descgenerator = OrderedDict()

        def s_choice(params):
            if s_option:
                params['s_option'] = s_option
            return params

        if extension:
            parsed = MBparser.parse_ext(extension)
            descgenerator['E'] = [partial(DescriptorsDict, **parsed)]
            s_option = parsed['s_option']
        else:
            s_option = None

        if fragments:
            descgenerator['F'] = [partial(Fragmentor, is_reaction=is_reaction, workpath=workpath, **s_choice(x))
                                  for x in MBparser.parse_fragmentor_opts(fragments)]

        if eed:
            descgenerator['D'] = [partial(Eed, is_reaction=is_reaction, workpath=workpath, **s_choice(x))
                                  for x in MBparser.parse_fragmentor_opts(eed)]

        if pka:
            descgenerator['P'] = [partial(Pkab, is_reaction=is_reaction, workpath=workpath, **s_choice(x))
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

        print('loaded %d descriptors generators' % len(descgens))
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
                model = g(x, list(y.values()), inp, in_structures=True, dispcoef=self.__disp_coef, fit=self.__fit,
                          scorers=self.__scorers, n_jobs=self.__n_jobs, nfold=self.__nfold,
                          rep_boost=self.__rep_boost, repetitions=self.__repetition,
                          normalize='scale' in y or self.__normalize)

                model.delete_work_path()
                models.add(model)
                if len(models) > self.__consensus:
                    models.pop()

        if 'tol' not in self.__description:
            self.__description['tol'] = models[0].get_model_stats()['dragostolerance']

        dump(dict(models=models, config=self.__description), gzip_open(self.__model, 'wb'))

    def __chk_est(self, est_params):
        if not est_params or 1 < len(est_params) < len(self.__generators) or len(est_params) > len(self.__generators):
            print('NUMBER of estimator params files SHOULD BE EQUAL to '
                  'number of descriptor generator params files or to 1')
            raise Exception('SVM estimators not configured')

        if len(est_params) == 1:
            tmp = []
            for i in range(len(self.__generators)):
                tmp.append(deepcopy(est_params[0]))
            est_params = tmp

        return est_params

    def __clean_desc_gens(self, select):
        tmp = [self.__generators[x] for x in select]

        for n in set(range(len(self.__generators))).difference(select):
            gen = self.__generators[n]
            if hasattr(gen, 'delete_work_path'):
                gen.delete_work_path()

        self.__generators = tmp

    def __gen_desc(self, input_file, output, fformat='svm', header=False):
        with open(input_file) as f:
            data = f.read()

        task_queue = Queue()
        done_queue = Queue()

        for i in range(self.__n_jobs):
            Process(target=worker, args=(task_queue, done_queue)).start()
        print('workers started')

        for n, dgen in enumerate(self.__generators, start=1):
            task_queue.put([n, dgen, data, '%s.%d' % (output, n), fformat, header, self.__is_reaction])
        print('generators sent to queue')

        print('unordered results of descriptors generation:')
        for i in range(len(self.__generators)):
            res = done_queue.get()
            print('\t%d: %s' % res)
            if not res[1]:
                return False

        # Tell child processes to stop
        for i in range(self.__n_jobs):
            task_queue.put('STOP')

        return True

    def __dragos_svm_fit(self, input_file, _type):
        """ files - basename for descriptors.
        """
        workpath = mkdtemp(prefix='gac_', dir=self.__workpath)
        files = join(workpath, 'drag')
        dragos_work = join(workpath, 'work')
        execparams = [GACONF, workpath, _type, str(self.__ga_maxconfigs), str(self.__repetition), str(self.__nfold),
                      str(self.__n_jobs)]

        print('descriptors generation for GAConf')
        if self.__gen_desc(input_file, files):
            print('GAConf exec:', execparams)
            if call(execparams) == 0:
                return self.__parse_dragos_results(dragos_work)

        raise Exception('GAConf failed')

    def __parse_dragos_results(self, dragos_work):
        cleared, svm = [], OrderedDict()
        with open(join(dragos_work, 'best_pop')) as f:
            for line in f:
                if len(cleared) == self.__consensus:
                    break
                dset, normal, *_, attempt, _, _ = line.split()
                parsed = list(self.get_svm_param([join(dragos_work, attempt, 'svm.pars')])[0].values())[0]
                if (parsed['kernel'], dset) not in svm:
                    svm[(parsed['kernel'], dset)] = {'scale' if normal == 'scaled' else 'orig': parsed}
                    cleared.append(int(dset[5:]) - 1)

        if not svm:
            raise Exception('Inalid best_pop file')

        print('GAConf results parsed')
        return list(svm.values()), cleared
