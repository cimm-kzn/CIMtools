#!/usr/bin/env python3.4
# -*- coding: utf-8 -*-
#
#  Copyright 2016 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of MODtools.
#
#  MODtools 
#  is free software; you can redistribute it and/or modify
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
import dill
import gzip
import os
import subprocess as sp
import tempfile
import threading
import time
from copy import deepcopy
from collections import OrderedDict
from functools import partial
from itertools import product, cycle
from sortedcontainers import SortedListWithKey
from CGRtools.files.RDFrw import RDFread
from CGRtools.files.SDFrw import SDFread
from .config import GACONF
from .descriptors.cxcalc import Pkab
from .descriptors.descriptoragregator import Descriptorsdict, Descriptorchain
from .descriptors.eed import Eed
from .descriptors.fragmentor import Fragmentor
from .estimators.svmodel import SVModel
from .parsers import MBparser


def descstarter(func, in_file, out_file, fformat, header, is_reaction):
    with open(in_file) as f:
        inp = list((RDFread(f) if is_reaction else SDFread(f)).read())

    dsc = func(structures=inp, parsesdf=True)
    if dsc:
        fformat(out_file, dsc['X'], dsc['Y'], header=header)
        return True

    print('BAD Descriptor generator params')
    return False


class Modelbuilder(MBparser):
    def __init__(self, **kwargs):
        self.__options = kwargs

        if not self.__options['output'] and (not os.path.exists(self.__options['description'])
                                             or os.path.isdir(self.__options['description'])):
            raise Exception('path to model description file invalid')

        description = self.parsemodeldescription(self.__options['description'])

        if any(x not in description for x in ('nlim', 'name', 'example', 'description')):
            raise Exception('Description Invalid')

        description['type'] = 2 if self.__options['isreaction'] else 1

        self.__description = description
        self.__descriptors_config()

    def run(self):
        if not self.__options['output']:
            if os.path.isdir(self.__options['model']) or \
               (os.path.exists(self.__options['model']) and not os.access(self.__options['model'], os.W_OK)) or \
               os.path.isdir(self.__options['model'] + '.save') or \
               (os.path.exists(self.__options['model'] + '.save') and
                    not (os.access(self.__options['model'] + '.save', os.W_OK) or
                         self.__options['model'] + '.save' == self.__options['reload'])) or \
               not os.access(os.path.dirname(self.__options['model']), os.W_OK):
                print('path for model saving not writable')
                return
            self.prepare_estimators()
            self.fit()
        else:
            self.__gendesc(self.__options['output'], fformat=self.__options['format'], header=True)

    def prepare_estimators(self):
        if self.__options['reload']:
            ests = dill.load(gzip.open(self.__options['reload'], 'rb'))
        else:
            ests = []

            svm = {'svr', 'svc'}.intersection(self.__options['estimator']).pop()
            if svm:
                estparams = self.__chkest(self.getsvmparam(self.__options['svm'])
                                          if self.__options['svm'] else self.__dragossvmfit(svm))

                if estparams:
                    ests.append((partial(SVModel, estimator=svm, probability=self.__options['probability'],
                                         max_iter=self.__options['max_iter']), estparams))
            # todo: implement RF

            dill.dump(ests, gzip.open(self.__options['model'] + '.save', 'wb'))
        if not ests:
            raise Exception('Estimators not configured')

        self.__estimators = ests

    def __descriptors_config(self):
        descgenerator = OrderedDict()
        if self.__options['fragments']:
            descgenerator['F'] = [partial(Fragmentor, is_reaction=self.__options['isreaction'], **x)
                                  for x in self.parsefragmentoropts(self.__options['fragments'])]

        if self.__options['extention']:
            descgenerator['E'] = [partial(Descriptorsdict, **self.parseext(self.__options['extention']))]

        if self.__options['eed']:
            descgenerator['D'] = [partial(Eed, is_reaction=self.__options['isreaction'], **x)
                                  for x in self.parsefragmentoropts(self.__options['eed'])]

        if self.__options['pka']:
            descgenerator['P'] = [partial(Pkab, is_reaction=self.__options['isreaction'], **x)
                                  for x in self.parsefragmentoropts(self.__options['pka'])]

        if self.__options['chains']:
            if self.__options['ad'] and len(self.__options['ad']) != len(self.__options['chains']):
                raise Exception('number of generators chains should be equal to number of ad modifiers')

            descgens = []
            for ch, ad in zip(self.__options['chains'], self.__options['ad'] or cycle([None])):
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
                descgens.extend([Descriptorchain(*[(g(), a) for gs in c
                                 for g, a in (gs if isinstance(gs, list) else [gs])]) for c in product(*combo)])
        else:
            descgens = [g() for x in descgenerator.values() for g in x]

        if not descgens:
            raise Exception('Descriptor generators not configured')

        self.__descgens = descgens

    def __order(self, model):
        s = (1 if self.__options['fit'] == 'rmse' else -1) * model.getmodelstats()[self.__options['fit']]
        v = self.__options['dispcoef'] * model.getmodelstats()['%s_var' % self.__options['fit']]
        return s + v

    def fit(self):
        models = SortedListWithKey(key=self.__order)
        with open(self.__options['input']) as f:
            inp = list((RDFread(f) if self.__options['isreaction'] else SDFread(f)).read())

        for g, e in self.__estimators:
            for x, y in zip(self.__descgens, e):
                models.add(g(x, list(y.values()), inp, parsesdf=True,
                           dispcoef=self.__options['dispcoef'], fit=self.__options['fit'],
                           scorers=self.__options['scorers'],
                           n_jobs=self.__options['n_jobs'], nfold=self.__options['nfold'],
                           rep_boost=self.__options['rep_boost'], repetitions=self.__options['repetition'],
                           normalize='scale' in y or self.__options['normalize']))

                if len(models) > self.__options['consensus']:
                    models.pop()

        if 'tol' not in self.__description:
            self.__description['tol'] = models[0].getmodelstats()['dragostolerance']

        print('name', self.__description['name'])
        print('description', self.__description['description'])
        print('tol', self.__description['tol'])
        print('nlim', self.__description.get('nlim'))
        dill.dump(dict(models=models, config=self.__description), gzip.open(self.__options['model'], 'wb'))

    def __chkest(self, estimatorparams):
        if not estimatorparams or 1 < len(estimatorparams) < len(self.__descgens) or \
                        len(estimatorparams) > len(self.__descgens):
            print('NUMBER of estimator params files SHOULD BE EQUAL to '
                  'number of descriptor generator params files or to 1')
            return False

        if len(estimatorparams) == 1:
            tmp = []
            for i in range(len(self.__descgens)):
                tmp.append(deepcopy(estimatorparams[0]))
            estimatorparams = tmp
        return estimatorparams

    def __gendesc(self, output, fformat='svm', header=False):
        queue = enumerate(self.__descgens, start=1)
        workpath = tempfile.mkdtemp(prefix='svm_', dir=self.__options['workpath'])
        while True:
            if threading.active_count() < self.__options['n_jobs']:
                tmp = next(queue, None)
                if tmp:
                    n, dgen = tmp
                    subworkpath = os.path.join(workpath, str(n))
                    os.mkdir(subworkpath)
                    dgen.setworkpath(subworkpath)
                    t = threading.Thread(target=descstarter,
                                         args=[dgen.get, self.__options['input'], '%s.%d' % (output, n),
                                               (self.savesvm if fformat == 'svm' else self.savecsv), header,
                                               self.__options['isreaction']])
                    t.start()
                else:
                    while threading.active_count() > 1:
                        time.sleep(2)
                    break
            time.sleep(2)

        return True

    def __dragossvmfit(self, tasktype):
        """ files - basename for descriptors.
        """
        workpath = tempfile.mkdtemp(prefix='gac_', dir=self.__options['workpath'])
        files = os.path.join(workpath, 'drag')
        dragos_work = os.path.join(workpath, 'work')

        execparams = [GACONF, workpath, tasktype]
        if self.__gendesc(files):
            if sp.call(execparams) == 0:
                best = {}
                with open(os.path.join(dragos_work, 'best_pop')) as f:
                    for n, line in enumerate(f):
                        if n == self.__options['best_pop']:
                            break  # get only first n's params.
                        dset, normal, *_, attempt, _, _ = line.split()
                        best.setdefault(int(dset[5:]), (normal, attempt))

                cleared, svmpar, scale = [], [], []
                for k, (nv, av) in best.items():
                    cleared.append(self.__descgens[k - 1])
                    svmpar.append(os.path.join(dragos_work, av, 'svm.pars'))
                    scale.append(nv)

                self.__descgens = cleared

                svm = []
                svmpar = self.getsvmparam(svmpar)
                if len(svmpar) == len(scale):
                    for x, y in zip(svmpar, scale):
                        svm.append({'scale' if y == 'scaled' else 'orig': list(x.values())[0]})
                    return svm
        return []
