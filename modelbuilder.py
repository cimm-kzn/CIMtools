# -*- coding: utf-8 -*-
#
# Copyright 2015, 2016 Ramil Nugmanov <stsouko@live.ru>
# This file is part of PREDICTOR.
#
# PREDICTOR is free software; you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
import argparse
import dill
import gzip
import os
import subprocess as sp
import sys
import tempfile
import threading
import time
from copy import deepcopy
from functools import partial
from itertools import product, cycle
from sortedcontainers import SortedListWithKey
from CGRtools.CGRcore import CGRcore
from CGRtools.FEAR import FEAR
from CGRtools.files.RDFrw import RDFread
from MODtools.config import GACONF
from MODtools.descriptors.cxcalc import Pkab
from MODtools.descriptors.descriptoragregator import Descriptorsdict, Descriptorchain
from MODtools.descriptors.eed import Eed
from MODtools.descriptors.fragmentor import Fragmentor
from MODtools.estimators.svmodel import SVModel
from MODtools.parsers import MBparser


class DefaultList(list):
    @staticmethod
    def __copy__(*_):
        return []


def descstarter(func, in_file, out_file, fformat, header):
    with open(in_file) as f:
        dsc = func(structures=f, parsesdf=True)
        if dsc:
            fformat(out_file, dsc['X'], dsc['Y'], header=header)
        else:
            print('BAD Descriptor generator params')
            return False
    return True


class Modelbuilder(MBparser):
    def __init__(self, **kwargs):
        self.__options = kwargs

        """ Descriptor generator Block
        """
        descgenerator = {}
        if self.__options['fragments']:
            descgenerator['F'] = [partial(Fragmentor, is_reaction=self.__options['isreaction'], **x)
                                  for x in self.parsefragmentoropts(self.__options['fragments'])]

        if self.__options['extention']:
            descgenerator['E'] = [partial(Descriptorsdict, is_reaction=self.__options['isreaction'],
                                          **self.parseext(self.__options['extention']))]

        if self.__options['eed']:
            descgenerator['D'] = [partial(Eed, is_reaction=self.__options['isreaction'], **x)
                                  for x in self.parsefragmentoropts(self.__options['eed'])]

        if self.__options['pka']:
            descgenerator['P'] = [partial(Pkab, is_reaction=self.__options['isreaction'], **x)
                                  for x in self.parsefragmentoropts(self.__options['pka'])]

        if self.__options['chains']:
            if self.__options['ad'] and len(self.__options['ad']) != len(self.__options['chains']):
                print('number of generators chains should be equal to number of ad modifiers')
                return

            self.__descgens = []
            for ch, ad in zip(self.__options['chains'], self.__options['ad'] or cycle([None])):
                gen_chain = [x for x in ch.split(':') if x in descgenerator]
                if ad:
                    ad_marks = [x in ('y', 'Y', '1', 'True', 'true') for x in ad.split(':')]
                    if len(ad_marks) != len(gen_chain):
                        print('length of generators chain should be equal to length of ad modifier')
                        return
                else:
                    ad_marks = cycle([True])

                ad_chain = {}
                for k, v in zip(gen_chain, ad_marks):
                    ad_chain.setdefault(k, []).append(v)

                combo = []
                for k, v in ad_chain.items():
                    if len(v) > 1:
                        if len(descgenerator[k]) != len(v):
                            print('length of same generators chain should be equal to number of same generators')
                            return
                        combo.append([list(zip(descgenerator[k], v))])
                    else:
                        combo.append(list(zip(descgenerator[k], cycle(v))))

                self.__descgens.extend(
                    [Descriptorchain(*[(g(), a) for gs in c
                                       for g, a in (gs if isinstance(gs, list) else [gs])]) for c in product(*combo)])
        else:
            self.__descgens = [g() for x in descgenerator.values() for g in x]

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

            if self.__options['reload']:
                ests, description, self.__descgens = dill.load(gzip.open(self.__options['reload'], 'rb'))
            else:
                if not os.path.exists(self.__options['description']) or os.path.isdir(self.__options['description']):
                    print('path to model description file invalid')
                    return

                description = self.parsemodeldescription(self.__options['description'])

                description['type'] = 2 if self.__options['isreaction'] else 1

                ests = []
                svm = {'svr', 'svc'}.intersection(self.__options['estimator']).pop()
                # rf = {'rf'}.intersection(self.__options['estimator']).pop()
                if svm:
                    if self.__options['svm']:
                        estparams = self.getsvmparam(self.__options['svm'])
                    else:
                        estparams = self.__dragossvmfit(svm)

                    estparams = self.__chkest(estparams)
                    if not estparams:
                        return
                    ests.append((partial(SVModel, estimator=svm, probability=self.__options['probability']),
                                 estparams))
                elif False:  # rf:  # todo: not implemented
                    if self.__options['rf']:
                        estparams = None
                        estparams = self.__chkest(estparams)
                        if not estparams:
                            ests.append((lambda *va, **kwa: None, estparams))
                    else:
                        return

                if not ests:
                    return

                dill.dump((ests, description, self.__descgens), gzip.open(self.__options['model'] + '.save', 'wb'))

            self.fit(ests, description)
        else:
            self.__gendesc(self.__options['output'], fformat=self.__options['format'], header=True)

    def __order(self, model):
        s = (1 if self.__options['fit'] == 'rmse' else -1) * model.getmodelstats()[self.__options['fit']]
        v = self.__options['dispcoef'] * model.getmodelstats()['%s_var' % self.__options['fit']]
        return s + v

    def fit(self, ests, description):
        models = SortedListWithKey(key=self.__order)
        for g, e in ests:
            for x, y in zip(self.__descgens, e):
                models.add(g(x, list(y.values()), open(self.__options['input']), parsesdf=True,
                           dispcoef=self.__options['dispcoef'], fit=self.__options['fit'],
                           scorers=self.__options['scorers'],
                           n_jobs=self.__options['n_jobs'], nfold=self.__options['nfold'],
                           rep_boost=self.__options['rep_boost'], repetitions=self.__options['repetition'],
                           normalize='scale' in y or self.__options['normalize']))

                if len(models) > self.__options['consensus']:
                    models.pop()

        if 'tol' not in description:
            description['tol'] = models[0].getmodelstats()['dragostolerance']
        print('name', description['name'])
        print('description', description['description'])
        print('tol', description['tol'])
        print('nlim', description.get('nlim'))
        dill.dump(dict(models=models, config=description), gzip.open(self.__options['model'], 'wb'))

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
                                               (self.savesvm if fformat == 'svm' else self.savecsv), header])
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
                    for line in f:
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


def argparser():
    rawopts = argparse.ArgumentParser(description="Model Builder",
                                      epilog="Copyright 2015, 2016 Ramil Nugmanov <stsouko@live.ru>",
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    rawopts.add_argument("--workpath", "-w", type=str, default='.', help="work path")

    rawopts.add_argument("--input", "-i", type=str, default='input.sdf', help="input SDF or RDF")
    rawopts.add_argument("--output", "-o", type=str, default=None, help="output SVM|CSV")
    rawopts.add_argument("--format", "-of", type=str, default='svm', choices=['svm', 'csv'], help="output format")

    rawopts.add_argument("--reload", type=str, default=None, help="saved state before fitting")

    rawopts.add_argument("--model", "-m", type=str, default='output.model', help="output model")

    rawopts.add_argument("--isreaction", "-ir", action='store_true', help="set as reaction model")

    rawopts.add_argument("--extention", "-e", action='append', type=str, default=None,
                         help="extention data files. -e extname:filename [-e extname2:filename2]")

    rawopts.add_argument("--fragments", "-f", type=str, default=None, help="ISIDA Fragmentor keys file")

    rawopts.add_argument("--eed", type=str, default=None, help="DRAGOS EED keys file")

    rawopts.add_argument("--pka", type=str, default=None, help="CXCALC pka keys file")

    rawopts.add_argument("--chains", "-c", action='append', type=str, default=None,
                         help="descriptors chains. where F-fragmentor, D-eed, E-extention, P-pka. "
                              "-c F:E [-c E:D:P]")
    rawopts.add_argument("-ad", action='append', type=str, default=None,
                         help="consider descriptor generator AD in descriptors chains. "
                              "example: -ad y:n:y [True = y Y True true 1] for -c F:E:P [ignore extension AD] "
                              "number of -ad should be equal to number of --chains or skipped")

    rawopts.add_argument("--description", "-ds", type=str, default='model.dsc', help="model description file")

    rawopts.add_argument("--svm", "-s", action='append', type=str, default=None,
                         help="SVM params. use Dragos Genetics if don't set."
                              "can be multiple [-s 1 -s 2 ...]"
                              "(number of files should be equal to number of configured descriptor generators) "
                              "or single for all")

    rawopts.add_argument("--nfold", "-n", type=int, default=5, help="number of folds")
    rawopts.add_argument("--repetition", "-r", type=int, default=1, help="number of repetitions")
    rawopts.add_argument("--rep_boost", "-R", type=int, default=25,
                         help="percentage of repetitions for use in greed search for optimization speedup")
    rawopts.add_argument("--n_jobs", "-j", type=int, default=2, help="number of parallel fit jobs")

    rawopts.add_argument("--estimator", "-E", action='append', type=str, default=DefaultList(['svr']),
                         choices=['svr', 'svc'],
                         help="estimator")
    rawopts.add_argument("--probability", "-P", action='store_true', help="estimates probability in SVC")

    rawopts.add_argument("--scorers", "-T", action='append', type=str, default=DefaultList(['rmse', 'r2']),
                         choices=['rmse', 'r2', 'ba', 'kappa', 'iap'],
                         help="needed scoring functions. -T rmse [-T r2]")
    rawopts.add_argument("--fit", "-t", type=str, default='rmse', choices=['rmse', 'r2', 'ba', 'kappa', 'iap'],
                         help="crossval score for parameters fit. (should be in selected scorers)")

    rawopts.add_argument("--dispcoef", "-p", type=float, default=0,
                         help="score parameter. mean(score) - dispcoef * sqrt(variance(score)). [-score for rmse]")

    rawopts.add_argument("--normalize", "-N", action='store_true', help="normalize X vector to range(0, 1)")

    rawopts.add_argument("--consensus", "-C", type=int, default=10, help="number of models for consensus")

    return vars(rawopts.parse_args())


if __name__ == '__main__':
    main = Modelbuilder(**argparser())
