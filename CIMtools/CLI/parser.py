# -*- coding: utf-8 -*-
#
#  Copyright 2017 Ramil Nugmanov <stsouko@live.ru>
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
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from ..version import version


class DefaultList(list):
    @staticmethod
    def __copy__(*_):
        return []


def argparser():
    rawopts = ArgumentParser(description="Model Builder",
                             epilog="Copyright 2015-2017 Ramil Nugmanov <stsouko@live.ru>",
                             formatter_class=ArgumentDefaultsHelpFormatter)
    rawopts.add_argument("--version", "-v", action="version", version=version(), default=False)
    rawopts.add_argument("--workpath", "-w", type=str, default='.', help="work path")

    rawopts.add_argument("--input", "-i", type=str, default='input.sdf', help="input SDF or RDF")
    rawopts.add_argument("--output", "-o", type=str, default=None, help="output SVM|CSV")
    rawopts.add_argument("--format", "-of", dest='out_format', type=str, default='svm', choices=['svm', 'csv'],
                         help="output format")

    rawopts.add_argument("--reload", type=str, default=None, help="reload saved state before fitting")
    rawopts.add_argument("--resume", type=str, default=None, help="resume from GAconf done")

    rawopts.add_argument("--model", "-m", type=str, default='output.model', help="output model")

    rawopts.add_argument("--isreaction", "-ir", dest='is_reaction', action='store_true', help="set as reaction model")

    rawopts.add_argument("--extension", "-e", action='append', type=str, default=None,
                         help="extension data files. -e extname:filename [-e extname2:filename2]")

    rawopts.add_argument("--fragments", "-f", type=str, default=None, help="ISIDA Fragmentor keys file")

    rawopts.add_argument("--eed", type=str, default=None, help="DRAGOS EED keys file")

    rawopts.add_argument("--pka", type=str, default=None, help="CXCALC pka keys file")

    rawopts.add_argument("--chains", "-c", action='append', type=str, default=None,
                         help="descriptors chains. where F-fragmentor, D-eed, E-extension, P-pka. "
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

    rawopts.add_argument("--ga_maxconfigs", "-gg", type=int, default=3000,
                         help="number of generations in Dragos Genetic SVM optimizer")
    rawopts.add_argument("--nfold", "-n", type=int, default=5, help="number of folds")
    rawopts.add_argument("--repetition", "-r", type=int, default=1, help="number of repetitions")
    rawopts.add_argument("--rep_boost", "-R", type=int, default=25,
                         help="percentage of repetitions for use in greed search for optimization speedup")
    rawopts.add_argument("--n_jobs", "-j", type=int, default=2, help="number of parallel fit jobs")

    rawopts.add_argument("--estimator", "-E", action='append', type=str, default=DefaultList(['svr']),
                         choices=['svr', 'svc'],
                         help="estimator")
    rawopts.add_argument("--max_iter", "-M", type=int, default=100000, help="number of iterations in SVM solver")
    rawopts.add_argument("--probability", "-P", action='store_true', help="estimates probability in SVC")

    rawopts.add_argument("--scorers", "-T", action='append', type=str, default=DefaultList(['rmse', 'r2']),
                         choices=['rmse', 'r2', 'acc', 'kappa', 'iap'],
                         help="needed scoring functions. -T rmse [-T r2]")
    rawopts.add_argument("--fit", "-t", type=str, default='rmse', choices=['rmse', 'r2', 'acc', 'kappa', 'iap'],
                         help="crossval score for parameters fit. (should be in selected scorers)")

    rawopts.add_argument("--dispcoef", "-p", type=float, default=0,
                         help="score parameter. mean(score) - dispcoef * sqrt(variance(score)). [-score for rmse]")

    rawopts.add_argument("--normalize", "-N", action='store_true', help="normalize X vector to range(0, 1)")

    rawopts.add_argument("--consensus", "-C", type=int, default=10, help="number of models for consensus")
    rawopts.add_argument("--best_pop", "-bp", type=int, default=20, help="number of models from bet_pop")

    return vars(rawopts.parse_args())
