# -*- coding: utf-8 -*-
#
#  Copyright 2016-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
#  This file is part of CIMtools.
#
#  CIMtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
from distutils.command.sdist import sdist
from distutils.util import get_platform
from importlib.util import find_spec
from pathlib import Path
from setuptools import setup, find_packages


class _sdist(sdist):
    def finalize_options(self):
        super().finalize_options()
        self.distribution.data_files.append(('bin', ['Fragmentor/fragmentor_win_2017.exe',
                                                     'Fragmentor/fragmentor_lin_2017',
                                                     'Fragmentor/fragmentor_mac_2017']))


cmd_class = {'sdist': _sdist}


if find_spec('wheel'):
    from wheel.bdist_wheel import bdist_wheel

    class _bdist_wheel(bdist_wheel):
        def finalize_options(self):
            super().finalize_options()
            self.root_is_pure = False
            platform = get_platform()
            if platform == 'win-amd64':
                self.distribution.data_files.append(('bin', ['Fragmentor/fragmentor_win_2017.exe']))
            elif platform == 'linux-x86_64':
                self.distribution.data_files.append(('bin', ['Fragmentor/fragmentor_lin_2017']))
            elif platform.startswith('macosx') and platform.endswith('x86_64'):
                self.distribution.data_files.append(('bin', ['Fragmentor/fragmentor_mac_2017']))

    cmd_class['bdist_wheel'] = _bdist_wheel


setup(
    name='CIMtools',
    version='4.0.15',
    packages=find_packages(),
    url='https://github.com/cimm-kzn/CIMtools',
    license='GPLv3',
    author='Dr. Ramil Nugmanov',
    author_email='nougmanoff@protonmail.com',
    python_requires='>=3.6.1',
    cmdclass=cmd_class,
    install_requires=['CGRtools>=4.0,<4.2', 'pandas>=0.22', 'scikit-learn>=0.24',
                      'pyparsing>=2.2', 'pyjnius>=1.3.0', 'StructureFingerprint'],
    extras_require={'gnnfp': ['tensorflow>=2.2.0']},
    package_data={'CIMtools.preprocessing.graph_encoder': ['weights.h5'],
                  'CIMtools.datasets': ['data/*.rdf', 'data/tautomer_database_release_3a.xlsx']},
    data_files=[('lib', ['RDtool/rdtool.jar'])],
    zip_safe=False,
    long_description=(Path(__file__).parent / 'README.rst').open().read(),
    classifiers=['Environment :: Plugins',
                 'Intended Audience :: Science/Research',
                 'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python',
                 'Programming Language :: Python :: 3 :: Only',
                 'Programming Language :: Python :: 3.6',
                 'Programming Language :: Python :: 3.7',
                 'Programming Language :: Python :: 3.8',
                 'Programming Language :: Python :: 3.9',
                 'Programming Language :: Python :: 3.10',
                 'Programming Language :: Python :: 3.11',
                 'Programming Language :: Python :: 3.12',
                 'Programming Language :: Python :: 3.13',                 
                 'Topic :: Scientific/Engineering',
                 'Topic :: Scientific/Engineering :: Chemistry',
                 'Topic :: Scientific/Engineering :: Information Analysis',
                 'Topic :: Software Development',
                 'Topic :: Software Development :: Libraries',
                 'Topic :: Software Development :: Libraries :: Python Modules'],
    command_options={'build_sphinx': {'source_dir': ('setup.py', 'doc'),
                                      'build_dir': ('setup.py', 'build/doc'),
                                      'all_files': ('setup.py', True)}}
)
