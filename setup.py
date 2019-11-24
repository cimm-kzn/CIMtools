# -*- coding: utf-8 -*-
#
#  Copyright 2016-2019 Ramil Nugmanov <stsouko@live.ru>
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
from distutils.util import get_platform
from pathlib import Path
from setuptools import setup, find_packages
from wheel.bdist_wheel import bdist_wheel


version = '4.0.1'

platform = get_platform()
if platform == 'win-amd64':
    fragmentor = ['Fragmentor/fragmentor_win_2017.exe']
elif platform == 'linux-x86_64':
    fragmentor = ['Fragmentor/fragmentor_lin_2017']
elif platform.startswith('macosx') and platform.endswith('x86_64'):
    fragmentor = ['Fragmentor/fragmentor_mac_2017']
else:
    fragmentor = []


class _bdist_wheel(bdist_wheel):
    def finalize_options(self):
        super().finalize_options()
        self.root_is_pure = False


setup(
    name='CIMtools',
    version=version,
    packages=find_packages(),
    url='https://github.com/stsouko/CIMtools',
    license='GPLv3',
    author='Dr. Ramil Nugmanov',
    author_email='stsouko@live.ru',
    python_requires='>=3.6.1',
    cmdclass={'bdist_wheel': _bdist_wheel},
    install_requires=['CGRtools[mrv]>=4.0,<4.1', 'pandas>=0.22.0,<0.26', 'scikit-learn>=0.20.1,<0.22',
                      'requests>=2.21,<2.23', 'pyparsing>=2.2.0,<2.5', 'numpy>=1.15.0,<1.18'],
    data_files=[('bin', fragmentor)],
    zip_safe=False,
    long_description=(Path(__file__).parent / 'README.md').open().read(),
    classifiers=['Environment :: Plugins',
                 'Intended Audience :: Science/Research',
                 'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python',
                 'Programming Language :: Python :: 3 :: Only',
                 'Programming Language :: Python :: 3.6',
                 'Programming Language :: Python :: 3.7',
                 'Topic :: Scientific/Engineering',
                 'Topic :: Scientific/Engineering :: Chemistry',
                 'Topic :: Scientific/Engineering :: Information Analysis',
                 'Topic :: Software Development',
                 'Topic :: Software Development :: Libraries',
                 'Topic :: Software Development :: Libraries :: Python Modules']
)
