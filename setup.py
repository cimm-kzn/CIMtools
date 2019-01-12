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
from pathlib import Path
from setuptools import setup, find_packages


version = '3.0.6'


setup(
    name='CIMtools',
    version=version,
    packages=find_packages(),
    zip_safe=True,
    url='https://github.com/stsouko/CIMtools',
    license='GPLv3',
    author='Dr. Ramil Nugmanov',
    author_email='stsouko@live.ru',
    python_requires='>=3.7.0',
    install_requires=['CGRtools>=3.0.10,<3.1', 'pandas>=0.22.0,<0.24', 'scikit-learn>=0.20.1,<0.21',
                      'requests>=2.21,<2.22', 'pyparsing>=2.2.0,<2.4'],
    long_description=(Path(__file__).parent / 'README.md').open().read(),
    classifiers=['Environment :: Plugins',
                 'Intended Audience :: Science/Research',
                 'Intended Audience :: Developers',
                 'Topic :: Scientific/Engineering :: Chemistry',
                 'Topic :: Software Development :: Libraries :: Python Modules',
                 'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python',
                 'Programming Language :: Python :: 3',
                 'Programming Language :: Python :: 3.7']
)
