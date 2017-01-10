#!/usr/bin/env python3
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
from setuptools import setup, find_packages
from MODtools.version import version

setup(
    name='modtools',
    version=version(),
    packages=find_packages(exclude=['CGRtools', 'CGRtools.cli', 'CGRtools.files', 'CGRtools.utils']),
    url='https://github.com/stsouko/MODtools',
    license='AGPLv3',
    author='Dr. Ramil Nugmanov',
    author_email='stsouko@live.ru',
    description='Modeler tools',
    scripts=['modelbuilder.py'],
    package_data={'': ['unwanted.elem', 'standardrules_dragos.rules']},
    requires=['CGRtools', 'networkx', 'periodictable', 'pandas', 'dill', 'sortedcontainers', 'sklearn', 'numpy',
              'requests'],
    long_description='Modeler tools distributive. include ISIDA Fragmentor and EED python wrappers',

    keywords="tools modeler cli ISIDA Framentor EED SVM IAP",
    classifiers=['Environment :: Console',
                 'Intended Audience :: End Users/Desktop',
                 'Intended Audience :: Developers',
                 ('License :: OSI Approved :: GNU Affero General Public License'
                  ' v3 or later (AGPLv3+)'),
                 'Operating System :: OS Independent',
                 'Programming Language :: Python',
                 'Programming Language :: Python :: 3',
                 'Programming Language :: Python :: 3.4',
                 'Programming Language :: Python :: 3.5',
                 ]
)
