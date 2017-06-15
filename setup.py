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
from CIMtools.version import version
from os.path import join, dirname
from setuptools import setup, find_packages

setup(
    name='CIMtools',
    version=version(),
    packages=find_packages(),
    url='https://github.com/stsouko/CIMtools',
    license='AGPLv3',
    author='Dr. Ramil Nugmanov',
    author_email='stsouko@live.ru',
    description='Modeler tools',
    entry_points={'console_scripts': ['cimtools=CIMtools.CLI:launcher']},
    scripts=['SETUP/colorstart.sh', 'SETUP/dragosgfstarter.sh', 'SETUP/eedstart.sh'],
    package_data={'CIMtools.preparers': ['unwanted.elem', 'standardrules_dragos.xml']},
    install_requires=['CGRtools>=2.7', 'pandas', 'scipy', 'sortedcontainers', 'sklearn', 'requests', 'multiprocess'],
    dependency_links=['git+https://github.com/stsouko/CGRtools.git@master#egg=CGRtools-2.7'],
    long_description=open(join(dirname(__file__), 'README.md')).read(),
    keywords="chemoinformatics tools modeler cli ISIDA Framentor EED SVM IAP",
    classifiers=['Environment :: Console',
                 'Intended Audience :: Science/Research',
                 'Intended Audience :: Developers',
                 'Topic :: Scientific/Engineering :: Chemistry',
                 'Topic :: Software Development :: Libraries :: Python Modules',
                 'License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python',
                 'Programming Language :: Python :: 3',
                 'Programming Language :: Python :: 3.4',
                 'Programming Language :: Python :: 3.5',
                 ]
)
