# -*- coding: utf-8 -*-
#
#  Copyright 2016-2018 Ramil Nugmanov <stsouko@live.ru>
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


version = '3.0.1'


setup(
    name='CIMtools',
    version=version,
    packages=find_packages(),
    zip_safe=True,
    url='https://github.com/stsouko/CIMtools',
    license='GPLv3',
    author='Dr. Ramil Nugmanov',
    author_email='stsouko@live.ru',
    package_data={'CIMtools.preprocessing.standardize': ['horvat.unwanted', 'horvat.xml']},
    install_requires=['CGRtools>=3.0.3,<3.1', 'pandas>=0.22.0,<0.24', 'scikit-learn>=0.20.1,<0.21',
                      'requests>=2.18.4,<2.20', 'scipy>=1.0.0,<1.2', 'pyparsing>=2.2.0,<2.3'],
    long_description=(Path(__file__).parent / 'README.md').open().read(),
    keywords="chemoinformatics tools modeler cli ISIDA Framentor EED SVM IAP",
    classifiers=['Environment :: Plugins',
                 'Intended Audience :: Science/Research',
                 'Intended Audience :: Developers',
                 'Topic :: Scientific/Engineering :: Chemistry',
                 'Topic :: Software Development :: Libraries :: Python Modules',
                 'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python',
                 'Programming Language :: Python :: 3',
                 'Programming Language :: Python :: 3.7',
                 ],
    command_options={'build_sphinx': {'project': ('setup.py', 'CIMtools'),
                                      'version': ('setup.py', version), 'source_dir': ('setup.py', 'doc'),
                                      'build_dir':  ('setup.py', 'build/doc'),
                                      'all_files': ('setup.py', True),
                                      'copyright': ('setup.py', 'Dr. Ramil Nugmanov <stsouko@live.ru>')}}
)
