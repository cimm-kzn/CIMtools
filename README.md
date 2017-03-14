MODTools
=========
classes for models preparation

usage
=====
modeler -h


INSTALL
=====

*step 1:*

    pip install -U git+https://github.com/stsouko/MODtools.git@master#egg=MODtools --process-dependency-links --allow-all-external

*step 2:*

type

    modeler -h

*step 3:*

add or edit ~/.MODtools.rc or /etc/.MODtools.rc file. local *.rc has higher priority.

    GACONF=~/GAconfig
    UTILS_DIR=~
    CHMXN_DIR=~/ChemAxon

replace default paths to actual.
* GACONF - Dragos's Genetic SVM optimizer
* CHMXN_DIR - dir with JChem distributive
* UTILS_DIR - dir with Dragos's utilities

edit `~/.MODtools.ini` or add `/etc/.MODtools.ini`
or add file into package `MODtools/.MODtools.ini` [useful for unpacked MODtools Lib]

priority: `MODtools/.MODtools.ini` >> `~/.MODtools.ini` >> `/etc/.MODtools.ini`
