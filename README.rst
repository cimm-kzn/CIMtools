CIMtools
=========
ChemoInformatics Modeling tools distributive.  
Include ISIDA Fragmentor python wrapper and RDtool atom-to-atom mapper python wrapper.

Documentation
=============

https://cimtools.readthedocs.io

DataSets
========
CIMtools include datasets of reactional rate constants (SN2, E2, Diels-Alder) and tautomerism equilibriums (Nicklaus dataset).

Has same API as sklearn. See tutorial for example.

    from CIMtools.datasets import load_sn2, load_e2, load_da, load_nicklaus_tautomers


INSTALL
=======

Linux Debian based
------------------

* Install python3.7, virtualenv and git::

    sudo apt install python3.7 python3.7-dev git python3-virtualenv
    
* Create new environment and activate it::

    virtualenv -p python3.7 venv
    source venv/bin/activate

Mac
---
* Install python3.7 and git using <https://brew.sh>::

    brew install git
    brew install python3

* Install virtualenv::

    pip install virtualenv

* Create new environment and activate it::

    virtualenv -p python3.7 venv
    source venv/bin/activate

Windows
-------

* Install python3.7 and git using <https://chocolatey.org/>::

    choco install git
    choco install python3
    
* Install virtualenv::

    pip install virtualenv

* Create new environment and activate it::

    virtualenv venv
    venv\Scripts\activate

General part
------------

* **stable version will be available through PyPI**::

    pip install CIMtools

* Install CGRtools library DEV version for features that are not well tested. Git lfs installation required <https://git-lfs.github.com/>::

    pip install -U git+https://github.com/cimm-kzn/CIMtools.git@master#egg=CIMtools

**If you still have questions, please open issue within github.**

SETUP
=====

For ChemAxon standardizer used pyjnius. First of all install JDK (not JRE) OpenJDK or Oracle.
Some times it can't to find java installation properly. Just set environment variables::

    JAVA_HOME = '/path/to/dir/which/contain/bin/dir'. for example /usr/lib/jvm/java-11-openjdk-amd64
    JVM_PATH = '/path/to/lib/server/libjvm.so'. For example '/usr/lib/jvm/java-11-openjdk-amd64/lib/server/libjvm.so' 

PACKAGING
=========

For wheel generation just type next command in source root::

    python setup.py bdist_wheel

COPYRIGHT
=========

2015-2021 Ramil Nugmanov <nougmanoff@protonmail.com> main developer

CONTRIBUTORS
============

* Assima Rakhimbekova <asima.astana@outlook.com>
* Tagir Akhmetshin <tagirshin@gmail.com>
* Zarina Ibragimova <zarinaIbr12@yandex.ru>
