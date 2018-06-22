CIMtools
=========
ChemoInformatics Modeler tools distributive. include ISIDA Fragmentor and EED python wrappers

usage
=====

cimtools -h


INSTALL
=======

*step 1:*

    pip install CIMtools

or latest repository version

    pip install -U git+https://github.com/stsouko/CIMtools.git@master#egg=CIMtools --process-dependency-links --allow-all-external

*step 2:*

type

    cimtools -h

*step 3:*

edit `~/.CIMtools.ini` or add `/etc/.CIMtools.ini`
or add file into package `CIMtools/.CIMtools.ini` [useful for manually downloaded CIMtools Lib]

    UTILS_DIR=~ # dir with Dragos's utilities (Utils dir should be in this dir)
    CHEMAXON=https://cimm.kpfu.ru/webservices # Chemaxon REST API
    JCHEM_DIR = ~/ChemAxon/JChem # dir with JChem distributive
    FRAGMENTOR = ~/fragmentor/fragmentor # path to fragmentor bin file (without -version_suffix [fragmentor-2015.22])


priority: `CIMtools/.CIMtools.ini` >> `~/.CIMtools.ini` >> `/etc/.CIMtools.ini`
