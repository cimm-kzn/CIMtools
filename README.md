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

edit `~/.MODtools.ini` or add `/etc/.MODtools.ini`
or add file into package `MODtools/.MODtools.ini` [useful for manually downloaded MODtools Lib]

    GACONF_PATH=~/GAconfig # GAconfig dir
    LIBSVM_PATH=~/GAconfig/libsvm-3.20 # LIBSVM binaries dir
    UTILS_DIR=~ # dir with Dragos's utilities (Utils dir should be in this dir)
    CHEMAXON=https://cimm.kpfu.ru/webservices # Chemaxon REST API
    JCHEM_DIR = ~/ChemAxon/JChem # dir with JChem distributive
    FRAGMENTOR = ~/fragmentor/fragmentor # path to fragmentor bin file (without -version_suffix [fragmentor-2015.22])


priority: `MODtools/.MODtools.ini` >> `~/.MODtools.ini` >> `/etc/.MODtools.ini`

*NOTE:*

`MODtools/.MODtools.ini` not visible for `colorstart.sh` `dragosgfstarter.sh` `eedstart.sh`

Edit this *.sh scripts manually and move it to one of the $PATH dirs
