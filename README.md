CIMtools
=========
ChemoInformatics Modeling tools distributive.  
include ISIDA Fragmentor python wrapper


INSTALL
=======

*step 1:*

    pip install CIMtools

or latest repository version

    pip install -U git+https://github.com/stsouko/CIMtools.git@master#egg=CIMtools --process-dependency-links

*step 2:* config ENVIRONMENT

* add CHEMAXON_REST (optional)

        export CHEMAXON_REST="url/to/webservices"

* add/edit for Colorize and EED using

        CLASSPATH should contains paths to lib/jchem.jar and infochim.u-strasbg Utils containing dir
        SETUP_DIR with path/to/infochim.u-strasbg/Utils dir
        FORCEFIELD= with path/to/infochim.u-strasbg/Utils/cvffTemplates.xml or another

* edit PATH for Fragmentor using. Add path to fragmentor bin files

* edit PATH for AtomMarkerPharmacophore using. add path to jchem bin files
