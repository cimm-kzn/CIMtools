CIMtools
=========
ChemoInformatics Modeling tools distributive.  
include ISIDA Fragmentor python wrapper


INSTALL
=======

Linux Debian based
------------------

* Install python3.7, virtualenv and git

    ```
    sudo apt install python3.7 python3.7-dev git python3-virtualenv
    ```
    
* Create new environment and activate it.

    ```
    virtualenv -p python3.7 venv
    source venv/bin/activate
    ```

Mac
---
* Install python3.7 and git using [brew](<https://brew.sh>)

    ```
    brew install git
    brew install python3
    ```
    
* Install virtualenv.

    ```
    pip install virtualenv
    ```

* Create new environment and activate it.

    ```
    virtualenv -p python3.7 venv
    source venv/bin/activate
    ```
    
Windows
-------

* Install python3.7 and git using [Chocolatey](<https://chocolatey.org/>)

    ```
    choco install git
    choco install python3
    ```
    
* Install virtualenv.

    ```
    pip install virtualenv
    ```

* Create new environment and activate it.

    ```
    virtualenv venv
    venv\Scripts\activate
    ```

General part
------------

* **stable version will be available through PyPI**

    ```
    pip install CIMtools
    ```    

* Install CGRtools library DEV version for features that are not well tested

    ```
    pip install -U git+https://github.com/stsouko/CIMtools.git@master#egg=CIMtools
    ```

**If you still have questions, please open issue within github.**

SETUP
=====

* add CHEMAXON_REST environment variable (optional)

        export CHEMAXON_REST="url/to/webservices" [in BASH]

PACKAGING
=========

For wheel generation just type next command in source root

    python setup.py bdist_wheel

On Linux additionally do repairing of package

    pip install auditwheel
    auditwheel repair dist/CIMtools-<version>-<python_version>-linux_x86_64.whl

COPYRIGHT
=========

2015-2020 Ramil Nugmanov <nougmanoff@protonmail.com> main developer   

CONTRIBUTORS
============

* Assima Rakhimbekova <asima.astana@outlook.com>
* Tagir Akhmetshin <tagirshin@gmail.com>
