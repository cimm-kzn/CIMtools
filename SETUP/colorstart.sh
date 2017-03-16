#!/usr/bin/env bash

UTILS_DIR=~
JCHEM_DIR=~/ChemAxon/JChem

if [ -f /etc/.CIMtools.ini ]; then
    . /etc/.CIMtools.ini
fi

if [ -f ~/.CIMtools.ini ]; then
    . ~/.CIMtools.ini
fi

export SETUP_DIR=${UTILS_DIR}/Utils
export FORCEFIELD=${SETUP_DIR}/cvffTemplates.xml
export CLASSPATH=${JCHEM_DIR}/lib/jchem.jar:${UTILS_DIR}:

java Utils.CA_Prop_Map2011 -f $1 -o $2 -stdoptions $3
