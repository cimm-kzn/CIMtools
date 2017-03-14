#!/usr/bin/env bash

UTILS_DIR=~
CHMXN_DIR=~/ChemAxon

if [ -f /etc/.MODtools.rc ]; then
    . /etc/.MODtools.rc
fi

if [ -f ~/.MODtools.rc ]; then
    . ~/.MODtools.rc
fi

export SETUP_DIR=${UTILS_DIR}/Utils
export FORCEFIELD=${SETUP_DIR}/cvffTemplates.xml
export CLASSPATH=${CHMXN_DIR}/JChem/lib/jchem.jar:${UTILS_DIR}:

java Utils/CA_Prop_Map2011 -f $1 -o $2 -stdoptions $3
