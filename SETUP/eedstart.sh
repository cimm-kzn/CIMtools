#!/usr/bin/env bash

[ $# -ge 1 -a -f "$1" ] && input="$1" || input="-"

UTILS_DIR=~
JCHEM_DIR=~/ChemAxon/JChem

if [ -f /etc/.CIMtools.ini ]; then
    . /etc/.CIMtools.ini
fi

if [ -f ~/.MODtools.ini ]; then
    . ~/.CIMtools.ini
fi

export CLASSPATH=${JCHEM_DIR}/lib/jchem.jar:${UTILS_DIR}:

cat ${input} | java Utils.react_desc -svm
