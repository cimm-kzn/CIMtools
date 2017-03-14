#!/usr/bin/env bash

[ $# -ge 1 -a -f "$1" ] && input="$1" || input="-"

UTILS_DIR=~
CHMXN_DIR=~/ChemAxon

if [ -f /etc/.MODtools.rc ]; then
    . /etc/.MODtools.rc
fi

if [ -f ~/.MODtools.rc ]; then
    . ~/.MODtools.rc
fi

export CLASSPATH=${CHMXN_DIR}/JChem/lib/jchem.jar:${UTILS_DIR}:

cat ${input} | java Utils.react_desc -svm
