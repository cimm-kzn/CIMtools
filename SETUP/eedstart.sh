#!/usr/bin/env bash

[ $# -ge 1 -a -f "$1" ] && input="$1" || input="-"

UTILS_DIR=/home/stsouko
CHMXN_DIR=/home/stsouko/ChemAxon

export CLASSPATH=${CHMXN_DIR}/JChem/lib/jchem.jar:${UTILS_DIR}:

cat ${input} | java Utils.react_desc -svm
