#!/usr/bin/env bash
# make work dir
# $1 path/to/file{.[0-9]+}
# $2 svr or svc.
# $3 maxconfig
# $4 repetitions
# $5 CV
# $workdir /path/to/worktmp.file{/}
GACONF=/path/to/ga

datadir=$1
workdir=${datadir}/work

if [[ $2 == 'svr' ]] ; then
    mod='SVMreg'
else
    mod='SVMclass'
fi

if [[ -z $3 ]] ; then
    maxconfigs=3000
else
    maxconfigs=$3
fi

if [[ -z $4 ]] ; then
    ntrials=5
else
    ntrials=$4
fi

if [[ -z $5 ]] ; then
    lo=5
else
    lo=$5
fi

# make SVMreg file
for i in `find ${datadir} -type f -name "*.svm"`; do
    awk '{print $1}' ${i} > ${datadir}/file.${mod};
    break;
done

${GACONF}/pilot_local.csh data_dir=${datadir} workdir=${workdir} maxconfigs=${maxconfigs} mode=${mod} ntrials=${ntrials} lo=${lo} > ${datadir}/GA.log 2>&1
rc=$?

killall pilot_local.csh local_SVMreg.csh svm-train svm-predict

if [[ ${rc} != 0 ]]; then exit ${rc}; fi  # return exitcode
