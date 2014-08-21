#!/bin/sh
#
# postprocess.sh
#
# Author: Florent Robinet
# robinet@lal.in2p3.fr
cd `dirname $0`
# GWOLLUM environment
source /home/detchar/opt/virgosoft/environment.v1r3.sh "" >> /dev/null

###### user parameters

hn=`hostname -d`
if [ "$hn" = "ligo-wa.caltech.edu" ]; then IFO="H1"
else  IFO="L1"; fi
run=${RUN_NAMES##* }
XMLDIR="/home/detchar/triggers/${run}/${IFO}"
types="std std2 low fine gw"

######################

# vars
now=`tconvert now`
logfile=`pwd`/logs/postprocess.${now}.txt
echo "Local time: `date`" > $logfile
echo "UTC time: `date -u`" >> $logfile

npostjob=`ps -efd | grep detchar | grep postprocess |wc -l`
echo "$npostjob postprocess jobs were detected" >> $logfile
if [ $npostjob -gt 4 ]; then
    echo "too many postprocess jobs -> exit" >> $logfile
    exit 0
fi

# clean log files
find ./logs -type f -mtime +4 -exec rm {} \; >> /dev/null 2>&1

# loop over production types
for t in $types; do
    if [ ! -d ./${t}/triggers ]; then continue; fi
    cd ./${t}/triggers

    # loop over channels
    for channel in ${IFO}:*; do
	if [ ! -d ./${channel} ]; then continue; fi
	if [ ! "$(ls -A ./${channel})" ]; then continue; fi
	mkdir -p ${OMICRON_ONLINE_TRIGGERS}/${channel}
	mkdir -p ${OMICRON_TRIGGERS}/${run}/${channel}
		
	# move root files
	echo "${channel}: moving root files" >> $logfile
	mv ${channel}/*.root ${OMICRON_ONLINE_TRIGGERS}/${channel} >> /dev/null 2>&1

	# XML outdir
	naked_chan=${channel#*:}
	naked_chan_onlyunderscore=${naked_chan//-/_}
	xmloutdir="${XMLDIR}/${naked_chan}_Omicron"

	# move xml files
	echo "${channel}: moving xml files" >> $logfile
	for file in ${channel}/*.xml; do
	    if [ ! -e $file ]; then continue; fi
	    gps=`echo $file | cut -d "-" -f 4`
	    gpsroot=$(( $gps / 100000 ))
	    mkdir -p ${xmloutdir}/$gpsroot
	    mv $file ${xmloutdir}/$gpsroot
	    gzip ${xmloutdir}/${gpsroot}/*.xml
	done
    done

    cd ../..
done

# archiving
rmdir ${OMICRON_ONLINE_TRIGGERS}/* >> /dev/null 2>&1
Online2Offline.sh -d 100000 -a >> $logfile 2>&1
rmdir ${OMICRON_ONLINE_TRIGGERS}/* >> /dev/null 2>&1


exit 0
