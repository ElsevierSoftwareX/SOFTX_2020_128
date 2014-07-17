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
types="std std2 low gw"

######################


# loop over production types
for t in $types; do
    if [ ! -d ./${t}/triggers ]; then continue; fi
    cd ./${t}/triggers

    # loop over channels
    for channel in ${IFO}:*; do
	if [ ! -d ./$channel ]; then continue; fi
	mkdir -p ${OMICRON_ONLINE_TRIGGERS}/${channel}
	mkdir -p ${OMICRON_TRIGGERS}/${run}/${channel}
	if [ ! -d ./$channel ]; then continue; fi

	# move root files
	mv ${channel}/*.root ${OMICRON_ONLINE_TRIGGERS}/${channel} >> /dev/null 2>&1

	# XML outdir
	naked_chan=${channel#*:}
	naked_chan_onlyunderscore=${naked_chan//-/_}
	xmloutdir="${XMLDIR}/${naked_chan}_Omicron"

	# move xml files
	for file in ${channel}/*.xml; do
	    if [ ! -e $file ]; then continue; fi
	    gps=`ligolw_print -t segment_summary -c start_time`
	    gpsroot=$(( $gps / 100000 ))
	    mkdir -p ${xmloutdir}/$gpsroot
	    mv $file ${xmloutdir}/$gpsroot
	    gzip ${xmloutdir}/${gpsroot}/*.xml
	done
    done

    cd ../..
done


# archiving
cd ${OMICRON_ONLINE_TRIGGERS}
rmdir * >> /dev/null 2>&1
for chan in ??:*; do
    if [ ! -d ${chan} ]; then continue; fi
    if [ "$(ls -A ./$chan)" ]; then
	Online2Offline.sh -c $chan -d 100000 -a
    else rmdir ./$chan
    fi
done

exit 0
