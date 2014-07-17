#!/bin/sh
#
# Online2Offline.sh
#
# merge trigger files and archive them
#
# Author: Florent Robinet
# florent.robinet@lal.in2p3.fr

printhelp(){
    echo ""
    echo "Usage:"
    echo "Online2Offline -c[CHANNEL_NAME]"
    echo ""
    echo "Example: Online2Offline -c V1:h_4096Hz"
    echo ""
    echo "TRIGGER SELECTION"
    echo "  -c  [CHANNEL_NAME]  triggers from channel [CHANNEL_NAME] only"
    echo ""
    echo "TIMING"
    echo "  -d  [DELAY]         time delay to mark files as old [s]"
    echo "                      Default = '100000'"
    echo ""
    echo "OUTPUT CONTROL"
    echo "  -a                  only perform archiving (no online merging)"
    echo ""
    echo "  -h                  print this help"
    echo ""
} 

##### Check the environment
if [ -z "$OMICRON_ONLINE_TRIGGERS" ]; then
    echo "Error: The Omicron environment is not set"
    exit 1
fi

##### default options
channel="${OMICRON_ONLINE_TRIGGERS}/??:*" # channel name
delay=100000
only_archive=0

##### read options
while getopts ":c:d:ah" opt; do
    case $opt in
	c)
	    channel="${OMICRON_ONLINE_TRIGGERS}/$OPTARG"
	    ;;
	d)
	    delay=`echo $OPTARG | awk '{print int($1)}'`
	    ;;
	a)
	    only_archive=1
	    ;;
	h)
	    printhelp
	    exit 0
	    ;;
	\?)
	    echo "Invalid option: -$OPTARG"
	    echo "type  'Online2Offline -h'  for help"
	    exit 1
	    ;;
    esac
done

##### timing
now=`tconvert now`
oldtime=$(( $now - $delay ))
now_base=$(( $now / $OMICRON_TRIGGERS_BASE ))
oldtime_base=$(( $oldtime / $OMICRON_TRIGGERS_BASE ))
now_base1000=$(( $now / 1000 ))
echo "Online2Offline: now=$now"
echo "Online2Offline: archive triggers before $(( $oldtime_base * $OMICRON_TRIGGERS_BASE ))"

##### select run
run="NONE"
for r in $RUN_NAMES; do
    r_s=${r}_START
    r_e=${r}_END
    if [[ $now -ge ${!r_s} && $now -lt ${!r_e} ]]; then
	run=$r
	break;
    fi
done

if [ $run = "NONE" ]; then
    echo "Invalid time: not configured for this time"
    echo "Possible runs = $RUN_NAMES"
    exit 1 
fi

# loop over channels
for chandir in ${channel}; do
    if [ ! -d ${chandir} ]; then continue; fi
    if [ ! "$(ls -A $chandir)" ]; then continue; fi
    echo "Online2Offline: archive triggers in ${chandir}..."

    # first online file
    first_file=`ls ${chandir}/ |head -1`
    first_start=`echo $first_file | awk -F_ '{print $((NF -1))}'`

    # merging online trigger files
    if [ $only_archive -eq 0 ]; then
	echo "                merging online triggers..."
	b1000=$(( $first_start / 1000 ))		
 
	# tmp place holder
	o2o_tmp=${TMP}/o2o_${RANDOM}-${now}
	mkdir -p $o2o_tmp
 	while [ $b1000 -lt $now_base1000 ]; do
	    nfiles=`ls ${chandir}/*_${b1000}*_*.root 2>&1 | wc -l`
	    if [ $nfiles -lt 2 ]; then
		b1000=$(( $b1000 + 1 ))
		continue;
	    fi
	    if triggermerge.exe $o2o_tmp "${chandir}/*_${b1000}*_*.root" 2>&1 | grep -q "no livetime"; then
		echo "no files -> skip"
	    else
		if [ "$(ls -A ${o2o_tmp})" ]; then
		    rm -f ${chandir}/*_${b1000}*_*.root
		    mv ${o2o_tmp}/*.root ${chandir}/
		fi
	    fi
	    b1000=$(( $b1000 + 1 ))
	done
	rm -fr $o2o_tmp
    fi

    # merge and archive files
    echo "                archiving online triggers..."
    b=$(( $first_start / $OMICRON_TRIGGERS_BASE ))
    channel_name=${channel##*/}
    mkdir -p ${OMICRON_TRIGGERS}/${run}/${channel_name}
    while [ $b -lt $oldtime_base ]; do
	triggermerge.exe ${OMICRON_TRIGGERS}/${run}/${channel_name} "${chandir}/*_${b}*_*.root"
	rm -f ${chandir}/*_${b}*_*.root
	b=$(( $b + 1 ))
    done

done

exit 0
