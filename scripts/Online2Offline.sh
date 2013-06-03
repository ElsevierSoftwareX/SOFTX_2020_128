#!/bin/bash
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
    echo "Example: Online2Offline -c h_4096Hz"
    echo ""
    echo "TRIGGER SELECTION"
    echo "  -c  [CHANNEL_NAME]  triggers from channel [CHANNEL_NAME]"
    echo ""
    echo "TIMING"
    echo "  -d  [DELAY]         time delay to mark files as old [s]"
    echo "                      Default = '100000'"
    echo ""
    echo "OUTPUT CONTROL"
    echo "  -h                  print this help"
    echo ""
} 

##### Check the environment
if [[ -z "$OMICRON_ONLINE_TRIGGERS" ]]; then
    echo "Error: The Omicron environment is not set"
    exit 1
fi

##### default options
channel="unknown" # channel name
delay=100000

##### read options
while getopts ":c:d:h" opt; do
    case $opt in
	c)
	    channel="$OPTARG"
	    ;;
	d)
	    delay=`echo $OPTARG | awk '{print int($1)}'`
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

##### Check channel
if [ ! -d ${OMICRON_ONLINE_TRIGGERS}/${channel} ]; then
    echo "Error: the online channel ${channel} does not exist"
    exit 1
fi

##### timing
now=`tconvert now`
oldtime=$(( $now - $delay + 0 ))
now_base=$(( $now / $OMICRON_TRIGGERS_BASE ))
oldtime_base=$(( $oldtime / $OMICRON_TRIGGERS_BASE ))
now_base1000=$(( $now / 1000 ))

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

##### get available channels
. GetOmicronChannels.sh -r $run > /dev/null 2>&1

##### check channel is available
if ! echo "$OMICRON_CHANNELS" | grep -q "$channel"; then
    echo "Invalid option: channel '${channel}' is not available (offline)"
    echo "type  'Online2Offline -h'  for help"
    exit 1
fi

# first online file
first_file=`ls ${OMICRON_ONLINE_TRIGGERS}/${channel}/ |head -1`
if [ "$first_file=" = "" ]; then
    echo "Error: there is no online files"
    exit 1
fi
first_start=`echo $first_file | awk -F_ '{print $((NF -1))}'`

# starting base
b=$(( $first_start / $OMICRON_TRIGGERS_BASE ))
echo "Starting base = ${b}"
b1000=$(( $first_start / 1000 ))
echo "Starting base1000 = ${b1000}"

# tmp place holder
mkdir -p ${TMP}/${channel}-${now}

# merge online files
while [ $b1000 -lt $now_base1000 ]; do
    echo "Merging ${OMICRON_ONLINE_TRIGGERS}/${channel}/${channel}_${b1000}*.root ..."
    if triggermerge.exe ${TMP}/${channel}-${now} ${channel} "${OMICRON_ONLINE_TRIGGERS}/${channel}/${channel}_${b1000}*.root" | grep -q "no livetime"; then
	echo "skip this files"
    else
	rm -f ${OMICRON_ONLINE_TRIGGERS}/${channel}/${channel}_${b1000}*.root
	mv ${TMP}/${channel}-${now}/*.root ${OMICRON_ONLINE_TRIGGERS}/${channel}/
    fi
    let "b1000+=1"
done
rm -fr ${TMP}/${channel}-${now}/


# merge and archive files
while [ $b -lt $oldtime_base ]; do
    echo "Merging and archiving ${OMICRON_ONLINE_TRIGGERS}/${channel}/${channel}_${b}*.root ..."
    triggermerge.exe ${OMICRON_TRIGGERS}/${run}/${channel} ${channel} "${OMICRON_ONLINE_TRIGGERS}/${channel}/${channel}_${b}*.root"
    rm -f ${OMICRON_ONLINE_TRIGGERS}/${channel}/${channel}_${b}*.root
    let "b+=1"
done

exit 0
