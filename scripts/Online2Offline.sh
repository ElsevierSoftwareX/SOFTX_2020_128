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
    echo "                      Default = '5000'"
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
delay=5000

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

##### timing
now=`tconvert now`
oldtime=$(( $now - $delay + 0 ))

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
    echo "Invalid option: channel '${channel}' is not available"
    echo "type  'Online2Offline -h'  for help"
    exit 1
fi

##### select online files to merge
cd $OMICRON_ONLINE_TRIGGERS
tomerge=""
# loop over online trigger files of a given channel
for file in ${channel}_*.root; do
    if [ ! -e $file ]; then continue; fi

    # timing
    s=`echo $file | awk -F_ '{print $((NF -1))}'`
    d=`echo $file | awk -F_ '{print $NF}' | awk 'BEGIN{FS=".root"} {print $1}'`
    e=$(($s+$d))
    
    # keep only old files
    if [ $e -gt $oldtime ]; then break; fi
	
    # update list of files
    tomerge="${tomerge} ${OMICRON_ONLINE_TRIGGERS}/${file}"

done
cd `dirname $0`

##### anything to merge?
if [ "$tomerge" = "" ]; then
    echo "Online2Offline: no files to merge"
    exit 1
fi

##### merge files
triggermerge.exe ${OMICRON_TRIGGERS}/${channel} ${channel} "${tomerge}"
rm -f ${tomerge}



exit 0