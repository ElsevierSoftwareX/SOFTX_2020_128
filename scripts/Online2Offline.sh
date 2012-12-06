#!/bin/bash
#
# Online2Offline.sh
#
# merge trigger files and archive them
#
# Author: Florent Robinet
# florent.robinet@lal.in2p3.fr

# move to the script directory
cd `dirname $0`
. ../environment.sh

printhelp(){
    echo ""
    echo "Usage:"
    echo "Online2Offline.sh -c[CHANNEL_NAME]"
    echo ""
    echo "Example: Online2Offline.sh -c h_4096Hz"
    echo ""
    echo "TRIGGER SELECTION"
    echo "  -c  [CHANNEL_NAME]  triggers from channel [CHANNEL_NAME]"
    echo ""
    echo "TIMING"
    echo "  -d  [DELAY]         time delay to mark files as old [s]"
    echo "                      Default = '5000'"
    echo ""
    echo "OUTPUT CONTROL"
    echo "By default 3 columns are printed: time - frequency - SNR"
    echo "  -r                  do not print frequency column"
    echo "  -u                  print duration column"
    echo "  -a                  print bandwidth column (not available with -C option)"
    echo "  -n                  do not print SNR column"
    echo "  -q                  print Q column (not available with -C option)"
    echo ""
    echo "  -l                  list available channels"
    echo "  -h                  print this help"
    echo ""
} 

##### default options
channel="unknwon" # channel name
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
	    echo "type  'Online2Offline.sh -h'  for help"
	    exit 1
	    ;;
    esac
done

##### gps interval
#shift $(($OPTIND - 1))
#tmin=`echo $1 | awk '{print int($1)}'`
#tmax=`echo $2 | awk '{print int($1)}'`


##### check channel
if [ ! -d ${OMICRON_TRIGGERS}/${channel} ]; then
    echo "Invalid option: channel '${channel}' is not available in ${OMICRON_TRIGGERS}"
    echo "create the directory or:"
    echo "type  'Online2Offline.sh -h'  for help"
    exit 1
fi


##### timing
now=`tconvert now`
oldtime=$(( $now - $delay + 0 ))

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
    echo "Online2Offline.sh: no files to merge"
    exit 1
fi

##### merge files
#mkdir -p ${OMICRON_TRIGGERS}/${channel}
triggermerge.exe ${OMICRON_TRIGGERS}/${channel} ${channel} "${tomerge}"
rm -f ${tomerge}



exit 0