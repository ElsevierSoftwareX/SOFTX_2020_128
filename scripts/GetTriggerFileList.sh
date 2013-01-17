#!/bin/bash
#
# GetTriggerFileList.sh
#
# list an optimized list of trigger files
#
# Author: Florent Robinet
# florent.robinet@lal.in2p3.fr

##### default options
channel="h_4096Hz" # channel name
outdir=`pwd` # output directory

printhelp(){
    echo ""
    echo "Usage:"
    echo "GetTriggerFileList -c[CHANNEL_NAME] [GPS_START] [GPS_STOP]"
    echo ""
    echo "Example: GetTriggerFileList -ch_4096Hz 934228815 934232415"
    echo ""
    echo "TRIGGER SELECTION OPTIONS"
    echo "  -c  [CHANNEL_NAME]  trigger files for channel [CHANNEL_NAME]"
    echo "                      Default = 'h_4096Hz'"
    echo ""
    echo "OUTPUT CONTROL"
    echo "  -h                  print this help"
    echo ""
} 

##### Check the environment
if [[ -z "$OMICRON_TRIGGERS" ]]; then
    echo "Error: The Omicron environment is not set"
    exit 1
fi

##### needs argument
if [ $# -lt 1 ]; then
    printhelp
    exit 1
fi

##### read options
while getopts ":c:h" opt; do
    case $opt in
	c)
	    channel="$OPTARG"
	    ;;
	h)
	    printhelp
	    exit 0
	    ;;
	\?)
	    echo "Invalid option: -$OPTARG"
	    echo "type  'GetTriggerFileList.sh -h'  for help"
	    exit 1
	    ;;
    esac
done

##### gps interval
shift $(($OPTIND - 1))
tmin=`echo $1 | awk '{print int($1)}'`
tmax=`echo $2 | awk '{print int($1)}'`
OPTIND=0

##### select run
run="NONE"
for r in $RUN_NAMES; do
    r_s=${r}_START
    r_e=${r}_END
    if [[ $tmin -ge ${!r_s} && $tmin -lt ${!r_e} ]]; then
	if [[ $tmax -gt ${!r_s} && $tmax -le ${!r_e} ]]; then
	  run=$r
	  break;
	fi
    fi
done

if [ $run = "NONE" ]; then
    echo "Invalid GPS times: the time interval must be entirely contained in a single run:"
    echo "Possible runs = $RUN_NAMES"
    exit 1 
fi
if [ ! -d ${OMICRON_TRIGGERS}/${run} ]; then
    echo "`basename $0`: the triggers for the run $run do not exist"
    exit 1
fi


##### get available channels
. GetOmicronChannels.sh -r $run > /dev/null 2>&1

##### check channel is available
if ! echo "$OMICRON_CHANNELS" | grep -q "$channel"; then
    echo "Invalid option: channel '${channel}' is not available"
    echo "type  'GetTriggerFileList -h'  for help"
    exit 1
fi

##### check timing
if [ $tmin -lt 700000000 ]; then
    echo "Invalid option: '$tmin' is not a reasonable starting time"
    echo "type  'GetTriggerFileList -h'  for help"
    exit 1
fi
if [ $tmax -lt 700000000 ]; then
    echo "Invalid option: '$tmax' is not a reasonable stop time"
    echo "type  'GetTriggerFileList -h'  for help"
    exit 1
fi
if [ $tmax -le $tmin ]; then
    echo "Invalid option: the time interval '$tmin-$tmax' is not reasonable"
    echo "type  'GetTriggerFileList -h'  for help"
    exit 1
fi

##### optimize timing
tmin_base=$(( $tmin / $OMICRON_TRIGGERS_BASE - 1 ))
tmax_base=$(( $tmax / $OMICRON_TRIGGERS_BASE + 1 ))

#### clean tmp dir with old trigger directories
for i in `find ${TMP} -maxdepth 1 -type d -mtime +2 -print | grep triggers`; do rm -fr $i > /dev/null 2>&1; done

#### tmp dir for online files
tag=$RANDOM
tmpdir=${TMP}/triggers.${RANDOM}
mkdir -p ${tmpdir}; rm -f ${tmpdir}/*.root > /dev/null 2>&1

##### map the trigger directory
triggers=""
first=1
while [ $tmin_base -le $tmax_base ]; do

    # start with offline triggers
    for file in ${OMICRON_TRIGGERS}/${run}/${channel}/${channel}_${tmin_base}*.root; do 
	if [ -e $file ]; then 
	    s=`echo $file | awk -F_ '{print $((NF -1))}'`
	    d=`echo $file | awk -F_ '{print $NF}' | awk 'BEGIN{FS=".root"} {print $1}'`
	    e=$(($s+$d))
	    if [[ $e -gt $tmin && $s -lt $tmax ]]; then
		triggers="$triggers $file"
	    fi	    
	fi 
    done

    # then online triggers
    for file in ${OMICRON_ONLINE_TRIGGERS}/${channel}/${channel}_${tmin_base}*.root; do 
	if [ -e $file ]; then 
	    s=`echo $file | awk -F_ '{print $((NF -1))}'`
	    d=`echo $file | awk -F_ '{print $NF}' | awk 'BEGIN{FS=".root"} {print $1}'`
	    e=$(($s+$d))
	    if [[ $e -gt $tmin && $s -lt $tmax ]]; then
		cp $file ${tmpdir}/
		if [ $first -eq 1 ]; then triggers="$triggers ${tmpdir}/*.root"; first=0; fi
	    fi	    
	fi 
    done

    let "tmin_base+=1"
done

echo "FILELIST $triggers"

exit 0