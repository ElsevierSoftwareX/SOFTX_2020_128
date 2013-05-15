\#!/bin/bash
#
# GetOmicronPlots.sh
#
# produce omicron-style plots
#
# Author: Florent Robinet
# florent.robinet@lal.in2p3.fr

##### default options
channel="required"  # channel name
outdir=`pwd`        # output directory
triggerfiles="NONE" # user trigger files

printhelp(){
    echo ""
    echo "Usage:"
    echo "GetOmicronPlots -c[CHANNEL_NAME] [GPS_START] [GPS_STOP]"
    echo ""
    echo "Example: GetOmicronPlots -ch_4096Hz 934228815 934232415"
    echo ""
    echo "TRIGGER SELECTION OPTIONS"
    echo "  -c  [CHANNEL_NAME]  triggers from channel [CHANNEL_NAME]"
    echo "                      this option is required"
    echo "  -t  [TRIGGER_FILES] the user provides his own trigger files"
    echo "                      The file pattern is [TRIGGER_FILES]"
    echo ""
    echo "OUTPUT CONTROL"
    echo "  -d  [OUTDIR]        output directory where to save plots"
    echo "                      Default = current directory"
    echo ""
    echo "  -h                  print this help"
    echo ""
} 

##### Check the Omicron environment
if [[ -z "$OMICRONROOT" ]]; then
    echo "Error: The Omicron environment is not set"
    exit 1
fi

##### needs argument
if [ $# -lt 1 ]; then
    printhelp
    exit 1
fi

##### read options
while getopts ":c:d:t:h" opt; do
    case $opt in
	c)
	    channel="$OPTARG"
	    ;;
	t)
	    triggerfiles="$OPTARG"
	    ;;
	d)
	    outdir="$OPTARG"
	    ;;
	h)
	    printhelp
	    exit 0
	    ;;
	\?)
	    echo "Invalid option: -$OPTARG"
	    echo "type  'GetOmicronTriggers -h'  for help"
	    exit 1
	    ;;
    esac
done

##### gps interval
shift $(($OPTIND - 1))
tmin=`echo $1 | awk '{print int($1)}'`
tmax=`echo $2 | awk '{print int($1)}'`
OPTIND=0

##### check timing
if [ $tmin -lt 700000000 ]; then
    echo "Invalid option: '$tmin' is not a reasonable starting time"
    echo "type  'GetOmicronPlots -h'  for help"
    exit 1
fi
if [ $tmax -le $tmin ]; then
    echo "Invalid option: the time interval '$tmin-$tmax' is not reasonable"
    echo "type  'GetOmicronPlots -h'  for help"
    exit 1
fi

##### channel is required
if [ "$channel" = "required" ]; then
    echo "GetOmicronPlots: the channel name is required with the -c option"
    echo "type  'GetOmicronPlots -h'  for help"
    exit 1
fi

##### check outdir
if [ ! -d $outdir ] ; then
    echo "Invalid option: the output directory $outdir cannot be found"
    echo "type  'GetOmicronPlots -h'  for help"
    exit 1
fi

##### case where the trigger files are provided
if [ ! "$triggerfiles" = "NONE" ]; then
    omicronplot.exe $channel $outdir "$triggerfiles" $tmin $tmax
    exit 0
fi

##### Check the trigger environment
if [[ -z "$OMICRON_TRIGGERS" ]]; then
    echo "Error: The Omicron trigger environment is not set"
    exit 1
fi

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

##### get available channels
. GetOmicronChannels.sh -r $run > /dev/null 2>&1

##### check channel is available
if ! echo "$OMICRON_CHANNELS" | grep -q "$channel"; then
    echo "Invalid option: channel '${channel}' is not available"
    echo "type  'GetOmicronPlots -h'  for help"
    exit 1
fi

##### map the trigger directory
triggers=`GetTriggerFileList.sh -c${channel} $tmin $tmax | grep "FILELIST" | sed 's|FILELIST ||g'`
if [ "$triggers" = "" ]; then
    echo "triggers are not available for $tmin-$tmax"
    exit 1
fi

##### no triggers?
if [ "$triggers" = "" ]; then
    echo "triggers are not available for $tmin-$tmax"
    exit 1
fi

##### retrieve triggers
omicronplot.exe $channel $outdir "$triggers" $tmin $tmax

exit 0