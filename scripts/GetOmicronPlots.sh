#!/bin/bash
#
# GetOmicronPlots.sh
#
# produce Omicron monitoring plots
#
# Author: Florent Robinet
# florent.robinet@lal.in2p3.fr

printhelp(){
    echo ""
    echo "Usage:"
    echo "GetOmicronPlots -c [CHANNEL_NAME] -s [GPS_START] -e [GPS_END]"
    echo " |__ plots 'standard' Omicron triggers of a given channel and between 2 GPS times"
    echo ""
    echo "GetOmicronPlots -t \"[FILE_PATTERN]\" -s [GPS_START] -e [GPS_END]"
    echo " |__ plots Omicron triggers saved in a list of root files"
    echo "     designated by a file pattern and between 2 GPS times"
    echo ""
    echo "For standard triggers, [GPS_START] and [GPS_END] should belong to the same run"
    echo ""
    echo "OPTIONS:"
    echo "  -c  [CHANNEL_NAME]    channel name. Default = V1:h_4096Hz"
    echo "  -t  \"[FILE_PATTERN]\"  the user provides his own trigger files"
    echo "                        The file pattern is [TRIGGER_PATTERN]"
    echo "  -s  [GPS_START]       starting GPS time. Required!"
    echo "  -e  [GPS_END]         stopping GPS time. Required!"
    echo "  -d  [OUTDIR]          output directory where to save plots"
    echo "                        Default = current directory"
    echo "  -f  [FILENAME]        this parameter forces the output file naming to:"
    echo "                        [FILENAME]_type.gif"
    echo ""
    echo "  -l                    list available channels"
    echo "  -h                    print this help"
    echo ""
    echo "Author: Florent Robinet (LAL - Orsay): robinet@lal.in2p3.fr"
    echo ""
} 

##### Check the Omicron environment
if [ -z "$OMICRONROOT" ]; then
    echo "Error: The Omicron environment is not set"
    exit 1
fi

##### default options
channel="V1:h_4096Hz"  # channel name
outdir=`pwd`        # output directory
triggerfiles="NONE" # user trigger files
tmin=0
tmax=0
ofilename=""

##### read options
while getopts ":c:t:s:e:d:f:lh" opt; do
    case $opt in
	c)
	    channel="$OPTARG"
	    ;;
	t)
	    triggerfiles="$OPTARG"
	    ;;
	s)
	    tmin=`echo $OPTARG | awk '{print int($1)}'`
	    ;;
	e)
	    tmax=`echo $OPTARG | awk '{print int($1)}'`
	    ;;
	d)
	    outdir="$OPTARG"
	    ;;
	f)
	    ofilename="$OPTARG"
	    ;;
	l)
	    for run in $RUN_NAMES; do
		echo ""
		GetOmicronChannels.sh -v -r $run
	    done
	    exit 0
	    ;;
	h)
	    printhelp
	    exit 0
	    ;;
	\?)
	    echo "`basename $0`: Invalid option: -$OPTARG"
	    echo "type  'GetOmicronPlots -h'  for help"
	    exit 1
	    ;;
    esac
done
OPTIND=1

##### check outdir
if [ ! -d $outdir ] ; then
    echo "`basename $0`: the output directory $outdir cannot be found"
    echo "type  'GetOmicronPlots -h'  for help"
    exit 1
fi

##### check timing
if [ $tmin -eq 0 ]; then
    echo "`basename $0`: A starting time is required"
    echo "type  'GetOmicronPlots -h'  for help"
    exit 1
fi
if [ $tmax -eq 0 ]; then
    echo "`basename $0`: A stopping time is required"
    echo "type  'GetOmicronPlots -h'  for help"
    exit 1
fi
if [ $tmax -le $tmin ]; then
    echo "`basename $0`: the time interval '$tmin-$tmax' is not reasonable"
    echo "type  'GetOmicronPlots -h'  for help"
    exit 1
fi

##### case where the trigger files are provided
if [ ! "$triggerfiles" = "NONE" ]; then
    omiplot.exe $outdir "$triggerfiles" $tmin $tmax $ofilename
    exit 0
fi

##### Check the trigger environment
if [ -z "$OMICRON_TRIGGERS" ]; then
    echo "`basename $0`: The Omicron trigger environment is not set"
    exit 1
fi

##### more checks
if [ $tmin -lt 700000000 ]; then
    echo "`basename $0`: '$tmin' is not a reasonable starting time"
    echo "type  'GetOmicronPlots -h'  for help"
    exit 1
fi
if [ $tmax -lt 700000000 ]; then
    echo "`basename $0`: '$tmax' is not a reasonable stop time"
    echo "type  'GetOmicronPlots -h'  for help"
    exit 1
fi

##### get file list
. ${OMICRON_SCRIPTS}/GetTriggerFileList.sh -c $channel -s $tmin -e $tmax

##### print triggers
if [ ! "$OMICRON_FILELIST" = "" ]; then
    omiplot.exe $outdir "$OMICRON_FILELIST" $tmin $tmax $ofilename
    exit 0
else 
    echo "`basename $0`: Omicron triggers are not available"
fi

exit 0
