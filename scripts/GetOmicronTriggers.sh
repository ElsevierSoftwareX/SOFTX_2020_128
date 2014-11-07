#!/bin/bash
#
# GetOmicronTriggers.sh
#
# print a selection of triggers
#
# Author: Florent Robinet
# florent.robinet@lal.in2p3.fr

printhelp(){
    echo ""
    echo "Usage:"
    echo "GetOmicronTriggers -c [CHANNEL_NAME] -s [GPS_START] -e [GPS_END]"
    echo " |__ prints the list of 'standard' Omicron triggers"
    echo "     of a given channel and between 2 GPS times"
    echo ""
    echo "GetOmicronTriggers -t \"[FILE_PATTERN]\" -s [GPS_START] -e [GPS_END]"
    echo " |__ prints the list of Omicron triggers saved in a list of root files"
    echo "     designated by a file pattern and between 2 GPS times"
    echo ""
    echo "By default 3 columns are printed: GPS time [s] - frequency [Hz] - SNR."
    echo "This can be changed with output options, see below."
    echo ""
    echo "[GPS_START] and [GPS_END] should belong to the same run"
    echo ""
    echo "TRIGGER OPTIONS:"
    echo "  -c  [CHANNEL_NAME]    channel name. Default = V1:h_4096Hz"
    echo "  -t  \"[FILE_PATTERN]\"  the user provides his own trigger files"
    echo "                        The file pattern is [TRIGGER_PATTERN]"
    echo "  -s  [GPS_START]       starting GPS time"
    echo "  -e  [GPS_END]         stopping GPS time"
    echo "  -x  [SNR_MIN]         minimum SNR value"
    echo "                        Default = '1'"
    echo "  -X  [SNR_MAX]         maximum SNR value"
    echo "                        Default = '10000'"
    echo "  -f  [FREQUENCY_MIN]   minimum Frequency value"
    echo "                        Default = '1'"
    echo "  -F  [FREQUENCY_MAX]   maximum Frequency value"
    echo ""
    echo ""
    echo "CLUSTER OPTIONS"
    echo "  -C  [TIME/TIMEFREQ]   [TIME] = print time-clustered triggers"
    echo "           /NONE        [TIMEFREQ]  = print time-frequency-clustered triggers"
    echo "                        [NONE]  = print unclustered triggers"
    echo "                        Default = 'TIME'"
    echo "  -T  [DELTAT]          clustering delta_t in seconds"
    echo "                        Default = 0.1s"
    echo "  -A  [AF]              clustering frequency asymmetry (<1)"
    echo "                        Default = 0.05"
    echo ""
    echo "OUTPUT CONTROL"
    echo "  -r                    do not print frequency column"
    echo "  -u                    print duration column"
    echo "  -a                    print bandwidth column (not available with -C option)"
    echo "  -o                    print tstart column"
    echo "  -j                    print tend column"
    echo "  -p                    print fstart column"
    echo "  -v                    print fend column"
    echo "  -m                    print amplitude column"
    echo "  -n                    do not print SNR column"
    echo "  -q                    print Q column (not available with -C option)"
    echo ""
    echo "  -l                    list available channels"
    echo "  -h                    print this help"
    echo ""
    echo "Author: Florent Robinet (LAL - Orsay): robinet@lal.in2p3.fr"
    echo ""
} 


##### Check the Omicron environment
if [ -z "$OMICRONROOT" ]; then
    echo "`basename $0`: The Omicron environment is not set"
    exit 1
fi

##### default options
channel="V1:h_4096Hz" # channel name
snrmin=1           # minimum SNR
snrmax=10000       # maximum SNR
freqmin=1          # minimum frequency
freqmax=10000      # maximum frequency
cluster="TIME"     # clusters
cluster_dt=0.1     # cluster deltat
cluster_af=0.05    # cluster Af
print_freq=1       # print frequency column
print_dur=0        # print duration column
print_bw=0         # print bandwidth column
print_ts=0         # print tstart column
print_te=0         # print tend column
print_fs=0         # print fstart column
print_fe=0         # print fend column
print_q=0          # print Q column
print_amp=0        # print Amplitude column
print_snr=1        # print SNR column
triggerfiles="NONE" # user trigger files
tmin=0             # starting time
tmax=0             # stopping time

##### read options
while getopts ":c:t:s:e:x:X:f:F:C:T:A:ruanojmpvqlh" opt; do
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
	x)
	    snrmin="$OPTARG"
	    ;;
	X)
	    snrmax="$OPTARG"
	    ;;
	f)
	    freqmin="$OPTARG"
	    ;;
	F)
	    freqmax="$OPTARG"
	    ;;
	C)
	    cluster="$OPTARG"
	    ;;
	T)
	    cluster_dt="$OPTARG"
	    ;;
	A)
	    cluster_af="$OPTARG"
	    ;;
	r)
	    print_freq=0
	    ;;
	u)
	    print_dur=1
	    ;;
	a)
	    print_bw=1
	    ;;
	o)
	    print_ts=1
	    ;;
	j)
	    print_te=1
	    ;;
	p)
	    print_fs=1
	    ;;
	v)
	    print_fe=1
	    ;;
	m)
	    print_amp=1
	    ;;
	n)
	    print_snr=0
	    ;;
	q)
	    print_q=1
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
	    echo "type  'GetOmicronTriggers -h'  for help"
	    exit 1
	    ;;
    esac
done
OPTIND=1

##### check timing
if [ $tmin -lt 700000000 ]; then
    echo "`basename $0`: '$tmin' is not a reasonable starting time"
    echo "type  'GetOmicronTriggers -h'  for help"
    exit 1
fi
if [ $tmax -lt 700000000 ]; then
    echo "`basename $0`: '$tmax' is not a reasonable stop time"
    echo "type  'GetOmicronTriggers -h'  for help"
    exit 1
fi
if [ $tmax -le $tmin ]; then
    echo "`basename $0`: the time interval '$tmin-$tmax' is not reasonable"
    echo "type  'GetOmicronTriggers -h'  for help"
    exit 1
fi

##### check SNR
snrmin_int=`echo $snrmin | awk '{print int($1)}'`
snrmax_int=`echo $snrmax | awk '{print int($1)}'`
if [ $snrmin_int -lt 1 ]; then
    echo "`basename $0`: SNR '$snrmin' is not a reasonable (<1)"
    echo "type  'GetOmicronTriggers.sh -h'  for help"
    exit 1
fi
if [ $snrmax_int -lt 1 ]; then
    echo "`basename $0`: SNR '$snrmax' is not a reasonable (<1)"
    echo "type  'GetOmicronTriggers.sh -h'  for help"
    exit 1
fi

##### check Frequency
freqmin_int=`echo $freqmin | awk '{print int($1)}'`
freqmax_int=`echo $freqmax | awk '{print int($1)}'`
if [ $freqmin_int -lt 1 ]; then
    echo "`basename $0`: Frequency '$freqmin' is not a reasonable (<1)"
    echo "type  'GetOmicronTriggers.sh -h'  for help"
    exit 1
fi
if [ $freqmax_int -lt 1 ]; then
    echo "`basename $0`: Frequency '$freqmax' is not a reasonable (<1)"
    echo "type  'GetOmicronTriggers.sh -h'  for help"
    exit 1
fi

##### case where the trigger files are provided
if [ ! "$triggerfiles" = "NONE" ]; then
    printtriggers.exe "$triggerfiles" $tmin $tmax $snrmin $snrmax $freqmin $freqmax $print_freq $print_dur $print_bw $print_ts $print_te $print_fs $print_fe $print_snr $print_amp $print_q $cluster $cluster_dt $cluster_af
    exit 0
fi

##### Check the Omicron environment
if [ -z "$OMICRON_TRIGGERS" ]; then
    echo "`basename $0`: The Omicron trigger environment is not set"
    exit 1
fi

##### get file list
. ${OMICRON_SCRIPTS}/GetTriggerFileList.sh -c $channel -s $tmin -e $tmax

##### print triggers
if [ ! "$OMICRON_FILELIST" = "" ]; then
    printtriggers.exe "$OMICRON_FILELIST" $tmin $tmax $snrmin $snrmax $freqmin $freqmax $print_freq $print_dur $print_bw $print_ts $print_te $print_fs $print_fe $print_snr $print_amp $print_q $cluster $cluster_dt $cluster_af
    exit 0
else 
    echo "`basename $0`: Omicron triggers are not available"
fi

exit 0
