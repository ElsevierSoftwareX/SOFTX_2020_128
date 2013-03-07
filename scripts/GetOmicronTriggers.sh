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
    echo "GetOmicronTriggers -c[CHANNEL_NAME] [GPS_START] [GPS_STOP]"
    echo ""
    echo "By default, standard Omicron triggers are printed unless the '-t' option is used"
    echo ""
    echo "Example: GetOmicronTriggers -ch_4096Hz 934228815 934232415"
    echo ""
    echo "TRIGGER SELECTION OPTIONS"
    echo "  -c  [CHANNEL_NAME]  triggers from channel [CHANNEL_NAME]"
    echo "                      Default = 'h_4096Hz'"
    echo "  -t  [TRIGGER_FILES] the user provides his own trigger files"
    echo "                      The file pattern is [TRIGGER_FILES]"
    echo "                      this option overides the channel option -c"
    echo "  -s  [SNR_MIN]       minimum SNR value"
    echo "                      Default = '1'"
    echo "  -S  [SNR_MAX]       maximum SNR value"
    echo "                      Default = '10000'"
    echo "  -f  [FREQUENCY_MIN] minimum Frequency value"
    echo "                      Default = '1'"
    echo "  -F  [FREQUENCY_MAX] maximum Frequency value"
    echo "                      Default = '10000'"
    echo ""
    echo "CLUSTERS"
    echo "  -C  [YES/NO]        [YES] = print clustered triggers"
    echo "                      [NO]  = print unclustered triggers"
    echo "                      Default = 'YES'"
    echo "  -T  [CLUSTERING_DT] Delta_t window for clustering (in s.)"
    echo "                      Default = '0.1'"
    echo ""
    echo "OUTPUT CONTROL"
    echo "By default 3 columns are printed: time - frequency - SNR"
    echo "  -r                  do not print frequency column"
    echo "  -u                  print duration column"
    echo "  -a                  print bandwidth column (not available with -C option)"
    echo "  -o                  print tstart column"
    echo "  -e                  print tend column"
    echo "  -p                  print fstart column"
    echo "  -v                  print fend column"
    echo "  -n                  do not print SNR column"
    echo "  -q                  print Q column (not available with -C option)"
    echo ""
    echo "  -l                  list available channels"
    echo "  -h                  print this help"
    echo ""
} 

##### needs argument
if [ $# -lt 1 ]; then
    printhelp
    exit 1
fi

##### Check the Omicron environment
if [[ -z "$OMICRONROOT" ]]; then
    echo "Error: The Omicron environment is not set"
    exit 1
fi

##### default options
channel="h_4096Hz" # channel name
snrmin=1           # minimum SNR
snrmax=10000       # maximum SNR
freqmin=1          # minimum frequency
freqmax=10000      # maximum frequency
cluster="YES"      # cluster or no
print_freq=1       # print frequency column
print_dur=0        # print duration column
print_bw=0         # print bandwidth column
print_ts=0         # print tstart column
print_te=0         # print tend column
print_fs=0         # print fstart column
print_fe=0         # print fend column
print_q=0          # print Q column
print_snr=1        # print SNR column
dt=0.1            # clustering dt
triggerfiles="NONE" # user trigger files

##### read options
while getopts ":c:t:s:S:f:F:C:T:ruanoepvqlh" opt; do
    case $opt in
	c)
	    channel="$OPTARG"
	    ;;
	t)
	    triggerfiles="$OPTARG"
	    ;;
	s)
	    snrmin="$OPTARG"
	    ;;
	S)
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
	    dt="$OPTARG"
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
	e)
	    print_te=1
	    ;;
	p)
	    print_fs=1
	    ;;
	v)
	    print_fe=1
	    ;;
	n)
	    print_snr=0
	    ;;
	q)
	    print_q=1
	    ;;
	l)
	    for run in $RUN_NAMES; do
		GetOmicronChannels.sh -r $run | grep -v OMICRON_CHANNELS | grep -v "own shell scripts" | grep -v GetOmicronChannels.sh
	    done
	    exit 0
	    ;;
	h)
	    printhelp
	    exit 0
	    ;;
	\?)
	    echo "Invalid option: -$OPTARG"
	    echo "type  'GetOmicronTriggers.sh -h'  for help"
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
    echo "type  'GetOmicronTriggers.sh -h'  for help"
    exit 1
fi
if [ $tmax -lt 700000000 ]; then
    echo "Invalid option: '$tmax' is not a reasonable stop time"
    echo "type  'GetOmicronTriggers.sh -h'  for help"
    exit 1
fi
if [ $tmax -le $tmin ]; then
    echo "Invalid option: the time interval '$tmin-$tmax' is not reasonable"
    echo "type  'GetOmicronTriggers.sh -h'  for help"
    exit 1
fi

##### check SNR
snrmin_int=`echo $snrmin | awk '{print int($1)}'`
snrmax_int=`echo $snrmax | awk '{print int($1)}'`
if [ $snrmin_int -lt 1 ]; then
    echo "Invalid option: SNR '$snrmin' is not a reasonable (<1)"
    echo "type  'GetOmicronTriggers.sh -h'  for help"
    exit 1
fi
if [ $snrmax_int -lt 1 ]; then
    echo "Invalid option: SNR '$snrmax' is not a reasonable (<1)"
    echo "type  'GetOmicronTriggers.sh -h'  for help"
    exit 1
fi

##### check Frequency
freqmin_int=`echo $freqmin | awk '{print int($1)}'`
freqmax_int=`echo $freqmax | awk '{print int($1)}'`
if [ $freqmin_int -lt 1 ]; then
    echo "Invalid option: Frequency '$freqmin' is not a reasonable (<1)"
    echo "type  'GetOmicronTriggers.sh -h'  for help"
    exit 1
fi
if [ $freqmax_int -lt 1 ]; then
    echo "Invalid option: Frequency '$freqmax' is not a reasonable (<1)"
    echo "type  'GetOmicronTriggers.sh -h'  for help"
    exit 1
fi

##### check cluster
if [ ! "$cluster" = "NO" ] ; then print_c=1; else print_c=0; fi

##### check dt
dt_int=`echo $dt | awk '{printf "%.4f\n", $1}'`
if [ "$dt" = "0.0000" ] ; then
    echo "Invalid option: the delta_t value '$dt' is not reasonable"
    echo "type  'GetOmicronTriggers.sh -h'  for help"
    exit 1
fi

##### case where the trigger files are provided
if [ ! "$triggerfiles" = "NONE" ]; then
    printtriggers.exe "$triggerfiles" $dt $tmin $tmax $snrmin $snrmax $freqmin $freqmax $print_freq $print_dur $print_bw $print_ts $print_te $print_fs $print_fe $print_snr $print_q $print_c
    exit 0
fi

##### Check the Omicron environment
if [[ -z "$OMICRON_TRIGGERS" ]]; then
    echo "Error: The Omicron environment is not set"
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
    echo "type  'GetOmicronTriggers.sh -h'  for help"
    exit 1
fi

##### map the trigger directory
triggers=`GetTriggerFileList.sh -c${channel} $tmin $tmax | grep "FILELIST" | sed 's|FILELIST ||g'`
if [ "$triggers" = "" ]; then
    echo "triggers are not available for $tmin-$tmax"
    exit 1
fi

##### retrieve triggers
printtriggers.exe "$triggers" $dt $tmin $tmax $snrmin $snrmax $freqmin $freqmax $print_freq $print_dur $print_bw $print_ts $print_te $print_fs $print_fe $print_snr $print_q $print_c

exit 0