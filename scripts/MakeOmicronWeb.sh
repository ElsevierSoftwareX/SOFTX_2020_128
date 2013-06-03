#!/bin/bash
#
# MakeOmicronWeb.sh
#
# generate omicron web pages
#
# Author: Florent Robinet
# florent.robinet@lal.in2p3.fr

##### default options
channel="h_4096Hz" # channel name
stopgps=0
dayoption=0
houroption=1

printhelp(){
    echo ""
    echo "Usage:"
    echo "MakeOmicronWeb -a -c[CHANNEL_NAME] [GPS]"
    echo ""
    echo "Examples:"
    echo "   MakeOmicronWeb -a -ch_4096Hz 934228815"
    echo "   MakeOmicronWeb -o -ch_4096Hz now"
    echo ""
    echo "TRIGGER SELECTION OPTIONS"
    echo "  -c  [CHANNEL_NAME]  triggers from channel [CHANNEL_NAME]"
    echo "                      Default = 'h_4096Hz'"
    echo ""
    echo "WEB FORMAT"
    echo "  -a                  Day page"
    echo "  -o                  Hour page"
    echo "  -s                  the time interval stops at [GPS]"
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
while getopts ":c:saoh" opt; do
    case $opt in
	c)
	    channel="$OPTARG"
	    ;;
	a)
	    dayoption=1
	    houroption=0
	    ;;
	o)
	    dayoption=0
	    houroption=1
	    ;;
	s)
	    stopgps=1
	    ;;
	h)
	    printhelp
	    exit 0
	    ;;
	\?)
	    echo "Invalid option: -$OPTARG"
	    echo "type  'MakeOmicronWeb.sh -h'  for help"
	    exit 1
	    ;;
    esac
done

##### gps interval
shift $(($OPTIND - 1))
gps=`echo $1 | awk '{print int($1)}'`
OPTIND=0

##### current time with a delay of 20s
if [ $gps -eq 0 ]; then
    gps_=`tconvert now`
    gps=$(( $gps_ - 20 ))
fi

##### select run
run="NONE"
for r in $RUN_NAMES; do
    r_s=${r}_START
    r_e=${r}_END
    if [[ $gps -ge ${!r_s} && $gps -lt ${!r_e} ]]; then
	run=$r
	break;
    fi
done

if [ $run = "NONE" ]; then
    echo "Invalid GPS time"
    exit 1 
fi

##### get available channels
. GetOmicronChannels.sh -r $run > /dev/null 2>&1

##### check channel is available
if ! echo "$OMICRON_CHANNELS" | grep -q "$channel"; then
    echo "Invalid option: channel '${channel}' is not available"
    echo "type  'MakeOmicronWeb -h'  for help"
    exit 1
fi

##### get date
year=`tconvert -f "%Y" ${gps}`
month=`tconvert -f "%m" ${gps}`
day=`tconvert -f "%d" ${gps}`
hour=`tconvert -f "%H" ${gps}`

# time of previous/next hour/day
if [ $dayoption -eq 0 ]; then
    previousdate=`tconvert "${year}-${month}-${day} ${hour}:00:00 UTC" - 1hour`
    nextdate=`tconvert "${year}-${month}-${day} ${hour}:00:00 UTC" + 1hour`
else
    previousdate=`tconvert "${year}-${month}-${day} 00:00:00 UTC" - 1day`
    nextdate=`tconvert "${year}-${month}-${day} 00:00:00 UTC" + 1day`
fi
previousyear=`tconvert -f "%Y" ${previousdate}`
previousmonth=`tconvert -f "%m" ${previousdate}`
previousday=`tconvert -f "%d" ${previousdate}`
previoushour=`tconvert -f "%H" ${previousdate}`
nextyear=`tconvert -f "%Y" ${nextdate}`
nextmonth=`tconvert -f "%m" ${nextdate}`
nextday=`tconvert -f "%d" ${nextdate}`
nexthour=`tconvert -f "%H" ${nextdate}`
previouslink="getomicronpage.php\?channel=${channel}\&timemode=date\&year=${previousyear}\&month=${previousmonth}\&day=${previousday}"
nextlink="getomicronpage.php\?channel=${channel}\&timemode=date\&year=${nextyear}\&month=${nextmonth}\&day=${nextday}"
if [ $dayoption -eq 0 ]; then
    previouslink="../../../../${previouslink}\&hour=${previoushour}"
    nextlink="../../../../${nextlink}\&hour=${nexthour}"
else
    previouslink="../../../${previouslink}\&hour=all"
    nextlink="../../../${nextlink}\&hour=all"
fi
uplink="../../../../getomicronpage.php\?channel=${channel}\&timemode=date\&year=${year}\&month=${month}\&day=${day}\&hour=all"
downlinkday="../../../getomicronpage.php\?channel=${channel}\&timemode=date\&year=${year}\&month=${month}\&day=${day}\&hour="
downlinkhour="../../../../getomicronpage.php\?channel=${channel}\&timemode=date\&year=${year}\&month=${month}\&day=${day}\&hour="


##### day production or...
if [ $dayoption -eq 1 ]; then # start of the day

    if [ $stopgps -eq 1 ]; then
	stop=$gps
	start=`expr ${stop} - 86400`
	outdir=${OMICRON_WEB}/latestday
    else
	startstring=`tconvert -f "%Y-%m-%d 00:00:00" ${gps}`
	start=`tconvert "${startstring}"`
	stop=`expr ${start} + 86400`
	outdir=${OMICRON_WEB}/${year}/${month}/${day}
    fi

##### ... or hour production.
else # start of the hour
    if [ $stopgps -eq 1 ]; then
	stop=$gps
	start=`expr ${stop} - 3600`
	outdir=${OMICRON_WEB}/latesthour
    else
	startstring=`tconvert -f "%Y-%m-%d %H:00:00" ${gps}`
	start=`tconvert "${startstring}"`
	stop=`expr ${start} + 3600`
	outdir=${OMICRON_WEB}/${year}/${month}/${day}/${hour}
    fi
fi

##### make web page
mkdir -p ${outdir}
if [ $stopgps -eq 1 ]; then
    GetOmicronWeb.sh -w -c${channel} -d${outdir} $start $stop
else
    if [ $dayoption -eq 1 ]; then
	GetOmicronWeb.sh -w -c${channel} -d${outdir} -p${previouslink} -n${nextlink} -t${downlinkday} $start $stop
    else
	GetOmicronWeb.sh -w -c${channel} -d${outdir} -p${previouslink} -n${nextlink} -u${uplink} -t${downlinkhour} $start $stop
    fi
fi

exit 0