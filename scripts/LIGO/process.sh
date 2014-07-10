#!/bin/sh
#
# process.sh
#
# Author: Florent Robinet
# robinet@lal.in2p3.fr
cd dirnname $0

###### user parameters

delay=300 # do not look at data after now-delay
IFO="H1"
data_type="${IFO}_C ${IFO}_R"
tmax=10000 # delay after which a new segment is started
proddir=`pwd`

######################


# GWOLLUM environment
source /home/detchar/opt/virgosoft/environment.sh "" >> /dev/null

# directories
mkdir -p ./logs
mkdir -p ./high
mkdir -p ./low
mkdir -p ./gw

# vars
now=`tconvert now`
logfile=./logs/process.${now}.txt
echo "Local time: `date`" > $logfile
echo "UTC time: `date`" >> $logfile

# for safety
condor_release -all >> /dev/null 2>&1

# check channels
if [ ! -e ./channels.${IFO} ]; then
    echo "`date -u`: ./channels.${IFO} is missing" >> $logfile
    exit 1
fi

# timing (high)
tstop=$(( ($now - $delay) / 2000 ))000
tstart=$(( $tstop - 2000 ))
let "tstop+=3"  # for overlap
let "tstart-=3" # for overlap
if [ ! -e ./high/segments.txt ]; then
    echo "`date -u`: no previous segment for high --> start a new one" >> $logfile
    echo "$tstart $tstop" > ./high/segments.txt
else
    prev_stop=`head -n 1 ./high/segments.txt | awk '{print $2}'`
    if [ $(( $tstop - $prev_stop )) -gt $tmax ]; then
	echo "`date -u:` the last processed segment is too old for high --> start a new one" >> $logfile
	echo "$tstart $tstop" > ${workdir}/segments.txt
    else
	tstart=$(( $prev_stop - 6 ))
	echo "`date -u:` create new segment for high" >> $logfile
        echo "$tstart $tstop" > ./high/segments.txt
    fi
fi


################################################################################
###########                    directory prep                        ###########
################################################################################

# channel lists
awk '$1>1024 {print}'  ./channels.${IFO} > ./high/channel.list
awk '$1<=1024 {print}' ./channels.${IFO} > ./low/channel.list

