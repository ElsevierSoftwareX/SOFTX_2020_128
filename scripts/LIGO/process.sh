#!/bin/bash
#
# process.sh
#
# Author: Florent Robinet
# robinet@lal.in2p3.fr
cd `dirname $0`

###### user parameters

delay=300 # do not look at data after now-delay
IFO="H1"
data_type="${IFO}_C ${IFO}_R"
tmax=16000 # delay after which a new segment is started to get back on track
proddir=`pwd`

######################


# GWOLLUM environment
#source /home/detchar/opt/virgosoft/environment.sh "" >> /dev/null

# directories
mkdir -p ./logs
mkdir -p ./high/triggers
mkdir -p ./low/triggers
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

################################################################################
###########                          input                           ###########
################################################################################

# timing
tstop=$(( ($now - $delay) / 2000 ))000
tstart=$(( $tstop - 2000 ))
let "tstop+=3"  # for overlap
let "tstart-=3" # for overlap
if [ ! -e ./segments.txt ]; then
    echo "`date -u`: no previous segment --> start a new one [$tstart; $tstop]" >> $logfile
    echo "$tstart $tstop" > ./segments.tmp
else
    prev_stop=`head -n 1 ./segments.txt | awk '{print $2}'`
    if [ $(( $tstop - $prev_stop )) -gt $tmax ]; then
	echo "`date -u`: the last processed segment is too old --> start a new one [$tstart; $tstop]" >> $logfile
	echo "$tstart $tstop" > ./segments.tmp
    else
	tstart=$(( $prev_stop - 6 ))
	echo "`date -u`: create new segment [$tstart; $tstop]" >> $logfile
        echo "$tstart $tstop" > ./segments.tmp
    fi
fi
if [ `segsum ./segments.tmp` -lt 2006 ]; then
    echo "`date -u`: the new segment is too short, start again later..." >> $logfile
    rm -f ./segments.tmp
    exit 0
fi
mv ./segments.tmp ./segments.txt


# data
echo "`date -u`: get LCF file" >> $logfile
for type in $data_type; do
    ligo_data_find -o ${IFO:0:1} -l -t $type -u file -s $tstart -e $tstop > ./frames.lcf
    if [ -s ./frames.lcf ]; then break; fi
done
if [ !-s ./frames.lcf ]; then
    echo "`date -u`: data are missing for this segment" >> $logfile
    rm -f ./frames.lcf
    exit 0
}

################################################################################
###########                    directory prep                        ###########
################################################################################

# channel lists
awk '$1>1024 {print}'  ./channels.${IFO} > ./high/channels.list
awk '$1<=1024 {print}' ./channels.${IFO} > ./low/channels.list

# segments
cp -f ./segments.txt ./high/segments.txt
cp -f ./segments.txt ./low/segments.txt

# LCF
lalcache2ffl ./frames.lcf > ./high/frames.ffl
cp -f ./high/frames.ffl ./low/frames.ffl
rm -f ./frames.lcf

# Omicron parameters
cd ./high
GetOmicronOptions -o ../parameters_high.txt -c ./channels.list
cd ../low
GetOmicronOptions -o ../parameters_low.txt -c ./channels.list
