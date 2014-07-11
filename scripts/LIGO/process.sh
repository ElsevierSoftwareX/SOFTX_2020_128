#!/bin/sh
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
source /home/detchar/opt/virgosoft/environment.v1r3.sh "" >> /dev/null

# directories
mkdir -p ./logs
mkdir -p ./std/triggers
mkdir -p ./low/triggers
mkdir -p ./gw

# vars
now=`tconvert now`
logfile=`pwd`/logs/process.${now}.txt
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
tstop=$(( ($now - $delay) / 1000 ))
if [ $(( $tstop % 2 )) -eq 1 ]; then let "tstop-=1"; fi
tstop="${tstop}000"
tstart=$(( $tstop - 2000 ))
if [ ! -e ./segments.txt ]; then
    echo "`date -u`: no previous segment --> start a new one [$tstart; $tstop]" >> $logfile
    echo "$tstart $tstop" > ./segments.tmp
else
    prev_stop=`head -n 1 ./segments.txt | awk '{print $2}'`
    if [ $(( $tstop - $prev_stop )) -gt $tmax ]; then
	echo "`date -u`: the last processed segment is too old --> start a new one [$tstart; $tstop]" >> $logfile
	echo "$tstart $tstop" > ./segments.tmp
    else
	tstart=$prev_stop
	echo "`date -u`: create new segment [$tstart; $tstop]" >> $logfile
        echo "$tstart $tstop" > ./segments.tmp
    fi
fi
if [ `segsum ./segments.tmp` -lt 2000 ]; then
    echo "`date -u`: the new segment is too short, start again later..." >> $logfile
    rm -f ./segments.tmp
    exit 0
fi
mv ./segments.tmp ./segments.txt


# data
echo "`date -u`: get LCF file" >> $logfile
for type in $data_type; do
    ligo_data_find -o ${IFO:0:1} -l -t $type -u file -s $(( $tstart - 100 )) -e $(( $tstop + 100 )) 1>./frames.lcf 2>> $logfile
    if [ -s ./frames.lcf ]; then break; fi
done
if [ ! -s ./frames.lcf ]; then
    echo "`date -u`: data are missing for this segment" >> $logfile
    rm -f ./frames.lcf
    exit 0
fi

################################################################################
###########                    directory prep                        ###########
################################################################################

# channel lists
awk '$2=="STD" {print}' ./channels.${IFO} > ./std/channels.list
awk '$2=="LOW" {print}' ./channels.${IFO} > ./low/channels.list

# segments
awk '{print $1-3,$2+3}' ./segments.txt > ./std/segments.txt
awk '{print $1-24,$2+24}' ./segments.txt > ./low/segments.txt

# LCF
cp ./frames.lcf ./std/frames.lcf
mv ./frames.lcf ./low/frames.lcf

# Omicron parameters
cd ./std
GetOmicronOptions -c ./channels.list -f ./frames.lcf -d ./triggers >> $logfile 2>&1
cd ../low
GetOmicronOptions -c ./channels.list -f ./frames.lcf -d ./triggers >> $logfile 2>&1
cd ..

# make dags
GetOmicronDAG -f -d `pwd`/std >> $logfile 2>&1
GetOmicronDAG -f -d `pwd`/low >> $logfile 2>&1

################################################################################
###########                         GO!                              ###########
################################################################################                 

cd ./std
condor_submit_dag omicron.dag
cd ../low
condor_submit_dag omicron.dag
cd ..

# cleaning
rm -f ./low/parameters_LOW_.txt

exit 0