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

################################################################################
###########                       check inputs                        ########### 
################################################################################                               

# GWOLLUM environment
source /home/detchar/opt/virgosoft/environment.v1r3.sh "" >> /dev/null

# directories
mkdir -p ./logs
mkdir -p ./std/triggers ./std/dags; rmdir ./std/triggers/* >> /dev/null 2>&1
mkdir -p ./std2/triggers ./std2/dags; rmdir ./std2/triggers/* >> /dev/null 2>&1
mkdir -p ./low/triggers ./low/dags; rmdir ./low/triggers/* >> /dev/null 2>&1
mkdir -p ./gw/triggers ./gw/dags; rmdir ./gw/triggers/* >> /dev/null 2>&1

# vars
now=`tconvert now`
logfile=`pwd`/logs/process.${now}.txt
echo "Local time: `date`" > $logfile
echo "UTC time: `date`" >> $logfile

# check channels
if [ ! -e ./channels.${IFO} ]; then
    echo "`date -u`: ./channels.${IFO} is missing" >> $logfile
    exit 1
fi

################################################################################
###########                       previous condor                     ##########
################################################################################

# for safety
condor_release -all >> /dev/null 2>&1

echo "`date -u`: Cleaning anything older than 4 days..." >> $logfile

# clean log files
find ./logs -type f -mtime +4 -exec rm {} \; >> /dev/null 2>&1

# clean omicron log files
find ./low/logs -type f -mtime +4 -exec rm {} \; >> /dev/null 2>&1
find ./std/logs -type f -mtime +4 -exec rm {} \; >> /dev/null 2>&1
find ./std2/logs -type f -mtime +4 -exec rm {} \; >> /dev/null 2>&1
find ./gw/logs -type f -mtime +4 -exec rm {} \; >> /dev/null 2>&1
find ./low/dags -type f -mtime +4 -exec rm {} \; >> /dev/null 2>&1
find ./std/dags -type f -mtime +4 -exec rm {} \; >> /dev/null 2>&1
find ./std2/dags -type f -mtime +4 -exec rm {} \; >> /dev/null 2>&1
find ./gw/dags -type f -mtime +4 -exec rm {} \; >> /dev/null 2>&1


# check if previous batch is still running
echo "`date -u`: check if previous batch is still running..." >> $logfile
if [ -e ./low/omicron.dag.lock ]; then
    echo "`date -u` condor (low) is already running" >> $logfile
    exit 2
fi
if [ -e ./std/omicron.dag.lock ]; then
    echo "`date -u` condor (std) is already running" >> $logfile
    exit 2
fi
if [ -e ./std2/omicron.dag.lock ]; then
    echo "`date -u` condor (std2) is already running" >> $logfile
    exit 2
fi
if [ -e ./gw/omicron.dag.lock ]; then
    echo "`date -u` condor (gw) is already running" >> $logfile
    exit 2
fi

# archive dagman
if [ -e ./low/omicron.dag.dagman.out ]; then
    mv ./low/omicron.dag.dagman.out ./low/dags/omicron.${now}.dag
    rm -f ./low/omicron.dag*
fi
if [ -e ./std/omicron.dag.dagman.out ]; then
    mv ./std/omicron.dag.dagman.out ./std/dags/omicron.${now}.dag
    rm -f ./std/omicron.dag*
fi
if [ -e ./std2/omicron.dag.dagman.out ]; then
    mv ./std2/omicron.dag.dagman.out ./std2/dags/omicron.${now}.dag
    rm -f ./std2/omicron.dag*
fi
if [ -e ./gw/omicron.dag.dagman.out ]; then
    mv ./gw/omicron.dag.dagman.out ./gw/dags/omicron.${now}.dag
    rm -f ./gw/omicron.dag*
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
awk '$2=="STD2" {print}' ./channels.${IFO} > ./std2/channels.list
awk '$2=="LOW" {print}' ./channels.${IFO} > ./low/channels.list
awk '$2=="GW" {print}' ./channels.${IFO} > ./gw/channels.list

# segments
awk '{print $1-3,$2+3}' ./segments.txt > ./std/segments.txt
awk '{print $1-3,$2+3}' ./segments.txt > ./std2/segments.txt
awk '{print $1-3,$2+3}' ./segments.txt > ./gw/segments.txt
awk '{print $1-24,$2+24}' ./segments.txt > ./low/segments.txt

# LCF
rm -f ./std/triggers/*.ffl
rm -f ./std2/triggers/*.ffl
rm -f ./low/triggers/*.ffl
rm -f ./gw/triggers/*.ffl
cp ./frames.lcf ./std/frames.lcf
cp ./frames.lcf ./std2/frames.lcf
cp ./frames.lcf ./gw/frames.lcf
mv ./frames.lcf ./low/frames.lcf

# Omicron parameters
cd ./std
GetOmicronOptions -c ./channels.list -f ./frames.lcf -d ./triggers >> $logfile 2>&1
cd ../std2
GetOmicronOptions -c ./channels.list -f ./frames.lcf -d ./triggers >> $logfile 2>&1
cd ../low
GetOmicronOptions -c ./channels.list -f ./frames.lcf -d ./triggers >> $logfile 2>&1
cd ../gw
GetOmicronOptions -c ./channels.list -f ./frames.lcf -d ./triggers >> $logfile 2>&1
cd ..

# make dags
GetOmicronDAG -f -d `pwd`/std >> $logfile 2>&1
GetOmicronDAG -f -d `pwd`/std2 >> $logfile 2>&1
GetOmicronDAG -f -d `pwd`/low >> $logfile 2>&1
GetOmicronDAG -f -d `pwd`/gw >> $logfile 2>&1


################################################################################
###########                         GO!                              ###########
################################################################################                 

cd ./std
if [ -e omicron.dag ]; then condor_submit_dag omicron.dag >> $logfile 2>&1; fi
cd ../std2
if [ -e omicron.dag ]; then condor_submit_dag omicron.dag >> $logfile 2>&1; fi
cd ../low
if [ -e omicron.dag ]; then condor_submit_dag omicron.dag >> $logfile 2>&1; fi
cd ../gw
if [ -e omicron.dag ]; then condor_submit_dag omicron.dag >> $logfile 2>&1; fi
cd ..

# cleaning
rm -f ./low/parameters_LOW_*.txt
rm -f ./gw/parameters_GW_*.txt
rm -f ./std/parameters_STD_*.txt
rm -f ./std2/parameters_STD2_*.txt

exit 0