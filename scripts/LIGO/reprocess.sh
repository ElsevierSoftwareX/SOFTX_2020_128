#!/bin/sh
#
# process.sh
#
# Author: Florent Robinet
# robinet@lal.in2p3.fr

if [ ! $# -eq 3 ]; then
    echo "./reprocess.sh [CHANNEL_LIST] [GPS_START] [GPS_END]"
    exit 1
fi
start=$2
end=$3
chanfile=$1
if [ -d ./reprocess ]; then
    echo "remove ./reprocess before running this script"
    exit 1
fi
if [ ! -e ${chanfile} ]; then
    echo "${chanfile} is missing"
    exit 1
fi

###### user parameters

hn=`hostname -d`
if [ "$hn" = "ligo-wa.caltech.edu" ]; then IFO="H1"
else  IFO="L1"; fi
data_type="${IFO}_C ${IFO}_R"
gw_type="${IFO}_ER_C00_L1"

######################

################################################################################
###########                       check inputs                        ########### 
################################################################################                               

# GWOLLUM environment
source /home/detchar/opt/virgosoft/environment.v1r3.sh "" >> /dev/null

# directories
mkdir -p ./reprocess/std/triggers
mkdir -p ./reprocess/std2/triggers
mkdir -p ./reprocess/low/triggers
mkdir -p ./reprocess/gw/triggers
mkdir -p ./reprocess/fine/triggers

################################################################################
###########                          input                           ###########
################################################################################

# timing
echo "get timing..."
tstart=$(( $start / 1000 ))
if [ $(( $tstart % 2 )) -eq 1 ]; then tstart=$(( $tstart - 1 )); fi
tstart="${tstart}000"
tstop=$(( $end / 1000 ))
if [ $(( $tstop % 2 )) -eq 1 ]; then tstop=$(( $tstop - 1 )); fi
tstop="${tstop}000"
if [ $tstart -ge $tstop ]; then
    echo "nothing to process"
    exit 1
fi
echo "$tstart $tstop" >  ./reprocess/segments.txt
echo "segment to process: $tstart $tstop"

# raw data
echo "get LCF file for raw data..."
for type in $data_type; do
    ligo_data_find -o ${IFO:0:1} -l -t $type -u file -s $(( $tstart - 100 )) -e $(( $tstop + 100 )) 1>./reprocess/frames.lcf
    if [ -s ./reprocess/frames.lcf ]; then break; fi
done
if [ ! -s ./reprocess/frames.lcf ]; then
    echo "data are missing for this segment"
    rm -f ./reprocess/frames.lcf
    exit 0
fi

# gw data
echo "get LCF file for gw data"
for type in $gw_type; do
    ligo_data_find -o ${IFO:0:1} -l -t $type -u file -s $(( $tstart - 100 )) -e $(( $tstop + 100 )) 1>./reprocess/frames.gw.lcf
    if [ -s ./reprocess/frames.gw.lcf ]; then break; fi
done

################################################################################
###########                    directory prep                        ###########
################################################################################

# channel lists
awk '$2=="STD" {print}' ${chanfile} > ./reprocess/std/channels.list
awk '$2=="STD2" {print}' ${chanfile} > ./reprocess/std2/channels.list
awk '$2=="LOW" {print}' ${chanfile} > ./reprocess/low/channels.list
awk '$2=="GW" {print}' ${chanfile} > ./reprocess/gw/channels.list
awk '$2=="FINE" {print}' ${chanfile} > ./reprocess/fine/channels.list

# segments (overlap from GetOmicronOptions)
awk '{print $1-3,$2+3}' ./reprocess/segments.txt > ./reprocess/std/segments.txt
awk '{print $1-3,$2+3}' ./reprocess/segments.txt > ./reprocess/std2/segments.txt
awk '{print $1-7,$2+7}' ./reprocess/segments.txt > ./reprocess/gw/segments.txt
awk '{print $1-24,$2+24}' ./reprocess/segments.txt > ./reprocess/low/segments.txt
awk '{print $1-7,$2+7}' ./reprocess/segments.txt > ./reprocess/fine/segments.txt

# LCF
cp -f ./reprocess/frames.lcf ./reprocess/std/frames.lcf
cp -f ./reprocess/frames.lcf ./reprocess/std2/frames.lcf
cp -f ./reprocess/frames.gw.lcf ./reprocess/gw/frames.lcf
cp -f ./reprocess/frames.lcf ./reprocess/fine/frames.lcf
cp -f ./reprocess/frames.lcf ./reprocess/low/frames.lcf

# Omicron parameters
if [ -s ./reprocess/frames.lcf ]; then
    cd ./reprocess/std
    GetOmicronOptions -c ./channels.list -f ./frames.lcf -d ./triggers -X
    cd ../../reprocess/std2
    GetOmicronOptions -c ./channels.list -f ./frames.lcf -d ./triggers -X
    cd ../../reprocess/low
    GetOmicronOptions -c ./channels.list -f ./frames.lcf -d ./triggers -X
    cd ../../reprocess/fine
    GetOmicronOptions -c ./channels.list -f ./frames.lcf -d ./triggers -X
    cd ../..
fi
if [ -s ./reprocess/frames.gw.lcf ]; then
    cd ./reprocess/gw
    GetOmicronOptions -c ./channels.list -f ./frames.lcf -d ./triggers -X
    cd ../..
fi

# make dags
GetOmicronDAG -f -d `pwd`/reprocess/std
GetOmicronDAG -f -d `pwd`/reprocess/std2
GetOmicronDAG -f -d `pwd`/reprocess/low
GetOmicronDAG -f -d `pwd`/reprocess/gw
GetOmicronDAG -f -d `pwd`/reprocess/fine


# cleaning
rm -f ./reprocess/low/parameters_LOW_*.txt
rm -f ./reprocess/gw/parameters_GW_*.txt
rm -f ./reprocess/gw/parameters_FINE_*.txt
rm -f ./reprocess/std/parameters_STD_*.txt
rm -f ./reprocess/std2/parameters_STD2_*.txt

exit 0
