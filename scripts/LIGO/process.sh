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

################################################################################
###########                    directory prep                        ###########
################################################################################

# channel lists
awk '$1>1024 {print}'  ./channels.${IFO} > ./high/channel.list
awk '$1<=1024 {print}' ./channels.${IFO} > ./low/channel.list

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
	if [ $(( $tstop -$tstart )) -lt $duration ]; then
	    echo "`date -u`: the new segment is too short to be processed. Wait longer..." >> $logfile
	    exit 0
	fi
	echo "`date -u` create new segment" >> ${worklog}
        echo "$tstart $tstop" > ${workdir}/segments.txt
    fi
fi


#############################################
# type-dependent parameters (for GetOmicronOption)
if [ "$type" = "high" ]; then
    duration=64
    overlap=4
    awk '$1>1024 {print}' ${proddir}/channel.list > ${workdir}/channel.list
elif [ "$type" = "low" ]; then
    duration=8192
    overlap=1280
    awk '$1<=1024 {print}' ${proddir}/channel.list > ${workdir}/channel.list
else
    echo "Type $type is not supported"
    exit 2
fi


mkdir -p ${workdir}/triggers

# Cleaning
echo "Cleaning anything loder than 4 days..." >> ${worklog}
now_4days=$(( $now - 345600 ))
for file in ${workdir}/log.1*.txt; do
    fileend="${file##*/log.}"
    gps="${fileend%.*}"
    if [ $gps -lt $now_4days ]; then rm -fr $file; fi
done
for file in ${workdir}/omicron.1*.dag; do
    fileend="${file##*/omicron.}"
    gps="${fileend%.*}"
    if [ $gps -lt $now_4days ]; then rm -fr $file; fi
done
find ${workdir}/logs -type f -mtime +4 -exec rm {} \; 

# check if previous batch is still running
if [ -e ${workdir}/omicron.dag.lock ]; then
    echo "`date -u` condor is already running" >> ${worklog}
    exit 2
fi
if [ -e ${workdir}/omicron.dag.dagman.out ]; then
    mv ${workdir}/omicron.dag.dagman.out ${workdir}/omicron.$now.dag
    rm -f ${workdir}/omicron.dag*
fi
rmdir ${workdir}/triggers/* >> /dev/null 2>&1

# timing
tstop=$(( $now - $delay ))
nseg=$(( ( $tmax - $overlap ) / ( $duration - $overlap ) ))
tstart=$(( $tstop - $nseg * ( $duration - $overlap ) - $overlap ))
if [ ! -e ${workdir}/segments.txt ]; then
    echo "`date -u` no previous segment --> start a new one" >> ${worklog}
    echo "$tstart $tstop" > ${workdir}/segments.txt
else
    prev_stop=`head -n 1 ${workdir}/segments.txt | awk '{print $2}'`
    if [ $(( $tstop - $prev_stop )) -gt $tmax ]; then
	echo "`date -u` the last processed segment is too old --> start a new one" >> ${worklog}
	echo "$tstart $tstop" > ${workdir}/segments.txt
    else
	tstart=$(( $prev_stop - $overlap ))
	nseg=$(( ( $tstop - $tstart - $overlap ) / ( $duration - $overlap ) ))
	tstop=$(( $tstart + $nseg * ( $duration - $overlap ) + $overlap ))
	if [ $(( $tstop -$tstart )) -lt $duration ]; then
	    echo "`date -u` the new segment is too short to be processed. Wait longer..." >> ${worklog}
	    exit 0
	fi
	echo "`date -u` create new segment" >> ${worklog}
        echo "$tstart $tstop" > ${workdir}/segments.txt
    fi
fi

# LCF
echo "`date -u` LCF file $tstart -> $tstop" >> ${worklog}
if [ "$IFO" = "L1" ]; then
    ligo_data_find -o L -l -t ${IFO}_C -u file -s $tstart -e $tstop > ${workdir}/frames.lcf
    ligo_data_find -o L -l -t ${IFO}_ER_C00_L1 -u file -s $tstart -e $tstop > ${workdir}/framesgw.lcf
    if [ ! -s ${workdir}/frames.lcf ]; then
	ligo_data_find -o L -l -t ${IFO}_R -u file -s $tstart -e $tstop > ${workdir}/frames.lcf
    fi
elif [ "$IFO" = "H1" ]; then
    ligo_data_find -o H -l -t ${IFO}_C -u file -s $tstart -e $tstop > ${workdir}/frames.lcf
    ligo_data_find -o H -l -t ${IFO}_ER_C00_L1 -u file -s $tstart -e $tstop > ${workdir}/framesgw.lcf
    if [ ! -s ${workdir}/frames.lcf ]; then
        ligo_data_find -o H -l -t ${IFO}_R -u file -s $tstart -e $tstop > ${workdir}/frames.lcf
    fi
else
    echo "`date -u` Unsupported IFO: $IFO" >> ${worklog}
    exit 1
fi

# FFL
if [ -s ${workdir}/frames.lcf ]; then
    echo "`date -u` LCF file to FFL file" >> ${worklog}
    lalcache2ffl ${workdir}/frames.lcf > ${workdir}/frames.ffl
else
    echo "`date -u` No frame cache files for the last day" >> ${worklog}
    exit 1
fi
if [ -s ${workdir}/framesgw.lcf ]; then
    echo "`date -u` LCF file to FFL file" >> ${worklog}
    lalcache2ffl ${workdir}/framesgw.lcf > ${workdir}/framesgw.ffl
fi

# Omicron option files
echo "`date -u` GetOmicronOptions" >> ${worklog}
cd ${workdir}
GetOmicronOptions -o ${workdir}/parameters.template -c ${workdir}/channel.list >> ${worklog}
cd `dirname $0`

# GW option file
if [ "$type" = "high" ]; then
    if [ -s ${workdir}/framesgw.ffl ]; then
        cp -f ${workdir}/parameters.gw ${workdir}/parameters_gw.txt
    else
        rm -f ${workdir}/parameters_gw.txt
    fi
fi

# Omicron DAG
echo "`date -u` GetOmicronDAG" >> ${worklog}
GetOmicronDAG -f -d ${workdir} >> ${worklog}

# cleaning
rm -f ${workdir}/frames.lcf
rm -f ${workdir}/parameters_*Hz_*.txt


# Submit DAG
cd ${workdir}
if [ -e ./omicron.dag ]; then
    echo "`date -u` Submit omicron DAG..." >> ${worklog}
    condor_submit_dag -f ./omicron.dag >> ${worklog} 2>&1
else
    echo "`date -u` No DAG to submit!" >> ${worklog}
    exit 1
fi

exit 0


