#!/bin/sh
#
# process.sh
#
# Author: Florent Robinet
# robinet@lal.in2p3.fr
cd `dirname $0`

###### user parameters
delay=300 # do not look at data after now-delay
tmax=16000 # delay after which a new segment is started to get back on track
PRODTYPES="STD1 STD2 LOW1 GW FINE"
data_type="${IFO}_C ${IFO}_R"
gw_type="${IFO}_ER_C00_L1"
######################

### THIS FUNCTION HAS TO MATCH GetOmicronOption
GetOverlap(){
    if [ "$1" = "STD1" ]; then overlap=4
    elif [ "$1" = "STD2" ]; then overlap=4
    elif [ "$1" = "LOW1" ]; then overlap=112
    elif [ "$1" = "GW" ]; then overlap=4
    elif [ "$1" = "FINE" ]; then overlap=4
    else overlap=0;
}
    

# LLO of LHO?
hn=`hostname -d`
if [ "$hn" = "ligo-wa.caltech.edu" ]; then IFO="H1"
else  IFO="L1"; fi

################################################################################
###########                         GENERAL                          ########### 
################################################################################                               

# GWOLLUM environment
source /home/detchar/opt/virgosoft/environment.v2r1.sh "" >> /dev/null

# for safety, release held jobs
condor_release -all >> /dev/null 2>&1


# start log file
mkdir -p ./logs
find ./logs -type f -mtime +4 -exec rm {} \; >> /dev/null 2>&1
now=`tconvert now`
logfile=`pwd`/logs/process.${now}.txt
echo "Local time: `date`" > $logfile
echo "UTC time: `date -u`" >> $logfile

# check channels
if [ ! -e ./channels.${IFO} ]; then
    echo "`date -u`: ./channels.${IFO} is missing" >> $logfile
    exit 1
fi

################################################################################
###########                       PROD TYPES                     ##########
################################################################################

for ptype in $PRODTYPES; do
    echo ""
    echo ""
    echo "`date -u` *******************   ${ptype} CONFIG   *******************" >> $logfile

    # directories
    mkdir -p ./${ptype}/triggers ./${ptype}/dags ./${ptype}/logs;

    # cleaning old files
    find ./${ptype}/logs -type f -mtime +4 -exec rm {} \; >> /dev/null 2>&1
    find ./${ptype}/dags -type f -mtime +4 -exec rm {} \; >> /dev/null 2>&1
    find ./${ptype} -type f -mtime +4 -exec rm {} \; >> /dev/null 2>&1

    # check if previous batch is still running
    echo "`date -u`: check if previous batch is still running..." >> $logfile
    if [ -e ./${ptype}/omicron.dag.lock ]; then
	echo "`date -u` condor (low) is already running" >> $logfile
	continue
    fi

    # archive dagman
    if [ -e ./${ptype}/omicron.dag.dagman.out ]; then
	mv ./${ptype}/omicron.dag.dagman.out ./${ptype}/dags/omicron.${now}.dag
    fi

    # more cleaning
    rm -f ./${ptype}/omicron.dag*
    rmdir ./${ptype}/triggers/* >> /dev/null 2>&1

    # channel list
    awk -v var="$ptype" '$2==var {print}' ./channels.${IFO} > ./$ptype/channels.list
    if [ ! -s ./$ptype/channels.list ]; then
	echo "`date -u`: no channels" >> $logfile
	continue
    fi

    # timing
    tstop=$(( ($now - $delay) / 1000 ))
    if [ $(( $tstop % 2 )) -eq 1 ]; then tstop=$(( $tstop - 1 )); fi
    tstop="${tstop}000"
    tstart=$(( $tstop - 2000 ))
    if [ ! -e ./${ptype}/segments.ref ]; then
	echo "`date -u`: no previous segment --> start a new one [$tstart; $tstop]" >> $logfile
	echo "$tstart $tstop" > ./${ptype}/segments.tmp
    else
	prev_stop=`head -n 1 ./${ptype}/segments.txt | awk '{print $2}'`
	if [ $(( $tstop - $prev_stop )) -gt $tmax ]; then
	    echo "`date -u`: the last processed segment is too old --> start a new one [$tstart; $tstop]" >> $logfile
	    echo "$tstart $tstop" > ./${ptype}/segments.tmp
	else
	    tstart=$prev_stop
	    echo "`date -u`: create new segment [$tstart; $tstop]" >> $logfile
            echo "$tstart $tstop" > ./${ptype}/segments.tmp
	fi
    fi
    if [ `segsum ./${ptype}/segments.tmp` -lt 2000 ]; then
	echo "`date -u`: the new segment is too short, start again later..." >> $logfile
	rm -f ./${ptype}/segments.tmp
	continue
    fi
    mv ./${ptype}/segments.tmp ./${ptype}/segments.ref

    # LCF
    echo "`date -u`: get LCF file" >> $logfile
    if [ "${ptype}" = "GW" ]; then
	for type in $gw_type; do
	    gw_data_find -o ${IFO:0:1} -l -t $type -u file -s $(( $tstart - 100 )) -e $(( $tstop + 100 )) 1>./${ptype}/frames.lcf 2>> $logfile
	    if [ -s ./${ptype}/frames.lcf ]; then break; fi
	done
    else
	for type in $data_type; do
	    gw_data_find -o ${IFO:0:1} -l -t $type -u file -s $(( $tstart - 200 )) -e $(( $tstop + 200 )) 1>./frames.lcf 2>> $logfile
	    if [ -s ./${ptype}/frames.lcf ]; then break; fi
	done
    fi
    
    if [ ! -s ./${ptype}/frames.lcf ]; then
	echo "`date -u`: data are missing for this segment" >> $logfile
	rm -f ./${ptype}/frames.lcf
	continue
    fi

    # segment to process
    GetOverlap
    awk -v var="$(( $overlap / 2 ))" '{print $1-var,$2+var}' ./${ptype}/segments.ref > ./${ptype}/segments.txt

    # generate option files
    cd ./${ptype}
    GetOmicronOptions -c ./channels.list -f ./frames.lcf -d ./triggers -X >> $logfile 2>&1
    cd ..

    # make dags
    GetOmicronDAG -d `pwd`/${ptype} >> $logfile 2>&1
    rm -f ./${ptype}/parameters_${ptype}_*.txt

    # GO!
    cd ./${ptype}
    #if [ -e omicron.dag ]; then condor_submit_dag omicron.dag >> $logfile 2>&1; fi
    cd ..
    
done

exit 0
