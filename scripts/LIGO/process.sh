#!/bin/sh
#
# process.sh
#
# Author: Florent Robinet
# robinet@lal.in2p3.fr
cd `dirname $0`

# LLO of LHO?
hn=`hostname -d`
if [ "$hn" = "ligo-wa.caltech.edu" ]; then 
    IFO="H1"
    LIGO_DATAFIND_SERVER="10.12.0.49:80"
else
    IFO="L1"; 
    LIGO_DATAFIND_SERVER="10.13.20.73:80"
fi
export LIGO_DATAFIND_SERVER

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
    fi
}

################################################################################
###########                         GENERAL                          ########### 
################################################################################                               

# environment
source /home/detchar/opt/virgosoft/environment.v2r1.sh "" >> /dev/null
source /home/detchar/opt/gwpysoft/etc/gwpy-user-env.sh

# for safety, release held jobs
condor_release -all >> /dev/null 2>&1

basedir=`pwd`
here=`pwd`

# special case for offline processing
if [ $# -eq 2 ]; then
    off_start=$1
    off_stop=$2
    basedir="`pwd`/${1}-${2}"
fi


# start log file
mkdir -p ${basedir}/logs
find ${basedir}/logs -type f -mtime +4 -exec rm {} \; >> /dev/null 2>&1
now=`tconvert now`
logfile=${basedir}/logs/process.${now}.txt
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
    echo "" >> $logfile
    echo "" >> $logfile
    echo "`date -u` *******************   ${ptype} CONFIG   *******************" >> $logfile

    # directories
    mkdir -p ${basedir}/${ptype}/triggers ${basedir}/${ptype}/dags ${basedir}/${ptype}/logs;

    # cleaning old files
    find ${basedir}/${ptype}/logs -type f -mtime +4 -exec rm {} \; >> /dev/null 2>&1
    find ${basedir}/${ptype}/dags -type f -mtime +4 -exec rm {} \; >> /dev/null 2>&1
    find ${basedir}/${ptype} -type f -mtime +2 -exec rm {} \; >> /dev/null 2>&1
    find ${basedir}/${ptype}/triggers -type f -mtime +2 -exec rm {} \; >> /dev/null 2>&1

    # check if previous batch is still running
    echo "`date -u`: check if previous batch is still running..." >> $logfile
    if [ -e ${basedir}/${ptype}/omicron.dag.lock ]; then
	echo "`date -u` condor is already running" >> $logfile
	continue
    fi

    # archive dagman
    echo "`date -u`: archive previous dagman..." >> $logfile
    if [ -e ${basedir}/${ptype}/omicron.dag.dagman.out ]; then
	mv ${basedir}/${ptype}/omicron.dag.dagman.out ${basedir}/${ptype}/dags/omicron.${now}.dag
    fi

    # more cleaning
    rm -f ${basedir}/${ptype}/omicron.dag*
    rmdir ${basedir}/${ptype}/triggers/* >> /dev/null 2>&1

    # channel list
    echo "`date -u`: make channel list..." >> $logfile
    grep -w "$ptype" ./channels.${IFO} > ${basedir}/$ptype/channels.list
    if [ ! -s ${basedir}/$ptype/channels.list ]; then
	echo "`date -u`: no channels" >> $logfile
	continue
    fi

    # reference timing
    echo "`date -u`: make reference timing..." >> $logfile
    if [ $# -eq 2 ]; then
	tstart=$1
	tstop=$2
	echo "`date -u`: create new segment [$tstart; $tstop]" >> $logfile
	echo "$tstart $tstop" > ${basedir}/${ptype}/segments.tmp

    else
	tstop=$(( ($now - $delay) / 1000 ))
	if [ $(( $tstop % 2 )) -eq 1 ]; then tstop=$(( $tstop - 1 )); fi
	tstop="${tstop}000"
	tstart=$(( $tstop - 2000 ))
	if [ ! -e ${basedir}/${ptype}/segments.ref ]; then
	    echo "`date -u`: no previous segment --> start a new one [$tstart; $tstop]" >> $logfile
	    echo "$tstart $tstop" > ${basedir}/${ptype}/segments.tmp
	else
	    prev_stop=`head -n 1 ${basedir}/${ptype}/segments.ref | awk '{print $2}'`
	    if [ $(( $tstop - $prev_stop )) -gt $tmax ]; then
		echo "`date -u`: the last processed segment is too old --> start a new one [$tstart; $tstop]" >> $logfile
		echo "$tstart $tstop" > ${basedir}/${ptype}/segments.tmp
	    else
		tstart=$prev_stop
		echo "`date -u`: create new segment [$tstart; $tstop]" >> $logfile
		echo "$tstart $tstop" > ${basedir}/${ptype}/segments.tmp
	    fi
	fi
    fi
    
    if [ `segsum ${basedir}/${ptype}/segments.tmp` -lt 2000 ]; then
	echo "`date -u`: the new segment is too short, start again later..." >> $logfile
	rm -f ${basedir}/${ptype}/segments.tmp
	continue
    fi
    mv ${basedir}/${ptype}/segments.tmp ${basedir}/${ptype}/segments.ref

    # LCF
    echo "`date -u`: get LCF file" >> $logfile
    if [ "${ptype}" = "GW" ]; then
	for type in $gw_type; do
	    gw_data_find -o ${IFO:0:1} -l -t $type -u file -s $(( $tstart - 100 )) -e $(( $tstop + 100 )) 1>${basedir}/${ptype}/frames.lcf 2>> $logfile
	    if [ -s ${basedir}/${ptype}/frames.lcf ]; then break; fi
	done
    else
	for type in $data_type; do
	    gw_data_find -o ${IFO:0:1} -l -t $type -u file -s $(( $tstart - 200 )) -e $(( $tstop + 200 )) 1>${basedir}/${ptype}/frames.lcf 2>> $logfile
	    if [ -s ${basedir}/${ptype}/frames.lcf ]; then break; fi
	done
    fi
    
    if [ ! -s ${basedir}/${ptype}/frames.lcf ]; then
	echo "`date -u`: data are missing for this segment" >> $logfile
	rm -f ${basedir}/${ptype}/frames.lcf
	continue
    fi

    # segments to process
    echo "`date -u`: make proc segments..." >> $logfile
    GetOverlap $ptype
    awk -v var="$(( $overlap / 2 ))" '{print $1-var,$2+var}' ${basedir}/${ptype}/segments.ref > ${basedir}/${ptype}/segments.tmp
    if [ "${ptype}" = "GW" ]; then
	if [ $# -eq 2 ]; then 
	    ligolw_segment_query_dqsegdb --segment-url=https://dqsegdb5.phy.syr.edu --query-segments --gps-start-time $(( $tstart - $overlap / 2 )) --gps-end-time $(( $tstop + $overlap / 2 )) --include-segments="${IFO}:DMT-CALIBRATED" | ligolw_print -t segment -c start_time -c end_time -d ' ' 1>${basedir}/${ptype}/segments.OK 2>> $logfile
	else
            /home/detchar/bin/gwdq-get-ligo-segments -c ${IFO}:GDS-CALIB_STATE_VECTOR -t ${IFO}_llhoft -b 2,3,4 $(( $tstart - $overlap / 2 )) $(( $tstop + $overlap / 2 )) 1>${basedir}/${ptype}/segments.OK 2>> $logfile
	fi
    else
        if [ $# -eq 2 ]; then 
	    ligolw_segment_query_dqsegdb --segment-url=https://dqsegdb5.phy.syr.edu --query-segments --gps-start-time $(( $tstart - $overlap / 2 )) --gps-end-time $(( $tstop + $overlap / 2 )) --include-segments="${IFO}:DMT-UP" | ligolw_print -t segment -c start_time -c end_time -d ' ' 1>${basedir}/${ptype}/segments.OK 2>> $logfile
        else
            /home/detchar/bin/gwdq-get-ligo-segments -c ${IFO}:GDS-CALIB_STATE_VECTOR -t ${IFO}_llhoft -b 2 $(( $tstart - $overlap / 2 )) $(( $tstop + $overlap / 2 )) 1>${basedir}/${ptype}/segments.OK 2>> $logfile
	fi
    fi
    segexpr 'intersection('${basedir}'/'${ptype}'/segments.tmp,'${basedir}'/'${ptype}'/segments.OK)' > ${basedir}/${ptype}/segments.txt
    rm -f ${basedir}/${ptype}/segments.tmp
    if [ ! -s ${basedir}/${ptype}/segments.txt ]; then
        echo "`date -u`: no segment to process" >> $logfile
        continue
    fi


    # generate option files
    cd ${basedir}/${ptype}
    GetOmicronOptions -c ./channels.list -f ./frames.lcf -d ./triggers -X >> $logfile 2>&1
    cd $here

    # make dags
    GetOmicronDAG -d ${basedir}/${ptype} >> $logfile 2>&1
    rm -f ${basedir}/${ptype}/parameters_${ptype}_*.txt

    # GO!
    cd ${basedir}/${ptype}
    if [ -e omicron.dag ]; then condor_submit_dag omicron.dag >> $logfile 2>&1; fi
    cd $here
done

exit 0
