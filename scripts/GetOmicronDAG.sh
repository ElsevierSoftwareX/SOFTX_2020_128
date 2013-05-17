#!/bin/bash
#
# GetOmicronDAG.sh
#
# produce a full condor DAG for Omicron processing
#
# Author: Florent Robinet
# florent.robinet@lal.in2p3.fr

##### default options
workdir=`pwd`    # working directory
usertag=""
merging=0

printhelp(){
    echo ""
    echo "Usage:"
    echo "GetOmicronDAG -d [WORKING_DIRECTORY]"
    echo ""
    echo ""
    echo "  -d  [WORKING_DIRECTORY]  path to the working directory"
    echo "                           if this option is not given, the current directory"
    echo "                           is considered"
    echo ""
    echo ""
    echo "The working directory must contain:"
    echo ""
    echo "1/ A segment file to process. The name of this file must be 'segments.txt'"
    echo "   This file must contain 2 columns exactly: [GPS start] [GPS end]"
    echo ""
    echo "2/ Omicron parameter files. The name of these files must be 'parameters_[ID].txt"
    echo "   where [ID] is a unique identifier (an integer for example)"
    echo "   These parameter files must be carefully checked by the user to be sure that"
    echo "   all the options are valid"
    echo "   This script will create specific output directories. As a consequence the"
    echo "   output directory option (OUTPUT/DIRECTORY) is not used. In other words,"
    echo "   the user is free to set this option to the value of his choice (it could"
    echo "   even be wrong)"
    echo ""

    echo ""
    echo "  -t [USER_TAG]       Set a user tag for this DAG [USER_TAG]"
    echo ""
    echo ""
    echo "  -m                  Flag to activate merging of trigger files"
    echo "                      This option is only available for ROOT file format"
    echo "                      --- NOT AVAILABLE OPTION ---"
    echo ""
    echo ""
    echo "  -h                  Print this help"
    echo ""
} 

##### Check the Omicron environment
if [[ -z "$OMICRONROOT" ]]; then
    echo "Error: The Omicron environment is not set"
    exit 1
fi
if [[ -z "$USER" ]]; then
    echo "Error: Please set the environment variable 'USER' with your user name"
    exit 1
fi

##### read options
while getopts ":d:t:mh" opt; do
    case $opt in
	d)
	    workdir="$OPTARG"
	    ;;
	t)
	    usertag="$OPTARG"
	    ;;
	m)
	    merging=1
	    ;;
	h)
	    printhelp
	    exit 0
	    ;;
	\?)
	    echo "Invalid option: -$OPTARG"
	    echo "type  'GetOmicronDAG -h'  for help"
	    exit 2
	    ;;
    esac
done

##### check workdir
if [ ! -d $workdir ] ; then
    echo "Invalid option: the working directory $workdir cannot be found"
    echo "type  'GetOmicronDAG -h'  for help"
    exit 3
fi

##### check segment file
if [ ! -e ${workdir}/segments.txt ] ; then
    echo "The segment file '${workdir}/segments.txt' cannot be found"
    echo "type  'GetOmicronDAG -h'  for help"
    exit 3
fi

##### check livetime
livetime=`segsum ${workdir}/segments.txt | awk '{print  int($1)}'`
if [ $livetime -eq 0 ] ; then
    echo "There is no livetime in your segment file '${workdir}/segments.txt'"
    echo "type  'GetOmicronDAG -h'  for help"
    exit 3
fi

##### create sub-directories
mkdir -p ${workdir}/triggers   # trigger directory
mkdir -p ${workdir}/segments   # segment directory
mkdir -p ${workdir}/parameters # parameter directory
mkdir -p ${workdir}/logs       # log directory

##### preparing parameter files
echo "*** preparing parameter files"
nproc=0;
for file in ${workdir}/parameters_*.txt; do
    if [ -e $file ]; then

	# remove output directory option
	sed '/OUTPUT[ \t]*DIRECTORY/d' $file > ${workdir}/parameters/parameters${usertag}_${nproc}.txt

	# add new output directory option
	echo "" >> ${workdir}/parameters/parameters${usertag}_${nproc}.txt
	echo "OUTPUT  DIRECTORY  ${workdir}/triggers" >> ${workdir}/parameters/parameters${usertag}_${nproc}.txt

	echo "parameters${usertag}_${nproc}.txt = $file"
	let "nproc+=1"
    fi
done

if [ $nproc -eq 0 ] ; then
    echo "There is no parameter files in '${workdir}/'"
    echo "type  'GetOmicronDAG -h'  for help"
    exit 3
fi

##### check number of channels
for file in ${workdir}/parameters/parameters${usertag}_*.txt; do
    channels=`grep DATA $file | grep -m 1 CHANNELS`
    nchannels=`echo $channels | wc -w`
    if [ $nchannels -le 2 ]; then 
	echo "There is no channel to process in $file"
	exit 3
    fi
    let "nchannels-=2"
    if [ $nchannels -gt 10 ]; then 
	echo ""
	echo "You have more than 10 channels to process in $file :"
	echo "$channels"
	echo ""
	echo "Are you sure you want to perform this search?"
	read -p "Press [ENTER] if yes, [CTRL-C] to cancel"
    fi

done

##### check merging
fileformat=`grep OUTPUT ${workdir}/parameters/parameters${usertag}_0.txt | grep -m1 FORMAT | awk '{print $3}'`
if [ ! "$fileformat" = "root" ]; then merging=0; fi

##### preparing segment files
echo "*** preparing segment files"
nseg=0
#FIXME : this is wrong if durations differ
chunkduration=`grep -m1 CHUNKDURATION ${workdir}/parameters/parameters${usertag}_0.txt | awk '{print $3}'`
overlapduration=`grep -m1 OVERLAPDURATION ${workdir}/parameters/parameters${usertag}_0.txt | awk '{print $3}'`
blockduration=`grep -m1 BLOCKDURATION ${workdir}/parameters/parameters${usertag}_0.txt | awk '{print $3}'`
dur=$(( $chunkduration - $overlapduration ))
ndur=$(( $OMICRON_TRIGGERS_BASE / $dur ))
duration=$(( $ndur * $dur + $overlapduration ))
seg_start=`head -1 ${workdir}/segments.txt | awk '{print $1}'`
seg_stop=`tail -1 ${workdir}/segments.txt | awk '{print $2}'`
seg_start_base=$(( $seg_start / $OMICRON_TRIGGERS_BASE ))
seg_stop_base=$(( $seg_stop / $OMICRON_TRIGGERS_BASE ))

# loop over segments
while read line; do
    if [ "$line" = "" ]; then continue; fi
    
    # timing
    dur=`echo $line | awk '{print int($2-$1)}'`
    ss=`echo $line | awk '{print int($1)}'`
    ee=`echo $line | awk '{print int($2)}'`

    # too short
    if [ $dur -lt $blockduration ]; then continue; fi

    # get number of sub-segments
    nsubseg=$(( $dur / $duration + 1 ))
    
    # loop over sub-segments
    seg=0
    while [ $seg -lt $nsubseg ]; do
	s=$(( $ss + $seg * $duration ))
	e=$(( $s + $duration ))
	if [ $e -gt $ee ]; then e=$ee; fi
    
	echo "$s $e" > ${workdir}/segments/segments${usertag}_${nseg}.txt
	let "seg+=1"
	let "nseg+=1"
 	ss=$(( $ss - $overlapduration ))
    done

done < ${workdir}/segments.txt

##### preparing segment files
sed -e "s|\[OMICRON_PATH\]|${OMICRONROOT}/${OMICRONCONFIG}|g" \
    -e "s|\[USER\]|${USER}|g" \
    ${OMICRON_SCRIPTS}/omicron.sub > ${workdir}/omicron.sub

##### make omicron jobs
echo "*** make omicron jobs"
n=0
rm -fr ${workdir}/omicron${usertag}.dag
while [ $n -lt $nseg ]; do # loop over segments
    p=0
    while [ $p -lt $nproc ]; do # loop over parameters
	echo "$nseg $n $p"
	echo "JOB omicron${usertag}_seg${n}_par${p} omicron.sub" >> ${workdir}/omicron${usertag}.dag
	echo "VARS omicron${usertag}_seg${n}_par${p} initialdir=\"${workdir}\" in_segments=\"./segments/segments${usertag}_${n}.txt\" in_parameters=\"./parameters/parameters${usertag}_${p}.txt\"" >> ${workdir}/omicron${usertag}.dag

	let "p+=1"
    done
    let "n+=1"
done

##### make merging jobs
# TO BE DONE




exit 0