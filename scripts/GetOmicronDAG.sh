#!/bin/bash
#
# GetOmicronDAG.sh
#
# produce a full condor DAG for Omicron processing
#
# Author: Florent Robinet
# florent.robinet@lal.in2p3.fr

echo "This tool is currently not available"
exit 1

printhelp(){
    echo ""
    echo "Usage:"
    echo "GetOmicronDAG -d [WORKING_DIRECTORY]"
    echo " |__ produce a CONDOR DAG in a working directory"
    echo ""
    echo "The working directory MUST contain:"
    echo ""
    echo "1/ A file with time segments to process. The name of this file must be 'segments.txt'"
    echo "   This file must contain 2 columns exactly: [GPS start] [GPS end]"
    echo ""
    echo "2/ Omicron parameter files. The name of these files must be 'parameters_[ID].txt"
    echo "   where [ID] is a unique identifier (an integer for example)"
    echo "   These parameter files must be carefully checked by the user to be sure that"
    echo "   all the options are valid. No check will be performed."
    echo "   GetOmicronDAG will create specific output directories. As a consequence the"
    echo "   output directory option (OUTPUT/DIRECTORY) is not used. In other words,"
    echo "   the user is free to set this option to the value of his choice (it could"
    echo "   even be wrong)"
    echo ""
    echo ""
    echo "OPTIONS:"
    echo "  -d  [WORKING_DIRECTORY]  Path to a working directory"
    echo "                           if this option is not given, the current directory"
    echo "                           is considered"
    echo "  -t [USER_TAG]            Set a user tag for this DAG [USER_TAG]"
    echo "                           This is useful when one wants to generates several"
    echo "                           dags in the same working directory"
    echo "  -f                       Flag to shortcut the user prompt"
    echo "  -h                       print this help"
    echo ""
    echo "Author: Florent Robinet (LAL - Orsay): robinet@lal.in2p3.fr"
} 

##### Check the Omicron environment
if [ -z "$OMICRONROOT" ]; then
    echo "Error: The Omicron environment is not set"
    exit 1
fi
if [ -z "$USER" ]; then
    echo "Error: Please set the environment variable 'USER' with your user name"
    exit 1
fi

##### default options
workdir=`pwd`    # working directory
usertag=""
forceprompt=0

##### read options
while getopts ":d:t:fh" opt; do
    case $opt in
	d)
	    workdir="$OPTARG"
	    ;;
	t)
	    usertag="$OPTARG"
	    ;;
	f)
	    forceprompt=1
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
OPTIND=1

##### check workdir
if [ ! -d $workdir ] ; then
    echo "`basename $0`: the working directory $workdir cannot be found"
    echo "type  'GetOmicronDAG -h'  for help"
    exit 3
fi

##### check segment file
if [ ! -e ${workdir}/segments.txt ] ; then
    echo "`basename $0`: The segment file '${workdir}/segments.txt' cannot be found"
    echo "type  'GetOmicronDAG -h'  for help"
    exit 3
fi

##### check livetime
livetime=`segsum ${workdir}/segments.txt | awk '{print  int($1)}'`
if [ $livetime -eq 0 ] ; then
    echo "`basename $0`: There is no livetime in your segment file '${workdir}/segments.txt'"
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
    echo "`basename $0`: There is no parameter files in '${workdir}/'"
    echo "type  'GetOmicronDAG -h'  for help"
    exit 3
fi

##### check number of channels
for file in ${workdir}/parameters/parameters${usertag}_*.txt; do
    if [ ! -e $file ]; then continue; fi
    channels=`grep DATA $file | grep -m 1 CHANNELS`
    nchannels=`echo $channels | wc -w`
    if [ $nchannels -le 2 ]; then 
	echo "`basename $0`: There is no channel to process in $file"
	exit 3
    fi
    let "nchannels-=2"
    if [ $forceprompt -eq 1 ]; then continue; fi
    if [ $nchannels -gt 10 ]; then 
	echo ""
	echo "You have more than 10 channels to process in $file :"
	echo "$channels"
	echo ""
	echo "Are you sure you want to perform this search?"
	read -p "Press [ENTER] if yes, [CTRL-C] to cancel"
    fi

done

##### preparing .sub
sed -e "s|\[OMICRON_PATH\]|${OMICRONROOT}/${OMICRONCONFIG}|g" \
    -e "s|\[USER\]|${USER}|g" \
    ${OMICRON_SCRIPTS}/omicron.sub > ${workdir}/omicron.sub
rm -fr ${workdir}/omicron${usertag}.dag

##### preparing segment files
echo "*** preparing segment files and jobs"
nseg=0

# loop over parameter files
p=0
while [ $p -lt $nproc ]; do # loop over parameters

    # timing of the proc
    chunkduration=`grep -m1 CHUNKDURATION ${workdir}/parameters/parameters${usertag}_${p}.txt | awk '{print $3}'`
    overlapduration=`grep -m1 OVERLAPDURATION ${workdir}/parameters/parameters${usertag}_${p}.txt | awk '{print $3}'`
    blockduration=`grep -m1 BLOCKDURATION ${workdir}/parameters/parameters${usertag}_${p}.txt | awk '{print $3}'`
    dur=$(( $chunkduration - $overlapduration ))
    ndur=$(( $OMICRON_TRIGGERS_BASE / $dur ))
    duration=$(( $ndur * $dur + $overlapduration ))
    seg_start=`head -1 ${workdir}/segments.txt | awk '{print int($1)}'`
    seg_stop=`awk '/./{line=$0} END{print line}' ${workdir}/segments.txt | awk '{print int($2)}'`
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
    
	    echo "    ** parameters #${p}, segment #${nseg}..."
	    echo "$s $e" > ${workdir}/segments/segments${usertag}_${p}_${nseg}.txt
	    echo "       $s $e ---> livetime = $(( $e - $s )) sec"

	    # Fill DAG with omicron job
	    echo "JOB omicron${usertag}_seg${nseg}_par${p} omicron.sub" >> ${workdir}/omicron${usertag}.dag
	    echo "VARS omicron${usertag}_seg${nseg}_par${p} initialdir=\"${workdir}\" in_segments=\"./segments/segments${usertag}_${p}_${nseg}.txt\" in_parameters=\"./parameters/parameters${usertag}_${p}.txt\"" >> ${workdir}/omicron${usertag}.dag

	    # increment
	    let "seg+=1"
	    let "nseg+=1"
 	    ss=$(( $ss - $overlapduration ))
	done

    done < ${workdir}/segments.txt

    let "p+=1"
done

exit 0
