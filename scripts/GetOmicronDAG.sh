#!/bin/bash
#
# GetOmicronDAG.sh
#
# produce a full condor DAG for Omicron processing
#
# Author: Florent Robinet
# florent.robinet@lal.in2p3.fr

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
    echo "   output directory option (OUTPUT/DIRECTORY) is irrelevant. In other words,"
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

##### read options
while getopts ":d:t:h" opt; do
    case $opt in
	d)
	    workdir="$OPTARG"
	    ;;
	t)
	    usertag="$OPTARG"
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

##### preparing .sub
sed -e "s|\[OMICRON_PATH\]|${OMICRONROOT}/${OMICRONCONFIG}|g" \
    -e "s|\[USER\]|${USER}|g" \
    ${OMICRON_SCRIPTS}/omicron.sub > ${workdir}/omicron.sub
rm -fr ${workdir}/omicron${usertag}.dag

##### preparing segment files
echo "*** preparing segment files and jobs"

# loop over parameter files
p=0
while [ $p -lt $nproc ]; do # loop over parameters

    # init
    nseg=0
    durcum=0
    rm -f ${workdir}/segments/segments${usertag}_${p}_*.txt
  
    # overlap
    overlapduration=`grep -m1 OVERLAPDURATION ${workdir}/parameters/parameters${usertag}_${p}.txt | awk '{print $3}'`

    # loop over segments
    while read line; do
	if [ "$line" = "" ]; then continue; fi
	
        # timing
	dur=`echo $line | awk '{print int($2-$1)}'`
	ss=`echo $line | awk '{print int($1)}'`
	ee=`echo $line | awk '{print int($2)}'`
	durcum=$(( $durcum + $dur ))

	# make jobs for this segment
	sss=$ss
	while [ $sss -lt $(( $ee - $overlapduration )) ]; do
	    eee=$(( $sss + $OMICRON_TRIGGERS_BASE ))
	    if [ $eee -gt $ee ]; then eee=$ee; fi
	    durcum=$(( $durcum + $eee - $sss ))
    	    echo "$sss $eee" >> ${workdir}/segments/segments${usertag}_${p}_${nseg}.txt

	    # make job
	    if [ $durcum -ge $OMICRON_TRIGGERS_BASE ]; then
	    	echo "JOB omicron${usertag}_seg${nseg}_par${p} omicron.sub" >> ${workdir}/omicron${usertag}.dag
		echo "VARS omicron${usertag}_seg${nseg}_par${p} initialdir=\"${workdir}\" in_segments=\"./segments/segments${usertag}_${p}_${nseg}.txt\" in_parameters=\"./parameters/parameters${usertag}_${p}.txt\"" >> ${workdir}/omicron${usertag}.dag
		echo "   JOB omicron${usertag}_seg${nseg}_par${p}:"
		echo "       ---> parameters = ${workdir}/parameters/parameters${usertag}_${p}.txt"
		echo "       ---> segments = ${workdir}/segments/segments${usertag}_${p}_${nseg}.txt"
		echo "       ---> livetime = "`segsum ${workdir}/segments/segments${usertag}_${p}_${nseg}.txt`" sec"
	    
		# increment
		nseg=$(( $nseg + 1 ))
		durcum=0
	    fi
	    
	    sss=$(( $eee - $overlapduration ))
	done

	# one last job
	if [ $durcum -gt 0 ]; then
	    echo "JOB omicron${usertag}_seg${nseg}_par${p} omicron.sub" >> ${workdir}/omicron${usertag}.dag
	    echo "VARS omicron${usertag}_seg${nseg}_par${p} initialdir=\"${workdir}\" in_segments=\"./segments/segments${usertag}_${p}_${nseg}.txt\" in_parameters=\"./parameters/parameters${usertag}_${p}.txt\"" >> ${workdir}/omicron${usertag}.dag
	    echo "   JOB omicron${usertag}_seg${nseg}_par${p}:"
	    echo "       ---> parameters = ${workdir}/parameters/parameters${usertag}_${p}.txt"
	    echo "       ---> segments = ${workdir}/segments/segments${usertag}_${p}_${nseg}.txt"
	    echo "       ---> livetime = "`segsum ${workdir}/segments/segments${usertag}_${p}_${nseg}.txt`" sec"
	    
	    # increment
	    nseg=$(( $nseg + 1 ))
	    durcum=0
	fi

    done < ${workdir}/segments.txt
    
    let "p+=1"
done

exit 0
