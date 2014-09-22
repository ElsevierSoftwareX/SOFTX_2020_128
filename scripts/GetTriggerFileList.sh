#!/bin/bash
#
# GetTriggerFileList.sh
#
# returns optimized list of trigger files
#
# Author: Florent Robinet
# florent.robinet@lal.in2p3.fr

printhelp(){
    echo ""
    echo "Usage:"
    echo "GetTriggerFileList -v -c [CHANNEL_NAME] -s [GPS_START] -e [GPS_END]"
    echo " |__ prints the list of files containing Omicron triggers"
    echo "     of a given channel and between 2 GPS times"
    echo ""
    echo "When sourcing this script, the file list is stored in the bash variable \$OMICRON_FILELIST"
    echo ""
    echo "When online trigger files are concerned, they are copied in a temporary directory."
    echo ""
    echo "[GPS_START] and [GPS_END] should belong to the same run"
    echo ""
    echo "OPTIONS:"
    echo "  -c  [CHANNEL_NAME]    channel name. Default = V1:h_4096Hz"
    echo "  -s  [GPS_START]       starting GPS time"
    echo "  -e  [GPS_END]         stopping GPS time"
    echo "  -v                    be verbose"
    echo "  -h                    print this help"
    echo ""
    echo "Author: Florent Robinet (LAL - Orsay): robinet@lal.in2p3.fr"
    echo ""
} 

##### Check the Omicron environment
if [ -z "$OMICRONROOT" ]; then
    echo "`basename $0`: The Omicron environment is not set"
    exit 1
fi

##### default options
gtf_channel="V1:h_4096Hz"
gtf_tmin=0
gtf_tmax=0
gtf_verbose=0
OMICRON_FILELIST=""
export OMICRON_FILELIST

##### read options
while getopts ":c:s:e:vh" gtf_opt; do
    case $gtf_opt in
	c)
	    gtf_channel="$OPTARG"
	    ;;
        s)
	    gtf_tmin=`echo $OPTARG | awk '{print int($1)}'`
            ;;
        e)
	    gtf_tmax=`echo $OPTARG | awk '{print int($1)}'`
            ;;
        v)
	    gtf_verbose=1
            ;;
	h)
	    printhelp
	    exit 0
	    ;;
	\?)
	    echo "`basename $0`: Invalid option: -$OPTARG"
	    exit 1
	    ;;
    esac
done
OPTIND=1

#### clean tmp dir with old trigger directories
gtf_now=`tconvert now`
for gtf_dir in ${TMP}/triggers.*.*; do
    if [ -d $gtf_dir ]; then
	gtf_g=`echo $gtf_dir | awk -F. '{print $((NF -1))}'`
	gtf_tdiff=$(( $gtf_now - $gtf_g ))
	if [ $gtf_tdiff -gt 10000 ]; then rm -fr $gtf_dir > /dev/null 2>&1; fi
    fi
done

##### get run
if [ $gtf_tmin -eq "0" ]; then
    gtf_tmin=$gtf_tmax
fi
if [ $gtf_tmax -eq "0" ]; then
    gtf_tmax=$gtf_tmin
fi
. ${GWOLLUM_SCRIPTS}/getrun.sh -g $gtf_tmax;

##### optimize timing
gtf_tmin_base=$(( $gtf_tmin / $OMICRON_TRIGGERS_BASE - 1 ))
gtf_tmax_base=$(( $gtf_tmax / $OMICRON_TRIGGERS_BASE + 1 ))

#### tmp dir for online files
gtf_tmpdir=${TMP}/triggers.${gtf_now}.${RANDOM}
mkdir -p ${gtf_tmpdir}; rm -f ${gtf_tmpdir}/*.root > /dev/null 2>&1

#### scan trigger directories
gtf_first=1
if [ $gtf_verbose -eq 1 ]; then 
    echo ""
    echo "List of triggers files for $gtf_channel between $gtf_tmin and ${gtf_tmax}:"; 
    echo ""
fi

while [ $gtf_tmin_base -le $gtf_tmax_base ]; do

    # start with offline triggers
    for gtf_file in ${OMICRON_TRIGGERS}/${RUN}/${gtf_channel}/${gtf_channel}_${gtf_tmin_base}*.root; do 
	if [ -e $gtf_file ]; then 
	    gtf_s=`echo $gtf_file | awk -F_ '{print $((NF -1))}'`
	    if [ $gtf_s -ge $gtf_tmax ]; then break; fi
	    gtf_d=`echo $gtf_file | awk -F_ '{print $NF}' | awk 'BEGIN{FS=".root"} {print $1}'`
	    gtf_e=$(($gtf_s+$gtf_d))
	    if [ $gtf_e -lt $gtf_tmin ]; then continue; fi

	    # keep this file
	    OMICRON_FILELIST="$OMICRON_FILELIST $gtf_file"
	    if [ $gtf_verbose -eq 1 ]; then echo "$gtf_file"; fi
	fi 
    done

    # then online triggers
    for gtf_file in ${OMICRON_ONLINE_TRIGGERS}/${gtf_channel}/${gtf_channel}_${gtf_tmin_base}*.root; do 
	if [ -e $gtf_file ]; then 
	    gtf_s=`echo $gtf_file | awk -F_ '{print $((NF -1))}'`
	    if [ $gtf_s -ge $gtf_tmax ]; then break; fi
	    gtf_d=`echo $gtf_file | awk -F_ '{print $NF}' | awk 'BEGIN{FS=".root"} {print $1}'`
	    gtf_e=$(($gtf_s+$gtf_d))
	    if [ $gtf_e -lt $gtf_tmin ]; then continue; fi

	    # keep this file
	    cp $gtf_file ${gtf_tmpdir}/
	    if [ $gtf_first -eq 1 ]; then 
		OMICRON_FILELIST="$OMICRON_FILELIST ${gtf_tmpdir}/*.root"; 
		if [ $gtf_verbose -eq 1 ]; then echo "${gtf_tmpdir}/*.root"; fi
		gtf_first=0; 
	    fi
	fi 
    done

    let "gtf_tmin_base+=1"

    gtf_newtmin=$(( $gtf_tmin_base * $OMICRON_TRIGGERS_BASE ))
    if [ $gtf_newtmin -gt $gtf_tmax ]; then break; fi
done

export OMICRON_FILELIST

#### clean up
unset gtf_channel
unset gtf_tmin
unset gtf_tmax
unset gtf_verbose
unset gtf_now
unset gtf_dir
unset gtf_g
unset gtf_tdiff
unset gtf_tmin_base
unset gtf_tmax_base
unset gtf_tmpdir
unset gtf_first
unset gtf_file
unset gtf_s
unset gtf_d
unset gtf_e

#### no exit
