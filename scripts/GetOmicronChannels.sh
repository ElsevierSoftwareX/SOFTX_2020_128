#!/bin/bash
#
# GetOmicronChannels.sh
#
# print available channels
#
# Author: Florent Robinet
# florent.robinet@lal.in2p3.fr

printhelp(){
    echo ""
    echo "Usage:"
    echo "GetOmicronChannels -v"
    echo " |__ prints the list of available channels available for the current run"
    echo ""
    echo "GetOmicronChannels -v -g [GPS]"
    echo " |__ prints the list of available channels available for the given GPS"
    echo ""
    echo "GetOmicronChannels -r [RUN_NAME]"
    echo " |__ prints the list of available channels available for the given run name"
    echo ""
    echo "When sourcing this script, the channel list is stored in the bash variable \$OMICRON_CHANNELS"
    echo ""
    echo "OPTIONS:"
    echo "  -g  [GPS]             for a given GPS time"
    echo "  -r  [RUN_NAME]        for a given run"
    echo "  -v                    be verbose"
    echo "  -h                    print this help"
    echo ""
    echo "Author: Florent Robinet (LAL - Orsay): robinet@lal.in2p3.fr"
    echo ""
} 

##### Check the environment
if [ -z "$OMICRON_TRIGGERS" ]; then
    echo "`basename $0`: The Omicron environment is not set"
    exit 1
fi

# default: take the latest run
goc_gps=`tconvert now`
RUN="NONE"
goc_verbose=0
OMICRON_CHANNELS=""
export OMICRON_CHANNELS

##### read options
while getopts ":r:g:vh" opt; do
    case $opt in
        r)
            RUN="$OPTARG"
            ;;
        g)
	    goc_gps=`echo $OPTARG | awk '{print int($1)}'`
            ;;
        v)
	    goc_verbose=1
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
goc_here=`pwd`

##### get run name if not requested
if [ "$RUN" = "NONE" ]; then . ${GWOLLUM_SCRIPTS}/getrun.sh -g $goc_gps; fi

##### scan available channels
if [ $goc_verbose -eq 1 ]; then echo "Available channels for run ${RUN}:"; fi
if echo $RUN_NAMES | grep -q -w $RUN; then # this run exists
    
    if [ -d ${OMICRON_TRIGGERS}/${RUN} ]; then # the trigger directory exists
	
	cd ${OMICRON_TRIGGERS}/${RUN}/
	OMICRON_CHANNELS=""
	for goc_chan in *; do
	    if [ ! -d $goc_chan ]; then continue; fi
	    if [ $goc_verbose -eq 1 ]; then echo ${goc_chan}; fi
	    OMICRON_CHANNELS="${goc_chan} ${OMICRON_CHANNELS}"
	done
	cd $goc_here

    fi

fi

export OMICRON_CHANNELS

##### clean 
unset goc_gps
unset goc_chan
unset goc_verbose
unset goc_here

##### no exit
