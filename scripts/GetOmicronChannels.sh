#!/bin/bash
#
# GetOmicronChannels.sh
#
# print available channels
#
# Author: Florent Robinet
# florent.robinet@lal.in2p3.fr

here=`pwd`

printhelp(){
    echo ""
    echo "Usage:"
    echo "GetOmicronChannels"
    echo ""
    echo "Example: GetOmicronChannels -r VSR2"
    echo ""
    echo "  -r  [RUN_NAME]  select run with name [RUN_NAME]"
    echo "                  Default = last run"
    echo "  -h              print this help"
    echo ""
} 

##### Check the environment
if [[ -z "$OMICRON_TRIGGERS" ]]; then
    echo "Error: The Omicron environment is not set"
    exit 1
fi

# default: take the latest run
for rr in $RUN_NAMES; do
    ru=$rr
done

##### read options
while getopts ":r:h" opt; do
    case $opt in
        r)
            ru="$OPTARG"
            ;;
	h)
	    printhelp
	    exit 0
	    ;;
	\?)
	    echo "Invalid option: -$OPTARG"
	    echo "type  'GetOmicronChannels -h'  for help"
	    exit 1
	    ;;
    esac
done
OPTIND=0

if [ ! -d $OMICRON_TRIGGERS ]; then
    echo "`basename $0`: the trigger directory does not exist"
    exit 1
fi

if [ ! -d ${OMICRON_TRIGGERS}/${ru} ]; then
    echo "`basename $0`: the triggers for the run $ru do not exist"
    exit 2
fi

echo "Available channels for run ${ru}:"
echo ""

## loop over available channels
cd ${OMICRON_TRIGGERS}/${ru}/
OMICRON_CHANNELS=""
for chan in *; do
    echo ${chan}
    OMICRON_CHANNELS="${chan} ${OMICRON_CHANNELS}"
done
cd $here

export OMICRON_CHANNELS

echo ""
echo "This list of channels is stored in the OMICRON_CHANNELS variable"
echo "You can use it in your own shell scripts"
