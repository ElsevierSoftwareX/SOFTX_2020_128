#!/bin/bash
#
# GetOmicronScan.sh
#
# generate an omicron scan
#
# Author: Florent Robinet
# florent.robinet@lal.in2p3.fr

printhelp(){
    echo ""
    echo "Usage:"
    echo "GetOmicronScan -g [GPS_CENTER]"
    echo " |__ scans 'standard' Omicron triggers around a given GPS time."
    echo ""
    echo ""
    echo "OPTIONS:"
    echo "  -g  [GPS_CENTER]      Central GPS time. required"
    echo "                        The GPS should be included in the input data"
    echo ""
    echo "SOURCE OPTIONS:"
    echo "  -i  [DATA_SOURCE]     Input data to use for the scan. Several inputs are possible:"
    echo "                        = \"TRIGGERS\"  : standard triggers are used as input"
    echo "                        = \"FRAMES\"    : standard raw frames are used as input"
    echo "                        = a file      : path to a ffl/lcf file"
    echo "                        = a directory : path a directory which contain sub-directories"
    echo "                                        where trigger files are located. The sub-directory"
    echo "                                        names should be the channel names."
    echo ""
    echo "CHANNEL OPTIONS:"
    echo "  -c  [CHANNEL_SEL]     A channel selection is used."
    echo "                        This option can be used in different ways:"
    echo "                        = \"ALL\"       : every possible channels will be scanned given"
    echo "                                        the source."
    echo "                        = \"STD\"       : a pre-defined standard selection will be used."
    echo "                        = a file      : this file should contain a single column with"
    echo "                                        the list of channels to scan"
    echo "                        = \"list\"      : where list is a list of channel to scan"
    echo "                                        separated by spaces"

    echo "  -m  [MAIN_CHANNEL]    Main channel: always plotted and plotted first"

    echo ""
    echo "CONTROL OPTIONS:"
    echo "  -d  [OUTDIR]          Path to output directory"
    echo "                        Default = current directory"
    echo "  -x  [SNRMIN]          Print a channel scan only if SNR > [SNRMIN]"
    echo "  -w  [WINDOWS]         List of time windows for the plots (between \"\")"
    
   
    echo ""
    echo "  -h                    print this help"
    echo ""
    echo "Author: Florent Robinet (LAL - Orsay): robinet@lal.in2p3.fr"
    echo ""
} 

##### Check the environment
if [ -z "$OMICRONROOT" ]; then
    echo "`basename $0`: The Omicron environment is not set"
    exit 1
fi

##### default options
tcenter=0
input="triggers"
chansel="std"
mainchan="none"
outdir=`pwd`
snrmin=8
windows="2 8 32"

here=`pwd`
triggerdir="none"
fflfile="none"
chanfile="none"

##### read options
while getopts ":g:i:c:m:d:x:w:h" opt; do
    case $opt in
	g)
	    tcenter=`echo "$OPTARG" | awk '{printf "%.3f", $1}'`
	    ;;
	i)
	    input="$OPTARG"
	    ;;
	c)
	    chansel="$OPTARG"
	    ;;
	m)
	    mainchan="$OPTARG"
	    ;;
	d)
	    outdir="$OPTARG"
	    ;;
	x)
	    snrmin=`echo $OPTARG | awk '{print int($1)}'`
	    ;;
	w)
	    windows="$OPTARG"
	    ;;
	h)
	    printhelp
	    exit 0
	    ;;
	\?)
	    echo "Invalid option: -$OPTARG"
	    echo "type  'GetOmicronScan -h'  for help"
	    exit 1
	    ;;
    esac
done
OPTIND=1

##### timing
tcenter_int=`echo $tcenter | awk '{print int($1)}'`
dateUTC=`tconvert -f "%A the %dth, %B %Y %H:%M:%S" ${tcenter_int}`
. ${GWOLLUM_SCRIPTS}/getrun.sh -g $tcenter_int;
tmin=$(( $tcenter_int - 1000 )) # to adjust
tmax=$(( $tcenter_int + 1000 )) # to adjust

##### check timing
if [ $tcenter_int -lt 700000000 ]; then
    echo "`basename $0`: '$tcenter' is not a reasonable central time"
    echo "type  'GetOmicronScan -h'  for help"
    exit 1
fi

##### check outdir
if [ ! -d $outdir ] ; then
    echo "`basename $0`: the output directory $outdir cannot be found"
    echo "type  'GetOmicronScan -h'  for help"
    exit 1
fi

##### check windows
nwin=0
winmax=0
for w in $windows; do
    if [ $w -lt 1 ]; then
	echo "`basename $0`: The time window must be an integer"
	echo "type  'GetOmicronScan -h'  for help"
	exit 1
    fi
    nwin=$(( nwin + 1 ))
    if [ $w -gt $winmax ]; then winmax=$w; fi
done
if [ $nwin -eq 0 ]; then
    echo "`basename $0`: There must be at least one time window"
    echo "type  'GetOmicronScan -h'  for help"
    exit 1
fi

##### standard triggers input
if [ "$input" = "TRIGGERS" ]; then

    if [ -z "$OMICRON_TRIGGERS" ]; then
	echo "`basename $0`: The Omicron trigger environment is not set"
	exit 1
    fi
    if [ ! -d  "${OMICRON_TRIGGERS}/${RUN}/${mainchan}" ]; then
	echo "`basename $0`: The main channel $mainchan does not exist in standard triggers"
	exit 1
    fi
    archive_stop=`printlivetime.exe "${OMICRON_TRIGGERS}/${RUN}/${mainchan}/*.root" | tail -n 1 | awk '{print int($2)}'`
    triggerdir="STD"
    if [ $tcenter_int -lt $archive_stop ]; then triggerdir=${OMICRON_TRIGGERS}/${RUN}
    else triggerdir=$OMICRON_ONLINE_TRIGGERS; fi

##### user triggers input
elif [ -d $input ]; then
    triggerdir=$input
    input="TRIGGERS"
    
##### std frames
elif [ "$input" = "FRAMES" ]; then
    fflfile="/virgoData/ffl/raw.ffl" # HARDCODED

##### user frames
elif [ -e $input ]; then
    fflfile=$input
    input="FRAMES"

##### WTF?
else
    echo "`basename $0`: the input data source option is not understood"
    echo "type  'GetOmicronScan -h'  for help"
    exit 1
fi

##### standard triggers selection
if [ "$chansel" = "STD" ]; then
    chanfile=${OMICRON_PARAMETERS}/scan.${RUN}.txt
    if [ ! -e $chanfile ]; then
	echo "`basename $0`: there is no standard channel selection for run=$RUN"
	echo "type  'GetOmicronScan -h'  for help"
	exit 1
    fi

##### all channels
elif [ "$chansel" = "ALL" ]; then
    echo "scan all"

##### user selection
elif [ -e $chansel ]; then
    chanfile=$chansel

##### list of channels
else
    chanfile=${TMP}/scan.${tcenter}.${RANDOM}
    rm -f $chanfile
    for chan in $chansel; do
	echo "$chan" >> $chanfile
    done
fi

#####################################################################################
##################                  PLOT TRIGGERS                  ##################
#####################################################################################
if [ "$input" = "TRIGGERS" ]; then

    # make parameter file
    parameterfile=${TMP}/parameters.${tcenter}.${RANDOM}
    echo "// Omiscan parameter file generated on `date`" > $parameterfile
    echo "DATA      TRIGGERS     $triggerdir" >> $parameterfile
    echo "PARAMETER WINDOWS      $windows"    >> $parameterfile
    echo "TRIGGER   SNRTHRESHOLD $snrmin"     >> $parameterfile
    echo "OUTPUT    FORMAT       gif"         >> $parameterfile
    echo "OUTPUT    VERBOSITY    2"           >> $parameterfile
    echo "OUTPUT    DIRECTORY    ${outdir}"   >> $parameterfile
    if [ -s $chanfile ]; then
	awk '{print "DATA      CHANNELS    ",$1}' $chanfile >> $parameterfile
    fi

    omiscan.exe $tcenter $parameterfile
    exit 0

fi


if [ "$input" = "FRAMES" ]; then

    # make parameter file
    parameterfile=${TMP}/parameters.${tcenter}.${RANDOM}
    echo "// Omiscan parameter file generated on `date`" > $parameterfile
    echo "DATA      FFL          $fflfile"    >> $parameterfile
    echo "PARAMETER WINDOWS      $windows"    >> $parameterfile
    echo "TRIGGER   SNRTHRESHOLD $snrmin"     >> $parameterfile
    echo "PARAMETER MISMATCHMAX  0.2"         >> $parameterfile
    echo "PARAMETER QRANGE       3.3166 141"  >> $parameterfile
    echo "OUTPUT    FORMAT       gif"         >> $parameterfile
    echo "OUTPUT    VERBOSITY    2"           >> $parameterfile
    echo "OUTPUT    DIRECTORY    ${outdir}"   >> $parameterfile
    if [ -s $chanfile ]; then
	scantype=
	awk '{print "DATA      CHANNELS    ",$1}' $chanfile >> $parameterfile
    fi

    omiscan.exe $tcenter $parameterfile
    exit 0

fi





exit 0
