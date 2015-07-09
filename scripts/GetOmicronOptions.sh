#!/bin/bash
#
# GetOmicronOptions.sh
#
# Author: Florent Robinet
# florent.robinet@lal.in2p3.fr

printhelp(){
    echo ""
    echo "Usage:"
    echo "GetOmicronOptions -c [CHANNEL_FILE] -f [FFL_FILE]"
    echo " |__ generates Omicron parameter files for an offline search."
    echo ""
    echo "GetOmicronOptions -c [CHANNEL_FILE] -o"
    echo " |__ generates Omicron parameter files for an online search."
    echo ""
    echo "The input channel file should contain a list of channels to process, one channel per line."
    echo "A keyword following the channel name is used to select a pre-defined parameter set"
    echo "listed below."
    echo ""
    echo "KEYWORD  FMIN[Hz]  FMAX[Hz]  SNR_THRESHOLD        OTHER"
    echo "-----------------------------------------------------------"
    echo " STD1        8      2048             7"
    echo " STD2        8      1024             7"
    echo " STD3        8       512             7"
    echo " LOW1      0.1        64             7         block=512s"
    echo " LOW2      0.1        64             7         block=8192s"
    echo " HIGH     1024      8192             7"
    echo " FINE        8      4096             6         fine tiling"
    echo " GW         32      4096             5         fine tiling"
    echo ""
    echo "The '-o' flag should be used for online analyses:"
    echo "Chunks are forced to 16s and overlaps to 2s (except for LOW)"
    echo ""
    echo "OPTIONS:"
    echo "  -c  [CHANNEL_FILE]    Channel file. Required option"
    echo "  -o                    Flag for online analyses"
    echo "  -n  [N_CHANNELS]      number of channels per option file"
   echo ""
    echo "PARAMETER OPTIONS:"
    echo "  -d  [TRIG_OUTDIR]     Trigger output directory"
    echo "                        Default: current directory"
    echo "  -f  [FFL_FILE]        Path to a FFL/LCF file to consider"
    echo "                        This option is required without the '-o' flag"
    echo "  -X                    Add XML output format with time clustering"
    echo ""
    echo "  -h                    Print this help"
    echo ""
    echo "Author: Florent Robinet (LAL - Orsay): robinet@lal.in2p3.fr"
    echo ""
} 

printoption(){
    if [ $1 = "LOW1" ]; then
	# online and offline configs are the same
	sampling=128
	freqmin=0.1
	freqmax=64
	chunk=512
	block=512
	overlap=112
	mmmax=0.35
	snrthr=7
    elif [ $1 = "LOW2" ]; then
	# online and offline configs are the same
	sampling=128
	freqmin=0.1
	freqmax=64
	chunk=8192
	block=8192
	overlap=192
	mmmax=0.35
	snrthr=7
    elif [ $1 = "HIGH" ]; then
	sampling=16384
	freqmin=1024
	freqmax=8192
	if [ $4 -eq 0 ]; then # offline
	    chunk=312
	    block=64
	    overlap=2
	else                  # online
	    chunk=16
	    block=16
	    overlap=2
	fi
	mmmax=0.35
	snrthr=7
    elif [ $1 = "GW" ]; then
	sampling=8192
	freqmin=32
	freqmax=4096
	if [ $4 -eq 0 ]; then # offline
	    chunk=544
	    block=64
	    overlap=4
	else                  # online
	    chunk=16
	    block=16
	    overlap=2
	fi
	mmmax=0.2
	snrthr=5	
    elif [ $1 = "FINE" ]; then
	sampling=8192
	freqmin=8
	freqmax=4096
	if [ $4 -eq 0 ]; then # offline
	    chunk=544
	    block=64
	    overlap=4
	else                  # online
	    chunk=16
	    block=16
	    overlap=2
	fi
	mmmax=0.2
	snrthr=6	
    elif [ $1 = "STD2" ]; then
	sampling=2048
	freqmin=8
	freqmax=1024
	if [ $4 -eq 0 ]; then # offline
	    chunk=544
	    block=64
	    overlap=4
	else                  # online
	    chunk=16
	    block=16
	    overlap=2
	fi
	mmmax=0.35
	snrthr=7
    elif [ $1 = "STD3" ]; then
        sampling=1024
        freqmin=8
        freqmax=512
        if [ $4 -eq 0 ]; then # offline
            chunk=544
            block=64
            overlap=4
        else                  # online
            chunk=16
            block=16
            overlap=2
        fi
        mmmax=0.35
        snrthr=7
    else
	sampling=4096
	freqmin=8
	freqmax=2048
	if [ $4 -eq 0 ]; then # offline
	    chunk=544
	    block=64
	    overlap=4
	else                  # online
	    chunk=16
	    block=16
	    overlap=2
	fi
	mmmax=0.35
	snrthr=7
    fi
    	
    # Static parameters
    echo "// ------------------------------------------------------------------"  > ./parameters_${1}_${2}.txt
    echo "// Omicron (v2r1) option file generated on `date`"                     >> ./parameters_${1}_${2}.txt
    echo "// Configuration type = $1"                                            >> ./parameters_${1}_${2}.txt
    echo "// ------------------------------------------------------------------" >> ./parameters_${1}_${2}.txt
    echo ""                                                                      >> ./parameters_${1}_${2}.txt

    if [ $4 -eq 0 ]; then 
	echo "DATA       FFL              $5"                                    >> ./parameters_${1}_${2}.txt
    else
	echo "// This is an online search: contain the trigger rate"             >> ./parameters_${1}_${2}.txt
	echo "OUTPUT     TRIGGERRATEMAX   500"                                   >> ./parameters_${1}_${2}.txt
	echo ""                                                                  >> ./parameters_${1}_${2}.txt
    fi
    for chan in $3; do
	echo "DATA       CHANNELS         $chan"                                 >> ./parameters_${1}_${2}.txt
    done
    echo "DATA       SAMPLEFREQUENCY  ${sampling}"                               >> ./parameters_${1}_${2}.txt
    echo ""                                                                      >> ./parameters_${1}_${2}.txt
    echo "PARAMETER  CHUNKDURATION    ${chunk}"                                  >> ./parameters_${1}_${2}.txt
    echo "PARAMETER  SEGMENTDURATION  ${block}"                                  >> ./parameters_${1}_${2}.txt
    echo "PARAMETER  OVERLAPDURATION  ${overlap}"                                >> ./parameters_${1}_${2}.txt
    echo "PARAMETER  QRANGE           3.3166  100.0"                             >> ./parameters_${1}_${2}.txt
    echo "PARAMETER  FREQUENCYRANGE   ${freqmin}  ${freqmax}"                    >> ./parameters_${1}_${2}.txt
    echo "PARAMETER  MISMATCHMAX      ${mmmax}"                                  >> ./parameters_${1}_${2}.txt
    echo "PARAMETER  SNRTHRESHOLD     ${snrthr}"                                 >> ./parameters_${1}_${2}.txt
    echo ""                                                                      >> ./parameters_${1}_${2}.txt
    echo "OUTPUT     DIRECTORY        $6"                                        >> ./parameters_${1}_${2}.txt
    echo "OUTPUT     PRODUCTS         triggers"                                  >> ./parameters_${1}_${2}.txt
    echo "OUTPUT     VERBOSITY        0"                                         >> ./parameters_${1}_${2}.txt
    if [ $7 -eq 1 ]; then
	echo "OUTPUT     FORMAT           rootxml"                               >> ./parameters_${1}_${2}.txt
	echo ""                                                                  >> ./parameters_${1}_${2}.txt
	echo "// clustering is only applied to XML"                              >> ./parameters_${1}_${2}.txt
	echo "PARAMETER  CLUSTERING       TIME"                                  >> ./parameters_${1}_${2}.txt
    else
	echo "OUTPUT     FORMAT           root"                                  >> ./parameters_${1}_${2}.txt
    fi
    echo ""                                                                      >> ./parameters_${1}_${2}.txt
} 

##### Check the Omicron environment
if [ -z "$OMICRONROOT" ]; then
    echo "Error: The Omicron environment is not set"
    exit 1
fi

##### default options
online=0
chanfile="./no.chance.this.file.exists"
fflfile="./no.chance.this.file.exists"
outdir=`pwd`
outxml=0

##### read options
while getopts ":oc:f:d:n:Xh" opt; do
    case $opt in
	o)
	    online=1
	    ;;
	c)
	    chanfile="$OPTARG"
	    if [ ! -e $chanfile ]; then
		echo "`basename $0`: the channel file $chanfile cannot be found"
		exit 2
	    fi
	    ;;
	f)
	    fflfile="$OPTARG"
	    if [ ! -e $fflfile ]; then
		echo "`basename $0`: the FFL file $fflfile cannot be found"
		exit 2
	    fi
	    ;;
	d)
	    outdir="$OPTARG"
	    if [ ! -d $outdir ]; then
		echo "`basename $0`: the output directory $outdir cannot be found"
		exit 2
	    fi
	    ;;
	n)
	    nmax="$OPTARG"
	    ;;
	X)
	    outxml=1
	    ;;
	h)
	    printhelp
	    exit 0
	    ;;
	\?)
	    echo "Invalid option: -$OPTARG"
	    echo "type  'GetOmicronOptions.sh -h'  for help"
	    exit 1
	    ;;
    esac
done
OPTIND=1


##### check option file
if [ ! -e $chanfile ]; then
    echo "`basename $0`: you must provide a channel file"
    echo "type  'GetOmicronOptions.sh -h'  for help"
    exit 3
fi

##### ffl file is required for offline analysis
if [ $online -eq 0 ]; then # offline
    if [ ! -e $fflfile ]; then
	echo "`basename $0`: a FFL file is required"
	exit 2
    fi
fi

##### user channel lists
for conf in STD1 STD2 STD3 LOW1 LOW2 HIGH GW FINE; do

    # select channels for this config
    if [ "$conf" = "STD1" ]; then
	awk '$2=="STD1" {print $1}' $chanfile | sort | uniq > ./channel.goo.$conf
    else
	awk -v var="$conf" '$2==var {print $1}' $chanfile | sort | uniq > ./channel.goo.$conf
    fi

    # cleaning
    rm -f ./parameters_${conf}_*.txt

    # skip if no channels
    if [ ! -s ./channel.goo.$conf ]; then continue; fi
    echo ""
    echo "GetOmicronOptions: Make $conf option files..."
    echo "   --> "`cat ./channel.goo.$conf | wc -l`" channels"

    # maximum number of channels per option file
    nchanmax=`grep "${conf}=" $chanfile | cut -d'=' -f2`
    if [ "$nchanmax" = "" ]; then nchanmax=5; fi # default
    echo "   --> $nchanmax channels / option file"

    # loop over channels for this conf
    n=0
    p=0
    channel_list=""
    while read channel; do

	# add this channel
	channel_list="$channel $channel_list"
	n=$(($n+1))
	
	# option file is full
	if [ $n -eq $nchanmax ]; then
	    #echo "${conf}_${p} = $channel_list"
	    printoption "${conf}" $p "$channel_list" $online $fflfile $outdir $outxml
	    channel_list=""
	    p=$(($p+1))
	    n=0
	fi
	
    done < ./channel.goo.${conf};


    #echo "${conf}_${p} = $channel_list"
    # left over --> last option file
    if [ $n -gt 0 ]; then printoption "${conf}" $p "$channel_list" $online $fflfile $outdir $outxml; fi
    echo "   --> $(( $p + 1)) option files were produced"
done


# cleaning
rm -f ./channel.goo.*

exit 0
