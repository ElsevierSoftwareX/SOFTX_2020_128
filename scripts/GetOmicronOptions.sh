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
    echo "A keyword following the channel name can be used to select a pre-defined parameter set"
    echo "listed below. By default, the STD parameter set is used."
    echo ""
    echo "KEYWORD  FMIN[Hz]  FMAX[Hz]  SNR_THRESHOLD        OTHER"
    echo "-----------------------------------------------------------"
    echo " STD         8      2048             8"
    echo " STD2        8      1024             8"
    echo " STD3        8       512             8"
    echo " LOW       0.1        64             8         block=512s"
    echo " LOW2      0.1        64             8         block=8192s"
    echo " LOW3      0.1        64             8         block=65536s"
    echo " HIGH     1024      8192             8"
    echo " FINE        8      4096             6         fine tiling"
    echo " GW         32      4096             5         fine tiling"
    echo ""
    echo "The '-o' flag should be used for online analyses:"
    echo "Chunks are forced to 16s and overlaps to 2s (except for LOW)"
    echo ""
    echo "OPTIONS:"
    echo "  -c  [CHANNEL_FILE]    Channel file. Required option"
    echo "  -o                    Flag for online analyses"
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
    if [ $1 = "LOW" ]; then
	# online and offline configs are the same
	sampling=128
	freqmin=0.1
	freqmax=64
	chunk=512
	block=512
	overlap=112
	mmmax=0.35
	snrthr=8
	ratemax=300
    elif [ $1 = "LOW2" ]; then
	# online and offline configs are the same
	sampling=128
	freqmin=0.1
	freqmax=64
	chunk=8192
	block=8192
	overlap=192
	mmmax=0.35
	snrthr=8
	ratemax=250
    elif [ $1 = "LOW3" ]; then
	# online and offline configs are the same
	sampling=128
	freqmin=0.1
	freqmax=64
	chunk=65536
	block=65536
	overlap=15536
	mmmax=0.35
	snrthr=8
	ratemax=250
    elif [ $1 = "HIGH" ]; then
	sampling=16384
	freqmin=1024
	freqmax=8192
	if [ $4 -eq 0 ]; then # offline
	    chunk=506
	    block=256
	    overlap=6
	    ratemax=300
	else                  # online
	    chunk=16
	    block=16
	    overlap=2
	    ratemax=500
	fi
	mmmax=0.35
	snrthr=8	
    elif [ $1 = "GW" ]; then # this was optimized for ligo scripts
	sampling=8192
	freqmin=32
	freqmax=4096
	if [ $4 -eq 0 ]; then # offline
	    chunk=314
	    block=64
	    overlap=14
	else                  # online
	    chunk=16
	    block=16
	    overlap=2
	fi
	mmmax=0.2
	snrthr=5	
	ratemax=2000
    elif [ $1 = "FINE" ]; then # this was optimized for ligo scripts
	sampling=8192
	freqmin=8
	freqmax=4096
	if [ $4 -eq 0 ]; then # offline
	    chunk=314
	    block=64
	    overlap=14
	else                  # online
	    chunk=16
	    block=16
	    overlap=2
	fi
	mmmax=0.2
	snrthr=6	
	ratemax=1000
    elif [ $1 = "STD2" ]; then
	sampling=2048
	freqmin=8
	freqmax=1024
	if [ $4 -eq 0 ]; then # offline
	    chunk=506
	    block=256
	    overlap=6
	    ratemax=300
	else                  # online
	    chunk=16
	    block=16
	    overlap=2
	    ratemax=500
	fi
	mmmax=0.35
	snrthr=8
    elif [ $1 = "STD3" ]; then
        sampling=1024
        freqmin=8
        freqmax=512
        if [ $4 -eq 0 ]; then # offline
            chunk=506
            block=256
            overlap=6
            ratemax=300
        else                  # online
            chunk=16
            block=16
            overlap=2
            ratemax=500
        fi
        mmmax=0.35
        snrthr=8
    else
	sampling=4096
	freqmin=8
	freqmax=2048
	if [ $4 -eq 0 ]; then # offline
	    chunk=506
	    block=256
	    overlap=6
	    ratemax=300
	else                  # online
	    chunk=16
	    block=16
	    overlap=2
	    ratemax=500
	fi
	mmmax=0.35
	snrthr=8
    fi
    	
    # Static parameters
    echo "// ------------------------------------------------------------------"  > ./parameters_${1}_${2}.txt
    echo "// Omicron option file generated on `date`"                            >> ./parameters_${1}_${2}.txt
    echo "// Configuration type = $1"                                            >> ./parameters_${1}_${2}.txt
    echo "// ------------------------------------------------------------------" >> ./parameters_${1}_${2}.txt
    echo ""                                                                      >> ./parameters_${1}_${2}.txt
    echo "PARAMETER  QRANGE           3.3166  100.0"                             >> ./parameters_${1}_${2}.txt
    echo "OUTPUT     VERBOSITY        0"                                         >> ./parameters_${1}_${2}.txt
    echo "OUTPUT     WRITEPSD         0"                                         >> ./parameters_${1}_${2}.txt
    echo "OUTPUT     WRITETIMESERIES  0"                                         >> ./parameters_${1}_${2}.txt
    echo ""                                                                      >> ./parameters_${1}_${2}.txt
    
    # Tunable parameters
    echo "DATA       CHANNELS         $3"                                        >> ./parameters_${1}_${2}.txt
    echo "DATA       SAMPLEFREQUENCY  ${sampling}"                               >> ./parameters_${1}_${2}.txt
    echo "PARAMETER  CHUNKDURATION    ${chunk}"                                  >> ./parameters_${1}_${2}.txt
    echo "PARAMETER  BLOCKDURATION    ${block}"                                  >> ./parameters_${1}_${2}.txt
    echo "PARAMETER  OVERLAPDURATION  ${overlap}"                                >> ./parameters_${1}_${2}.txt
    echo "PARAMETER  FREQUENCYRANGE   ${freqmin}  ${freqmax}"                    >> ./parameters_${1}_${2}.txt
    echo "PARAMETER  MISMATCHMAX      ${mmmax}"                                  >> ./parameters_${1}_${2}.txt
    echo "TRIGGER    SNRTHRESHOLD     ${snrthr}"                                 >> ./parameters_${1}_${2}.txt
    echo "TRIGGER    RATEMAX          ${ratemax}"                                >> ./parameters_${1}_${2}.txt
    echo ""                                                                      >> ./parameters_${1}_${2}.txt

    # ffl
    if [ -e $5 ]; then 
	echo "DATA       FFL              $5"                                    >> ./parameters_${1}_${2}.txt
    fi

    # outdir
    echo "OUTPUT     DIRECTORY        $6"                                        >> ./parameters_${1}_${2}.txt
    
    # XML output
    if [ $7 -eq 1 ]; then 
	echo "OUTPUT     FORMAT           rootxml"                               >> ./parameters_${1}_${2}.txt
	echo "TRIGGER    CLUSTERING       TIME noroot"                           >> ./parameters_${1}_${2}.txt
    else
	echo "OUTPUT     FORMAT           root"                                  >> ./parameters_${1}_${2}.txt
    fi
    
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
while getopts ":oc:f:d:Xh" opt; do
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
rm -f ./channel.goo.*
awk '$2=="" || $2=="STD" {print $1}' $chanfile | sort | uniq > ./channel.goo.STD
awk '$2=="STD2" {print $1}' $chanfile | sort | uniq          > ./channel.goo.STD2
awk '$2=="STD3" {print $1}' $chanfile | sort | uniq          > ./channel.goo.STD3
awk '$2=="LOW" {print $1}' $chanfile | sort | uniq           > ./channel.goo.LOW
awk '$2=="LOW2" {print $1}' $chanfile | sort | uniq          > ./channel.goo.LOW2
awk '$2=="LOW3" {print $1}' $chanfile | sort | uniq          > ./channel.goo.LOW3
awk '$2=="HIGH" {print $1}' $chanfile | sort | uniq          > ./channel.goo.HIGH
awk '$2=="GW" {print $1}' $chanfile | sort | uniq            > ./channel.goo.GW
awk '$2=="FINE" {print $1}' $chanfile | sort | uniq          > ./channel.goo.FINE

rm -f ./parameters_STD_*.txt \
    ./parameters_STD2_*.txt \
    ./parameters_STD3_*.txt \
    ./parameters_FINE_*.txt \
    ./parameters_HIGH_*.txt \
    ./parameters_LOW_*.txt \
    ./parameters_LOW2_*.txt \
    ./parameters_LOW3_*.txt \
    ./parameters_GW_*.txt

for conf in STD STD2 STD3 LOW LOW2 LOW3 HIGH GW FINE; do
    if [ ! -s ./channel.goo.$conf ]; then continue; fi
    n=0
    p=0
    channel_list=""
    
    # maximum number of channels per option file
    if [ "$conf" = "STD" ];    then nchanmax=20
    elif [ "$conf" = "STD2" ]; then nchanmax=15
    elif [ "$conf" = "STD3" ]; then nchanmax=40
    elif [ "$conf" = "LOW" ];  then nchanmax=20
    elif [ "$conf" = "LOW2" ]; then nchanmax=10
    elif [ "$conf" = "LOW3" ]; then nchanmax=4
    elif [ "$conf" = "HIGH" ]; then nchanmax=4
    elif [ "$conf" = "FINE" ]; then nchanmax=5
    else nchanmax=5 #GW
    fi
    
    while read channel; do
	channel_list="$channel $channel_list"
	n=$(($n+1))
	# this is one option file
	if [ $n -eq $nchanmax ]; then
	    echo "${conf}_${p} = $channel_list"
	    printoption "${conf}" $p "$channel_list" $online $fflfile $outdir $outxml
	    channel_list=""
	    p=$(($p+1))
	    n=0
	fi
	
    done < ./channel.goo.${conf};
    echo "${conf}_${p} = $channel_list"
    if [ $n -gt 0 ]; then printoption "${conf}" $p "$channel_list" $online $fflfile $outdir $outxml; fi
done


# cleaning
rm -f ./channel.goo.*

exit 0
