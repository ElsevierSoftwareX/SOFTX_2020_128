#!/bin/sh
#
# GetOmicronOptions.sh
#
# Author: Florent Robinet
# florent.robinet@lal.in2p3.fr

printhelp(){
    echo ""
    echo "Usage:"
    echo "GetOmicronOptions -c [CHANNEL_FILE]"
    echo ""
    echo "This script generates Omicron parameter files."
    echo ""
    echo "  -c  [CHANNEL_FILE]  path to a text file listing all the channels to process"
    echo "                      Different settings are available and are specified by"    
    echo "                      a keyword given after the channel name:"
    echo "                         STD:  8    to 2048Hz, SNR>8 (default)"
    echo "                         STD2: 8    to 1024Hz, SNR>8"
    echo "                         LOW:  0.1  to 64Hz,   SNR>8"
    echo "                         HIGH: 1024 to 8192Hz, SNR>8"
    echo "                         FINE: 8    to 4096Hz, SNR>6, fine tiling"
    echo "                         GW:   32   to 4096Hz, SNR>5, fine tiling"
    echo ""
    echo "  -o                  flag for online analyses"
    echo "                      chunks are forced ro 16s and overlaps to 2s"
    echo "                      (except for LOW)"
    echo ""

    echo "  -d  [TRIG_OUTDIR]   trigger output directory"
    echo "  -f  [FFL_FILE]      path to a FFL/LCF file to consider"
    echo ""
    echo "  -X                  add XML output format with time clustering"
    echo ""
    echo "  -h                  prints this help message"
    echo ""
} 

printoption(){
    if [ $1 = "LOW" ]; then
	# online and offline configs are the same
	sampling=128
	freqmin=0.1
	freqmax=64
	chunk=2048
	block=2048
	overlap=48
	chunk=2048
	block=2048
	overlap=48
	mmmax=0.35
	snrthr=8
    elif [ $1 = "HIGH" ]; then
	sampling=16384
	freqmin=1024
	freqmax=8192
	if [ $4 -eq 0 ]; then # offline
	    chunk=506
	    block=256
	    overlap=6
	else                  # online
	    chunk=16
	    block=16
	    overlap=2
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
    elif [ $1 = "STD2" ]; then
	sampling=2048
	freqmin=8
	freqmax=1024
	if [ $4 -eq 0 ]; then # offline
	    chunk=506
	    block=256
	    overlap=6
	else                  # online
	    chunk=16
	    block=16
	    overlap=2
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
	else                  # online
	    chunk=16
	    block=16
	    overlap=2
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
    echo "TRIGGER    RATEMAX          1000"                                      >> ./parameters_${1}_${2}.txt
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
outdir="."
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
		echo "GetOmicronOptions: the channel file $chanfile cannot be found"
		exit 2
	    fi
	    ;;
	f)
	    fflfile="$OPTARG"
	    if [ ! -e $fflfile ]; then
		echo "GetOmicronOptions: the FFL file $fflfile cannot be found"
		exit 2
	    fi
	    ;;
	d)
	    outdir="$OPTARG"
	    if [ ! -d $outdir ]; then
		echo "GetOmicronOptions: the output directory $outdir cannot be found"
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
OPTIND=0


##### check option file
if [ ! -e $chanfile ]; then
    echo "GetOmicronOptions: you must provide a channel file"
    echo "type  'GetOmicronOptions.sh -h'  for help"
    exit 3
fi

##### ffl file is required for offline analysis
if [ $online -eq 0 ]; then # offline
    if [ ! -e $fflfile ]; then
	echo "GetOmicronOptions: a FFL file is required"
	exit 2
    fi
fi

##### user channel lists
awk '$2=="" || $2=="STD" {print $1}' $chanfile | sort | uniq > ./channel.getomicronoptions.std
awk '$2=="STD2" {print $1}' $chanfile | sort | uniq          > ./channel.getomicronoptions.std2
awk '$2=="LOW" {print $1}' $chanfile | sort | uniq           > ./channel.getomicronoptions.low
awk '$2=="HIGH" {print $1}' $chanfile | sort | uniq          > ./channel.getomicronoptions.high
awk '$2=="GW" {print $1}' $chanfile | sort | uniq            > ./channel.getomicronoptions.gw
awk '$2=="FINE" {print $1}' $chanfile | sort | uniq          > ./channel.getomicronoptions.fine

rm -f ./parameters_STD_*.txt ./parameters_STD2_*.txt ./parameters_FINE_*.txt ./parameters_HIGH_*.txt ./parameters_LOW_*.txt./parameters_GW_*.txt

############## STD SETTING
if [ -s ./channel.getomicronoptions.std ]; then
    n=0
    p=0
    channel_list=""
    nchanmax=6 # maximum number of channels per option file
    while read channel; do
	channel_list="$channel $channel_list"
	n=$(($n+1))
	# this is one option file
	if [ $n -eq $nchanmax ]; then
	    echo "STD_${p} = $channel_list"
	    printoption "STD" $p "$channel_list" $online $fflfile $outdir $outxml
	    channel_list=""
	    p=$(($p+1))
	    n=0
	fi
	
    done < ./channel.getomicronoptions.std;
    echo "STD_${p} = $channel_list"
    if [ $n -gt 0 ]; then printoption "STD" $p "$channel_list" $online $fflfile $outdir $outxml; fi
fi

############## STD2 SETTING
if [ -s ./channel.getomicronoptions.std2 ]; then
    n=0
    p=0
    channel_list=""
    nchanmax=15 # maximum number of channels per option file
    while read channel; do
	channel_list="$channel $channel_list"
	n=$(($n+1))
	# this is one option file
	if [ $n -eq $nchanmax ]; then
	    echo "STD2_${p} = $channel_list"
	    printoption "STD2" $p "$channel_list" $online $fflfile $outdir $outxml
	    channel_list=""
	    p=$(($p+1))
	    n=0
	fi
	
    done < ./channel.getomicronoptions.std2;
    echo "STD2_${p} = $channel_list"
    printoption "STD2" $p"$channel_list" $online $fflfile $outdir $outxml
fi

############## LOW SETTING
if [ -s ./channel.getomicronoptions.low ]; then
    n=0
    p=0
    channel_list=""
    nchanmax=20 # maximum number of channels per option file
    while read channel; do
	channel_list="$channel $channel_list"
	n=$(($n+1))
	# this is one option file
	if [ $n -eq $nchanmax ]; then
	    echo "LOW_${p} = $channel_list"
	    printoption "LOW" $p "$channel_list" $online $fflfile $outdir $outxml
	    channel_list=""
	    p=$(($p+1))
	    n=0
	fi
	
    done < ./channel.getomicronoptions.low;
    echo "LOW_${p} = $channel_list"
    if [ $n -gt 0 ]; then printoption "LOW" $p "$channel_list" $online $fflfile $outdir $outxml; fi
fi
############## HIGH SETTING
if [ -s ./channel.getomicronoptions.high ]; then
    n=0
    p=0
    channel_list=""
    nchanmax=4 # maximum number of channels per option file
    while read channel; do
	channel_list="$channel $channel_list"
	n=$(($n+1))
	# this is one option file
	if [ $n -eq $nchanmax ]; then
	    echo "HIGH_${p} = $channel_list"
	    printoption "HIGH" $p "$channel_list" $online $fflfile $outdir $outxml
	    channel_list=""
	    p=$(($p+1))
	    n=0
	fi
	
    done < ./channel.getomicronoptions.high;
    echo "HIGH_${p} = $channel_list"
    if [ $n -gt 0 ]; then printoption "HIGH" $p "$channel_list" $online $fflfile $outdir $outxml; fi
fi
############## FINE SETTING
if [ -s ./channel.getomicronoptions.fine ]; then
    n=0
    p=0
    channel_list=""
    nchanmax=5 # maximum number of channels per option file
    while read channel; do
	channel_list="$channel $channel_list"
	n=$(($n+1))
	# this is one option file
	if [ $n -eq $nchanmax ]; then
	    echo "FINE_${p} = $channel_list"
	    printoption "FINE" $p "$channel_list" $online $fflfile $outdir $outxml
	    channel_list=""
	    p=$(($p+1))
	    n=0
	fi
	
    done < ./channel.getomicronoptions.fine;
    echo "FINE_${p} = $channel_list"
    if [ $n -gt 0 ]; then printoption "FINE" $p "$channel_list" $online $fflfile $outdir $outxml; fi
fi
############## GW SETTING
if [ -s ./channel.getomicronoptions.gw ]; then
    n=0
    p=0
    channel_list=""
    nchanmax=5 # maximum number of channels per option file
    while read channel; do
	channel_list="$channel $channel_list"
	n=$(($n+1))
	# this is one option file
	if [ $n -eq $nchanmax ]; then
	    echo "GW_${p} = $channel_list"
	    printoption "GW" $p "$channel_list" $online $fflfile $outdir $outxml
	    channel_list=""
	    p=$(($p+1))
	    n=0
	fi
	
    done < ./channel.getomicronoptions.gw;
    echo "GW_${p} = $channel_list"
    if [ $n -gt 0 ]; then printoption "GW" $p "$channel_list" $online $fflfile $outdir $outxml; fi
fi

# cleaning
rm -f ./channel.getomicronoptions.*

exit 0
