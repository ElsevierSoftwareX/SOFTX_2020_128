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
    echo "This script produces Omicron parameter files."
    echo ""
    echo "  -o  [OPTION FILE]   path to an Omicron parameter file used as a template"
    echo "                      Default: use a standardized parameter file"
    echo "  -c  [CHANNEL_FILE]  path to a text file listing all the channels to process"
    echo "                      Different settings are available and are specified by"    
    echo "                      a keyword given after the channel name:"
    echo "                         STD:  8    to 1024Hz, SNR>8 (default)"
    echo "                         STD2: 8    to 2048Hz, SNR>8 (default)"
    echo "                         LOW:  0.1  to 64Hz,   SNR>8"
    echo "                         HIGH: 1024 to 8192Hz, SNR>8"
    echo "                         GW:   8    to 4096Hz, SNR>5, fine tiling"
    echo ""
    echo ""

    echo "  -d  [TRIG_OUTDIR]   trigger output directory (inactive if option file given)"
    echo "  -f  [FFL_FILE]      path to a FFL/LCF file to consider (inactive if option file given)"
    echo ""
    echo "  -h                  prints this help message"
    echo ""
} 

printoption(){
    if [ $1 = "LOW" ]; then
	sampling=128
	freqmin=0.1
	freqmax=64
	chunk=2048
	block=2048
	overlap=48
	mmmax=0.35
	snrthr=8
    elif [ $1 = "HIGH" ]; then
	sampling=16384
	freqmin=1024
	freqmax=8192
	chunk=506
	block=256
	overlap=6
	mmmax=0.35
	snrthr=8	
    elif [ $1 = "GW" ]; then
	sampling=8192
	freqmin=32
	freqmax=4096
	chunk=506
	block=256
	overlap=6
	mmmax=0.2
	snrthr=5	
    elif [ $1 = "STD2" ]; then
	sampling=2048
	freqmin=8
	freqmax=1024
	chunk=506
	block=256
	overlap=6
	mmmax=0.35
	snrthr=8
    else
	sampling=4096
	freqmin=8
	freqmax=2048
	chunk=506
	block=256
	overlap=6
	mmmax=0.35
	snrthr=8
    fi
    	
    sed -e '/DATA[ \t]*CHANNELS/d' \
	-e '/DATA[ \t]*SAMPLEFREQUENCY/d' \
	-e '/PARAMETER[ \t]*CHUNKDURATION/d' \
	-e '/PARAMETER[ \t]*BLOCKDURATION/d' \
	-e '/PARAMETER[ \t]*OVERLAPDURATION/d' \
	-e '/PARAMETER[ \t]*FREQUENCYRANGE/d' \
	-e '/PARAMETER[ \t]*MISMATCHMAX/d' \
	-e '/TRIGGER[ \t]*SNRTHRESHOLD/d' \
	$2 > ./parameters_${1}_${3}.txt
    
    echo "" >> ./parameters_${1}_${3}.txt
    echo "DATA  CHANNELS  $4" >> ./parameters_${1}_${3}.txt
    echo "DATA  SAMPLEFREQUENCY      ${sampling}" >> ./parameters_${1}_${3}.txt
    echo "PARAMETER CHUNKDURATION    ${chunk}" >> ./parameters_${1}_${3}.txt
    echo "PARAMETER BLOCKDURATION    ${block}" >> ./parameters_${1}_${3}.txt
    echo "PARAMETER OVERLAPDURATION  ${overlap}" >> ./parameters_${1}_${3}.txt
    echo "PARAMETER FREQUENCYRANGE   ${freqmin}  ${freqmax}" >> ./parameters_${1}_${3}.txt
    echo "PARAMETER MISMATCHMAX      ${mmmax}" >> ./parameters_${1}_${3}.txt
    echo "TRIGGER SNRTHRESHOLD       ${snrthr}" >> ./parameters_${1}_${3}.txt
    
} 

##### Check the Omicron environment
if [ -z "$OMICRONROOT" ]; then
    echo "Error: The Omicron environment is not set"
    exit 1
fi

##### default options
optfile="./no.chance.this.file.exists"
chanfile="./no.chance.this.file.exists"
framefile="./no.chance.this.file.exists"
fflfile="./no.chance.this.file.exists"
outdir="./no.chance.this.dir.exists"

##### read options
while getopts ":o:c:f:d:h" opt; do
    case $opt in
	o)
	    optfile="$OPTARG"
	    if [ ! -e $optfile ]; then
		echo "GetOmicronOptions: the option file $optfile cannot be found"
		exit 2
	    fi
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


##### check option file
if [ ! -e $chanfile ]; then
    echo "GetOmicronOptions: you must provide a channel file"
    echo "type  'GetOmicronOptions.sh -h'  for help"
    exit 3
fi

##### check option file
if [ ! -e $optfile ]; then
    echo ""
    echo "******************************************"
    echo "**  No template option file was given!  **"
    echo "**  A default option file will be used  **"
    echo "******************************************"
    echo ""
    if [ ! -d $outdir ]; then
	echo "Please enter a path to an output directory: "
	read outdir
	outdir="${outdir}/triggers"
    fi
    if [ ! -e $fflfile ]; then
	echo "Please enter a path to an FFL/LCF file: "
	read fflfile
    fi
    sed -e "s|\[FFL]|$fflfile|g" \
	-e "s|\[OUTDIR]|$outdir|g" \
	$OMICRON_SCRIPTS/parameters.template > ./parameters.template
    optfile="./parameters.template"
fi

##### user channel lists
awk '$2=="" || $2=="STD" {print $1}' $chanfile | sort | uniq > ./channel.getomicronoptions.std
awk '$2=="STD2" {print $1}' $chanfile | sort | uniq           > ./channel.getomicronoptions.std2
awk '$2=="LOW" {print $1}' $chanfile | sort | uniq           > ./channel.getomicronoptions.low
awk '$2=="HIGH" {print $1}' $chanfile | sort | uniq          > ./channel.getomicronoptions.high
awk '$2=="GW" {print $1}' $chanfile | sort | uniq            > ./channel.getomicronoptions.gw

##### break into freq
rm -f ./parameters_STD_*.txt ./parameters_STD2_*.txt ./parameters_HIGH_*.txt ./parameters_LOW_*.txt./parameters_GW_*.txt

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
	    printoption "STD" $optfile $p "$channel_list"
	    channel_list=""
	    p=$(($p+1))
	    n=0
	fi
	
    done < ./channel.getomicronoptions.std;
    echo "STD_${p} = $channel_list"
    printoption "STD" $optfile $p "$channel_list"
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
	    printoption "STD2" $optfile $p "$channel_list"
	    channel_list=""
	    p=$(($p+1))
	    n=0
	fi
	
    done < ./channel.getomicronoptions.std2;
    echo "STD2_${p} = $channel_list"
    printoption "STD2" $optfile $p "$channel_list"
fi

############## LOW SETTING
if [ -s ./channel.getomicronoptions.low ]; then
    n=0
    p=0
    channel_list=""
    nchanmax=50 # maximum number of channels per option file
    while read channel; do
	channel_list="$channel $channel_list"
	n=$(($n+1))
	# this is one option file
	if [ $n -eq $nchanmax ]; then
	    echo "LOW_${p} = $channel_list"
	    printoption "LOW" $optfile $p "$channel_list"
	    channel_list=""
	    p=$(($p+1))
	    n=0
	fi
	
    done < ./channel.getomicronoptions.low;
    echo "LOW_${p} = $channel_list"
    printoption "LOW" $optfile $p "$channel_list"
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
	    printoption "HIGH" $optfile $p "$channel_list"
	    channel_list=""
	    p=$(($p+1))
	    n=0
	fi
	
    done < ./channel.getomicronoptions.high;
    echo "HIGH_${p} = $channel_list"
    printoption "HIGH" $optfile $p "$channel_list"
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
	    printoption "GW" $optfile $p "$channel_list"
	    channel_list=""
	    p=$(($p+1))
	    n=0
	fi
	
    done < ./channel.getomicronoptions.gw;
    echo "GW_${p} = $channel_list"
    printoption "GW" $optfile $p "$channel_list"
fi

# cleaning
rm -f ./parameters.template channel.getomicronoptions.*

exit 0
