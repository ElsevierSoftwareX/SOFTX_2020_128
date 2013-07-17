#!/bin/bash
#
# GetOmicronOptions.sh
#
# Author: Florent Robinet
# florent.robinet@lal.in2p3.fr

printhelp(){
    echo ""
    echo "Usage:"
    echo "GetOmicronOptions -o [OPTION FILE]"
    echo ""
} 

##### Check the Omicron environment
if [[ -z "$OMICRONROOT" ]]; then
    echo "Error: The Omicron environment is not set"
    exit 1
fi

##### default options
optfile="./no.chance.this.file.exists"
chanfile="./no.chance.this.file.exists"

##### read options
while getopts ":o:c:h" opt; do
    case $opt in
	o)
	    optfile="$OPTARG"
	    ;;
	c)
	    chanfile="$OPTARG"
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
if [ ! -e $optfile ]; then
    echo ""
    echo "******************************************"
    echo "**  No template option file was given!  **"
    echo "**  A default option file will be used  **"
    echo "******************************************"
    echo ""
    cp -f $OMICRON_SCRIPTS/parameters.template ./parameters.template
    optfile="./parameters.template"
fi

##### user channel list
if [ -e $chanfile ]; then
    cp $chanfile ./channel.list

else
##### get reference frame file
    fflfile=`grep DATA $optfile | grep -m1 FFL | awk '{print $3}'`
    if [ -e $fflfile ]; then
	framefile=`head -1 $fflfile | awk '{print $1}'`
	gps=`head -1 $fflfile | awk '{print $2}'`
    else # let's try LCF 
	lcffile=`grep DATA $optfile | grep -m1 LCF | awk '{print $3}'`
	if [ -e $lcffile ]; then
	    framefile=`head -1 $lcffile | awk '{print $5}'`
	    gps=`head -1 $lcffile | awk '{print $3}'`
	else
	    echo "No FFL/LCF file in the option file"
	    exit 2
	fi
    fi
    if [ ! -e $framefile ]; then
	echo "No reference frame file"
	exit 2
    else
	echo "Reference frame file: $framefile"
    fi

##### get channel list
    ${FRROOT}/${FRCONFIG}/FrDump.exe -i $framefile -d 4 -f $gps -l $gps | grep Vector: | grep -v -w Auxiliary | sed 's|Vector:||g' | sed 's|dx=||g' | awk '$6<0.0039{print int(1.0/$6+0.5),$1}' | sort -n | uniq> ./channel.list
    if [ ! -s ./channel.list ]; then
	echo "No channel"
	exit 2
    fi
fi

##### break into freq
rm -f ./channel.*Hz
rm -f ./parameters_*Hz_*.txt
while read line; do

    freq=`echo $line | awk '{print $1}'`
    
    # keep only f0 > 16Hz
    if [ $freq -le 16 ]; then continue; fi

    echo $line >> ./channel.${freq}Hz
done < ./channel.list


##### loop over channels
for file in ./channel.*Hz; do
    freq=`head -1 $file | awk '{print $1}'`
        
    echo ""
    echo "Making parameter files for channels at $freq Hz..."
    echo ""

    # maximum number of channels to process
    if [ $freq -le 1024 ]; then nmax=30;
    elif [ $freq -le 2048 ]; then nmax=20;
    elif [ $freq -le 4096 ]; then nmax=10;
    else nmax=5; fi
    
    # sampling frequency
    freqsample=$freq
    if [ $freqsample -gt 4096 ]; then freqsample=4096; fi
    
    # frequency range
    freqmax=$(( $freqsample / 2 ))
    if [ $freq -le 1024 ]; then 
	freqmin=0.1; 
	overlap=160;
	chunk=864;
	block=512;
    else
	freqmin=16; 
	overlap=4;
	chunk=484;
	block=64;
    fi

    n=0
    p=0
    channel_list=""
    freq_list=""
    while read line; do

	channel=`echo $line | awk '{print $2}'`
	channel_list="$channel $channel_list"
	freq_list="$freq $freq_list"
	let "n+=1"

	if [ $n -eq $nmax ]; then
	    echo "#${p}"
	    sed -e '/DATA[ \t]*CHANNELS/d' \
		-e '/DATA[ \t]*NATIVEFREQUENCY/d' \
		-e '/DATA[ \t]*SAMPLEFREQUENCY/d' \
		-e '/PARAMETER[ \t]*CHUNKDURATION/d' \
		-e '/PARAMETER[ \t]*BLOCKDURATION/d' \
		-e '/PARAMETER[ \t]*OVERLAPDURATION/d' \
		-e '/PARAMETER[ \t]*FREQUENCYRANGE/d' \
		$optfile > ./parameters_${freq}Hz_${p}.txt

	    echo "" >> ./parameters_${freq}Hz_${p}.txt
	    echo "DATA  CHANNELS  ${channel_list}" >> ./parameters_${freq}Hz_${p}.txt
	    echo "DATA  NATIVEFREQUENCY  ${freq_list}" >> ./parameters_${freq}Hz_${p}.txt
	    echo "DATA  SAMPLEFREQUENCY  ${freqsample}" >> ./parameters_${freq}Hz_${p}.txt
	    echo "PARAMETER CHUNKDURATION  ${chunk}" >> ./parameters_${freq}Hz_${p}.txt
	    echo "PARAMETER BLOCKDURATION  ${block}" >> ./parameters_${freq}Hz_${p}.txt
	    echo "PARAMETER OVERLAPDURATION  ${overlap}" >> ./parameters_${freq}Hz_${p}.txt
	    echo "PARAMETER FREQUENCYRANGE  ${freqmin}  ${freqmax}" >> ./parameters_${freq}Hz_${p}.txt
	    channel_list=""
	    freq_list=""
	    let "p+=1"
	    n=0
	fi
    done < $file

    # one more file
    if [ $n -gt 0 ]; then
	echo "#${p}"
	sed -e '/DATA[ \t]*CHANNELS/d' \
	    -e '/DATA[ \t]*NATIVEFREQUENCY/d' \
	    -e '/DATA[ \t]*SAMPLEFREQUENCY/d' \
	    -e '/PARAMETER[ \t]*CHUNKDURATION/d' \
	    -e '/PARAMETER[ \t]*BLOCKDURATION/d' \
	    -e '/PARAMETER[ \t]*OVERLAPDURATION/d' \
	    -e '/PARAMETER[ \t]*FREQUENCYRANGE/d' \
	    $optfile > ./parameters_${freq}Hz_${p}.txt
	
	echo "" >> ./parameters_${freq}Hz_${p}.txt
	echo "DATA  CHANNELS  ${channel_list}" >> ./parameters_${freq}Hz_${p}.txt
	echo "DATA  NATIVEFREQUENCY  ${freq_list}" >> ./parameters_${freq}Hz_${p}.txt
	echo "DATA  SAMPLEFREQUENCY  ${freqsample}" >> ./parameters_${freq}Hz_${p}.txt
	echo "PARAMETER CHUNKDURATION  ${chunk}" >> ./parameters_${freq}Hz_${p}.txt
	echo "PARAMETER BLOCKDURATION  ${block}" >> ./parameters_${freq}Hz_${p}.txt
	echo "PARAMETER OVERLAPDURATION  ${overlap}" >> ./parameters_${freq}Hz_${p}.txt
	echo "PARAMETER FREQUENCYRANGE  ${freqmin}  ${freqmax}" >> ./parameters_${freq}Hz_${p}.txt
    fi

done

# cleaning
rm -f ./channel.*Hz

exit 0