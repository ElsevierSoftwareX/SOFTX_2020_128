#!/bin/sh
#
# MakeOmicronScan.sh
#
# generate omicron scan from scratch
#
# Author: Florent Robinet
# florent.robinet@lal.in2p3.fr

##### default options
here=`pwd`
outdir=`pwd` # output directory
channelfile="none"
omegafile="none"
fflfile="none"
gps=0
windows="2 8 32"
snrthr=8
writeroot=""
verbose=0
style="GWOLLUM"

printhelp(){
    echo ""
    echo "Usage:"
    echo "MakeOmicronScan -c [CHANNEL_LIST] -f [FFL_FILE] -g [GPS_CENTER]"
    echo ""
    echo "*** MAIN OPTIONS"
    echo "  -c  [CHANNEL_LIST]    path to a text file listing the channels to be scanned"
    echo "                        listed in a single column."
    echo ""
    echo "  -f  [FFL_FILE]        path to a FFL file pointing to the data to be used."
    echo "                        the user must make sure the data are available"
    echo "                        for the requested time stretch"
    echo "                        Several formats are supported:"
    echo "                        - Virgo FFL"
    echo "                        - LIGO lalcache"
    echo "                        - LIGO framecache"
    echo ""
    echo "  -g  [GPS_CENTER]      central GPS time of of the scan"
    echo ""
    echo "  -O  [OMEGA_CONFIG]    path to an omega-scan configuration file."
    echo "                        this option can be given instead of a channel list file."
    echo "                        the channel list will be extracted from the config file"
    echo ""
    echo "*** PARAMETERS"
    echo "  -w  [TIME_WINDOWS]    set of time windows"
    echo "                        By default [TIME_WINDOWS] = \"1 8 64\""
    echo ""
    echo "  -s  [SNR_THRESHOLD]   SNR threshold below which the channel is not plotted"
    echo "                        By default [SNR_THRESHOLD] = 8"
    echo ""
    echo "*** OUTPUT"
    echo "  -d  [OUTDIR]          path to output directory"
    echo "                        Default = current directory"
    echo ""
    echo "  -v  [VERBOSE_LEVEL]   verbosity level"
    echo "                        0 --> minimum printing"
    echo "                        1 --> progress printing"
    echo "                        2 --> parameter printing"
    echo "                        3 --> full printing"
    echo "                        Default = 0"
    echo ""
    echo "  -R                    with this flag, all the plots are saved in a ROOT file"
    echo ""
    echo "  -S                    this flag commands the plotting style"
    echo "                        the standard ROOT style will be used"
    echo "                        instead of the GWOLLUM style"
    echo ""
    echo "  -h                    print this help"
    echo ""
} 


##### Check the environment
if [ -z "$OMICRONROOT" ]; then
    echo "Error: The Omicron environment is not set"
    exit 1
fi

##### read options
while getopts ":c:O:f:l:W:g:d:w:s:v:SRh" opt; do
    case $opt in
	c)
	    channelfile="$OPTARG"
	    ;;
	O)
	    omegafile="$OPTARG"
	    ;;
	f)
	    fflfile="$OPTARG"
	    ;;
	g)
	    gps="$OPTARG"
	    ;;
	w)
	    windows="$OPTARG"
	    ;;
	s)
	    snrthr="$OPTARG"
	    ;;
	d)
	    outdir="$OPTARG"
	    ;;
	v)
	    verbose="$OPTARG"
	    ;;
	S)
	    style="STANDARD"
	    ;;
	R)
	    writeroot="root"
	    ;;
	h)
	    printhelp
	    exit 0
	    ;;
	\?)
	    echo "Invalid option: -$OPTARG"
	    echo "type  'MakeOmicronScan -h'  for help"
	    exit 1
	    ;;
    esac
done

##### format gps
gps=`echo $gps | awk '{printf "%.3f", $1}'`

##### check gps
if [ "$gps" = "0.000" ]; then
    echo "ERROR: A GPS time must be provided with the '-g' option"
    echo "type  'MakeOmicronScan -h'  for help"
    exit 1
fi
dateUTC=`tconvert $gps`

##### check outdir
if [ ! -d $outdir ]; then
    echo "ERROR: The output directory $outdir does not exist"
    echo "type  'MakeOmicronScan -h'  for help"
    exit 1
fi

##### prepare outdir
echo "MakeOmicronScan: Output directory = ${outdir}/${gps}"
echo "                 You can check the progress of your scan when loading this address in your web browser"
outdir="${outdir}/${gps}"
mkdir -p ${outdir};
template=${outdir}/index.template
cp -f ${OMICRON_HTML}/template/template.omicronscan.html $template
cp -f ${OMICRON_HTML}/template/comparison_mode.html ${outdir}/
cp -f ${GWOLLUM_DOC}/style.css ${outdir}/style.css
cp -f ${GWOLLUM_DOC}/Pics/gwollum_logo_min_trans.gif ${outdir}/icon.gif
cp -f ${OMICRON_HTML}/pics/omicronlogo_xxl.gif ${outdir}/omicronlogo_xxl.gif

##### check channel list
if [ "$channelfile" = "none" ]; then
    if [ "$omegafile" = "none" ]; then
	echo "ERROR: A channel file must be provided with the '-c' option"
	echo "type  'MakeOmicronScan -h'  for help"
	exit 1
    fi
    if [ ! -e $omegafile ]; then
	echo "ERROR: The omega config file $omegafile does not exist"
	echo "type  'MakeOmicronScan -h'  for help"
	exit 1
    fi
    grep "channelName" $omegafile | awk '{print $2}' | sed "s|'||g" | awk '{ if (!h[$0]) { print $0; h[$0]=1 } }' > ${outdir}/channels.list
else
    if [ ! -e $channelfile ]; then
	echo "ERROR: The channel file $channelfile does not exist"
	echo "type  'MakeOmicronScan -h'  for help"
	exit 1
    else
	cat $channelfile | awk '{ if (!h[$0]) { print $0; h[$0]=1 } }' > ${outdir}/channels.list
    fi
fi
if [ ! -s ${outdir}/channels.list ]; then
    echo "ERROR: The channel file $channelfile is empty"
    echo "type  'MakeOmicronScan -h'  for help"
    exit 1
fi

##### check FFL/LAL/CACHE
if [ ! -e $fflfile ]; then
    echo "ERROR: The ffl file $fflfile is missing"
    echo "type  'MakeOmicronScan -h'  for help"
    exit 1
fi

# detect framecache format
ncol=`awk -F'|' '{print NF; exit}' $fflfile`
if [ $ncol -eq 6 ]; then # framecache
    ${GWOLLUM_SCRIPTS}/framecache2ffl.sh $fflfile > ${outdir}/ffl_omiscan_${gps}.ffl
    fflfile=${outdir}/ffl_omiscan_${gps}.ffl
fi

# extract channels
printchannels.exe $fflfile > ${outdir}/channels.total

##### check windows
winids="["
nwin=0
winmax=0
for w in $windows; do
    if [ $w -lt 1 ]; then
	echo "ERROR: The time window must be an integer"
	echo "type  'MakeOmicronScan -h'  for help"
	exit 1
    fi
    winids="${winids}'dt${w}',"
    nwin=$(( nwin + 1 ))
    if [ $w -gt $winmax ]; then winmax=$w; fi
done
winids="${winids%?}]"
if [ $nwin -eq 0 ]; then
    echo "ERROR: There must be at least one time window"
    echo "type  'MakeOmicronScan -h'  for help"
    exit 1
fi

currentdate=`date -u`

##### waiting index
title="Omiscan of $gps"
message="Your Omiscan is currently running. The web report will appear here when it's done.<br>In the meantime you can check the <a href=\"./log.txt\">log file</a> to monitor the progress of the scan.<br>-- <i>$USER - $currentdate</i> --"
sed -e "s|\[TITLE\]|${title}|g" \
    -e "s|\[MESSAGE\]|${message}|g" \
    $GWOLLUM_HTML/template/template.waiting.html > ${outdir}/index.html


##### prepare option file
echo "MakeOmicronScan: Make option file..."
paramfile=${outdir}/parameters.txt
paramfile_low=${outdir}/parameters.low.txt
paramfile_std=${outdir}/parameters.std.txt
echo "//***************** OmiScan option file *****************" > $paramfile
echo "DATA        FFL           ${fflfile}"   >> $paramfile
echo "PARAMETER   WINDOW        ${windows}"   >> $paramfile
echo "TRIGGER     SNRTHRESHOLD  ${snrthr}"    >> $paramfile
echo "OUTPUT      DIRECTORY     ${outdir}"    >> $paramfile
echo "OUTPUT      LEVEL         ${verbose}"   >> $paramfile
echo "OUTPUT      PLOTSTYLE     ${style}"     >> $paramfile
echo "OUTPUT      FORMAT        gif ${root}"  >> $paramfile

cp -f $paramfile $paramfile_low
cp -f $paramfile $paramfile_std

echo "PARAMETER   CHUNKDURATION   $(( 2 * winmax ))"  >> $paramfile_std
echo "PARAMETER   BLOCKDURATION   $(( 2 * winmax ))"  >> $paramfile_std
echo "PARAMETER   OVERLAPDURATION 2"                  >> $paramfile_std
echo "DATA	  SAMPLEFREQUENCY 2048"               >> $paramfile_std
echo "PARAMETER   FREQUENCYRANGE  16  1024"           >> $paramfile_std

echo "PARAMETER   CHUNKDURATION   $(( 2 * winmax ))"  >> $paramfile_low
echo "PARAMETER   BLOCKDURATION   $(( 2 * winmax ))"  >> $paramfile_low
echo "PARAMETER   OVERLAPDURATION 2"                  >> $paramfile_low
echo "DATA	  SAMPLEFREQUENCY 256"                >> $paramfile_low
echo "PARAMETER   FREQUENCYRANGE  0.1  128"           >> $paramfile_low

while read line; do
    channelname=`echo $line | awk '{print $1}'`
    sampling=`grep -w $channelname ${outdir}/channels.total | awk '{print $2}'`
    if [ $sampling -gt 1024 ]; then 
	echo "DATA CHANNELS  $channelname"            >> $paramfile_std
    else
	echo "DATA CHANNELS  $channelname"            >> $paramfile_low
    fi
done < ${outdir}/channels.list

##### run omiscan
echo "MakeOmicronScan: Run OmiScan (std)..."
omiscan.exe $gps $paramfile_std > ${outdir}/log.txt 2>&1
omiscan.exe $gps $paramfile_low >> ${outdir}/log.txt 2>&1
exit

if [ ! -e ${outdir}/summary.txt ]; then
    echo "                 omiscan failed"
    exit 2
fi
if ! grep -qw "omiscan done" ${outdir}/summary.txt; then
    echo "                 omiscan not ended"
    exit 2
fi

##### make thunbnails
echo "MakeOmicronScan: Make thumbnails..."
while read line; do
    channel=`echo $line | awk '{print $1}'`
    for file in ${outdir}/plots/${channel}*.gif; do     
	if [ -e $file ]; then 	
	    convert -density 100 -thumbnail 320  ${file} ${outdir}/plots/th_${file##*/}     
	fi 
    done
done < ${outdir}/summary.txt

##### fill html report
echo "MakeOmicronScan: Make web report..."
while read channel; do

    if ! grep -qw "$channel" ${outdir}/summary.txt; then
	echo "                 channel ${channel} was not scanned --> skip"
	continue
    fi
    echo "                 ${channel}..."

    default="<tr>"

    cd ${outdir}/plots
    for win in $windows; do
	ln -sf ./th_${channel}_asd.gif ./th_${channel}_asd_dt${win}.gif
        ln -sf ./${channel}_asd.gif ./${channel}_asd_dt${win}.gif
	default="${default}<td><a id=\"a_${channel}_dt${win}\" href=\"./plots/${channel}_map_dt${win}.gif\"><img id=\"img_${channel}_dt${win}\" src=\"./plots/th_${channel}_map_dt${win}.gif\" alt=\"${channel}_map\" /></a></td>"
    done
    cd $here
    default="${default}</tr>"

    maps="<a href=\"javascript:showImage('./plots', '$channel', 'map', ${winids});\">full</a>"
    tseries="<a href=\"javascript:showImage('./plots', '$channel', 'raw', ${winids});\">raw</a>"
    other="<a href=\"javascript:showImage('./plots', '$channel', 'asd', ${winids});\">asd</a> <a href=\"javascript:showImage('./plots', '$channel', 'projt', ${winids});\">SNR vs time</a> <a href=\"javascript:showImage('./plots', '$channel', 'projf', ${winids});\">SNR vs frequency</a>"

    if [ $writeroot -eq 1 ]; then
	other="${other} <a href=\"./plots/${channel}.root\">Get ROOT file</a>"
    fi

    # add q map links
    q=0;
    for Qval in `grep -w "$channel" ${outdir}/summary.txt | awk '{for(i=1;i<9;i++) $i="";print}'`; do
	maps="${maps} <a href=\"javascript:showImage('./plots', '$channel', 'map_Q${q}', ${winids});\">Q=${Qval}</a>"
	let "q+=1"
    done

    # loudest event
    loudest_gps=`grep -w "$channel" ${outdir}/summary.txt | awk '{printf "%.3f",$4}'`
    loudest_freq=`grep -w "$channel" ${outdir}/summary.txt | awk '{printf "%.2f",$5}'`
    loudest_snr=`grep -w "$channel" ${outdir}/summary.txt | awk '{printf "%.2f",$6}'`
    loudest_snr=`grep -w "$channel" ${outdir}/summary.txt | awk '{printf "%.2f",$6}'`
    loudest_q=`grep -w "$channel" ${outdir}/summary.txt | awk '{printf "%.2f",$7}'`

    # get channel description
    if [ -e ${OMICRON_PARAMETERS}/channels.txt ]; then
	chan_des=` grep -w ${channel} ${OMICRON_PARAMETERS}/channels.txt`
	if [ "$chan_des" = "" ]; then
	    description="No description for this channel"
	else
	    description=`echo $chan_des | awk '{$1="";print}'`
	fi
    else
	description="No description for this channel"
    fi

    sed -i "/<\!-- channels -->/i<tr><td colspan=\"${q}\"><h2>${channel}</h2><h3>${description}</h3>Loudest tile: GPS=${loudest_gps}, f=${loudest_freq}Hz, Q=${loudest_q}, <b>SNR=${loudest_snr}</b></td></tr>" $template
    sed -i "/<\!-- channels -->/i<tr><td colspan=\"${q}\">Maps: $maps \| Time series: $tseries \| Other: $other</td></tr>" $template
    sed -i "/<\!-- channels -->/i${default}" $template
done < ${outdir}/channels.list

##### fill static variables
sed -e "s|\[GPS_CENTER\]|${gps}|g" \
    -e "s|\[DATEUTC\]|${dateUTC}|g" \
    -e "s|\[USERNAME\]|${USER}|g" \
    $template > ${outdir}/index.html

##### cleaning
rm -f $template

##### finish log
currentdate=`date -u`

exit 0