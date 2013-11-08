#!/bin/bash
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
lalfile="none"
fflfile="none"
cachefile="none"
gps=0
windows="2 8 32"
snrthr=8
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
    echo ""
    echo "  -l  [LAL_CACHE_FILE]  path to a lal-cache file pointing to the data to be used."
    echo "                        the user must make sure the data are available"
    echo "                        for the requested time stretch"
    echo ""
    echo "  -W  [FRAMECACHE_FILE] path to a frame-cache file pointing to the data to be used."
    echo "                        the user must make sure the data are available"
    echo "                        for the requested time stretch"
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
    echo "                        By default [SNR_THRESHOLD] = 5"
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
    echo "  -S                    this flag commands the plotting style"
    echo "                        the standard ROOT style will be used"
    echo "                        instead of the GWOLLUM style"
    echo ""
    echo "  -h                    print this help"
    echo ""
} 


##### Check the environment
if [[ -z "$OMICRONROOT" ]]; then
    echo "Error: The Omicron environment is not set"
    exit 1
fi

##### read options
while getopts ":c:O:f:l:W:g:d:w:s:v:Sh" opt; do
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
	l)
	    lalfile="$OPTARG"
	    ;;
	W)
	    cachefile="$OPTARG"
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
mkdir -p ${outdir}; #rm -fr ${outdir}/*
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
    grep "channelName" $omegafile | awk '{print $2}' | sed "s|'||g" | uniq -u > ${outdir}/channels.list
else
    if [ ! -e $channelfile ]; then
	echo "ERROR: The channel file $channelfile does not exist"
	echo "type  'MakeOmicronScan -h'  for help"
	exit 1
    else
	uniq -u $channelfile > ${outdir}/channels.list
    fi
fi

##### check FFL/LAL/CACHE
if [ "$fflfile" = "none" ]; then
    if [ "$lalfile" = "none" ]; then
	if [ "$cachefile" = "none" ]; then
	    echo "ERROR: A FFL/lal-cache/framecache file must be provided with the '-f'/'-l'/'-W' option"
	    echo "type  'MakeOmicronScan -h'  for help"
	    exit 1
	else
	    if [ ! -e $cachefile ]; then
		echo "ERROR: The lal-cache file $cachefile does not exist"
		echo "type  'MakeOmicronScan -h'  for help"
		exit 1
	    fi
	    fflfile=$outdir/ffl_omiscan_${gps}.ffl
	    $GWOLLUM_SCRIPTS/framecache2ffl.sh $cachefile > $fflfile
	fi
    else
	if [ ! -e $lalfile ]; then
	    echo "ERROR: The lal-cache file $lalfile does not exist"
	    echo "type  'MakeOmicronScan -h'  for help"
	    exit 1
	fi
	fflfile=$outdir/ffl_omiscan_${gps}.ffl
	$GWOLLUM_SCRIPTS/lalcache2ffl.sh $lalfile > $fflfile
    fi
fi
if [ ! -e $fflfile ]; then
    echo "ERROR: The FFL file $fflfile does not exist"
    echo "type  'MakeOmicronScan -h'  for help"
    exit 1
fi

##### check outdir
if [ ! -d $outdir ]; then
    echo "ERROR: The output directory $outdir does not exist"
    echo "type  'MakeOmicronScan -h'  for help"
    exit 1
fi

##### check windows
winids="["
nwin=0
for w in $windows; do
    if [ $w -lt 1 ]; then
	echo "ERROR: The time window must be an integer"
	echo "type  'MakeOmicronScan -h'  for help"
	exit 1
    fi
    winids="${winids}'dt${w}',"
    let "nwin+=1"
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
echo "//***************** OmiScan option file *****************" > ${outdir}/options.txt
echo "VERBOSITY   LEVEL       ${verbose}"  >> ${outdir}/options.txt
echo "DATA        FFL         ${fflfile}"  >> ${outdir}/options.txt
echo "PARAMETER   WINDOW      ${windows}"  >> ${outdir}/options.txt
echo "PARAMETER   SNR_THRESHOLD ${snrthr}" >> ${outdir}/options.txt
echo "OUTPUT      DIRECTORY   ${outdir}"   >> ${outdir}/options.txt
echo "OUTPUT      STYLE       ${style}"    >> ${outdir}/options.txt

##### run omiscan
echo "MakeOmicronScan: Run OmiScan..."
omiscan.exe ${outdir}/channels.list ${outdir}/options.txt $gps > ${outdir}/log.txt 2>&1
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
    other="<a href=\"javascript:showImage('./plots', '$channel', 'asd', ${winids});\">asd</a>"
    
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