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
fflfile="none"
gps="none"
windows="2 8 32"
snrthr=5

printhelp(){
    echo ""
    echo "Usage:"
    echo "MakeOmicronScan -c [CHANNEL_LIST] -f [FFL_FILE] -g [GPS_CENTER]"
    echo ""
    echo "*** MAIN OPTIONS"
    echo "  -c  [CHANNEL_LIST]    path to a text file listing the channels to be scanned"
    echo "                        listed in a single column."
    echo ""
    echo "  -f  [FFL_FILE]        path to FFL file pointing to the data to be used."
    echo "                        data must be available at least within a time window"
    echo "                        twice as large as the largest requested time window"
    echo "                        centered on [GPS_CENTER]"
    echo ""
    echo "  -g  [GPS_CENTER]      central GPS time of of the scan"
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
    echo "  -h                    print this help"
    echo ""
} 


##### Check the environment
if [[ -z "$OMICRONROOT" ]]; then
    echo "Error: The Omicron environment is not set"
    exit 1
fi

##### read options
while getopts ":c:f:g:d:w:s:h" opt; do
    case $opt in
	c)
	    channelfile="$OPTARG"
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
dateUTC=`tconvert $gps`

##### check gps
if [ "$gps" = "none" ]; then
    echo "ERROR: A GPS time must be provided with the '-g' option"
    echo "type  'MakeOmicronScan -h'  for help"
    exit 1
fi

##### check FFL
if [ "$fflfile" = "none" ]; then
    echo "ERROR: A FFL file must be provided with the '-f' option"
    echo "type  'MakeOmicronScan -h'  for help"
    exit 1
fi
if [ ! -e $fflfile ]; then
    echo "ERROR: The FFL file $fflfile does not exist"
    echo "type  'MakeOmicronScan -h'  for help"
    exit 1
fi

##### check channel list
if [ "$channelfile" = "none" ]; then
    echo "ERROR: A channel file must be provided with the '-c' option"
    echo "type  'MakeOmicronScan -h'  for help"
    exit 1
fi
if [ ! -e $channelfile ]; then
    echo "ERROR: The channel file $channelfile does not exist"
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

##### prepare outdir
echo "MakeOmicronScan: Output directory = ${outdir}/${gps}"
outdir="${outdir}/${gps}"
mkdir -p ${outdir}; #rm -fr ${outdir}/*
uniq -u $channelfile > ${outdir}/channels.list

##### prepare option file
echo "MakeOmicronScan: Make option file..."
echo "//***************** OmiScan option file *****************" > ${outdir}/options.txt
echo "VERBOSITY   LEVEL       1"           >> ${outdir}/options.txt
echo "DATA        FFL         ${fflfile}"  >> ${outdir}/options.txt
echo "PARAMETER   WINDOW      ${windows}"  >> ${outdir}/options.txt
echo "PARAMETER   SNR_THRESHOLD ${snrthr}" >> ${outdir}/options.txt
echo "OUTPUT      DIRECTORY   ${outdir}"   >> ${outdir}/options.txt

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


##### html template and web material
echo "MakeOmicronScan: Make web report..."
template=${outdir}/index.template
cp -f ${OMICRON_HTML}/template/template.omicronscan.html $template
cp -f ${OMICRON_HTML}/template/comparison_mode.html ${outdir}/
cp -f ${GWOLLUM_DOC}/style.css ${outdir}/style.css
cp -f ${GWOLLUM_DOC}/Pics/gwollum_logo_min_trans.gif ${outdir}/icon.gif
cp -f ${OMICRON_HTML}/pics/omicronlogo_xxl.gif ${outdir}/omicronlogo_xxl.gif

##### fill html report
while read channel; do

    if ! grep -qw "$channel" ${outdir}/summary.txt; then
	echo "                 channel ${channel} was not scanned --> skip"
	continue
    fi
    echo "                 ${channel}..."

    default="<tr>"

    cd ${outdir}/plots
    for win in $windows; do
	ln -sf ./th_${channel}_psd.gif ./th_${channel}_psd_dt${win}.gif
        ln -sf ./${channel}_psd.gif ./${channel}_psd_dt${win}.gif
	default="${default}<td><a id=\"a_${channel}_dt${win}\" href=\"./plots/${channel}_map_dt${win}.gif\"><img id=\"img_${channel}_dt${win}\" src=\"./plots/th_${channel}_map_dt${win}.gif\" alt=\"${channel}_map\" /></a></td>"
    done
    cd $here
    default="${default}</tr>"

    maps="<a href=\"javascript:showImage('./plots', '$channel', 'map', ${winids});\">full</a>"
    tseries="<a href=\"javascript:showImage('./plots', '$channel', 'raw', ${winids});\">raw</a>"
    other="<a href=\"javascript:showImage('./plots', '$channel', 'psd', ${winids});\">psd</a>"
    
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