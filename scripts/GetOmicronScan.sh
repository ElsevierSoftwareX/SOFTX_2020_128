#!/bin/bash
#
# GetOmicronScan.sh
#
# generate omicron scan
#
# Author: Florent Robinet
# florent.robinet@lal.in2p3.fr

##### default options
here=`pwd`
main="h_4096Hz" # channel name
outdir=`pwd` # output directory
snrmin=8
parameterfile=""
tagname="NONE"
chandir="NONE"

printhelp(){
    echo ""
    echo "Usage:"
    echo "GetOmicronScan -m [MAIN_CHANNEL] [GPS_CENTER]"
    echo ""
    echo "TRIGGER SELECTION OPTIONS"
    echo "  -s  [SNRMIN]         print aux. scan only if SNR > [SNRMIN]"
    echo "                       Default = 8"
    echo "  -t  [TAGNAME]        Predefined channel selection (Virgo only):"
    echo "                       [TAGNAME] = std  -->  standard setting"
    echo "                       [TAGNAME] = min  -->  minimal setting with only the GW channels"
    echo "                       [TAGNAME] = sel  -->  only a \"happy few\" channels are selected"
    echo "  -p  [PARAMETER_FILE] if given, this option allows you to scan a user-defined selection"
    echo "                       of channels"
    echo "                       [PARAMETER_FILE] is the path to a user parameter file"
    echo "                       This file should contain a single column with the list"
    echo "                       of channels to scan"
    echo "                       When this option is given the option '-t' is ignored"
    echo ""
    echo "  -r  [TRIGGER_DIR]    The user can provide his own directory where the triggers are located"
    echo "                       [TRIGGER_DIR] must contain subdirectories for each channels:"
    echo "                       [TRIGGER_DIR]/channel_1/"
    echo "                       [TRIGGER_DIR]/channel_2/"
    echo "                       ..."
    echo "                       [TRIGGER_DIR]/channel_N/"
    echo ""
    echo "OUTPUT CONTROL"
    echo "  -m  [MAIN_CHANNEL]   main channel: always plotted"
    echo "                       Default = h_4096Hz"
    echo "  -d  [OUTDIR]         path to output directory"
    echo "                       Default = current directory"
    echo ""
    echo "  -h                   print this help"
    echo ""
} 


##### Check the environment
if [[ -z "$OMICRONROOT" ]]; then
    echo "Error: The Omicron environment is not set"
    exit 1
fi

##### needs argument
if [ $# -lt 1 ]; then
    printhelp
    exit 1
fi

##### read options
while getopts ":m:s:p:d:t:r:h" opt; do
    case $opt in
	m)
	    main="$OPTARG"
	    ;;
	s)
	    snrmin=`echo $OPTARG | awk '{print int($1)}'`
	    ;;
	p)
	    parameterfile="$OPTARG"
	    ;;
	t)
	    tagname="$OPTARG"
	    ;;
	r)
	    chandir="$OPTARG"
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
	    echo "type  'GetOmicronScan -h'  for help"
	    exit 1
	    ;;
    esac
done

##### gps interval
shift $(($OPTIND - 1))
tcenter=`echo $1 | awk '{printf "%.3f", $1}'`
tcenter_=`echo $tcenter | awk '{print int($1)}'`
OPTIND=0

##### check timing
if [ $tcenter_ -lt 700000000 ]; then
    echo "Invalid option: '$tcenter' is not a reasonable central time"
    echo "type  'GetOmicronScan -h'  for help"
    exit 1
fi

##### scan user channel directory
if [ ! $chandir = "NONE" ]; then
    if [ ! -d $chandir ] ; then
	echo "Invalid option: the channel directory $chandir cannot be found"
	echo "type  'GetOmicronScan -h'  for help"
	exit 1
    fi
    cd $chandir
    for dir in *; do
	if [ ! -d ${dir} ]; then continue; fi
	OMICRON_CHANNELS="${OMICRON_CHANNELS} ${dir}"
    done
    echo "channel candidates: $OMICRON_CHANNELS"
    cd $here

else
    ##### select run
    run="NONE"
    for r in $RUN_NAMES; do
	r_s=${r}_START
	r_e=${r}_END
	if [[ $tcenter_ -ge ${!r_s} && $tcenter_ -lt ${!r_e} ]]; then
	    run=$r
	    break;
	fi
    done
    
    if [ $run = "NONE" ]; then
	echo "Invalid GPS time"
	exit 1 
    fi

    ##### get available channels
    . GetOmicronChannels.sh -r $run > /dev/null 2>&1
fi
##### check main channel is available
if ! echo "$OMICRON_CHANNELS" | grep -q "$main"; then
    echo "Invalid option: channel '${main}' is not available"
    echo "type  'GetOmicronScan -h'  for help"
    exit 1
fi

##### elaborate timing
dateUTC=`tconvert -f "%A the %dth, %B %Y %H:%M:%S" ${tcenter_}`
tstart=$(( $tcenter_ - 100 ))
tstop=$(( $tcenter_ + 100 ))

##### check outdir
if [ ! -d $outdir ] ; then
    echo "Invalid option: the output directory $outdir cannot be found"
    echo "type  'GetOmicronScan -h'  for help"
    exit 1
fi

##### check parameter file if given
if [ ! $parameterfile = "" ] ; then
    if [ ! -e $parameterfile ] ; then
	echo "Invalid option: the parameter file $parameterfile cannot be found"
	echo "type  'GetOmicronScan -h'  for help"
	exit 1
    fi
    if [ ! -s $parameterfile ] ; then
	echo "Invalid option: the parameter file $parameterfile is empty"
	echo "type  'GetOmicronScan -h'  for help"
	exit 1
    fi
fi

##### get parameter file from tag
if [ "$parameterfile" = "" ]; then
    if [ "$tagname" = "std" ]; then
	if [ -e ${OMICRON_PARAMETERS}/scan.parameters.${run}.std.txt ]; then
	    parameterfile=${OMICRON_PARAMETERS}/scan.parameters.${run}.std.txt
	else
	    echo "Invalid option: no 'std' parameter file for ${run}"
	    echo "type  'GetOmicronScan -h'  for help"
	    exit 1
	fi

    elif [ "$tagname" = "min" ]; then
	if [ -e ${OMICRON_PARAMETERS}/scan.parameters.${run}.min.txt ]; then
	    parameterfile=${OMICRON_PARAMETERS}/scan.parameters.${run}.min.txt
	else
	    echo "Invalid option: no 'min' parameter file for ${run}"
	    echo "type  'GetOmicronScan -h'  for help"
	    exit 1
	fi

    elif [ "$tagname" = "sel" ]; then
	if [ -e ${OMICRON_PARAMETERS}/scan.parameters.${run}.sel.txt ]; then
	    parameterfile=${OMICRON_PARAMETERS}/scan.parameters.${run}.sel.txt
	else
	    echo "Invalid option: no 'sel' parameter file for ${run}"
	    echo "type  'GetOmicronScan -h'  for help"
	    exit 1
	fi
	
    else # scan everything
	tagname=$tagname
    fi
fi

##### update channel list if parameter file is given
if [ ! $parameterfile = "" ] ; then
    OMICRON_CHANNELS_NEW=""
    while read line; do
	OMICRON_CHANNELS_NEW="$OMICRON_CHANNELS_NEW $line"
    done < $parameterfile
    OMICRON_CHANNELS=$OMICRON_CHANNELS_NEW
fi

##### build directory
outdir="${outdir}/${tcenter}"
mkdir -p ${outdir}; rm -fr ${outdir}/*

##### html template and web material
template=${outdir}/index.template
cp -f ${OMICRON_HTML}/template/template.omicronscan.html $template
cp -f ${OMICRON_HTML}/template/comparison_mode.html ${outdir}/
cp -f ${GWOLLUM_DOC}/style.css ${outdir}/style.css
cp -f ${OMICRON_HTML}/pics/omicronlogo_xxl.gif ${outdir}/omicronlogo_xxl.gif
if [ ! $parameterfile = "" ] ; then cp -f $parameterfile ${outdir}/parameters.txt; fi
currentdate=`date -u`

##### start log
echo "**********************************" > ${outdir}/log.txt
echo "**  Start Omiscan of $tcenter" >> ${outdir}/log.txt
echo "**********************************" >> ${outdir}/log.txt
echo "" >> ${outdir}/log.txt
echo "Current date: $currentdate" >> ${outdir}/log.txt
echo "" >> ${outdir}/log.txt


##########################################################################
##########################     MAIN CHANNEL     ##########################
##########################################################################
echo "Scanning ${main}..." >> ${outdir}/log.txt
# map the main channel triggers
if [ $chandir = "NONE" ]; then triggers=`GetTriggerFileList.sh -c${main} $tstart $tstop | grep "FILELIST" | sed 's|FILELIST ||g'`
else triggers="${chandir}/${main}/*.root"
fi

if [ "$triggers" = "" ]; then
    echo "triggers are not available for $main at $tcenter"
    exit 1
fi
eventmap.exe ${outdir} "${triggers}" $tcenter >> ${outdir}/log.txt 2>&1

# generic names for plots and thumbnails
mv ${outdir}/map_${tcenter}_dt4.gif ${outdir}/${main}_${tcenter}_dt4.gif
mv ${outdir}/map_${tcenter}_dt16.gif ${outdir}/${main}_${tcenter}_dt16.gif
mv ${outdir}/map_${tcenter}_dt64.gif ${outdir}/${main}_${tcenter}_dt64.gif
mv ${outdir}/info_${tcenter}.txt ${outdir}/${main}_${tcenter}.txt 
convert -density 100 -thumbnail 320  ${outdir}/${main}_${tcenter}_dt4.gif ${outdir}/${main}_${tcenter}_dt4_th.gif
convert -density 100 -thumbnail 320  ${outdir}/${main}_${tcenter}_dt16.gif ${outdir}/${main}_${tcenter}_dt16_th.gif
convert -density 100 -thumbnail 320  ${outdir}/${main}_${tcenter}_dt64.gif ${outdir}/${main}_${tcenter}_dt64_th.gif
    
# get plot info
loudgps=`grep -m 1 -w LOUDESTGPS_dt4 ${outdir}/${main}_${tcenter}.txt  | head -1 | awk '{printf "%.3f", $2}'`
loudf=`grep -m 1 -w LOUDESTFREQ_dt4 ${outdir}/${main}_${tcenter}.txt  | head -1 | awk '{printf "%.1f", $2}'`
loudsnr=`grep -m 1 -w LOUDESTSNR_dt4 ${outdir}/${main}_${tcenter}.txt  | head -1 | awk '{printf "%.2f", $2}'`

# get channel description
if [ -e ${OMICRON_PARAMETERS}/channels.txt ]; then

    chan_des=` grep -w ${main} ${OMICRON_PARAMETERS}/channels.txt`
    if [ "$chan_des" = "" ]; then
	description="No description for this channel"
    else
	description=`echo $chan_des | sed -e 's/^\w*\ *//'`
    fi
else
    description="No description for this channel"
fi

# fill html report
sed -i "/<\!-- channels -->/i\<tr><td><table><tr><td colspan=\"3\"><h2>${main}</h2><h3>${description}</h3>Loudest tile: GPS=${loudgps}, f=${loudf}Hz, SNR=${loudsnr}</td></tr>" $template
sed -i "/<\!-- channels -->/i\<tr><td><a href=\"./${main}_${tcenter}_dt4.gif\"><img src=\"./${main}_${tcenter}_dt4_th.gif\"></a></td><td><a href=\"./${main}_${tcenter}_dt16.gif\"><img src=\"./${main}_${tcenter}_dt16_th.gif\"></a></td><td><a href=\"./${main}_${tcenter}_dt64.gif\"><img src=\"./${main}_${tcenter}_dt64_th.gif\"></a></td></tr></table></td></tr>" $template

##########################################################################
##########################     AUX. CHANNEL     ##########################
##########################################################################

naux=0
naux_with_triggers=0
naux_above_thr=0
for channel in $OMICRON_CHANNELS; do

    # exclude main channel
    if [ "$channel" = "$main" ]; then continue; fi
    echo "Scanning ${channel}..." >> ${outdir}/log.txt
    let "naux+=1"

    # map the channel triggers
    if [ $chandir = "NONE" ]; then triggers=`GetTriggerFileList.sh -c${channel} $tstart $tstop | grep "FILELIST" | sed 's|FILELIST ||g'`
    else triggers="${chandir}/${channel}/*.root"
    fi

    if [ "$triggers" = "" ]; then 
	echo "  no trigger ---> Do not plot" >> ${outdir}/log.txt
	continue
    fi

    let "naux_with_triggers+=1"
    eventmap.exe ${outdir} "${triggers}" $tcenter >> ${outdir}/log.txt 2>&1

    if [ ! -e ${outdir}/map_${tcenter}_dt4.gif ]; then
	echo "  no trigger files ---> Do not plot" >> ${outdir}/log.txt
	continue
    fi

    # generic names for plots
    mv ${outdir}/map_${tcenter}_dt4.gif ${outdir}/${channel}_${tcenter}_dt4.gif
    mv ${outdir}/map_${tcenter}_dt16.gif ${outdir}/${channel}_${tcenter}_dt16.gif
    mv ${outdir}/map_${tcenter}_dt64.gif ${outdir}/${channel}_${tcenter}_dt64.gif
    mv ${outdir}/info_${tcenter}.txt ${outdir}/${channel}_${tcenter}.txt 

    #get loudest SNR within 16s
    loudsnr=`grep -m 1 -w LOUDESTSNR_dt4 ${outdir}/${channel}_${tcenter}.txt  | head -1 | awk '{printf "%.2f", $2}'`
    loudsnr_=`echo $loudsnr | awk '{print int($1)}'`

    # apply snr threshold
    if [ $loudsnr_ -lt $snrmin ]; then
	rm -f ${outdir}/${channel}_${tcenter}*
	echo "  SNR = ${loudsnr} < $snrmin ---> Do not plot" >> ${outdir}/log.txt
	continue;
    fi
    let "naux_above_thr+=1"

    # get plot info
    loudgps=`grep -m 1 -w LOUDESTGPS_dt4 ${outdir}/${channel}_${tcenter}.txt  | head -1 | awk '{printf "%.3f", $2}'`
    loudf=`grep -m 1 -w LOUDESTFREQ_dt4 ${outdir}/${channel}_${tcenter}.txt  | head -1 | awk '{printf "%.1f", $2}'`

    # get channel description
    if [ -e ${OMICRON_PARAMETERS}/channels.txt ]; then
	
	chan_des=` grep -w ${channel} ${OMICRON_PARAMETERS}/channels.txt`
	if [ "$chan_des" = "" ]; then
	    description="No description for this channel"
	else
	    description=`echo $chan_des | sed -e 's/^\w*\ *//'`
	fi
    else
	description="No description for this channel"
    fi

    # thumbnails
    convert -density 100 -thumbnail 320  ${outdir}/${channel}_${tcenter}_dt4.gif ${outdir}/${channel}_${tcenter}_dt4_th.gif
    convert -density 100 -thumbnail 320  ${outdir}/${channel}_${tcenter}_dt16.gif ${outdir}/${channel}_${tcenter}_dt16_th.gif
    convert -density 100 -thumbnail 320  ${outdir}/${channel}_${tcenter}_dt64.gif ${outdir}/${channel}_${tcenter}_dt64_th.gif

    # fill html report
    sed -i "/<\!-- channels -->/i\<tr><td><table><tr><td colspan=\"3\"><h2>${channel}</h2><h3>${description}</h3>Loudest tile: GPS=${loudgps}, f=${loudf}Hz, SNR=${loudsnr}</td></tr>" $template
    sed -i "/<\!-- channels -->/i\<tr><td><a href=\"./${channel}_${tcenter}_dt4.gif\"><img src=\"./${channel}_${tcenter}_dt4_th.gif\"></a></td><td><a href=\"./${channel}_${tcenter}_dt16.gif\"><img src=\"./${channel}_${tcenter}_dt16_th.gif\"></a></td><td><a href=\"./${channel}_${tcenter}_dt64.gif\"><img src=\"./${channel}_${tcenter}_dt64_th.gif\"></a></td></tr></table></td></tr>" $template

done


##### fill static variables
sed -e "s|\[GPS_CENTER\]|${tcenter}|g" \
    -e "s|\[DATEUTC\]|${dateUTC}|g" \
    $template > ${outdir}/index.html

##### cleaning
rm -f $template

##### finish log
currentdate=`date -u`

##### start log
echo "" >> ${outdir}/log.txt
echo "Number of scanned channels:                               $naux" >> ${outdir}/log.txt
echo "Number of scanned channels with triggers:                 $naux_with_triggers" >> ${outdir}/log.txt
echo "Number of scanned channels with triggers above threshold: $naux_above_thr" >> ${outdir}/log.txt
echo "Current date: $currentdate" >> ${outdir}/log.txt
echo "**********************************" >> ${outdir}/log.txt
echo "**  End Omiscan of $tcenter" >> ${outdir}/log.txt
echo "**********************************" >> ${outdir}/log.txt


exit 0