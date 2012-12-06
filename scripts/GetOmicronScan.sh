#!/bin/bash
#
# GetOmicronWeb.sh
#
# generate omicron web page
#
# Author: Florent Robinet
# florent.robinet@lal.in2p3.fr

##### default options
main="h_4096Hz" # channel name
outdir=`pwd` # output directory
snrmin=8

if [ ! test -z "$OMICRON_HTML" ]; then
  echo "The Omicron environment is not sourced"
fi


# move to the script directory
cd `dirname $0`
. $GWOLLUMROOT/local/environment.sh

printhelp(){
    echo ""
    echo "Usage:"
    echo "GetOmicronScan.sh -m [MAIN_CHANNEL] [GPS_CENTER]"
    echo ""
    echo "TRIGGER SELECTION OPTIONS"
    echo "  -s  [SNRMIN]        print aux. scan only if SNR > [SNRMIN]"
    echo "                      Default = 8"
    echo ""
    echo "OUTPUT CONTROL"
    echo "  -m  [MAIN_CHANNEL]  main channel: always plotted"
    echo "                      Default = h_4096Hz"
    echo "  -d  [OUTDIR]        _full_ pOMICRON_HTMLath to output directory"
    echo "                      Default = current directory"
    echo ""
    echo "  -h                  print this help"
    echo ""
} 


##### needs argument
if [ $# -lt 1 ]; then
    printhelp
    exit 1
fi

##### read options
while getopts ":m:s:d:h" opt; do
    case $opt in
	m)
	    main="$OPTARG"
	    ;;
	s)
	    snrmin=`echo $OPTARG | awk '{print int($1)}'`
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
	    echo "type  'GetOmicronScan.sh -h'  for help"
	    exit 1
	    ;;
    esac
done

##### gps interval
shift $(($OPTIND - 1))
tcenter=`echo $1 | awk '{printf "%.3f", $1}'`
tcenter_=`echo $tcenter | awk '{print int($1)}'`

##### check main channel is available
if ! echo "$OMICRON_ONLINE_CHANNELS $OMICRON_CHANNELS" | grep -q "$main"; then
    echo "Invalid option: channel '${main}' is not available"
    echo "type  'GetOmicronScan.sh -h'  for help"
    exit 1
fi

##### check timing
if [ $tcenter_ -lt 700000000 ]; then
    echo "Invalid option: '$tcenter' is not a reasonable central time"
    echo "type  'GetOmicronScan.sh -h'  for help"
    exit 1
fi

##### elaborate timing
dateUTC=`tconvert -f "%A the %dth, %B %Y %H:%M:%S" ${tcenter_}`
tstart=$(( $tcenter_ - 100 ))
tstop=$(( $tcenter_ + 100 ))

##### check outdir
if [ ! -d $outdir ] ; then
    echo "Invalid option: the output directory $outdir cannot be found"
    echo "type  'GetOmicronScan.sh -h'  for help"
    exit 1
fi

##### build directory
outdir="${outdir}/${tcenter}"
mkdir -p ${outdir}; rm -fr ${outdir}/*

##### html template and web material
template=${outdir}/index.template
cp ${OMICRON_HTML}/template.omicronscan.html $template
cp ${OMICRON_HTML}/comparison_mode.html ${outdir}/
ln -sf ${GWOLLUM_DOC}/style.css ${outdir}/style.css
ln -sf ${GWOLLUM_DOC}/Pics/omicronlogo_xxl.gif ${outdir}/omicronlogo_xxl.gif
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
triggers=`GetTriggerFileList.sh -c${main} $tstart $tstop | grep "FILELIST" | sed 's|FILELIST ||g'`
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
if [ -e ${OMICRON_TRIGGERS}/${main}/description.txt ]; then
    description=`cat ${OMICRON_TRIGGERS}/${main}/description.txt`
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
    triggers=`GetTriggerFileList.sh -c${channel} $tstart $tstop | grep "FILELIST" | sed 's|FILELIST ||g'`
    if [ "$triggers" = "" ]; then 
	echo "  no trigger ---> Do not plot" >> ${outdir}/log.txt
	continue
    fi
    let "naux_with_triggers+=1"
    eventmap.exe ${outdir} "${triggers}" $tcenter >> ${outdir}/log.txt 2>&1

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
    if [ -e ${OMICRON_TRIGGERS}/${channel}/description.txt ]; then
	description=`cat ${OMICRON_TRIGGERS}/${channel}/description.txt`
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