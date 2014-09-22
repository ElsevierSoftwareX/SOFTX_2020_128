#!/bin/bash
#
# GetOmicronWeb.sh
#
# generate omicron web page
#
# Author: Florent Robinet
# florent.robinet@lal.in2p3.fr

here=`pwd`

echo "This tool is not available for now."
exit 0

printhelp(){
    echo ""
    echo "Usage:"
    echo "GetOmicronWeb -c [CHANNEL_NAME] -s [GPS_START] -e [GPS_END]"
    echo " |__ web report for 'standard' Omicron triggers of a given channel and between 2 GPS times"
    echo ""
    echo "GetOmicronWeb -t \"[FILE_PATTERN]\" -s [GPS_START] -e [GPS_END]"
    echo " |__ web report for triggers saved in a list of root files"
    echo "     designated by a file pattern and between 2 GPS times"
    echo ""
    echo "For standard triggers, [GPS_START] and [GPS_END] should belong to the same run"
    echo ""
    echo "OPTIONS:"
    echo "  -c  [CHANNEL_NAME]    channel name. Default = V1:h_4096Hz"
    echo "  -t  \"[FILE_PATTERN]\"  the user provides his own trigger files"
    echo "                        The file pattern is [TRIGGER_PATTERN]"
    echo "  -s  [GPS_START]       starting GPS time. Required!"
    echo "  -e  [GPS_END]         stopping GPS time. Required!"
    echo ""
    echo "OUTPUT OPTIONS"
    echo "  -d  [OUTDIR]          output directory where to save plots"
    echo "                        Default = current directory"
    echo "  -n  [HTMLPATH]        add -next- button pointing to [HTMLPATH]"
    echo "  -p  [HTMLPATH]        add -previous- button pointing to [HTMLPATH]"
    echo "  -u  [HTMLPATH]        add -up- button pointing to [HTMLPATH]"
    echo "  -o  [DAYPATH]/        add hour links found in [DAYPATH]/"
    echo "  -w                    full web architecture flag"
    echo ""
    echo "  -l                    list available channels"
    echo "  -h                    print this help"
    echo ""
    echo "Author: Florent Robinet (LAL - Orsay): robinet@lal.in2p3.fr"
    echo ""
} 

##### Check the environment
if [ -z "$OMICRONROOT" ]; then
    echo "Error: The Omicron environment is not set"
    exit 1
fi

##### default options
channel="V1:h_4096Hz"  # channel name
outdir=`pwd`        # output directory
triggerfiles="NONE" # user trigger files
tmin=0
tmax=0
next=""
previous=""
up=""
hlinks=""
fullweb="no"
triggerfiles="NONE" # user trigger files

##### read options
while getopts ":c:t:s:e:d:n:p:u:o:wlh" opt; do
    case $opt in
	c)
	    channel="$OPTARG"
	    ;;
	t)
	    triggerfiles="$OPTARG"
	    ;;
	s)
	    tmin=`echo $OPTARG | awk '{print int($1)}'`
	    ;;
	e)
	    tmax=`echo $OPTARG | awk '{print int($1)}'`
	    ;;
	d)
	    outdir="$OPTARG"
	    ;;
	n)
	    next="$OPTARG"
	    ;;
	p)
	    previous="$OPTARG"
	    ;;
	u)
	    up="$OPTARG"
	    ;;
	o)
	    hlinks="$OPTARG"
	    ;;
	w)
	    fullweb="yes"
	    ;;
	l)
	    for run in $RUN_NAMES; do
		echo ""
		GetOmicronChannels.sh -v -r $run
	    done
	    exit 0
	    ;;
	h)
	    printhelp
	    exit 0
	    ;;
	\?)
	    echo "`basename $0`: Invalid option: -$OPTARG"
	    echo "type  'GetOmicronWeb -h'  for help"
	    exit 1
	    ;;
    esac
done
OPTIND=1

##### check outdir
if [ ! -d $outdir ] ; then
    echo "`basename $0`: the output directory $outdir cannot be found"
    echo "type  'GetOmicronWeb -h'  for help"
    exit 1
fi

##### check timing
if [ $tmin -eq 0 ]; then
    echo "`basename $0`: A starting time is required"
    echo "type  'GetOmicronWeb -h'  for help"
    exit 1
fi
if [ $tmax -eq 0 ]; then
    echo "`basename $0`: A stopping time is required"
    echo "type  'GetOmicronWeb -h'  for help"
    exit 1
fi
if [ $tmax -le $tmin ]; then
    echo "`basename $0`: the time interval '$tmin-$tmax' is not reasonable"
    echo "type  'GetOmicronWeb -h'  for help"
    exit 1
fi

##### build directory
mkdir -p ${outdir}/${channel}

##### check lock
if [ -e ${outdir}/${channel}/lock ]; then
    echo "GetOmicronWeb: the directory ${outdir}/${channel} is locked"
    exit 1
fi
touch ${outdir}/${channel}/lock

##### starting date
date -u > ${outdir}/${channel}/plot.log.txt 2>&1

##### make plots
if [ "$triggerfiles" = "NONE" ]; then
    GetOmicronPlots.sh -c ${channel} -d ${outdir}/${channel} -s $tmin -e $tmax >> ${outdir}/${channel}/plot.log.txt 2>&1
else
    GetOmicronPlots.sh -t "$triggerfiles" -d ${outdir}/${channel} -s $tmin -e $tmax >> ${outdir}/${channel}/plot.log.txt 2>&1
fi

##### ending date
date -u >> ${outdir}/${channel}/plot.log.txt 2>&1

##### clean
cd ${outdir}/${channel}
rm -f ./map.gif ./seg.gif ./rate.gif ./snrfreq.gif ./loudest.gif ./snr.gif ./freq.gif ./snrtime.gif ./info.txt

##### no data
if grep -q "triggers are not available" ./plot.log.txt; then

    if [ "$fullweb" = "no" ]; then
	cp -f ${GWOLLUM_DOC}/Pics/WebReport/nodata1.gif ./map.gif
	cp -f ${GWOLLUM_DOC}/Pics/WebReport/nodata1.gif ./rate.gif
	cp -f ${GWOLLUM_DOC}/Pics/WebReport/nodata1.gif ./snrfreq.gif
	cp -f ${GWOLLUM_DOC}/Pics/WebReport/nodata1.gif ./snr.gif
	cp -f ${GWOLLUM_DOC}/Pics/WebReport/nodata1.gif ./freq.gif
	cp -f ${GWOLLUM_DOC}/Pics/WebReport/nodata1.gif ./snrtime.gif
	cp -f ${GWOLLUM_DOC}/Pics/WebReport/nodata1.gif ./loudest.gif
	cp -f ${GWOLLUM_DOC}/Pics/WebReport/nodata1.gif ./info.txt
    else
	if [ -e ${outdir}/../../../../nodata.gif ]; then
	    ln -s ../../../../nodata.gif ./map.gif
	    ln -s ../../../../nodata.gif ./rate.gif
	    ln -s ../../../../nodata.gif ./snrfreq.gif
	    ln -s ../../../../nodata.gif ./snr.gif
	    ln -s ../../../../nodata.gif ./freq.gif
	    ln -s ../../../../nodata.gif ./snrtime.gif
	    ln -s ../../../../nodata.gif ./loudest.gif
	    ln -s ../../../../nodata.gif ./info.txt
	elif [ -e ${outdir}/../../../../../nodata.gif ]; then
	    ln -s ../../../../../nodata.gif ./map.gif
	    ln -s ../../../../../nodata.gif ./rate.gif
	    ln -s ../../../../../nodata.gif ./snrfreq.gif
	    ln -s ../../../../../nodata.gif ./snr.gif
	    ln -s ../../../../../nodata.gif ./freq.gif
	    ln -s ../../../../../nodata.gif ./snrtime.gif
	    ln -s ../../../../../nodata.gif ./loudest.gif
	    ln -s ../../../../../nodata.gif ./info.txt
	else
	    cp -f ${GWOLLUM_DOC}/Pics/WebReport/nodata1.gif ./map.gif
	    cp -f ${GWOLLUM_DOC}/Pics/WebReport/nodata1.gif ./rate.gif
	    cp -f ${GWOLLUM_DOC}/Pics/WebReport/nodata1.gif ./snrfreq.gif
	    cp -f ${GWOLLUM_DOC}/Pics/WebReport/nodata1.gif ./snr.gif
	    cp -f ${GWOLLUM_DOC}/Pics/WebReport/nodata1.gif ./freq.gif
	    cp -f ${GWOLLUM_DOC}/Pics/WebReport/nodata1.gif ./snrtime.gif
	    cp -f ${GWOLLUM_DOC}/Pics/WebReport/nodata1.gif ./loudest.gif
	    cp -f ${GWOLLUM_DOC}/Pics/WebReport/nodata1.gif ./info.txt
	fi
    fi

##### generic names for plots
else
    for file in ./*.gif ./*.txt; do
	if echo $file | grep -q "map"; then mv $file ./map.gif;
	elif echo $file | grep -q "segments"; then mv $file ./seg.gif;
	elif echo $file | grep -q "rate"; then mv $file ./rate.gif;
	elif echo $file | grep -q "snrfreq"; then mv $file ./snrfreq.gif;
	elif echo $file | grep -q "snrtime"; then mv $file ./snrtime.gif;
	elif echo $file | grep -q "snr"; then mv $file ./snr.gif;
	elif echo $file | grep -q "freq"; then mv $file ./freq.gif;
	elif echo $file | grep -q "info"; then mv $file ./info.txt;
	elif echo $file | grep -q "loudest"; then mv $file ./loudest.gif;
	else continue; fi
    done
fi
cd $here

##### elaborate timing
datestart=`tconvert -f "%Y-%m-%d %H:%M" ${tmin}`
dateend=`tconvert -f "%Y-%m-%d %H:%M" ${tmax}`

##### trigger information
if grep -q "triggers are not available" ${outdir}/${channel}/plot.log.txt; then
    ntrigger=0; ncluster=0; ncluster0=0; ncluster1=0; ncluster1=0; ncluster2=0; ncluster3=0;
    freqmin=0; freqmax=0; snrmin=0; snrmax=0; deltat=0; livetime=0;
else
    ntrigger=`grep -m 1 -w NTRIGGERS ${outdir}/${channel}/info.txt | head -1 | awk '{print int($2)}'`
    ncluster=`grep -m 1 -w NCLUSTERS ${outdir}/${channel}/info.txt | head -1 | awk '{print int($2)}'`
    ncluster0=`grep -m 1 -w NCLUSTERS0 ${outdir}/${channel}/info.txt | head -1 | awk '{print int($2)}'`
    ncluster1=`grep -m 1 -w NCLUSTERS1 ${outdir}/${channel}/info.txt | head -1 | awk '{print int($2)}'`
    ncluster2=`grep -m 1 -w NCLUSTERS2 ${outdir}/${channel}/info.txt | head -1 | awk '{print int($2)}'`
    ncluster3=`grep -m 1 -w NCLUSTERS3 ${outdir}/${channel}/info.txt | head -1 | awk '{print int($2)}'`
    rate0=`grep -m 1 -w RATE0 ${outdir}/${channel}/info.txt | head -1 | awk '{printf "%.3f", $2}'`
    rate1=`grep -m 1 -w RATE1 ${outdir}/${channel}/info.txt | head -1 | awk '{printf "%.3f", $2}'`
    rate2=`grep -m 1 -w RATE2 ${outdir}/${channel}/info.txt | head -1 | awk '{printf "%.3f", $2}'`
    rate3=`grep -m 1 -w RATE3 ${outdir}/${channel}/info.txt | head -1 | awk '{printf "%.3f", $2}'`
    freqmin=`grep -m 1 -w FREQMIN ${outdir}/${channel}/info.txt | head -1 | awk '{print int($2)}'`
    freqmax=`grep -m 1 -w FREQMAX ${outdir}/${channel}/info.txt | head -1 | awk '{print int($2)}'`
    snrmin=`grep -m 1 -w SNRMIN ${outdir}/${channel}/info.txt | head -1 | awk '{printf "%.1f", $2}'`
    snrmax=`grep -m 1 -w SNRMAX ${outdir}/${channel}/info.txt | head -1 | awk '{printf "%.1f", $2}'`
    timesnrmax=`grep -m 1 -w TIMESNRMAX ${outdir}/${channel}/info.txt | head -1 | awk '{printf "%.3f", $2}'`
    freqsnrmax=`grep -m 1 -w FREQSNRMAX ${outdir}/${channel}/info.txt | head -1 | awk '{printf "%.1f", $2}'`
    deltat=`grep -m 1 -w DELTAT ${outdir}/${channel}/info.txt | head -1 | awk '{printf "%.3f", $2}'`
    livetime=`grep -m 1 -w LIVETIME ${outdir}/${channel}/info.txt | head -1 | awk '{print int($2)}'`
fi

##### web page materials
template=${outdir}/${channel}.template
cp ${OMICRON_HTML}/template/template.omicronmonitor.html $template
if [ "$fullweb" = "no" ]; then
    cp -f ${GWOLLUM_DOC}/style.css ${outdir}
    cp -f ${GWOLLUM_DOC}/Pics/gwollum_logo_min_trans.gif ${outdir}/icon.gif
    cp -f ${GWOLLUM_DOC}/Pics/gwollum_logo.gif ${outdir}/logo.gif
    cp -f ${OMICRON_HTML}/pics/omicronlogo_l.gif ${outdir}/omicronlogo.gif
else
    if [ -e ${outdir}/../../../style.css ]; then
	ln -sf ../../../style.css ${outdir}/style.css
	ln -sf ../../../omicronlogo_l.gif ${outdir}/omicronlogo.gif
	ln -sf ../../../icon.gif ${outdir}/icon.gif
	ln -sf ../../../omicronlogo_s.gif ${outdir}/logo.gif
    elif [ -e ${outdir}/../../../../style.css ]; then
	ln -sf ../../../../style.css ${outdir}/style.css
	ln -sf ../../../../omicronlogo_l.gif ${outdir}/omicronlogo.gif
	ln -sf ../../../../icon.gif ${outdir}/icon.gif
	ln -sf ../../../../omicronlogo_s.gif ${outdir}/logo.gif
    else
	cp -f ${GWOLLUM_DOC}/style.css ${outdir}
	cp -f ${GWOLLUM_DOC}/Pics/gwollum_logo_min_trans.gif ${outdir}/icon.gif
	cp -f ${OMICRON_HTML}/pics/omicronlogo_s.gif ${outdir}/logo.gif
	cp -f ${OMICRON_HTML}/pics/omicronlogo_l.gif ${outdir}/omicronlogo.gif
    fi
fi
currentdate=`date -u`

##### time navigation
if [ ! "${previous}${next}${up}${hlinks}" = "" ]; then
    sed -i "/<\!-- time navigation -->/i\<table class=\"omicrontiming\">" $template

    # previous/up/next
    if [ ! "${previous}${next}${up}" = "" ]; then
	sed -i "/<\!-- time navigation -->/i\<tr><td>" $template
	if [ ! "${previous}" = "" ]; then
    	    sed -i "/<\!-- time navigation -->/i\<a href=\"${previous}\">previous</a>" $template
	fi
	sed -i "/<\!-- time navigation -->/i\</td><td>" $template
	if [ ! "${up}" = "" ]; then
    	    sed -i "/<\!-- time navigation -->/i\<a href=\"${up}\">up</a>" $template
	fi
	sed -i "/<\!-- time navigation -->/i\</td><td>" $template
	if [ ! "${next}" = "" ]; then
    	    sed -i "/<\!-- time navigation -->/i\<a href=\"${next}\">next</a>" $template
	fi
	sed -i "/<\!-- time navigation -->/i\</td></tr>" $template
    fi

    # previous/up/next
    if [ ! "${hlinks}" = "" ]; then
	sed -i "/<\!-- time navigation -->/i\<tr><td><a href=\"${hlinks}00\">00h</a></td><td><a href=\"${hlinks}01\">01h</a></td><td><a href=\"${hlinks}02\">02h</a></td></tr>" $template
	sed -i "/<\!-- time navigation -->/i\<tr><td><a href=\"${hlinks}03\">03h</a></td><td><a href=\"${hlinks}04\">04h</a></td><td><a href=\"${hlinks}05\">05h</a></td></tr>" $template
	sed -i "/<\!-- time navigation -->/i\<tr><td><a href=\"${hlinks}06\">06h</a></td><td><a href=\"${hlinks}07\">07h</a></td><td><a href=\"${hlinks}08\">08h</a></td></tr>" $template
	sed -i "/<\!-- time navigation -->/i\<tr><td><a href=\"${hlinks}09\">09h</a></td><td><a href=\"${hlinks}10\">10h</a></td><td><a href=\"${hlinks}11\">11h</a></td></tr>" $template
	sed -i "/<\!-- time navigation -->/i\<tr><td><a href=\"${hlinks}12\">12h</a></td><td><a href=\"${hlinks}13\">13h</a></td><td><a href=\"${hlinks}14\">14h</a></td></tr>" $template
	sed -i "/<\!-- time navigation -->/i\<tr><td><a href=\"${hlinks}15\">15h</a></td><td><a href=\"${hlinks}16\">16h</a></td><td><a href=\"${hlinks}17\">17h</a></td></tr>" $template
	sed -i "/<\!-- time navigation -->/i\<tr><td><a href=\"${hlinks}18\">18h</a></td><td><a href=\"${hlinks}19\">19h</a></td><td><a href=\"${hlinks}20\">20h</a></td></tr>" $template
	sed -i "/<\!-- time navigation -->/i\<tr><td><a href=\"${hlinks}21\">21h</a></td><td><a href=\"${hlinks}22\">22h</a></td><td><a href=\"${hlinks}23\">23h</a></td></tr>" $template
    fi
    sed -i "/<\!-- time navigation -->/i\</table>"  $template
fi

##### remove the segment plot if no data
if grep -q "triggers are not available" ${outdir}/${channel}/plot.log.txt; then
    sed -i '/<!-- to_remove -->/d' $template
fi

##### fill static variables
sed -e "s|\[DATESTART\]|${datestart}|g" \
    -e "s|\[DATEEND\]|${dateend}|g" \
    -e "s|\[CHANNEL\]|${channel}|g" \
    -e "s|\[NTRIGGERS\]|${ntrigger}|g" \
    -e "s|\[NCLUSTERS\]|${ncluster}|g" \
    -e "s|\[NCLUSTERS0\]|${ncluster0}|g" \
    -e "s|\[NCLUSTERS1\]|${ncluster1}|g" \
    -e "s|\[NCLUSTERS2\]|${ncluster2}|g" \
    -e "s|\[NCLUSTERS3\]|${ncluster3}|g" \
    -e "s|\[RATE0\]|${rate0}|g" \
    -e "s|\[RATE1\]|${rate1}|g" \
    -e "s|\[RATE2\]|${rate2}|g" \
    -e "s|\[RATE3\]|${rate3}|g" \
    -e "s|\[SNR0\]|5.5|g" \
    -e "s|\[SNR1\]|8|g" \
    -e "s|\[SNR2\]|10|g" \
    -e "s|\[SNR3\]|20|g" \
    -e "s|\[FREQMIN\]|${freqmin}|g" \
    -e "s|\[FREQMAX\]|${freqmax}|g" \
    -e "s|\[SNRMIN\]|${snrmin}|g" \
    -e "s|\[SNRMAX\]|${snrmax}|g" \
    -e "s|\[TIMESNRMAX\]|${timesnrmax}|g" \
    -e "s|\[FREQSNRMAX\]|${freqsnrmax}|g" \
    -e "s|\[DELTAT\]|${deltat}|g" \
    -e "s|\[LIVETIME\]|${livetime}|g" \
    -e "s|\[CREATIONDATE\]|${currentdate}|g" \
    $template > ${outdir}/${channel}.html

##### cleaning
rm -f $template
rm -f ${outdir}/${channel}/lock


exit 0
