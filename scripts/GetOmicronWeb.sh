#!/bin/bash
#
# GetOmicronWeb.sh
#
# generate omicron web page
#
# Author: Florent Robinet
# florent.robinet@lal.in2p3.fr

##### default options
channel="h_4096Hz" # channel name
outdir=`pwd` # output directory
next=""
previous=""
up=""
hlinks=""

# move to the script directory
cd `dirname $0`
. $GWOLLUMROOT/local/environment.sh

printhelp(){
    echo ""
    echo "Usage:"
    echo "GetOmicronWeb.sh -c[CHANNEL_NAME] [GPS_START] [GPS_STOP]"
    echo ""
    echo "Example: GetOmicronWeb.sh -ch_4096Hz 934228815 934232415"
    echo ""
    echo "TRIGGER SELECTION OPTIONS"
    echo "  -c  [CHANNEL_NAME]  triggers from channel [CHANNEL_NAME]"
    echo "                      Default = 'h_4096Hz'"
    echo ""
    echo "OUTPUT CONTROL"
    echo "  -d  [OUTDIR]        _full_ path to output directory"
    echo "                      Default = current directory"
    echo "  -n  [HTMLPATH]      add -next- button pointing to [HTMLPATH]"
    echo "  -p  [HTMLPATH]      add -previous- button pointing to [HTMLPATH]"
    echo "  -u  [HTMLPATH]      add -up- button pointing to [HTMLPATH]"
    echo "  -t  [DAYPATH]/       add hour links found in [DAYPATH]/"
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
while getopts ":c:d:n:p:u:t:h" opt; do
    case $opt in
	c)
	    channel="$OPTARG"
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
	t)
	    hlinks="$OPTARG"
	    ;;
	h)
	    printhelp
	    exit 0
	    ;;
	\?)
	    echo "Invalid option: -$OPTARG"
	    echo "type  'GetOmicronWeb.sh -h'  for help"
	    exit 1
	    ;;
    esac
done

##### gps interval
shift $(($OPTIND - 1))
tmin=`echo $1 | awk '{print int($1)}'`
tmax=`echo $2 | awk '{print int($1)}'`


##### check channel is available
if ! echo "$OMICRON_ONLINE_CHANNELS $OMICRON_CHANNELS" | grep -q "$channel"; then
    echo "Invalid option: channel '${channel}' is not available"
    echo "type  'GetOmicronWeb.sh -h'  for help"
    exit 1
fi

##### check timing
if [ $tmin -lt 700000000 ]; then
    echo "Invalid option: '$tmin' is not a reasonable starting time"
    echo "type  'GetOmicronWeb.sh -h'  for help"
    exit 1
fi
if [ $tmax -lt 700000000 ]; then
    echo "Invalid option: '$tmax' is not a reasonable stop time"
    echo "type  'GetOmicronWeb.sh -h'  for help"
    exit 1
fi
if [ $tmax -le $tmin ]; then
    echo "Invalid option: the time interval '$tmin-$tmax' is not reasonable"
    echo "type  'GetOmicronWeb.sh -h'  for help"
    exit 1
fi

##### check outdir
if [ ! -d $outdir ] ; then
    echo "Invalid option: the output directory $outdir cannot be found"
    echo "type  'GetOmicronWeb.sh -h'  for help"
    exit 1
fi

##### build directory
mkdir -p ${outdir}/${channel}

##### make plots
./GetOmicronPlots.sh -c${channel} -d${outdir}/${channel} $tmin $tmax > ${outdir}/${channel}/plot.log.txt 2>&1

##### clean
cd ${outdir}/${channel}
rm -f ./map.gif ./rate.gif ./snrfreq.gif ./loudest.gif ./snr.gif ./freq.gif ./snrtime.gif ./info.txt

##### no data
if grep -q "triggers are not available" ./plot.log.txt; then
    ln -s ${GWOLLUM_DOC}/Pics/WebReport/nodata1.gif ./map.gif
    ln -s ${GWOLLUM_DOC}/Pics/WebReport/nodata1.gif ./rate.gif
    ln -s ${GWOLLUM_DOC}/Pics/WebReport/nodata1.gif ./snrfreq.gif
    ln -s ${GWOLLUM_DOC}/Pics/WebReport/nodata1.gif ./snr.gif
    ln -s ${GWOLLUM_DOC}/Pics/WebReport/nodata1.gif ./freq.gif
    ln -s ${GWOLLUM_DOC}/Pics/WebReport/nodata1.gif ./snrtime.gif
    ln -s ${GWOLLUM_DOC}/Pics/WebReport/nodata1.gif ./loudest.gif
    ln -s ${GWOLLUM_DOC}/Pics/WebReport/nodata1.gif ./info.txt

##### generic names for plots
else
    for file in ./*.gif ./*.txt; do
	if echo $file | grep -q "map"; then mv $file ./map.gif;
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
cd `dirname $0`

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
cp ${GWOLLUM_HTML}/template.omicronmonitor.html $template
ln -sf ${GWOLLUM_DOC}/style.css ${outdir}/style.css
ln -sf ${GWOLLUM_DOC}/Pics/gwollum_logo_min_trans.gif ${outdir}/icon.gif
ln -sf ${GWOLLUM_DOC}/Pics/gwollum_logo.gif ${outdir}/logo.gif
ln -sf ${GWOLLUM_DOC}/Pics/omicronlogo_l.gif ${outdir}/omicronlogo.gif
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



exit 0