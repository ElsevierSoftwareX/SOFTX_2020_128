#!/bin/bash
#
# PrepareOmicronWeb.sh
#
# prepare omicron web directory
#
# Author: Florent Robinet
# florent.robinet@lal.in2p3.fr

webdir="./no.chance.this.exists"

printhelp(){
    echo ""
    echo "Usage:"
    echo "PrepareOmicronWeb -d [WEB_DIRECTORY]"
    echo ""
    echo "WEB FORMAT"
    echo "  -d    [WEB_DIRECTORY]   path to web directory"
    echo ""
    echo "OUTPUT CONTROL"
    echo "  -h                      print this help"
    echo ""
} 

##### Check the environment
if [[ -z "$OMICRONROOT" ]]; then
    echo "Error: The Omicron environment is not set"
    exit 1
fi

##### read options
while getopts ":c:d:h" opt; do
    case $opt in
	d)
	    webdir="$OPTARG"
	    ;;
	h)
	    printhelp
	    exit 0
	    ;;
	\?)
	    echo "Invalid option: -$OPTARG"
	    echo "type  'PrepareOmicronWeb -h'  for help"
	    exit 1
	    ;;
    esac
done


##### check web directory
if [ ! -d $webdir ]; then
    echo "Invalid option: directory '$webdir' does not exist"
    echo "type  'PrepareOmicronWeb -h'  for help"
    exit 1
fi

cd $webdir

mkdir -p ./latestday
mkdir -p ./latesthour

# soft links
ln -sf $OMICRON_HTML/monitoringweb/getomicronpage.php
ln -sf $OMICRON_HTML/monitoringweb/getomiscan.php
ln -sf $OMICRON_HTML/monitoringweb/logo.html
ln -sf $OMICRON_HTML/monitoringweb/omiscan.html
ln -sf $OMICRON_HTML/monitoringweb/index.html
ln -sf $GWOLLUM_DOC/Pics/gwollum_logo_min_trans.gif icon.gif
ln -sf $GWOLLUM_DOC/Pics/WebReport/nodata1.gif nodata.gif
ln -sf $OMICRON_HTML/pics/omicronlogo_l.gif
ln -sf $OMICRON_HTML/pics/omicronlogo_s.gif
ln -sf $OMICRON_HTML/pics/omicronlogo_xl.gif
ln -sf $GWOLLUM_DOC/style.css
ln -sf $GWOLLUM_HTML/php/tconvert.php

cp -f $OMICRON_HTML/monitoringweb/header.html .
cp -f $OMICRON_HTML/monitoringweb/calendar.html .

exit 0