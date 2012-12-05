#!/bin/sh
#
# configure.sh
#
# Define the Omicron environment
#
# Author: Florent Robinet
# florent.robinet@lal.in2p3.fr

echo ""
echo ""
echo "   !!! WELCOME TO OMICRON's WORLD !!!"
echo ""
echo "Let's check your environment..."
echo ""

# omicron
omicron=`pwd`


######################################################################
##############################  GWOLLUM  #############################
######################################################################
echo "Do you have GWOLLUM installed in your system?"
if test $GWOLLUM_LOCATION; then
    echo "YES!  " $GWOLLUM_LOCATION
    echo ""
    gwollumpath=$GWOLLUM_LOCATION
else
    echo "NO!"
    echo "You need to install GWOLLUM in your system:"
    echo ""
    echo "https://wwwcascina.virgo.infn.it/DataAnalysis/GWOLLUM/doc/index.html"
    echo ""
    echo "Or, give the full path to your GWOLLUM installation:"
    read gwollumpath
    echo ""
fi

# check GWOLLUM library
if [ ! -e ${gwollumpath}/lib/libDataQuality.so ]; then
    echo "${gwollumpath}/lib/libDataQuality.so cannot be found, try again..."
    exit 1
fi

# CSH
echo "Create CSH environment..."
echo '# Source this file to access Omicron' > ./etc/omicron-user-env.csh
echo 'setenv OMICRON_LOCATION "'${omicron}'"' >> ./etc/omicron-user-env.csh
echo 'source '${gwollumpath}'/etc/gwollum-user-env.csh' >> ./etc/omicron-user-env.csh
echo 'setenv OMICRON_HTML     "${OMICRON_LOCATION}/html"' >> ./etc/omicron-user-env.csh
echo 'setenv OMICRON_DOC      "${OMICRON_LOCATION}/doc"' >> ./etc/omicron-user-env.csh
echo 'setenv OMICRON_SCRIPTS  "${OMICRON_LOCATION}/scripts"' >> ./etc/omicron-user-env.csh
echo 'setenv OMICRON_BIN      "${OMICRON_LOCATION}/bin"' >> ./etc/omicron-user-env.csh
echo 'setenv OMICRON_LIB      "${OMICRON_LOCATION}/lib"' >> ./etc/omicron-user-env.csh
echo 'setenv PATH "${OMICRON_BIN}:${OMICRON_SCRIPTS}:${PATH}"' >> ./etc/omicron-user-env.csh
echo 'if ( $?LD_LIBRARY_PATH ) then' >> ./etc/omicron-user-env.csh
echo '	setenv LD_LIBRARY_PATH "${OMICRON_LIB}:${LD_LIBRARY_PATH}"' >> ./etc/omicron-user-env.csh
echo 'else' >> ./etc/omicron-user-env.csh
echo '	setenv LD_LIBRARY_PATH "${OMICRON_LIB}"' >> ./etc/omicron-user-env.csh
echo 'endif' >> ./etc/omicron-user-env.csh
echo '' >> ./etc/omicron-user-env.csh

# SH
echo "Create SH environment..."
echo '# Source this file to access OMICRON' > ./etc/omicron-user-env.sh
echo 'OMICRON_LOCATION="'${omicron}'"' >> ./etc/omicron-user-env.sh
echo 'source '${gwollumpath}'/etc/gwollum-user-env.sh' >> ./etc/omicron-user-env.sh
echo 'OMICRON_HTML="${OMICRON_LOCATION}/html"' >> ./etc/omicron-user-env.sh
echo 'OMICRON_DOC="${OMICRON_LOCATION}/doc"' >> ./etc/omicron-user-env.sh
echo 'OMICRON_SCRIPTS="${OMICRON_LOCATION}/scripts"' >> ./etc/omicron-user-env.sh
echo 'OMICRON_BIN="${OMICRON_LOCATION}/bin"' >> ./etc/omicron-user-env.sh
echo 'OMICRON_LIB="${OMICRON_LOCATION}/lib"' >> ./etc/omicron-user-env.sh
echo 'PATH="${OMICRON_BIN}:${OMICRON_SCRIPTS}:${PATH}"' >> ./etc/omicron-user-env.sh
echo 'LD_LIBRARY_PATH="${OMICRON_LIB}:${LD_LIBRARY_PATH}"' >> ./etc/omicron-user-env.sh
echo 'export OMICRON_LOCATION' >> ./etc/omicron-user-env.sh
echo 'export OMICRON_HTML' >> ./etc/omicron-user-env.sh
echo 'export OMICRON_DOC' >> ./etc/omicron-user-env.sh
echo 'export OMICRON_SCRIPTS' >> ./etc/omicron-user-env.sh
echo 'export OMICRON_BIN' >> ./etc/omicron-user-env.sh
echo 'export OMICRON_LIB' >> ./etc/omicron-user-env.sh
echo 'export PATH' >> ./etc/omicron-user-env.sh
echo 'export LD_LIBRARY_PATH' >> ./etc/omicron-user-env.sh

echo ""
echo "******************************************************"
echo ""
echo "Now you can update your environment by sourcing"
echo "one of the files in ${omicron}/etc"
echo ""
echo "Then you can type 'make'"
echo ""
echo "******************************************************"


exit 0


