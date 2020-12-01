# source this file to access Omicron

##### PKG_CONFIG_PATH
PKG_CONFIG_PATH=`echo "$PKG_CONFIG_PATH" | /bin/sed -e 's|/home/florent/Soft/Omicron/lib/pkgconfig:||g;'`
PKG_CONFIG_PATH="/home/florent/Soft/Omicron/lib/pkgconfig:$PKG_CONFIG_PATH"
export PKG_CONFIG_PATH

##### PATH
PATH=`echo "$PATH" | /bin/sed -e 's|/home/florent/Soft/Omicron/bin:||g;'`
PATH="/home/florent/Soft/Omicron/bin:$PATH"
export PATH

##### LD_LIBRARY_PATH
LD_LIBRARY_PATH=`echo "$LD_LIBRARY_PATH " | /bin/sed -e 's|/home/florent/Soft/Omicron/lib:||g;'`
LD_LIBRARY_PATH="/home/florent/Soft/Omicron/lib:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH



##### Omicron environment
nodename=`uname -n`

OMICRON_HTML="/home/florent/Soft/Omicron/share/omicron/html"
export OMICRON_HTML

# Lyon
if [[ $nodename == *"cca"* ]]; then
    OMICRON_TRIGGERS="/hpss/in2p3.fr/group/virgo/DETCHAR/triggers/Omicron"

# Florent
elif [[ $nodename == "florent" ]]; then
    OMICRON_TRIGGERS="/home/florent/Analysis/Omicron/triggers"

# Cascina
elif [[ $nodename == *"virgo"* ]]; then
    OMICRON_TRIGGERS="/data/omicron"

# LIGO
elif [[ $nodename == *"ldas"* ]]; then
    OMICRON_TRIGGERS="/home/detchar/triggers"

# unknown
else
    OMICRON_TRIGGERS=""
fi
export OMICRON_TRIGGERS
