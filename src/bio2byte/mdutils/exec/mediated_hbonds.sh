#!/usr/bin/env bash
# Author:   David Bickel
# Date:     08/02/2023
#
# usage: ./analyze_mediated_hbonds.sh <MD-DIRECTORY>
#

# FUNCTION DEFINITONS
usage() {
    cat <<HELPTXT

usage: $0 <MD-DIRECTORY>

This script is a shell wrapper for the identically named python
script. It takes an <MD-DIRECTORY>, automatically processes the
trajectories in their - retaining their solvent molecules - and
then calculates hydrogen bonds within the solute as well as 
hydrogen bonds that are mediated by one solvent molecule.

HELPTXT
}

check_success() {
    ERR=$?
    if [ "$ERR" -ne 0 ]; then
        echo "[ERR:${ERR}] $*"
        exit $ERR
    fi
}

# CHECK IF GROMACS IS AVAILABLE
which gmx &> /dev/null
check_success "Cannot find GROMACS executable."

# CHECK IF HBOND ANALYSIS SCRIPT IS AVAILABLE
HBOND_ANALYSIS_EXE="mdutils-hbonds"
which "$HBOND_ANALYSIS_EXE" &> /dev/null
check_success "Cannot find hydrogen bond analysis script. Is bio2byte-mdutils properly installed?"

# SOME GLOBAL VARIABLE DEFINITIONS
unalias grep ls &> /dev/null
CENTER_GROUP=Protein
OUTPUT_GROUP=System
PBC_OPTION=mol
UR_OPTION=compact
OUTPUTFILE=hbonds.csv

# PARSE COMMAND LINE PARAMETERS
if [ -z "$1" ]; then
    mdDir=$(pwd)
elif [ -d "$1" ]; then
    mdDir=$(readlink -f "$1")
else
    usage
    exit 1
fi

# DIRECTORIES TO PROCESS
prodDirs=( $(find "$mdDir" -type d -name "3_prod" | sort | xargs) )
[ ${#prodDirs[@]} -gt 0 ]; check_success "Cannot find an directories '3_prod' in: '$mdDir'"

# INITIALIZE OUTPUT DIRECTORY
analdir="$mddir/analy"
mkdir -p "$analdir"
cd "$analdir" && echo "Goto $analdir" || exit 1
# [ -f "$OUTPUTFILE" ] && exit 0

_log=$(mktemp logfile_XXXXXXXXX.log)
for proddir in ${prodDirs[@]}; do

    [ ! -e residuetypes.dat ] &&  ln -s "$proddir/../1_prepi/residuetypes.dat"
    
    echo "Getting trajectories in '$proddir'..."
    runN=$(echo "$proddir" | grep -oP 'run\d+')
    topFiles=( $(find "$proddir" -maxdepth 1 -name "*.tpr" | sort | xargs ) )
    trjFiles=( $(find "$proddir" -maxdepth 1 -name "*.trr" -o -name "*.xtc" | sort | xargs) )
    echo "  Topology:     $(basename $topFiles)"
    echo -n "  Trajectories:"
    for t in ${trjFiles[@]}; do echo -n " $(basename $t)"; done; echo
    
    _xtc0=$(mktemp --dry-run ${runN}_XXXXXXXXX.xtc)
    if [ "${#trjFiles[@]}" -gt 1 ]; then
        echo "Concatenating trajectories..."
        gmx trjcat -f ${trjFiles[@]} -o $_xtc0 &>> $_log
        check_success "Concatenating trajectories failed ($runN)."
    else
        ln -sv "${trjFiles}" "$_xtc0"
    fi
    
    echo "Imaging trajectories..."
    _xtc1=$(mktemp --dry-run ${runN}_XXXXXXXXX.xtc)
    printf "$CENTER_GROUP\n$OUTPUT_GROUP\n" | \
    gmx trjconv -s "$topFiles" \
                -f "$_xtc0" -b 1 \
                -o "$_xtc1" \
                -pbc "$PBC_OPTION" -ur "$UR_OPTION" \
                -center &>> "$_log"
    check_success "Wrapping trajectory failed ($runN)."
    echo " => $_xtc1"
    rm "$_xtc0"
    
    echo "Analyze hydrogen bonds..."
    [ -f hbonds.csv ] && append_flag='--append'        
    $HBOND_ANALYSIS_EXE $append_flag \
        -s "$topFiles" -f "$_xtc1" -o "$OUTPUTFILE" \
        --sele0 "not resname SOL" --sele1 "not resname SOL" \
        --solvent "resname SOL"
    rm "$_xtc1"
done

rm "$_log"
