#!/usr/bin/env bash
# Author:   David Bickel
# Date:     23/03/2022
#
# usage: ./process_trajectories.sh <MD-DIRECOTRY>
#

# FUNCTION DEFINITONS
usage() {
    cat <<HELPTXT

usage: $0 <MD-DIRECTORY>

This script automatically processes a MD trajectory. Automation 
is achieved through a strict folder structure. Inside the <MD-DIRECTORY>
The script is looking for one ore more directories named '3_prod'.
They should contain run files (*.tpr) as well as trajectories (*.xtc).

Finally all found trajectories centered on the protein, imaged
inside the PBC box, and solvent molecules are stripped from the 
files. The resulting trajectories are concatenated and written
out as:
    * stripped.tpr
    * stripped.xtc

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

# SET SOME GLOBAL VARIABLES
unalias grep ls &> /dev/null
CENTER_GRP=Protein
OUTPUT_GRP=Protein
PBC_OPTION=mol
UR_OPTION=compact
OUTPUT_TPR="stripped.tpr"
OUTPUT_XTC="stripped.xtc"

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
analyDir="$mdDir/analy"
mkdir -v "$analyDir"; check_success "Cannot create analy directory: $analyDir"
cd "$analyDir" && echo "Goto $analyDir"

# INITIALIZE A LOGFILE (only retained in case of error)
_log=$(mktemp logfile_XXXXXXXXX.log)    

# PROCESS REPLICA TRAJECTORIES INDIVIDUALLY
runTrajectories=( )
for prodDir in ${prodDirs[@]}; do

    if [ ! -f residuetypes.dat ]; then
        if [ -f "$prodDir/../1_prepi/residuetypes.dat" ]; then
            cp -p "$prodDir/../1_prepi/residuetypes.dat" residuetypes.dat
        elif [ -f "$VSC_HOME/residuetypes.dat" ]; then
            cp -p "$VSC_HOME/residuetypes.dat" residuetypes.dat
        else
            echo "[ERROR] Cannto find file: residuetypes.dat" 1>&2
            exit 1
        fi
    fi

    echo "Processing trajectory information in: $prodDir"
    runN=$(echo "$prodDir" | grep -oP 'run\d+')
    topFiles=( $(find "$prodDir" -maxdepth 1 -name "*.tpr" | sort | xargs ) )
    trjFiles=( $(find "$prodDir" -maxdepth 1 -name "*.trr" -o -name "*.xtc" | sort | xargs) )
    echo "  Topology:     $(basename $topFiles)"
    echo -n "  Trajectories:"; 
    for t in ${trjFiles[@]}; do echo -n " $(basename $t)"; done; echo
    
    
    _xtc0=$(mktemp --dry-run ${runN}_XXXXXXXXX.xtc)
    if [ "${#trjFiles[@]}" -gt 1 ]; then
        echo "  Concatenating trajectories..."
        gmx trjcat -f ${trjFiles[@]} -o $_xtc0 &>> $_log
        check_success "Concatenating trajectories failed ($runN)."
    else
        ln -sv "${trjFiles}" "$_xtc0"
    fi

    echo -n "  Stripping and wrapping trajectory..."
    _xtc1=$(mktemp --dry-run ${runN}_XXXXXXXXX.xtc)
    printf "$CENTER_GRP\n$OUTPUT_GRP" | \
    gmx trjconv -s "$topFiles" \
                -f "$_xtc0" \
                -o "$_xtc1" \
                -pbc "$PBC_OPTION" -ur "$UR_OPTION" \
                -center &>> "$_log"
    check_success "Wrapping trajectory failed ($runN)."
    echo " => $_xtc1"
    
    if [ ! -f "$OUTPUT_TPR" ]; then
        echo -n "  Write stripped topology file..."
        _ndx=$(mktemp --dry-run ${runN}_XXXXXXXXX.ndx)
        printf "q\n" | \
            gmx make_ndx -f "$topFiles" -o "$_ndx" &>> "$_log"
        printf "$OUTPUT_GRP" | \
            gmx convert-tpr -s "$topFiles" \
                            -n "$_ndx" \
                            -o "$OUTPUT_TPR" &>> "$_log"
        check_success "Writing stripped topology file failed ($runN)."
        echo " => $OUTPUT_TPR"
        rm "$_ndx"
    fi
    
    runTrajectories+=( "$_xtc1" )
    rm "$_xtc0"
done

echo -n "Concatenating all..."
for t in ${runTrajectories[@]}; do echo c; done | \
    gmx trjcat -b 1 -settime -f ${runTrajectories[@]} -o $OUTPUT_XTC &>> "$_log"
echo " => $OUTPUT_XTC"
rm -v ${runTrajectories[@]} "$_log"
