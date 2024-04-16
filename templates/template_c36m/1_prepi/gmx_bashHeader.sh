#
# This file contains a number of bash functions for running
# MD simulations with gromacs.
# It is intended to be sourced by other scripts and thus 
# needs no "executable" flags.
#
# Index:
#
#   1. Check if gromacs is sourced
#
#   2. Helper functions
#       - check_success()
#       - get_last()
#       - snapshotPDB()
#

#
# +++ 1.  Generic check, if Gromacs is sourced +++
#
if $(which gmx &> /dev/null); then
    :
else
    echo "[ERROR] Gromacs not found! Please, ensure the software is sourced and available in \$PATH."
    exit 1
fi


#
# +++ 2.  Helper functions +++
#

# Check if the prior command executed successfully
check_success() {
    ERRCODE=$?
    if [ $ERRCODE -ne 0 ]; then
        echo "[ERROR:${ERRCODE}] $*"
        exit $ERRCODE
    fi
}


# Returns the last one out of a list of files
get_last() {
    lastArg=${!#}
    [ -f "$lastArg" ] && echo $lastArg
    unset lastArg
}


# Updates the tinit in the template *.mdp ($1) with the last
# timestep from a prior *.trr/*.cpt ($2)
update_tinit() {
    mdp_="$1"
    trr_="$2"
    t_=$(gmx check -f "$trr_" 2>&1 | grep -oP 'Last frame +-?\d+ +time +\K\d+\.?\d*')

    if $(grep -E '^tinit *=' "$mdp_" &> /dev/null); then
        sed -i "s/^\(tinit *=\).*/\1 $t_/" "$mdp_"
    else
        echo "tinit = $t_" >> "$mdp_"
    fi

    unset trr_ mdp_ t_
}


# Generates a wrapped PDB representation of the input files
#
# usage:  snapshotPDB <GRO> <TOP>
#
snapshotPDB() {
    crdFile="$1"
    topFile="$2"

    echo Generating snapshot \"$(basename $crdFile).pdb\" from:
    echo "  Topology: $topFile"
    echo "  Coordinates: $crdFile"

    [ -n ${GMX_MAXBACKUP+x} ] && gmxBakBak=$GMX_MAXBACKUP
    export GMX_MAXBACKUP=-1
    
    mdpFile=$(mktemp -p . --suffix=.mdp)
    tprFile=$(mktemp -p . --suffix=.tpr)
    mdoFile=$(mktemp -p . --suffix=.mdp)
    cat <<_EOF > $mdpFile 
; Minimal *.mdp for generation of *.tpr only
integrator = steep
nsteps = 1
_EOF
    gmx grompp \
        -f  $mdpFile  \
        -c  $crdFile \
        -p  $topFile \
        -po $mdoFile \
        -o  $tprFile &> /dev/null
    
    echo -e '1\n0' | \
    gmx trjconv \
        -f $crdFile \
        -s $tprFile \
        -o "$(basename $crdFile).pdb" \
        -center -pbc mol -ur compact &> /dev/null
    
    rm "$mdpFile" "$tprFile" "$mdoFile"
    unset GMX_MAXBACKUP
    [ -n ${gmxBakBak+x} ] && export GMX_MAXBACKUP=$gmxBakBak
    unset gmxBakBak
}
