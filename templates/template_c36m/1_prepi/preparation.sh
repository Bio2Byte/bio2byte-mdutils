#!/usr/bin/env bash
#
# If mathematical analysis should ever hold a prominent place in chemistry - 
# an aberration which is happily almost impossible - it would occasion a rapid 
# and widespread degeneration of that science. (Aguste Comte, 1830)
#
# usage:  preparation.sh <INPUT_PDB>
#

# Load helper functions
check_success() {
    ERRCODE=$?
    if [ $ERRCODE -ne 0 ]; then
        echo "[ERROR:${ERRCODE}] $*"
        exit $ERRCODE
    else
        return $ERRCODE
    fi
}

snapshotPDB() {
    topFile="$1"
    crdFile="$2"

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

INPUT_PDB="$1"
SYSTNAME="${INPUT_PDB%.*}"

# Check if the INPUT_PDB is provided
test -f $INPUT_PDB
check_success "Input structure [.pdb] not found!"
# Check if GROMACS is installed and available
which gmx &> /dev/null
check_success "Gromacs not found! Please, ensure the software is sourced and available in \$PATH."

# ::: PREPARE FOR GMX

gmx pdb2gmx \
    -f "$INPUT_PDB" \
    -o "${SYSTNAME}_4gmx.gro" \
    -p "${SYSTNAME}.top" \
    -ter -ignh
check_success

# ::: GENERATE PBC BOX

gmx editconf \
    -f "${SYSTNAME}_4gmx.gro" \
    -o "${SYSTNAME}_boxed.gro" \
    -c -d 1.2 -bt dodecahedron
check_success

# ::: SOLVATE
if $( grep -iq 'tip3p' "${SYSTNAME}.top"); then
    SOLCRD=spc216.gro
else
    SOLCRD=tip4p.gro
fi
gmx solvate \
    -cp "${SYSTNAME}_boxed.gro" \
    -cs "$SOLCRD" \
    -o "${SYSTNAME}_solv.gro" \
    -p "${SYSTNAME}.top"
check_success

# :: ADD IONS

cat <<_EOF > add_ions.mdp
; add_ions.mdp
;   Used to create a *.tpr for adding ions only.

; Energy minimization parameters
integrator = steep
emtol      = 100
nsteps     = 10000

; Neighbor list
cutoff-scheme           = Verlet
rlist                   = 1.2
pbc                     = xyz
verlet-buffer-tolerance = -1

; Electrostatics
coulombtype  = Cut-off
rcoulomb     = 1.2

; VdW
vdwtype      = Cut-off
vdw-modifier = Potential-shift
rvdw         = 1.2
rvdw-switch  = 1.0
DispCorr     = no

continuation = yes
_EOF

gmx grompp \
    -f add_ions.mdp \
    -c "${SYSTNAME}_solv.gro" \
    -p "${SYSTNAME}.top" \
    -po add_ions.mdpout \
    -o add_ions.tpr
check_success

echo SOL | gmx genion \
    -s add_ions.tpr \
    -p "${SYSTNAME}.top" \
    -o "${SYSTNAME}_solv_ions.gro" \
    -np 0 -pname NA -pq 1 \
    -nn 0 -nname CL -nq -1 \
    -conc 0.10 -neutral
check_success

# ::: Create links for final files 
ln -sf "${SYSTNAME}.top" TOPOL.top
ln -sf "${SYSTNAME}_solv_ions.gro" COORD.gro

snapshotPDB TOPOL.top COORD.gro


