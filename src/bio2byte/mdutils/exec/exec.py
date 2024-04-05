"""
A python-wrapper for the bash scripts that were used, so that they become
invocable within the package.
"""

import os, sys
from shutil import which
import subprocess

BASH_EXECUTABLE = which("bash")

PROCESS_TRAJECTORIES_SCRIPT = os.path.join(os.path.dirname(__file__), "process_trajectories.sh")

ANALYZE_HBONDS_SCRIPT = os.path.join(os.path.dirname(__file__), "mediated_hbonds.sh")

def process_trajectories():
    cmd = " ".join([BASH_EXECUTABLE, PROCESS_TRAJECTORIES_SCRIPT, *sys.argv[1:]])
    p = subprocess.run(cmd, shell=True)
    sys.exit(p.returncode)

def hbonds_wrapper():
    cmd = " ".join([BASH_EXECUTABLE, ANALYZE_HBONDS_SCRIPT, *sys.argv[1:]])
    p = subprocess.run(cmd, shell=True)
    sys.exit(p.returncode)