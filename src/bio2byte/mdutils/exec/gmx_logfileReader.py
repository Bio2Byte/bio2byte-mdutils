#!/usr/bin/env python
import argparse
import os
import sys
import re
import csv


def parse_cmdline_arguments(arguments):

    aparser = argparse.ArgumentParser()
    aparser.add_argument("logfiles", metavar="<logfile>", nargs="*",)
    aparser.add_argument("--stdin", action='store_true',
        help=("Reads logfiles from STDIN"))
    aparser.add_argument("-o", "--outfile", metavar="<outfile>",
                         default="log_report.csv")
    
    args = aparser.parse_args(arguments)


    return args


class gmxLogfileParser:

    def __init__(self, additional_fields=[]):
        self.__entries = []
        self.__fieldnames = list(additional_fields) + [
            "n_atoms",
            "n_vsites",
            "n_particles",          # n_atoms + n_visites
            "n_particles_solutes",  # n_atoms - 3 * n_visites
            "n_particles_solvent",  # 4 * n_visites
            "nsteps",
            "threads_MPI",
            "threads_oMP",
            "wall_time",
            "cpu_info",
            "gpu_info",
        ] 

    def write_csv(self, outfile):

        with open(outfile, "w") as ofhandle:
            writer = csv.DictWriter(ofhandle, fieldnames=self.__fieldnames)
            writer.writeheader()
            for row in self.__entries:
                writer.writerow(row)
        
    def parse(self, logfile, additional_fields):

        new_entry = additional_fields.copy()

        with open(logfile, "r") as fhandle:
            # Read hardware setup
            for ln in fhandle:
                if ln.startswith("Hardware detected"): break
            new_entry.update(self._parse_hardware(fhandle))

            # Read nsteps
            for ln in fhandle:
                if ln.startswith("   nsteps"):
                    break
            new_entry["nsteps"] = int(ln.split()[-1])

            # Read threading info
            ptn = re.compile("Using (\d+) MPI ")
            for ln in fhandle:
                try:
                    m = ptn.match(ln)
                    new_entry["threads_MPI"] = int(m.group(1))
                    break
                except AttributeError:
                    pass
            ptn = re.compile("Using (\d+) OpenMP ")
            for ln in fhandle:
                try:
                    m = ptn.match(ln)
                    new_entry["threads_oMP"] = int(m.group(1))
                    break
                except AttributeError:
                    pass

            # Read atoms
            n_atoms, n_vsites = 0, 0
            for ln in fhandle:
                if ln.startswith("There are:"):
                    n_atoms = int(ln.split()[2])
                    ln = next(fhandle)
                    n_vsites = int(ln.split()[2])
                    break
            new_entry["n_atoms"] = n_atoms
            new_entry["n_vsites"] = n_vsites
            new_entry["n_particles"] = n_atoms + n_vsites
            new_entry["n_particles_solutes"] = n_atoms - 3 * n_vsites
            new_entry["n_particles_solvent"] = 4 * n_vsites
            
            # Read performance table
            for ln in fhandle:
                if ln.startswith("     R E A L   C Y C L E   A N D   T I M E   A C C O U N T I N G"):
                    break

            # Read performance table
            for ln in fhandle:
                if ln.startswith(" Total"):
                    lns = ln.split()
                    new_entry["wall_time"] = float(lns[1])
                    break
        self.__entries.append(new_entry)
    
    def _parse_hardware(self, fhandle):
        hardware_info = dict(cpu_info="", gpu_info="")
        cpu_pattern = re.compile("\s*Brand:\s+(\w.*)$")
        gpu_pattern = re.compile("\s*#\d:\s+(\w.*?),")
        
        ln = next(fhandle).strip()
        while ln:
            if cpu_pattern.match(ln):
                m = cpu_pattern.match(ln)
                hardware_info["cpu_info"] = m.group(1)
            elif gpu_pattern.match(ln):
                m = gpu_pattern.match(ln)
                hardware_info["gpu_info"] = m.group(1)
            ln = next(fhandle).strip()
        
        return hardware_info
    

def main():
    args = parse_cmdline_arguments(sys.argv[1:])

    LogfileParser = gmxLogfileParser(additional_fields=[
        "system_name", "xtc_size", "log_file"])

    if args.stdin:
        stdin_logfiles = list(map(lambda s: s.strip(), sys.stdin.readlines()))
    else:
        stdin_logfiles = []

    for logfile in stdin_logfiles + args.logfiles:

        print("Processing %s..." % logfile)

        try:
            systname = re.search("/(gg\w+gg)/", os.path.abspath(logfile)).group(1)
            xtcsize = os.path.getsize(os.path.splitext(logfile)[0] + '.xtc')
        except AttributeError:
            systname = "unknown"
        except OSError:
            xtcsize = -1

        LogfileParser.parse(logfile, additional_fields={
            "system_name": systname,
            "xtc_size": xtcsize,
            "log_file": logfile
        })

    LogfileParser.write_csv(args.outfile)

    return 0


if __name__ == "__main__":
    sys.exit(main())

    
