import glob, os, re
import numpy as np
import pandas as pd
from typing import Iterable, Tuple

aminoacids1to3 = dict(
	A="ALA", C="CYS", D="ASP", E="GLU", F="PHE",
	G="GLY", H="HIS", I="ILE", K="LYS", L="LEU",
	M="MET", N="ASN", P="PRO", Q="GLN", R="ARG",
	S="SER", T="THR", V="VAL", W="TRP", Y="TYR",
)

aminoacids3to1 = {j: i for i, j in aminoacids1to3.items()}


class BmrbComparisonReader:
    
    def __init__(self, apopattern=".+", modpattern=".+", identity_margin=(8,8)):
        self.apopattern = re.compile(apopattern)
        self.modpattern = re.compile(modpattern)
        self.identity_margin = identity_margin
        self.dataframe = pd.DataFrame()
        
    def read_file(self, filepath, *, append=True):
        """ Read BMRB data from a comparison file """
        rawdata = pd.read_csv(filepath)
        nterm, cterm = self.identity_margin
        new_data = []
        for i, (aporesn, aporesi, modresn, modresi)  in rawdata[[
                "apo.resn", "apo.resi", "mod.resn", "mod.resi"]].iterrows():
            try:
                assert self.modpattern.match(modresn) and self.apopattern.match(aporesn)
                df = rawdata.loc[i-nterm:i+cterm].copy()                
                seqident = (df["apo.resn"] == df["mod.resn"])
                assert seqident.sum() == nterm + cterm and not seqident.loc[i]
                assert (df["apo.rci"] < .8).all()
            except (AssertionError, TypeError):
                continue
            else:
                df.reset_index(inplace=True)
                df["apo.bmrbId"] = df["apo.bmrbId"].astype(int)
                df["mod.bmrbId"] = df["mod.bmrbId"].astype(int)
                df["aporesn"] = aporesn
                df["modresn"] = modresn
                df["idx"] = np.arange(-nterm, cterm+1, dtype=int)
                df["seq"] = self.seq2str(df["apo.resn"],  df["mod.resn"])
                new_data.append(df)
        if append:
            self.dataframe = pd.concat([self.dataframe, *new_data], ignore_index=True)
        else:
            self.dataframe = pd.concat(new_data, ignore_index=True)
    
    def get_delta_dataframe(self):
        """Returns a delta dataframe between apo and mod (mod - apo)"""
        df = self.dataframe.loc[:, 
            ["seq", "idx", 'apo.resn', 'mod.resn', 'apo.bmrbId', 'mod.bmrbId']]
        df["delta.rci"]   = self.dataframe["mod.rci"]   - self.dataframe["apo.rci"]  
        df["delta.helix"] = self.dataframe["mod.helix"] - self.dataframe["apo.helix"]
        df["delta.beta"]  = self.dataframe["mod.beta"]  - self.dataframe["apo.beta"] 
        df["delta.coil"]  = self.dataframe["mod.coil"]  - self.dataframe["apo.coil"] 
        df["delta.ppII"]  = self.dataframe["mod.ppII"]  - self.dataframe["apo.ppII"]
        return df
    
    @staticmethod
    def seq2str(aposeq, modseq):
        """Returns a sequence-like string indikating any non-standard aa or substitutions"""
        seqstr = ""
        for ares, mres in zip(aposeq, modseq):
            if ares == mres:
                seqstr += aminoacids3to1.get(ares.upper(), ares.capitalize())
            else:
                seqstr += "[{0}|{1}]".format(
                    aminoacids3to1.get(ares.upper(), ares.capitalize()),
                    aminoacids3to1.get(mres.upper(), mres.capitalize()))
        return seqstr

        
class MdDataParser:
    
    def __init__(self, apopatterns: Iterable[str], modpatterns: Iterable[str]):
        self.apopattern = [apopatterns] if isinstance(apopatterns, str) else apopatterns
        self.modpattern = [modpatterns] if isinstance(modpatterns, str) else modpatterns
        self.df = pd.DataFrame()
        self.segid = 0
        
    def get_corresponding_apodir(self, mod_dirpath, *, return_count=False):
        dpath, dname = os.path.split(mod_dirpath)
        count = 0
        for apo, mod in zip(self.apopattern, self.modpattern):
            for m in re.finditer(f"g{mod}|{mod}g", dname):
                count += 1
                i,j = m.span()
                dname = dname[:i] + dname[i:j].replace(mod, apo) + dname[j:]
        if return_count:
            return os.path.join(dpath, dname), count
        else:
            return os.path.join(dpath, dname)
    
    def pair_apo_mod(self, dirpaths):
        for modpath in dirpaths:
            apopath, count = self.get_corresponding_apodir(modpath, return_count=True)
            if count == 1 and apopath in dirpaths:
                yield apopath, modpath
    
    def read_dirs(self, apo_dir, mod_dir, *, method="window/3/", append=True):
        
        apores = re.findall("g|[A-Z1]{3}", os.path.basename(apo_dir))
        modres = re.findall("g|[A-Z1]{3}", os.path.basename(mod_dir))
        modidx = next(i for i,(a,m) in enumerate(zip(apores, modres)) if a != m)

        apoConstava = pd.read_csv(os.path.join(apo_dir, "analy", "constava.csv"))
        apoConstava = apoConstava.loc[apoConstava["Method"] == method]
        apoCircvar  = pd.read_csv(os.path.join(apo_dir, "analy", "circvar.csv"))
        apoDssp     = pd.read_csv(os.path.join(apo_dir, "analy", "dssp_summary.dat"))
        apoDf = apoConstava.merge(apoCircvar, on=["#ResIndex", "ResName"])
        apoDf = apoDf.merge(apoDssp, on=["#ResIndex", "ResName"])
        apoDf = apoDf.rename(columns={c: f"apo.{c}" for c in apoDf.columns if not c.startswith("#")})
        del apoConstava, apoCircvar, apoDssp

        modConstava = pd.read_csv(os.path.join(mod_dir, "analy", "constava.csv"))
        modConstava = modConstava.loc[modConstava["Method"] == method]
        modCircvar  = pd.read_csv(os.path.join(mod_dir, "analy", "circvar.csv"))
        modDssp     = pd.read_csv(os.path.join(mod_dir, "analy", "dssp_summary.dat"))
        modDf = modConstava.merge(modCircvar, on=["#ResIndex", "ResName"])
        modDf = modDf.merge(modDssp, on=["#ResIndex", "ResName"])
        modDf = modDf.rename(columns={c: f"mod.{c}" for c in modDf.columns if not c.startswith("#")})
        del modConstava, modCircvar, modDssp

        merged = apoDf.merge(modDf, on="#ResIndex")
        merged["idx"] = np.arange(-modidx, len(apores) - modidx, 1, dtype=int)
        merged["segid"] = self.segid
        self.segid += 1
        if append:
            self.df = pd.concat([self.df, merged], axis=0, ignore_index=True)
        else:
            self.df = merged
            
    def get_delta_df(self):
        
        df = self.df.loc[:, ["segid", "idx", "apo.ResName", "mod.ResName"]]
        
        df["d.Var"] = self.df["mod.ConStaVa"] - self.df["apo.ConStaVa"]
        df["d.coreHelix"] = self.df["mod.coreHelix"] - self.df["apo.coreHelix"]
        df["d.surrHelix"] = self.df["mod.surrHelix"] - self.df["apo.surrHelix"]
        df["d.coreSheet"] = self.df["mod.coreSheet"] - self.df["apo.coreSheet"]
        df["d.surrSheet"] = self.df["mod.surrSheet"] - self.df["apo.surrSheet"]
        df["d.Turn"] = self.df["mod.Turn"] - self.df["apo.Turn"]
        df["d.Other"] = self.df["mod.Other"] - self.df["apo.Other"]
        
        df["d.dssp[H]"] = self.df["mod.dssp[H]"] - self.df["apo.dssp[H]"]
        df["d.dssp[B]"] = self.df["mod.dssp[B]"] - self.df["apo.dssp[B]"]
        df["d.dssp[E]"] = self.df["mod.dssp[E]"] - self.df["apo.dssp[E]"]
        df["d.dssp[G]"] = self.df["mod.dssp[G]"] - self.df["apo.dssp[G]"]
        df["d.dssp[I]"] = self.df["mod.dssp[I]"] - self.df["apo.dssp[I]"]
        df["d.dssp[P]"] = self.df["mod.dssp[P]"] - self.df["apo.dssp[P]"]
        df["d.dssp[T]"] = self.df["mod.dssp[T]"] - self.df["apo.dssp[T]"]
        df["d.dssp[S]"] = self.df["mod.dssp[S]"] - self.df["apo.dssp[S]"]
        df["d.dssp[C]"] = self.df["mod.dssp[C]"] - self.df["apo.dssp[C]"]
        
        df["d.CircVar"] = self.df["mod.CV"] - self.df["apo.CV"]
        
        df["d.helix"] = df["d.coreHelix"] + df["d.surrHelix"]
        df["d.beta"] = df["d.coreSheet"] + df["d.surrSheet"]
        df["d.coil"] = df["d.Other"]
        
        return df