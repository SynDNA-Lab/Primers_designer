import re
import os
import pandas as pd 

import config
import task
from primer import Candidates, Regions
from selection import Selection
from dataclasses import dataclass, field



class Primer3Interface:
    selection: Selection
    result_path: str
    primer_sites: pd.DataFrame = field(init=False)


    def run(self):
        self.write_settings()
        self.run_primer3()
        self.parse_results()
        self.generate_fasta()


    def write_settings(self):
        pass


    def run_primer3(self):
        pass

    
    def parse_results(self):
        with open(self.result_path) as file:
            data = file.read()

        pattern = r"(?s)(?=PRIMER_PAIR_([0-9]+)_PENALTY)(.*?)(.*?PRIMER_PAIR_[0-9]+_PRODUCT_TM=[0-9]+.[0-9]+)"
        re_result = re.findall(pattern , data)

        dfdat = []
        for r in re_result:
            d = r[-1].split("\n")
            e = [x for x in d if x]
            dfdat.append([x.split("=")[1] for x in e])

        cols = ["primer_pair_penalty",
                "primer_left_penalty",
                "primer_right_penalty",
                "primer_left_sequence",
                "primer_right_sequence",
                "primer_left",
                "primer_right",
                "primer_left_tm",
                "primer_right_tm",
                "primer_left_gc_percent",
                "primer_right_gc_percent",
                "primer_left_self_any_th",
                "primer_right_self_any_th",
                "primer_left_self_end_th",
                "primer_right_self_end_th",
                "primer_left_hairpin_th",
                "primer_right_hairpin_th",
                "primer_left_end_stability",
                "primer_right_end_stability",
                "primer_pair_compl_any_th",
                "primer_pair_compl_end_th",
                "primer_pair_product_size",
                "primer_pair_product_tm"]

        df2 = pd.DataFrame(data=dfdat, columns=cols)
        
        bla = []
        seqids = re.findall(r"SEQUENCE_ID=.+?(?=\n)", data)
        res =  re.findall(r"PRIMER_PAIR_NUM_RETURNED=.+?(?=\n)", data)

        for sid, rs in zip(seqids, res):
            s = sid.split("=")[1]
            for ids in range(int(rs.split("=")[1])):
                bla.append(f"{s}_{ids}")
        df2.insert(0, "name", bla)
        # self.primer_sites


    def generate_fasta(self) -> str:
        buffer = ""
        for _, col in self.primer_sites.iterrows():
            buffer += f">{col['name']}_fwd\n{col['primer_left_sequence']}\n{col['name']}_rev\n{col['primer_right_sequence']}\n"

        with open("generated_primer.fasta", "w") as file:
            file.write(buffer)


# Aufgaben:
# - erstellt die notwendige primer3 settings datei
# - startet primer3 
# - parsed die primer3 results
# - stored die primer3 results
# - generates all the fasta results

# Deleted make_settings(), replace later

def process(result_path:str) -> None:
    '''
    Parent function to process everything Primer3 related
    '''

    task.run(f"{config.primer3_path} settings > {result_path}")
    parse(result_path)
    generate_fasta()



def parse(result_path:str) -> None:
    with open(result_path) as file:
        result_buffer = file.read()
        output = result_buffer.split("=\n")

    for out in output[1:]:
        for i in re.split("PRIMER_PAIR_[0-9]*_PENALTY", out)[1:]:
            name = [re.split("PRIMER_PAIR_[0-9]*_PENALTY", out)[0].split("=")[1].split("\n")[0]]
            buffer = [n.split("=")[1] for n in i.split("\n") if "=" in n]
            name.extend(buffer)
            Candidates(name) #generate candidates


def generate_fasta():
    with open("primercheck.fasta", "w") as file:
        for idx, cand in enumerate(Candidates.all):
            file.write(f">{cand.target}_{idx}_fwd\n")
            file.write(f"{cand.left_sequence}\n")
            file.write(f">{cand.target}_{idx}_rev\n")
            file.write(f"{cand.right_sequence}\n")

