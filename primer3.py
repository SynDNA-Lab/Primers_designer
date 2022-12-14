import re
import shutil
import pandas as pd 
from subprocess import Popen, PIPE

from selection import Selection
from dataclasses import dataclass, field



@dataclass
class Primer3Interface:
    selection: Selection
    primer3_path: str 
    result_path: str = field(default="primer3_result")
    primer_sites: pd.DataFrame = field(init=False)


    def run(self):
        self.write_settings()
        self.run_primer3()
        self.parse_results()
        self.generate_fasta()


    def write_settings(self):
        # create a primer3 settings file
        shutil.copyfile("settings.bak", "settings")
        with open("settings", "a") as file:
            for site in targets:
                start = min(site[0], 400)
                seq = target_sequence[site[0]-start:site[1]+400]
                file.write(f"SEQUENCE_ID=recomb_site_{site[0]}\n")
                file.write(f"SEQUENCE_TEMPLATE={seq}\n")
                file.write(f"SEQUENCE_TARGET={start},{site[1]-site[0]}\n")
                file.write("=\n")


    def run_primer3(self) -> None:
        cmd = f"{self.primer3_path} settings > {self.result_path}"
        process = Popen(args=cmd, stdout=PIPE, stderr=PIPE, shell=True)
        _, stderr = process.communicate()
        if stderr:
            print(stderr.decode("ascii"))


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
