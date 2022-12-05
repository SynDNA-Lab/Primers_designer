import re
import os

import config
import task
from primer import Candidates, Regions


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

