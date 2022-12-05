import task
import config
from primer import Candidates



def create_index(ref_file:str, orgname:str)->None:
    task.run(f"cd {config.homepath}/bowtie/index/ && bowtie-build -f {ref_file} {orgname} && cd {config.homepath}/input/script")


def write_primers()->None:
    with open("primercheck1.fasta", "w") as file:
        for idx, cand in enumerate(Candidates.all):
            file.write(f">{cand.target}_{idx}_fwd\n")
            file.write(f"{cand.left_sequence}\n")
            file.write(f">{cand.target}_{idx}_rev\n")
            file.write(f"{cand.right_sequence}\n")


def run(inp_file:str, out_file:str, index:str)->None:
    # idx example : s_cerevisiae
    task.run(f"bowtie -x {index} -a -f {inp_file} -v 3 > {out_file}")


def parse(bt_out_file:str)->list[str]:
    with open(bt_out_file) as file:
        data = file.read()

    bowtie_data = []
    for line in data.split("\n"):
        if line:
            bowtie_data.append(line.split("\t"))

    return bowtie_data