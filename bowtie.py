import pandas as pd 
from subprocess import Popen, PIPE
from dataclasses import dataclass, field

from config import Config
import os


BASENAME = "target"

@dataclass
class BowtieResult:
    result_path: str
    result: pd.DataFrame = field(init=False)

    def __post_init__(self) -> None:
        self.parse_data()


    def split_data(self, lst:list[str], instruction:str) -> list[str]:
        if instruction == "orientation":
            return [l.split("_")[-1] for l in lst]
        elif instruction == "id":
            return ["_".join([l.split("_")[-3], l.split("_")[-2]]) for l in lst]
        else:
            raise ValueError(f"Unknown bowtie parser instruction {instruction}")

    def split_data_corrected(self, lst: list[str], instruction: str) -> list[str]:
        if instruction == "orientation":
            return [l.split("_")[-1] if len(l.split("_")) >= 2 else None for l in lst]
        elif instruction == "id":
            return [
                "_".join([l.split("_")[-3], l.split("_")[-2]])
                if len(l.split("_")) >= 3 else None
                for l in lst
            ]
        else:
            raise ValueError(f"Unknown bowtie parser instruction {instruction}")


    def parse_data(self) -> None:
        expected_columns = ["name", "strand", "reference", "start", "sequence", "quality", "instances", "mismatch_descriptor"]
        
        if not os.path.exists(self.result_path) or os.stat(self.result_path).st_size == 0:
            df = pd.DataFrame(columns=expected_columns + ["id", "orientation"])
        else:
            df = pd.read_csv(self.result_path, sep="\t", header=None)
            df.columns = expected_columns
            df = df.drop_duplicates()
            df["id"] = self.split_data(df["name"], "id")
            df["orientation"] = self.split_data(df["name"], "orientation")

        self.result = df




@dataclass
class BowtieInterface:
    config: Config
    output_target: str = field(default="bt_target.csv")
    output_genome: str = field(default="bt_genome.csv")
    output_host : str = field(default ="bt_host.csv")
    result_target: BowtieResult = field(init=False)
    result_genome: BowtieResult = field(init=False)
    result_host: BowtieResult  = field(init=False)


    def __post_init__(self) -> None:
        self.create_index()
        self.run_bowtie(index=f"{self.config.bowtie_path}/target/{BASENAME}", fasta_path="potential_primers.fasta", output_path=self.output_target)
        self.run_bowtie(index=f"{self.config.bowtie_path}/s_cerevisiae/s_cerevisiae", fasta_path="potential_primers.fasta", output_path=self.output_genome)
        if os.path.exists("bowtie_index/host/host.1.ebwt") :
            self.run_bowtie(index=f"{self.config.bowtie_path}/host/host", fasta_path="potential_primers.fasta", output_path=self.output_host)
            self.result_host = BowtieResult(result_path=self.output_host)
        self.result_target = BowtieResult(result_path=self.output_target)
        self.result_genome = BowtieResult(result_path=self.output_genome)
        
    

    def run_command(self, cmd: str) -> None:
        process = Popen(args=cmd, stdout=PIPE, stderr=PIPE, shell=True)
        _, stderr = process.communicate()
        if stderr:
            print(stderr.decode("ascii"))
    

    def create_index(self) -> None:
        homepath = self.config.home_path
        btpath = self.config.bowtie_path 
        targetpath = "../../target.fasta" 

        cmd = f"cd {btpath}/target && bowtie-build -f {targetpath} {BASENAME} && cd ../.."
        self.run_command(cmd=cmd)

        if os.path.exists('bowtie_index/s_cerevisiae/s_cerevisiae.1.ebwt') == False : 
            print("Scerevisiae Bowtie_index creation")
            targetpath = "Scerevisia_genome.fasta"
            Basename = "s_cerevisiae"
            cmd = f"cd {btpath}/s_cerevisiae && bowtie-build -f {targetpath} {Basename} && cd ../.."
            self.run_command(cmd=cmd)
        
        if os.path.exists("bowtie_index/host/host.fasta") and (os.path.exists("bowtie_index/host/host.1.ebwt") == False) : 
            targetpath = "host.fasta"
            Basename = "host"
            cmd = f"cd {btpath}/host && bowtie-build -f {targetpath} {Basename} && cd ../.."
            self.run_command(cmd=cmd)
        
    def run_bowtie(self, index:str, fasta_path:str, output_path:str) -> None:
        cmd = f"bowtie -x {index} -a -f {fasta_path} -v 3 > {output_path}"
        self.run_command(cmd=cmd)
