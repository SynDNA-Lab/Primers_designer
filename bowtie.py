import config
import pandas as pd 
from dataclasses import dataclass, field
from subprocess import Popen, PIPE



@dataclass
class BowtieResult:
    result_path: str
    result: pd.DataFrame = field(init=False)

    def __init__(self) -> None:
        self.parse_data(self.result_path)


    def split_data(lst:list[str], index:int) -> list[str]:
        return [l.split("_")[index] for l in lst]


    def parse_data(self) -> None:
        df = pd.read_csv(self.result_path, sep="\t", header=None)
        df.columns =  ["name", "strand", "reference", "start", "sequence", "quality", "instances", "mismatch_descriptor"]
        df = df.drop_duplicates()
        df["id"] = self.split_name(df.name, -2)
        df["orientation"] = self.split_name(df.name, -1)
        
        self.result = df



@dataclass
class BowtieInterface:
    target_path: str
    output_target: str = field(default="bt_target.csv")
    output_genome: str = field(default="bt_genome.csv")
    result_target: BowtieResult = field(init=False)
    result_genome: BowtieResult = field(init=False)


    def __init__(self) -> None:
        self.create_index(ref_file=, org_name=)
        self.run_bowtie(input_path=, output_path=)
        self.run_bowtie(input_path=, output_path=)
        self.result_target = BowtieResult(result_path=)
        self.result_genome = BowtieResult(result_path=)


    def run_command(self, cmd: str) -> None:
        process = Popen(args=cmd, stdout=PIPE, stderr=PIPE, shell=True)
        _, stderr = process.communicate()
        if stderr:
            print(stderr.decode("ascii"))
    

    def create_index(self, ref_file: str, org_name:str) -> None:
        cmd = f"cd {config.homepath}/bowtie/index/ && bowtie-build -f {ref_file} {org_name} && cd {config.homepath}/input/script"
        self.run_command(cmd=cmd)


    def run_bowtie(self, input_path:str, output_path:str) -> None:
        cmd = f"bowtie -x {index} -a -f {input_path} -v 3 > {output_path}"
        self.run_command(cmd=cmd)
