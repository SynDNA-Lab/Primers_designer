import yaml
import logging
from dataclasses import dataclass, field



@dataclass
class Config:
    config_path: str 

    # General Settings:
    jobname: str = field(init=False)
    target_path: str = field(init=False)
    top: str = field(init=False)

    # Primer Settings:
    offtarget_size_cutoff: int = field(init=False, default=10_000) 
    sponge_value: int = field(init=False, default=5) # number of maximal offtarget binding sites across loci

    # Paths:
    home_path: str = field(init=False, default="/home/Primer/PrimerDesigner")
    primer3_path: str = field(init=False, default="/home/Primer/PrimerDesigner/primer3/src/primer3_core")
    bowtie_path: str = field(init=False, default="/home/Primer/PrimerDesigner/bowtie/index")


    def __post_init__(self) -> None:
        self.load()


    def load(self) -> None:
        logging.info("Loading the configuraiton file ...")
        with open(self.config_path) as file:
            config = yaml.load(file, Loader=yaml.FullLoader)

        self.jobname = config["jobname"]
        self.target_path = config["target"]
        self.top = config["top"]
        
        self.offtarget_size_cutoff = config["offtarget_size_cutoff"] 
        self.sponge_value = config["sponge_value"]

        self.home_path = config["home"]
        self.primer3_path = config["primer3_path"]
        self.bowtie_path = config["bowtie_path"]
