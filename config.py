import yaml
import logging
from dataclasses import dataclass, field



@dataclass
class Config:
    config_path: str 

    # General Settings:
    jobname: str = field(init=False)
    target_path: str = field(init=False)
    chromosome_id: int = field(init=False)
    output_path: str = field(init=False)

    # Primer Settings:
    max_loci_offtargets: int = field(init=False, default=3)
    offtarget_size_cutoff: int = field(init=False, default=10_000) #prev. max_pcr_offtarget
    sp_value: int = field(init=False, default=5) # number of maximal offtarget binding sites across loci

    # Paths:
    home_path: str = field(init=False, default="/home/Primer/PrimerDesigner")
    primer3_path: str = field(init=False, default="/home/Primer/PrimerDesigner/primer3/src/primer3_core")
    bowtie_path: str = field(init=False, default="/home/Primer/PrimerDesigner/bowtie/index")


    def __post_init__(self) -> None:
        self.load()


    def load(self) -> None:
        logging.debug("Started config.py load procedure")
        with open(self.config_path) as file:
            config = yaml.load(file, Loader=yaml.FullLoader)

        self.jobname = config["jobname"]
        self.target_path = config["target"]
        self.subject_path = config["chromosome"]
        self.output_path = config["pcr_size_target"] # fix this here
        # check if these even exist
        self.max_loci_offtargets = ...
        self.offtarget_size_cutoff = ...
        self.sp_value = ...
        self.home_path = ...
        self.primer3_path = ...
        self.bowtie_path = ...
