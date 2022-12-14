import argparse
import logging

from config import Config
from target import Target
from selection import RoxP
from primer3 import Primer3Interface
from bowtie import BowtieInterface
from offtarget import OfftargetChecker



logging.basicConfig(filename='report.log', encoding='utf-8', level=logging.DEBUG)


def main(args:argparse.Namespace) -> None:
    config = Config(config_path=args.config)
    target = Target(target_path=config.target_path)
    selection = RoxP(target=target)
    p3interface = Primer3Interface(selection=selection, primer3_path=config.primer3_path)
    p3interface.run()
    #BowtieInterface()
    bowtie = BowtieInterface(config=config)

    check_offtargets = OfftargetChecker(
        primer_candidates = p3interface.primer_sites,
        bowtie_target = bowtie.result_target,
        bowtie_genome = bowtie.result_genome,
        sponge_value = config.sponge_value,
        max_pcr_size = config.offtarget_size_cutoff
    )

    '''
    target_path: str
    config: Config
    output_target: str = field(default="bt_target.csv")
    output_genome: str = field(default="bt_genome.csv")
    '''
    
    # Cleanup routine


def cleanup() -> None:
    '''
    Create a report folder
    move all the generated files into this folder
    save report.txt into this folder 
    '''
    raise NotImplementedError


if __name__ == "__main__":
    # Parser
    parser = argparse.ArgumentParser(description="qTagGer - An Automtic Amplicon Primer Designer")
    parser.add_argument("-c", "--config", required=True, type=str, help="Path to the config file (YAML)")
    args = parser.parse_args()
    main(args)
