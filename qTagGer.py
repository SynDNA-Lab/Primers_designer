import argparse
import logging
import sys
import shutil
import os
from datetime import datetime

from config import Config
from target import Target
from selection import RoxP
from primer3 import Primer3Interface
from bowtie import BowtieInterface
from offtarget import OfftargetChecker



logging.basicConfig(stream=sys.stdout, encoding='utf-8', level=logging.INFO)


def main(args:argparse.Namespace) -> None:
    config = Config(config_path=args.config)
    logging.info("Successfully loaded config file")
    target = Target(target_path=config.target_path)
    logging.info("Successfully loaded target record")
    selection = RoxP(target=target)
    logging.info("Successfully performed selection")
    logging.info("Running Primer3...")
    p3interface = Primer3Interface(selection=selection, primer3_path=config.primer3_path)
    p3interface.run()
    logging.info("Primer3 finished")
    logging.info("Running Bowtie")
    bowtie = BowtieInterface(config=config)
    logging.info("Successfully performed bowtie alignments")
    logging.info("Offtarget checking")
    check_offtargets = OfftargetChecker(
        primer_candidates = p3interface.primer_sites,
        bowtie_target = bowtie.result_target,
        bowtie_genome = bowtie.result_genome,
        config=config
    )
    cleanup(config)



def cleanup(config:Config) -> None:
    now = datetime.now() 
    folder_name = f"{config.jobname}_{now.strftime('%Y%m%d_%H%M')}"
    os.makedirs(folder_name, exist_ok=True)
    shutil.move("bt_genome.csv", f"{folder_name}/bt_genome.csv")
    shutil.move("bt_target.csv", f"{folder_name}/bt_target.csv")
    shutil.move("potential_primers.fasta", f"{folder_name}/potential_primers.fasta")
    shutil.move("primer3_result", f"{folder_name}/primer3_result")
    shutil.move("settings", f"{folder_name}/settings")
    shutil.move("qTagGer_Output.csv", f"{folder_name}/qTagGer_Output.csv")
    os.remove("target.fasta")

    

if __name__ == "__main__":
    # Parser
    parser = argparse.ArgumentParser(description="qTagGer - An Automtic Amplicon Primer Designer")
    parser.add_argument("-c", "--config", required=True, type=str, help="Path to the config file (YAML)")
    args = parser.parse_args()
    main(args)
