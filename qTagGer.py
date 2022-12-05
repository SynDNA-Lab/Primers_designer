# Main tool to generate primer designs 
import os
import sys
import argparse
import logging

import config
import primer3
import target
import offtarget
from primer import Candidates


logging.basicConfig(filename='report.log', encoding='utf-8', level=logging.DEBUG)

def exists(filepath:str) -> None:
    try:
        os.path.exists(filepath)
    except Exception as e:
        logging.error(f"Unable to find {filepath}")
        logging.error(e)
        sys.exit()


if __name__ == "__main__":
    # Parser
    parser = argparse.ArgumentParser(description="PrimerDesigner - An Automtic Amplicon Primer Designer")
    parser.add_argument("-c", "--config", required=True, type=str, help="Path to the config file (YAML)")
    args = parser.parse_args()

    exists(args.config)
    config.load(args.config)

    exists(args.target)
    logging.info(f"Working on the file {args.target}")
    target.load(args.target)
    target.select_amplicons(args.target)
    primer3.run(result_path="primer3_result.txt")
    primer3.parse(result_path="primer3_result.txt")
    primer3.generate_fasta()
    offtarget.filter() #not complete