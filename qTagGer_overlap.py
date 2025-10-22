
#TAACTTTAAATAATTGGCATTATTTAAAGTTA sequence test
import os
import sys
import shutil
import argparse
import logging
from datetime import datetime

from config import Config
from target import Target
from overlap_annotator import OverlapSequenceAnnotator 
from primer3 import Primer3Interface
from bowtie import BowtieInterface
from offtarget import OfftargetChecker
from selection import  Overlap
from Vizualisation import Visualization
from Bio import SeqIO
from selection import GenbankfileHandling

logging.basicConfig(stream=sys.stdout, encoding='utf-8', level=logging.INFO)

def main(args:argparse.Namespace) -> None:
    genome_path = "test_fragments.gb"
    config = Config()
    logging.info("Successfully loaded config file")
    if "fa" in genome_path or "fasta" in genome_path : 
        annotator = OverlapSequenceAnnotator(search_seq="TGAATATAGGTTATTTTTATTACTACATGCTTGAG",genome_path=genome_path)
        annotator.run()
        logging.info("Successfully performed genome annotation")
        target = Target(target_path="annotated_genome.gb")
        logging.info("Successfully loaded target record")
        selection = Overlap(target=target)
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
        visualization = Visualization(
            sequence = str(selection.regions[0]).upper(), 
            primers_df=check_offtargets.final_primers,
            overlap_seq=annotator.search_seq)
        cleanup(config)
    elif "gb" in genome_path or "genbank" in genome_path : 
        target = Target(target_path=genome_path)
        selection = GenbankfileHandling(genbank_path=genome_path, target=target)
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
        visualization = Visualization(
            sequence = str(selection.regions[0]).upper(), 
            primers_df=check_offtargets.final_primers,
            overlap_seq=selection.overlap)
        cleanup(config)
    else : 
        print("Error")

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
    parser.add_argument("-c", "--config", required=False, type=str, help="Path to the config file (YAML)")
    args = parser.parse_args()
    main(args)
