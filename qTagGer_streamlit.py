import os
import sys
import shutil
import argparse
import logging
from datetime import datetime
import streamlit as st 
import tempfile 
from config import Config
from target import Target
from selection import Overlap
from primer3 import Primer3Interface
from bowtie import BowtieInterface
from offtarget import OfftargetChecker
from overlap_annotator import PCROverlapAnnotator


logging.basicConfig(stream=sys.stdout, encoding='utf-8', level=logging.INFO)

def main() -> None:
    st.title("Starting the design of primers for Overlapping sequence")
    st.header("Downloading the configuration")
    config_path = st.text_input("Enter configuration path")
    # Use session state to track submit click
    if 'config_submitted' not in st.session_state:
        st.session_state['config_submitted'] = False
    # Handle button click
    if st.button('Submit config'):
        if config_path:  # Path is not empty
            st.session_state['config_submitted'] = True
            st.success(f"Config file path submitted: {config_path}")
        else:
            st.error("Please enter a valid configuration path before submitting.")
            st.stop()  # Stop execution here

    # Only run config loading after successful submit
    if st.session_state['config_submitted']:
        try:
            config = Config(config_path=config_path)
            st.write("Config loaded successfully.")
            logging.info("Successfully loaded config file")
            st.header("Genome and PCR fragment loading")
            st.header("Downloading the genome file")
            uploaded_genome = st.file_uploader("Select a genome fasta file...", type = ['fa','fasta'])
            if uploaded_genome is not None :
                genome_path = uploaded_genome.name
                current_dir = os.getcwd()
                save_path = os.path.join(current_dir, uploaded_genome.name)
                with open(save_path,"wb") as f : 
                    f.write(uploaded_genome.read()) ##read might be the issue 
                    st.success(f'Genome successfully uploaded')
            
            #genome_path = st.text_input("Enter genome fasta file")
            seq1 = st.text_input("Enter PCR fragment 1 sequence")
            seq2 = st.text_input("Enter PCR fragment 2 sequence")
            # Use session state to track submit click
            if 'sequence_submitted' not in st.session_state:
                st.session_state['sequence_submitted'] = False
            # Handle button click
            if st.button('Submit sequence'):
                if genome_path:  # Path is not empty
                    st.session_state['sequence_submitted'] = True
                    st.success(f"Genome fasta file path submitted: {genome_path}")
                else:
                    st.error("Please enter a valid genome path before submitting.")
                    st.stop()  # Stop execution here

            # Only run config loading after successful submit
            if st.session_state['sequence_submitted']:
                try:
                    annotator = PCROverlapAnnotator(genome_path=genome_path,seq1=seq1,seq2=seq2)
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
                    st.dataframe(check_offtargets.final_primers)
                    cleanup(config)
                    try:
                        if os.path.exists(save_path):
                            os.remove(save_path)
                    except Exception as e:
                        st.error(f"Error during cleanup: {e}")
                except Exception as e:
                    st.error(f"Failed to load config: {e}")
        except Exception as e:
            st.error(f"Failed to load config: {e}")
    
        
    


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
    main()
