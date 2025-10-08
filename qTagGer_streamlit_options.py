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
from overlap_annotator import PCROverlapAnnotator, OverlapSequenceAnnotator


logging.basicConfig(stream=sys.stdout, encoding='utf-8', level=logging.INFO)
st.title("Primer Design for Overlap :dna:")
option = st.sidebar.radio(
    "Choose the type of sequence :",
    ("Overlap sequence", "PCR product sequence")
)
def main_overlap() -> None : 
    st.header("Design from :red[_PCR sequences_] ")
    st.subheader("Downloading the configuration")
    uploaded_config = st.file_uploader("Select configuration file")
    if uploaded_config is not None :
        config_path = uploaded_config.name
        current_dir = os.getcwd()
        save_path = os.path.join(current_dir, uploaded_config.name)
        with open(save_path,"wb") as f : 
                f.write(uploaded_config.read()) ##read might be the issue 
                st.success(f'Genome successfully uploaded')

        config = Config(config_path=config_path)
        st.write("Config loaded successfully.")
        logging.info("Successfully loaded config file")
        st.header("Genome and overlap sequence loading")
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
        seq = st.text_input("Enter Overlap sequence")
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
                annotator = OverlapSequenceAnnotator(genome_path=genome_path,search_seq=seq)
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
                # try:
                #     if os.path.exists(config_path):
                #         os.remove(config_path)
                # except Exception as e:
                #         st.error(f"Error during cleanup: {e}")
                # try:
                #     if os.path.exists(genome_path):
                #         os.remove(genome_path)
                # except Exception as e:
                #         st.error(f"Error during cleanup: {e}")
            except Exception as e:
                st.error(f"Failed to load config: {e}")


def main_pcr() -> None:
    st.header("Design of primers from :red[_Overlap sequence_]")
    st.subheader("Downloading the configuration")
    uploaded_config = st.file_uploader("Select configuration file")
    if uploaded_config is not None :
        config_path = uploaded_config.name
        current_dir = os.getcwd()
        save_path = os.path.join(current_dir, uploaded_config.name)
        with open(save_path,"wb") as f : 
                f.write(uploaded_config.read()) ##read might be the issue 
                st.success(f'Genome successfully uploaded')

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
    


if option == "Overlap sequence" :
    main_overlap()
elif option == 'PCR product sequence' :
    main_pcr()
