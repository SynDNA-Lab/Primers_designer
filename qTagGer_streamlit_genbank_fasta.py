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
from selection import Overlap,GenbankfileHandling
from primer3 import Primer3Interface
from bowtie import BowtieInterface
from offtarget import OfftargetChecker
from overlap_annotator import PCROverlapAnnotator, OverlapSequenceAnnotator
from Visualization import Visualization

def handle_file_upload(file, filename: str) -> str:
    """Save uploaded file to disk and return the saved path."""
    path = os.path.join(os.getcwd(), filename)
    with open(path, "wb") as f:
        f.write(file.read())
    return path


def main_overlap(uploaded_genome) -> None:
    st.title("Design from a :red[_pasted overlap sequence_]")
    config = Config()

    genome_path = None
    if uploaded_genome:
        genome_path = handle_file_upload(uploaded_genome, uploaded_genome.name)

    st.subheader("Modification of :green[_primer3 settings_] ?")
    show_uploader = st.checkbox("Set a specific primer3 settings file")

    uploaded_settings = None
    save_path_settings = None
    specific_file = False

    if show_uploader:
        uploaded_settings = st.file_uploader("Drop your file here", type=["bak"])
        st.write("Primer3 webpage for designing: https://primer3.org/manual.html")

        if uploaded_settings:
            save_path_settings = handle_file_upload(uploaded_settings, "settings_spec.bak")
            st.success("Primer3 settings file uploaded")
            specific_file = True
        else:
            st.warning("No settings file uploaded yet.")

    st.subheader("Paste the overlap sequence")
    seq = st.text_input("Enter Overlap sequence")

    if 'sequence_submitted' not in st.session_state:
        st.session_state['sequence_submitted'] = False

    if st.button("Submit sequence"):
        if genome_path:
            st.session_state['sequence_submitted'] = True
            st.success(f"Overlap sequence submitted: {genome_path}")
        else:
            st.error("Please upload a valid genome FASTA file before submitting.")
            st.stop()

    if st.session_state['sequence_submitted']:
        try:
            annotator = OverlapSequenceAnnotator(genome_path=genome_path, search_seq=seq)
            annotator.run()
            logging.info("Annotation completed")

            target = Target(target_path="annotated_genome.gb")
            logging.info("Target loaded")

            selection = Overlap(target=target)
            logging.info("Overlap selection done")

            logging.info("Running Primer3...")
            p3interface = Primer3Interface(
                selection=selection,
                primer3_path=config.primer3_path,
                specific_file=specific_file
            )
            p3interface.run()
            logging.info("Primer3 finished")

            logging.info("Running Bowtie...")
            bowtie = BowtieInterface(config=config)
            logging.info("Bowtie alignment done")

            logging.info("Checking off-targets...")
            check_offtargets = OfftargetChecker(
                primer_candidates=p3interface.primer_sites,
                bowtie_target=bowtie.result_target,
                bowtie_genome=bowtie.result_genome,
                config=config
            )

            st.dataframe(check_offtargets.final_primers)

            visualization = Visualization(
                sequence=str(selection.regions[0]).upper(),
                primers_df=check_offtargets.final_primers,
                overlap_seq=annotator.search_seq
            )
            st.image("primers_plot.png")

            cleanup(config)
            for path in [save_path_settings, genome_path]:
                if path and os.path.exists(path):
                    os.remove(path)

        except Exception as e:
            st.error(f"Failed to process DNA fragments: {e}")



def main_fragments(uploaded_genome) -> None:
    st.header("Design of primers from :red[_pasted DNA fragments_]")
    config = Config()

    genome_path = None
    if uploaded_genome : 
        genome_path = handle_file_upload(uploaded_genome, uploaded_genome.name)
    
    st.subheader("Modification of :green[_primer3 settings_] ?")
    show_uploader = st.checkbox("Set a specific primer3 settings file")

    uploaded_settings = None
    save_path_settings = None
    specific_file = False

    if show_uploader:
        uploaded_settings = st.file_uploader("Drop your file here", type=["bak"])
        st.write("Primer3 webpage for designing: https://primer3.org/manual.html")

        if uploaded_settings:
            save_path_settings = handle_file_upload(uploaded_settings, "settings_spec.bak")
            st.success("Primer3 settings file uploaded")
            specific_file = True
        else:
            st.warning("No settings file uploaded yet.")
    
    st.subheader("Paste the fragments's sequence")
    seq1 = st.text_input("Enter fragment 1 sequence")
    seq2 = st.text_input("Enter fragment 2 sequence")

    if 'sequence_submitted' not in st.session_state:
        st.session_state['sequence_submitted'] = False

    if st.button('Submit sequence'):

        if genome_path:  
            st.session_state['sequence_submitted'] = True
            st.success(f"Genome fasta file path submitted: {genome_path}")
        else:
            st.error("Please enter a valid genome path before submitting.")
            st.stop() 

    if st.session_state['sequence_submitted']:
        try:
            annotator = PCROverlapAnnotator(genome_path=genome_path,seq1=seq1,seq2=seq2)
            annotator.run()
            logging.info("Successfully performed expected DNA assembly annotation")

            target = Target(target_path="annotated_genome.gb")
            logging.info("Successfully loaded target record")

            selection = Overlap(target=target)
            logging.info("Successfully performed selection")

            logging.info("Running Primer3...")
            p3interface = Primer3Interface(selection=selection, primer3_path=config.primer3_path,specific_file=specific_file)
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

            visualization = Visualization(
                sequence = str(selection.regions[0]).upper(), 
                primers_df=check_offtargets.final_primers,
                overlap_seq=annotator.overlap_seq)
            st.image('primers_plot.png')

            cleanup(config)
            for path in [save_path_settings, genome_path]:
                if path and os.path.exists(path):
                    os.remove(path)
        except Exception as e:
            st.error(f"Failed to load config: {e}")


def main_annotated_genbank(uploaded_genome) -> None:
    st.header("Design of primers from an :red[_annotated_genbank file_]")
    config = Config()

    genome_path = None
    if uploaded_genome : 
        genome_path = handle_file_upload(uploaded_genome, uploaded_genome.name)
    
    st.subheader("Modification of :green[_primer3 settings_] ?")
    show_uploader = st.selectbox("Set a specific primer3 settings file ?",("Yes", "No")
        )

    uploaded_settings = None
    save_path_settings = None
    specific_file = False
    settings_file_uploaded = False
    if show_uploader == "Yes":
        uploaded_settings = st.file_uploader("Drop your file here", type=["bak"])
        st.write("Primer3 webpage for designing: https://primer3.org/manual.html")

        if uploaded_settings:
            save_path_settings = handle_file_upload(uploaded_settings, "settings_spec.bak")
            st.success("Primer3 settings file uploaded")
            specific_file = True
            settings_file_uploaded = True
        else:
            st.warning("No settings file uploaded yet.")
    if show_uploader == "No" or settings_file_uploaded == True : 
        target = Target(target_path=genome_path)
        logging.info("Successfully loaded target record")

        selection = GenbankfileHandling(genbank_path=genome_path, target=target)
        logging.info("Successfully performed selection")

        logging.info("Running Primer3...")
        p3interface = Primer3Interface(selection=selection, primer3_path=config.primer3_path,specific_file=specific_file)
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

        visualization = Visualization(
            sequence = str(selection.regions[0]).upper(), 
            primers_df=check_offtargets.final_primers,
            overlap_seq=selection.overlap)
        st.image('primers_plot.png')

        cleanup(config)
        for path in [save_path_settings, genome_path]:
            if path and os.path.exists(path):
                os.remove(path)
        


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


logging.basicConfig(stream=sys.stdout, encoding='utf-8', level=logging.INFO)


st.title("Primer Design for Overlap :dna:")

uploaded_genome = st.file_uploader(
    "Select the expected DNA construct sequence...", 
    type=['fa', 'fasta', 'genbank', 'gb']
)
if uploaded_genome : 
    st.success("Expected DNA assembly successfully uploaded")
file_type = None
sequence_type = None

if 'Continue' not in st.session_state:
    st.session_state['Continue'] = False
if 'Option' not in st.session_state:
    st.session_state['Option'] = None

if uploaded_genome:
    file_name = uploaded_genome.name.lower()
    
    if any(ext in file_name for ext in ['fa', 'fasta']):
        file_type = "fasta"
        sequence_type = st.selectbox(
            "Choose the type of sequence:", 
            ("Overlap", "Multiple DNA fragments")
        )

    elif any(ext in file_name for ext in ['gb', 'genbank']):
        file_type = "genbank"
        sequence_type = st.selectbox(
            "Choose the type of sequence:", 
            ("Overlap", "Multiple DNA fragments", "Annotated genbank file")
        )

    if st.button("Continue") and sequence_type:
        st.session_state['Continue'] = True
        st.session_state['Option'] = sequence_type

if st.session_state.get('Continue') and st.session_state.get('Option'):
    if file_type == "fasta":
        if st.session_state['Option'] == "Overlap":
            main_overlap(uploaded_genome)
        elif st.session_state['Option'] == "Multiple DNA fragments":
            main_fragments(uploaded_genome)

    elif file_type == "genbank":
        if st.session_state['Option'] == "Annotated genbank file":
            main_annotated_genbank(uploaded_genome)
        elif st.session_state['Option'] == "Overlap":
            main_overlap(uploaded_genome)
        elif st.session_state['Option'] == "Multiple DNA fragments":
            main_fragments(uploaded_genome)
