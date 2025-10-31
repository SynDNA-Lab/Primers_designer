import os
import sys
import shutil
import argparse
import logging
from datetime import datetime
import streamlit as st 
import tempfile 

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from config import Config
from target import Target
from selection import SequenceOfInterest,GenbankfileHandling
from primer3 import Primer3Interface
from bowtie import BowtieInterface
from offtarget import OfftargetChecker
from sequence_annotator import OverlappingFragmentAnnotator,SequenceAnnotator
from Visualization import Visualization

def handle_file_upload(file, filename: str) -> str:
    """Save uploaded file to disk and return the saved path."""
    path = os.path.join(os.getcwd(), filename)
    with open(path, "wb") as f:
        f.write(file.read())
    return path


def main_sequence_of_interest(uploaded_expected_structure) -> None:
    st.title("Design from :red[_one pasted sequence_]")
    config = Config()

    expected_structure_path = None
    if type(uploaded_expected_structure) == str :
        expected_structure_path = uploaded_expected_structure
    else : 
        expected_structure_path = handle_file_upload(uploaded_expected_structure, uploaded_expected_structure.name)

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

    st.subheader("Paste the sequence of interest")
    seq = st.text_input("Enter Sequence of interest")
    st.write("Enter minimum and maximum size for the PCR product")
    min_max = st.select_slider("Select a range of size for the PCR product",
                                              options=list(range(100, 801)), value=(300, 500),)
    st.write("You selected size between", min_max[0], "and", min_max[1])

    if 'sequence_submitted' not in st.session_state:
        st.session_state['sequence_submitted'] = False

    if st.button("Submit sequence"):
        if expected_structure_path:
            st.session_state['sequence_submitted'] = True
            st.success(f"Sequence of interest sequence submitted: {expected_structure_path}")
        else:
            st.error("Please upload a FASTA sequence file before submitting.")
            st.stop()

    if st.session_state['sequence_submitted']:
        try:

            annotator = SequenceAnnotator(expected_structure_path=expected_structure_path, search_seq=seq)
            annotator.run()
            logging.info("Annotation completed")

            target = Target(target_path="annotated_sequence_expected_structure.gb")
            logging.info("Target loaded")

            selection = SequenceOfInterest(target=target)
            logging.info("Sequence of interest selection done")

            logging.info("Running Primer3...")
            p3interface = Primer3Interface(
                regions=selection.regions[0],
                primer3_path=config.primer3_path,
                specific_file=specific_file,
                size_pcr_product = min_max
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
            name = check_offtargets.final_primers.loc[0,'name']
            name = "_".join(name.split("_")[:-1])

            st.dataframe(check_offtargets.final_primers[[
                                "name",
                                "primer_pair_product_size",
                                "primer_pair_penalty",
                                "primer_left_sequence",
                                "primer_right_sequence",
                                "primer_left_tm",
                                "primer_right_tm",
                                "primer_left_gc_percent",
                                "primer_right_gc_percent"
                                                        ]])

            visualization = Visualization(
                sequence=str(selection.target.record.seq).upper(),
                primers_df=check_offtargets.final_primers,
                sequence_of_interest_seq=annotator.search_seq
            )
            st.image("primers_plot.png")

            cleanup(config)
            for path in [save_path_settings, expected_structure_path]:
                if path and os.path.exists(path):
                    os.remove(path)

        except Exception as e:
            st.error(f"Failed to process DNA fragments: {e}")



def main_fragments(uploaded_expected_structure) -> None:
    st.header("Design of primers from :red[_pasted DNA fragments_]")
    config = Config()

    expected_structure_path = None
    if type(uploaded_expected_structure) == str :
        expected_structure_path = uploaded_expected_structure
    else : 
        expected_structure_path = handle_file_upload(uploaded_expected_structure, uploaded_expected_structure.name)
    
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
    st.write("Enter minimum and maximum size for the PCR product")
    min_max = st.select_slider("Select a range of size for the PCR product", options=list(range(100, 801)), value=(300, 500),)
    st.write("You selected size between", min_max[0], "and", min_max[1])

    if 'sequence_submitted' not in st.session_state:
        st.session_state['sequence_submitted'] = False

    if st.button('Submit sequence'):

        if expected_structure_path:  
            st.session_state['sequence_submitted'] = True
            st.success(f"Expected DNA structure fasta file path submitted: {expected_structure_path}")
        else:
            st.error("Please enter a valid expected DNA structure path before submitting.")
            st.stop() 

    if st.session_state['sequence_submitted']:
        try:
            annotator = OverlappingFragmentAnnotator(expected_structure_path=expected_structure_path,seq_frag1=seq1,seq_frag2=seq2)
            annotator.run()
            logging.info("Successfully performed expected DNA assembly annotation")

            target = Target(target_path="annotated_sequence_expected_structure.gb")
            logging.info("Successfully loaded target record")

            selection = SequenceOfInterest(target=target)
            logging.info("Successfully performed selection")

            logging.info("Running Primer3...")
            p3interface = Primer3Interface(
                regions=selection.regions[0],
                primer3_path=config.primer3_path,
                specific_file=specific_file,
                size_pcr_product = min_max
            )
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

            st.dataframe(check_offtargets.final_primers[[
                                "name",
                                "primer_pair_product_size",
                                "primer_pair_penalty",
                                "primer_left_sequence",
                                "primer_right_sequence",
                                "primer_left_tm",
                                "primer_right_tm",
                                "primer_left_gc_percent",
                                "primer_right_gc_percent"
                                                        ]])
            visualization = Visualization(
                sequence = str(selection.target.record.seq).upper(), 
                primers_df=check_offtargets.final_primers,
               sequence_of_interest_seq=annotator.sequence_of_interest_seq)
            st.image('primers_plot.png')

            cleanup(config)
            for path in [save_path_settings, expected_structure_path]:
                if path and os.path.exists(path):
                    os.remove(path)
        except Exception as e:
            st.error(f"Failed to process DNA fragments : {e}")


def main_annotated_genbank(uploaded_expected_structure) -> None:
    st.header("Design of primers from an :red[_annotated_genbank file_]")
    config = Config()

    expected_structure_path = None
    if uploaded_expected_structure : 
        expected_structure_path = handle_file_upload(uploaded_expected_structure, uploaded_expected_structure.name)
    
    st.subheader("Modification of :green[_primer3 settings_] ?")

    show_uploader = st.checkbox("Set a specific primer3 settings file ?")

    uploaded_settings = None
    save_path_settings = None
    specific_file = False

    st.write("Enter minimum and maximum size for the PCR product")
    min_max = st.select_slider("Select a range of size for the PCR product", options=list(range(100, 801)), value=(300, 500),)
    st.write("You selected size between", min_max[0], "and", min_max[1])
    if show_uploader :
        uploaded_settings = st.file_uploader("Drop your file here", type=["bak"])
        st.write("Primer3 webpage for designing: https://primer3.org/manual.html")

        if uploaded_settings:
            save_path_settings = handle_file_upload(uploaded_settings, "settings_spec.bak")
            st.success("Primer3 settings file uploaded")
            specific_file = True
        else:
            st.warning("No settings file uploaded yet.")

    target = Target(target_path=expected_structure_path)
    logging.info("Successfully loaded target record")

    selection = GenbankfileHandling(genbank_path=expected_structure_path, target=target)
    logging.info("Successfully performed selection")
    sequence_labels = [region.name for region in selection.regions]

    selected_indices = st.multiselect(
        "Select one or more sequences of interest to design primers for:",
        options=list(range(len(sequence_labels))),
        format_func=lambda x: sequence_labels[x],
    )

    if selected_indices:
        st.markdown("### Preview of selected sequences")
        for idx in selected_indices:
            with st.expander(f"{sequence_labels[idx]}"):
                st.code(selection.list_sequence_of_interest[idx], language='text')
                
    if st.button("Run Primer Design on Selected Sequences"):
        if not selected_indices:
            st.warning("Please select at least one sequence.")
        else:
            progress_bar = st.progress(0)
            total = len(selected_indices)

            for count, i in enumerate(selected_indices, start=1):
                logging.info(f"Running Primer3")

                # Run Primer3
                p3interface = Primer3Interface(
                    regions=selection.regions[i],
                    primer3_path=config.primer3_path,
                    specific_file=specific_file,
                    size_pcr_product=min_max
                )
                p3interface.run()

                # Run Bowtie
                logging.info(f"Running Bowtie")
                bowtie = BowtieInterface(config=config)

                # Off-target checking
                logging.info(f"Off-target checking")
                check_offtargets = OfftargetChecker(
                    primer_candidates=p3interface.primer_sites,
                    bowtie_target=bowtie.result_target,
                    bowtie_genome=bowtie.result_genome,
                    config=config
                )

                name = check_offtargets.final_primers.loc[0, 'name']
                name = "_".join(name.split("_")[:-1])

                st.markdown("---")
                st.subheader(f"Results for sequence of interest :blue[_{name}_]")

                st.dataframe(check_offtargets.final_primers[[
                    "name",
                    "primer_pair_product_size",
                    "primer_pair_penalty",
                    "primer_left_sequence",
                    "primer_right_sequence",
                    "primer_left_tm",
                    "primer_right_tm",
                    "primer_left_gc_percent",
                    "primer_right_gc_percent"
                ]])

                # Visualization
                visualization = Visualization(
                    sequence=str(selection.target.record.seq).upper(),
                    primers_df=check_offtargets.final_primers,
                    sequence_of_interest_seq=selection.list_sequence_of_interest[i]
                )

                st.image('primers_plot.png', caption=f"Primers Visualization for {name}")

                # Update progress
                progress_bar.progress(count / total)
                st.markdown("---")

            st.success("Primer design completed for all selected sequences!")


    cleanup(config)
    for path in [save_path_settings, expected_structure_path]:
        if path and os.path.exists(path):
            os.remove(path)
        
def get_sequence_id(region_str: str) -> str:
    for line in region_str.splitlines():
        if line.startswith("SEQUENCE_ID="):
            return line.split("=")[1].strip()
    return "Unknown"

def cleanup(config:Config) -> None:
    now = datetime.now()
    folder_name = f"{config.jobname}_{now.strftime('%Y%m%d_%H%M')}"
    os.makedirs(folder_name, exist_ok=True)

    files_to_move = [
        "bt_genome.csv",
        "bt_target.csv",
        "potential_primers.fasta",
        "primer3_result",
        "settings",
        "qTagGer_Output.csv"
    ]

    for file in files_to_move:
        if os.path.exists(file):
            try:
                shutil.move(file, os.path.join(folder_name, file))
            except Exception as e:
                logging.warning(f"Could not move {file}: {e}")
        else:
            logging.info(f"File {file} not found, skipping cleanup.")

    # Remove temporary files if they exist
    for temp_file in ["annotated_sequence_expected_structure.gb", "target.fasta"]:
        if os.path.exists(temp_file):
            try:
                os.remove(temp_file)
            except Exception as e:
                logging.warning(f"Could not remove {temp_file}: {e}")


logging.basicConfig(stream=sys.stdout, encoding='utf-8', level=logging.INFO)


st.title("Primer Design for Sequence of Interest :dna:")

uploaded_expected_structure = st.file_uploader(
    "Select the expected DNA construct sequence...", 
    type=['fa', 'fasta', 'genbank', 'gb']
)
if uploaded_expected_structure : 
    st.success("Expected DNA assembly successfully uploaded")
file_type = None
sequence_type = None

if 'Continue' not in st.session_state:
    st.session_state['Continue'] = False
if 'Option' not in st.session_state:
    st.session_state['Option'] = None

pasted_expected_structure_construct = st.text_area("Or paste your expected DNA structure sequence here:")

if pasted_expected_structure_construct:
    sequence_type = st.selectbox(
        "Choose the type of sequence to paste :",
        ("Single sequence", "Overlapping DNA fragments")
    )
    if st.button("Continue"):
        st.session_state['Continue'] = True
        st.session_state['Option'] = sequence_type

if st.session_state.get('Continue') and st.session_state.get('Option') and pasted_expected_structure_construct:
    pasted_expected_structure_construct = pasted_expected_structure_construct.replace("\n", "").replace("\r", "").replace(" ", "")
    record = SeqRecord(Seq(pasted_expected_structure_construct), id="pasted_sequence", description="")
    
    with open("pasted_sequence.fasta", "w") as output_handle:
        SeqIO.write(record, output_handle, "fasta")

    if st.session_state['Option'] == "Single sequence":
        main_sequence_of_interest("pasted_sequence.fasta")
    elif st.session_state['Option'] == "Overlapping DNA fragments":
        main_fragments("pasted_sequence.fasta")

elif uploaded_expected_structure:
    file_name = uploaded_expected_structure.name.lower()
    if any(ext in file_name for ext in ['fa', 'fasta']):
        file_type = "fasta"
        sequence_type = st.selectbox(
            "Choose the type of sequence:", 
            ("Single sequence", "Overlapping DNA fragments")
        )

    elif any(ext in file_name for ext in ['gb', 'genbank']):
        file_type = "genbank"
        sequence_type = st.selectbox(
            "Choose the type of sequence:", 
            ("Single sequence", "Overlapping DNA fragments", "Annotated genbank file")
        )

    if st.button("Continue") and sequence_type:
        st.session_state['Continue'] = True
        st.session_state['Option'] = sequence_type

if st.session_state.get('Continue') and st.session_state.get('Option'):
    if file_type == "fasta":
        if st.session_state['Option'] == "Single sequence":
            main_sequence_of_interest(uploaded_expected_structure)
        elif st.session_state['Option'] == "Overlapping DNA fragments":
            main_fragments(uploaded_expected_structure)

    elif file_type == "genbank":
        if st.session_state['Option'] == "Annotated genbank file":
            main_annotated_genbank(uploaded_expected_structure)
        elif st.session_state['Option'] == "Single sequence":
            main_sequence_of_interest(uploaded_expected_structure)
        elif st.session_state['Option'] == "Overlapping DNA fragments":
            main_fragments(uploaded_expected_structure)
