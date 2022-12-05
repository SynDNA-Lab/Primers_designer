import re
import sys
import shutil
import logging
import pathlib
import warnings

from Bio import SeqIO
from Bio import SeqRecord
from Bio import BiopythonWarning
from primer import Regions

warnings.simplefilter('ignore', BiopythonWarning)



# To be implemented:
#   - The 400 in the section below should be replaced by the optimal product size in primer3
#   - list of gene names to design amplicons for
#   - regex expression for string sequence

def select_amplicons():
    target_sequence = str(record.seq)
    pattern = r"TAACTTTAAATAATTGGCATTATTTAAAGTTA"
    #pattern = r"ATAACTTCGTATA[a-z,A-Z]{8}TATACGAAGTTAT"
    it = re.finditer(pattern, target_sequence.upper())
    targets = [(m.start(0), m.end(0)) for m in it]
    
    shutil.copyfile("settings.bak", "settings")
    with open("settings", "a") as file:
        for site in targets:
            start = min(site[0], 400)
            seq = target_sequence[site[0]-start:site[1]+400]
            file.write(f"SEQUENCE_ID=recomb_site_{site[0]}\n")
            file.write(f"SEQUENCE_TEMPLATE={seq}\n")
            file.write(f"SEQUENCE_TARGET={start},{site[1]-site[0]}\n")
            file.write("=\n")


def load(target_path:str) -> SeqRecord.SeqRecord:
    logging.info("Startet the target.load() procedure")
    global record 

    # Check format and load the record
    file_extension = pathlib.Path(target_path).suffix
    if file_extension in ["gb", "gbk"]:
        logging.info("Genbank file detected")
        records = []
        for seq_record in SeqIO.parse(target_path, "genbank"):
            records.append(seq_record)

    elif file_extension == "dna":
        logging.info("Snapgene file detected")
        records = []
        for seq_record in SeqIO.parse(target_path, "snapgene"):
            records.append(seq_record)

    elif file_extension in ["fa", "mpfa", "fna", "fsa", "fasta"]:
        logging.info("Fasta file detected")
        records = []
        for seq_record in SeqIO.parse(target_path, "fasta"):
            records.append(seq_record)

    else:
        logging.error(f"Unable to load the file {target_path}, unknown format {file_extension}") 
        sys.exit()
    
    logging.info("Successfully loaded the target sequence record")
    record = records[0]
    return record