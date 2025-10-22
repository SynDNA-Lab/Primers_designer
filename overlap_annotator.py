#python3 genome_overlap_annotator.py genome.fasta ACGTACGTAC GTACTGATCG

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord


class PCROverlapAnnotator:
    def __init__(self,genome_path,seq1,seq2):
        self.genome_path = genome_path
        self.seq1 = seq1
        self.seq2 = seq2
        self.output_path = "annotated_genome.gb"
        self.genome_record = None
        self.overlap_seq = None

    def __post_init__(self): 
        self.run()

    def load_genome(self): 
        file_name = self.genome_path.lower()
    
        if any(ext in file_name for ext in ['fa', 'fasta']):
            with open(self.genome_path, "r") as handle: 
                self.genome_record = SeqIO.read(handle, "fasta") 
        elif any(ext in file_name for ext in ['gb', 'genbank']):
            with open(self.genome_path, "r") as handle: 
                self.genome_record = SeqIO.read(handle, "genbank") 

    def find_sequence(self, sequence):
        genome_seq = str(self.genome_record.seq).upper()
        start = genome_seq.find(sequence)
        if start == -1:
            print(f"Sequence not found: {sequence}")
            return None
        end = start + len(sequence)
        print(f"Sequence found: {sequence} at positions {start} - {end}")
        return (start, end)

    def annotate_overlap(self):
        pos1 = self.find_sequence(self.seq1)
        pos2 = self.find_sequence(self.seq2)

        if not pos1 or not pos2:
            print("One or both sequences were not found in the genome.")
            return

        # Calculate overlap
        overlap_start = max(pos1[0], pos2[0])
        overlap_end = min(pos1[1], pos2[1])

        if overlap_start >= overlap_end:
            print("No overlap found between the sequences.")
            return

        print(f"Overlapping region: {overlap_start} - {overlap_end}")
        self.overlap_seq = str(self.genome_record.seq)[overlap_start:overlap_end]
        # Create a feature for the overlapping region
        overlap_feature = SeqFeature(
            FeatureLocation(overlap_start, overlap_end),
            type="misc_feature",
            qualifiers={
                "note": "Overlap between input sequences",
                "label": "OverlapRegion"
            }
        )
        self.genome_record.features.append(overlap_feature)

    def write_genbank(self):
        gb_record = SeqRecord(
            Seq(str(self.genome_record.seq)),
            id=self.genome_record.id,
            name=self.genome_record.name,
            description=self.genome_record.description,
            annotations={"molecule_type": "DNA"},
            features=self.genome_record.features
        )
        with open(self.output_path, "w") as handle:
            SeqIO.write(gb_record, handle, "genbank")
        print(f"Annotated GenBank written to {self.output_path}")

    def run(self):
        self.load_genome()
        self.annotate_overlap()
        self.write_genbank()
        
class OverlapSequenceAnnotator: 

    def __init__(self, genome_path, search_seq, output_path="annotated_genome.gb"): 
        self.genome_path = genome_path 
        self.search_seq = search_seq.upper() 
        self.output_path = output_path 
        self.genome_record = None 
    
    def __post_init__(self): 
        self.run()

    def load_genome(self): 
        file_name = self.genome_path.lower()
    
        if any(ext in file_name for ext in ['fa', 'fasta']):
            with open(self.genome_path, "r") as handle: 
                self.genome_record = SeqIO.read(handle, "fasta") 
        elif any(ext in file_name for ext in ['gb', 'genbank']):
            with open(self.genome_path, "r") as handle: 
                self.genome_record = SeqIO.read(handle, "genbank") 

    def find_sequence(self): 
        genome_seq = str(self.genome_record.seq).upper() 
        index = genome_seq.find(self.search_seq) 
        return index if index != -1 else None 

    def annotate_sequence(self, start_index): 
        end_index = start_index + len(self.search_seq) 
        feature = SeqFeature(FeatureLocation(start_index, end_index),type="misc_feature",qualifiers={"note": "Matched search sequence", "label": "FoundSequence"} ) 
        self.genome_record.features.append(feature) 

    def write_genbank(self): 
        gb_record = SeqRecord(
            Seq(str(self.genome_record.seq)),id=self.genome_record.id,name=self.genome_record.name, 
            description=self.genome_record.description, 
            annotations={"molecule_type": "DNA"}, 
            features=self.genome_record.features) 

        with open(self.output_path, "w") as output_handle: 
            SeqIO.write(gb_record, output_handle, "genbank") 

    def run(self): 
        self.load_genome() 
        pos = self.find_sequence() 
        if pos is not None: 
            self.annotate_sequence(pos) 
            self.write_genbank() 
            return True, pos 
        else: 
            return False, None 

 