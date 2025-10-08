#exemple of call python3 genome_annotator.py genome.fasta GATTACA output.gb
import argparse
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os


class GenomeAnnotator:
    def __init__(self, genome_path, search_seq):
        self.genome_path = genome_path
        self.search_seq = search_seq.upper()
        self.output_path = "annotated_genome.gb"
        self.genome_record = None

    def load_genome(self):
        # Load genome assuming FASTA format
        with open(self.genome_path, "r") as handle:
            self.genome_record = SeqIO.read(handle, "fasta")
        print(f"Loaded genome: {self.genome_record.id} with length {len(self.genome_record.seq)}")

    def find_sequence(self):
        genome_seq = str(self.genome_record.seq).upper()
        index = genome_seq.find(self.search_seq)
        if index == -1:
            print("Sequence not found in genome.")
            return None
        print(f"Sequence found at position: {index}")
        return index

    def annotate_sequence(self, start_index):
        end_index = start_index + len(self.search_seq)
        feature = SeqFeature(
            FeatureLocation(start_index, end_index),
            type="misc_feature",
            qualifiers={
                "note": "Matched search sequence",
                "label": "FoundSequence"
            }
        )
        self.genome_record.features.append(feature)

    def write_genbank(self):
        gb_record = SeqRecord(
            Seq(str(self.genome_record.seq)),
            id=self.genome_record.id,
            name=self.genome_record.name,
            description=self.genome_record.description,
            annotations={"molecule_type": "DNA"},
            features=self.genome_record.features
        )

        with open(self.output_path, "w") as output_handle:
            SeqIO.write(gb_record, output_handle, "genbank")
        print(f"GenBank file written to {self.output_path}")

    def run(self):
        self.load_genome()
        pos = self.find_sequence()
        if pos is not None:
            self.annotate_sequence(pos)
            self.write_genbank()


def main():
    parser = argparse.ArgumentParser(description="Annotate a DNA sequence in a genome and export as GenBank.")
    parser.add_argument("genome_file", help="Path to input genome file (FASTA format)")
    parser.add_argument("sequence", help="DNA sequence to search for")
    args = parser.parse_args()

    annotator = GenomeAnnotator(args.genome_file, args.sequence)
    annotator.run()


if __name__ == "__main__":
    main()
