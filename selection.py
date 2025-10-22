import re
from abc import ABC, abstractmethod
from dataclasses import dataclass, field

from target import Target
from roi import ROI
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord


@dataclass
class Selection(ABC):
    target: Target
    regions: list[ROI] = field(init=False)

    @abstractmethod
    def selection(self) -> list[ROI]:
        pass


@dataclass
class RoxP(Selection):
    target: Target
    regions: list[ROI] = field(init=False)

    def __post_init__(self) -> None:
        self.selection()

    def selection(self) -> list[ROI]:
        # roxP recognition site
        pattern = r"TAACTTTAAATAATTGGCATTATTTAAAGTTA"
        it = re.finditer(pattern, self.target.seq)
        targets = [(m.start(0), m.end(0)) for m in it]
        rois = []
        for idx, trg in enumerate(targets):
            #breakpoint()
            rois.append(ROI(
                name=f"target_{trg[0]}",
                target_start = trg[0],
                target_end = trg[1],
                genome_sequence = self.target.seq
            ))
        print(rois[0])
        self.regions = rois

@dataclass
class Overlap(Selection):
    target: Target
    regions: list[ROI] = field(init=False)
    def __post_init__(self) -> None:
        self.selection()

    def selection(self) -> list[ROI]:
        targets = []
        for features in self.target.record.features : 
            if features.type == 'misc_feature' :
                targets.append([features.location.start,features.location.end])
        rois = []
        for idx, trg in enumerate(targets):
            #breakpoint()
            rois.append(ROI(
                name=f"target_{trg[0]}",
                target_start = trg[0],
                target_end = trg[1],
                genome_sequence = self.target.seq
            ))
        print(rois[0])
        self.regions = rois

@dataclass
class GenbankfileHandling(Selection) : 
    target : Target
    genbank_path : str
    genome_record : SeqRecord = field(init=False)
    regions: list[ROI] = field(init=False)
    overlap :str = field(init=False)
    def __post_init__(self) : 
        self.load_genome()
        self.selection()
    
    def load_genome(self):
        with open(self.genbank_path, "r") as handle:
            self.genome_record = SeqIO.read(handle, "genbank")
        print(f"Genome loaded: {self.genome_record.id} ({len(self.genome_record.seq)} bp)")
    
    def selection(self) : 
        overlaps = []
        fragments = []

        for features in self.genome_record.features : 
            if features.type == 'misc_feature' :
                overlaps.append([int(features.location.start),int(features.location.end)])
        
            elif features.type == 'misc_binding' : 
                fragments.append([int(features.location.start),int(features.location.end)])

        rois = []

        if fragments and len(fragments) % 2 != 0 :
            raise ValueError("Number of fragments is odd and it should be even")
        if overlaps :
            for bound in overlaps:
                rois.append(ROI(
                    name=f"target_{bound[0]}",
                    target_start = bound[0],
                    target_end = bound[1],
                    genome_sequence = str(self.genome_record.seq)
                ))
            self.overlap = str(self.genome_record.seq)[bound[0]:bound[1]]
        elif fragments :
            for fragment in range (0,len(fragments),2) :
                overlap_start = max(fragments[fragment][0], fragments[fragment + 1][0])
                overlap_end = min(fragments[fragment][1], fragments[fragment + 1][1])

                if overlap_start >= overlap_end:
                    print("No overlap found between the sequences.")

                rois.append(ROI(
                    name=f"target_{overlap_start}",
                    target_start = overlap_start,
                    target_end = overlap_end,
                    genome_sequence = str(self.genome_record.seq)
                ))
            self.overlap = str(self.genome_record.seq)[overlap_start:overlap_end]
        self.regions = rois
        print(rois)
        print('rois')