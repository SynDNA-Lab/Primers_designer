# classes for both regions of interest as well as primer3 outputs
import attr

from typing import List
from dataclasses import dataclass, field

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation



# Regions of iterest used for generatign the primer3 settings file
@attr.s
class Regions:
    all = []
    refgenome: SeqRecord
    feature: SeqFeature
    name: str = field(init=False)
    sequence: str = field(init=False)
    genstart: int = field(init=False)
    genend: int = field(init=False)
    genstrand: int = field(init=False)
    roistart: int = field(init=False)
    roiend: int = field(init=False)
    excludesites: List[int] = field(default_factory=list)
    offset:int = 400

    def generate_feature(self):
        return SeqFeature(
            FeatureLocation(
                self.feature.location.start - self.offset,
                self.feature.location.end + self.offset,
            )
        )

    def get_sequence(self):
        f = self.generate_feature()
        self.sequence = f.extract(self.refgenome.seq)

    # Replace that with something more modular
    def region_oi(self):
        f = self.generate_feature()
        regionrange = set(range(f.location.start, f.location.end))

        overlaps = []
        for feat in self.refgenome.features[1:]:
            if feat.type == "misc_feature":
                if "rox" in feat.qualifiers["label"][0]:
                    continue
            elif (feat.type == "tRNA") and (len(feat) < 300):
                featurerange = set(range(feat.location.start, feat.location.end))
                overlap = regionrange.intersection(featurerange)
                if overlap:
                    overlaps.append([min(overlap), max(overlap)])

        return overlaps

    def generate(self, filename):
        #generate 10 instead of 5
        with open(filename, "a") as file:
            file.write(f"SEQUENCE_ID={self.name}\n")
            file.write(f"SEQUENCE_TEMPLATE={str(self.sequence)}\n")
            #file.write(f"SEQUENCE_TARGET={self.genstart},{self.genend-self.genstart}\n")
            file.write(f"SEQUENCE_TARGET={self.roistart},{self.roiend-self.roistart}\n")

            if self.excludesites:
                exclstr = " ".join([f"{x[0]},{x[1]}" for x in self.excludesites])
                file.write(f"SEQUENCE_EXCLUDED_REGION={exclstr}\n")
            file.write("=\n")

    def __post_init__(self):
        self.get_sequence()
        self.genstart = self.feature.location.start
        self.genend = self.feature.location.end
        self.genstrand = self.feature.location.strand
        self.roistart = self.offset
        self.roiend = self.roistart + len(self.feature)
        self.name = f"primer_region_{self.genstart}"

        excludebuffer = self.region_oi()
        for excld in excludebuffer:
            self.excludesites.append(
                [excld[0] - (self.genstart - 400), excld[1] - excld[0]]
                #[excld[0], excld[1]-excld[0]]
            )
        self.all.append(self)


# Primer3 Output
@dataclass
class Candidates:
    all = []
    target:str 
    pair_penalty:float 
    left_penalty:float 
    right_penalty:float 
    left_sequence:str 
    right_sequence:str 
    left:str 
    right:str 
    left_tm:float 
    right_tm:float 
    left_gc_percent:float 
    right_gc_percent:float 
    left_self_any_th:float 
    right_self_any_th:float 
    left_self_end_th:float 
    right_self_end_th:float 
    left_hairpin_th:float 
    right_hairpin_th:float 
    left_end_stability:float 
    right_end_stability:float 
    pair_compl_any_th:float 
    pair_compl_end_th:float 
    pair_product_size:int 
    pair_product_tm:float 

    def __init__(self, data):
        self.target = data[0]  
        self.pair_penalty = data[1]  
        self.left_penalty = data[2]  
        self.right_penalty = data[3]  
        self.left_sequence = data[4]  
        self.right_sequence = data[5]  
        self.left = data[6]  
        self.right = data[7]  
        self.left_tm = data[8]  
        self.right_tm = data[9]  
        self.left_gc_percent = data[10]  
        self.right_gc_percent = data[11] 
        self.left_self_any_th = data[12] 
        self.right_self_any_th = data[13] 
        self.left_self_end_th = data[14] 
        self.right_self_end_th = data[15] 
        self.left_hairpin_th = data[16]  
        self.right_hairpin_th = data[17] 
        self.left_end_stability = data[18] 
        self.right_end_stability = data[19] 
        self.pair_compl_any_th = data[20] 
        self.pair_compl_end_th = data[21] 
        self.pair_product_size = data[22]  
        self.pair_product_tm = data[23]  

        self.all.append(self)