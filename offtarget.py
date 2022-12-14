import pandas as pd 
import itertools
from dataclasses import dataclass, field

from bowtie import BowtieResult



@dataclass
class OfftargetChecker:
    primer_candidates: pd.DataFrame
    bowtie_target: BowtieResult 
    bowtie_genome: BowtieResult 
    sponge_value: int
    max_pcr_size: int
    offtarget_ids: list[int] = field(init=False)


    def __init__(self) -> None:
        pass


    def get_sponges(self, data_frame:pd.DataFrame) -> list[str]:
        grouped = data_frame.groupby(["id", "orientation"])
        sponge_df = grouped.filter(lambda x: len(x)>= self.sponge_value)
        return sponge_df.id.unique()


    def filter_offtarget(self, data_frame:pd.DataFrame) -> list[str]:
        grouped = data_frame.groupby(["id", "reference"])
        df_filtered = grouped.filter(lambda x: len(x)>=2)
        return_list = []
        for name, group in df_filtered.groupby(["id", "reference"]):
            if len(group.strand.unique())<2:
                continue
            fwd = group[group.strand=="-"].start.unique()
            rev = group[group.strand=="+"].start.unique()
            
            for f, r in itertools.product(fwd, rev):
                difference = r-f
                if (difference > 0) and (difference <= self.max_pcr_size):
                    return_list.append(group.id.values[0])

        return return_list