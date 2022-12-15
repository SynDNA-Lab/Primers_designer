import re
from abc import ABC, abstractmethod
from dataclasses import dataclass, field

from target import Target
from roi import ROI



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
        self.regions = rois