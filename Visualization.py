import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow
import matplotlib.colors as mcolors
import matplotlib.patches as patches


class Visualization:
    def __init__(self, sequence, primers_df, overlap_seq):
        self.sequence = str(sequence.upper())
        self.primers = primers_df
        self.overlap = overlap_seq.upper()
        self.processed_features = []

        self.primer_colors = list(mcolors.TABLEAU_COLORS.values()) + list(mcolors.CSS4_COLORS.values())

        self.draw_diagram()
        

    def _find_primer_positions(self):
        """Process primers and store start/end/strand info."""
        for i, row in self.primers.iterrows():
            name = row['name']
            fw_seq = row['primer_left_sequence'].upper()
            rv_seq = row['primer_right_sequence'].upper()

            # Forward primer
            fw_start = self.sequence.find(fw_seq)
            if fw_start == -1:
                print(f"⚠️ Forward primer '{name}' not found.")
                continue
            fw_end = fw_start + len(fw_seq)

            # Reverse primer (search using reverse complement)
            rv_seq_rc = str(Seq(rv_seq).reverse_complement())
            rv_end = self.sequence.find(rv_seq_rc)
            if rv_end == -1:
                print(f"⚠️ Reverse primer '{name}' not found.")
                continue
            rv_start = rv_end + len(rv_seq_rc) 

            color = self.primer_colors[i % len(self.primer_colors)]

            self.processed_features.append({
                'name': name + '.fw',
                'start': fw_start,
                'end': fw_end,
                'strand': '+',
                'color': color,
                'index': i
            })
            self.processed_features.append({
                'name': name + '.rv',
                'start': rv_start,
                'end': rv_end,
                'strand': '-',
                'color': color,
                'index': i
            })

    def _find_overlap_position(self):
        """Find and store overlap position"""
        start = self.sequence.find(self.overlap)
        if start == -1:
            print("⚠️ Overlap not found.")
            return None
        end = start + len(self.overlap)
        return (start, end)

    def draw_diagram(self, output_path="primers_plot.png"):
        self._find_primer_positions()
        overlap_pos = self._find_overlap_position()
        overlap_center = (overlap_pos[0] + overlap_pos[1]) // 2

        # Set up plot
        fig, ax = plt.subplots(figsize=(12, 4))
        seq_length = len(self.sequence)

        # Find overlap center for coordinate shift
        if overlap_pos:
            overlap_center = (overlap_pos[0] + overlap_pos[1]) // 2
        else:
            overlap_center = seq_length // 2  # fallback if overlap not found

        # Genome baseline (centered around 0)
        ax.hlines(y=0, xmin=-overlap_center, xmax=seq_length - overlap_center, color='black', linewidth=2)


        # Draw primers
        for feature in self.processed_features:
            y = 1 + feature['index'] * 0.6  # Stack primers vertically
            dx = feature['end'] - feature['start']
            arrow_props = dict(
                width=0.2,
                head_width=0.3,
                head_length=3,
                length_includes_head=True,
                color=feature['color']
            )

            ax.add_patch(FancyArrow(
                feature['start'] - overlap_center, y,
                dx if feature['strand'] == '+' else -abs(dx),
                0,
                **arrow_props
            ))

            # Label
            ax.text(feature['start'] + dx / 2 - overlap_center, y + 0.25, feature['name'],
                    ha='center', fontsize=8)

        # Draw overlap
        if overlap_pos:
            max_height = 1 + len(self.primers) * 0.6
            label_offset = 0.4
            overlap_x = overlap_pos[0]
            overlap_width = overlap_pos[1] - overlap_pos[0]
            rect = patches.Rectangle(
                (overlap_x - overlap_center, 0),        # (x, y)
                overlap_width,         # width
                max_height,            # height
                linewidth=0,
                edgecolor=None,
                facecolor='red',
                alpha=0.5,
                label="Overlap"
            )
            ax.add_patch(rect)
            # Optional: add label centered inside the box
            ax.text(
                overlap_x + overlap_width / 2 - overlap_center,
                max_height +label_offset,
                "Overlap",
                ha='center',
                va='top',
                fontsize=9,
                color='red'
            )
        # Format
        if overlap_pos:
            overlap_center = (overlap_pos[0] + overlap_pos[1]) // 2
            window_size = seq_length // 2  # or any other value like 500
            x_start = max(0, overlap_center - window_size)
            x_end = min(seq_length, overlap_center + window_size)
            ax.set_xlim(-window_size, window_size)
            ax.set_xlabel("Position relative to overlap center")
        else:
            ax.set_xlim(-window_size, window_size)
            ax.set_xlabel("Position relative to overlap center")

        ax.set_ylim(-1, 1 + len(self.primers) * 0.8)
        ax.set_yticks([])
        ax.set_title("Primer Binding Sites and Overlap Region")
        plt.tight_layout()

        # Save and/or show
        plt.savefig(output_path, dpi=300)
        plt.show()

        return output_path