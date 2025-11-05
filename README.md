# PrimEase
This project is a variation of the qTagGer tool originally developed by RGSchindler.
It enables the automatic design of primers for regions of interest, making it easy to integrate into cloning or primer design pipelines. 
The tool simplifies the process of verifying successful recombination of different DNA fragments through PCR 
by allowing users to quickly generate the necessary primers.

## Installation 

### From Github repository
```bash
# Clone the repository
git clone https://github.com/SynDNA-Lab/Primers_designer.git
cd Primers_designer

#installation of the virtual environment and activation
conda env create -f environment.yml --name primers_design
conda activate primers_design

#run the main script
streamlit run main_streamlit.py
```
### Use
This script opens the tool in the browser window.

It offers two main options:
  - Option 1: Drop (or paste) a FASTA sequence, then paste the sequence of interest or two overlapping sequences. The overlapping region will be treated as the sequence of interest.
  - Option 2: Drop an already annotated GenBank file with the desired sequences annotated either as type = "fragment" or as type = "misc_feature" with locus_tag = "fragment".

It is also possible to modify the settings of the Primer3 tool by dropping a modified settings file.
