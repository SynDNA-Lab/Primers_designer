# Primer Designer Region of Interest
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
streamlit run qTagGer_streamlit_genbank_fasta.py
```
