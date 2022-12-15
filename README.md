# qTagGer
## Overview
qTagGer creates a set of diagnostic PCR primers for a user defined motif (i.e. loxPsym) within a user defined input sequence and utilizes [Primer3](https://doi.org/10.1385/1-59259-192-2:365) for primer design. Off-target detection is performed by applying [bowtie](https://doi.org/10.1186/gb-2009-10-3-r25) to detect all sequences with multiple matches in the genome with a defined mismatch value (default value ≤ 4 bp). All off-target and Primer3 config parameters can be customized in the corresponding `config.yaml` and `settings.bak` file. 

## Getting Started
### Dependencies
qTagGer was tested on Ubuntu 22.04 LTS and Python 3.10.4 using the following programs/libraries:
### Programs
| Name | Version |
|------|---------|
| [Bowtie](https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.3.1/) | 1.3.1 |
| [Primer3](https://github.com/primer3-org/primer3) | 2.6.1 |
### Python Libraries
| Name | Version |
|------|---------|
| [pandas](https://pandas.pydata.org/) | ~= 1.4.3 |
| [PyYAML](https://pyyaml.org/) | ~= 6.0 |
| [Biopython](https://biopython.org/) | ~= 1.79 |

## Usage
> python3 qTagGer.py -c config.yml

## Docker

## Configuration
| Name | Type | Example | Description |
|------|------|---------|-------------|
| jobname | str | "NeoChrIII" | Jobtitle for the current run|
| top | int | 5 | Number of primer pairs per site |
| target | str | "neochrIII.gb" | Path to target record |
| offtarget_size_cutoff | int | 10000 | Threshold for max considered PCR offtarget product |
| sponge_value | int | 5 | Number of max allowed primer binding sites|
| home | str | "/home/qTagGer/qTagGer" | Path to parent folder |
| primer3_path | str | "/home/qTagGer/qTagGer/primer3/src/primer3_core" | Path to primer3_core file |
| bowtie_path | str | "/home/qTagGer/qTagGer/bowtie/index" | Path to bowtie index files |



## Publication
Lindeboom, T. A., Sanchez Olmos, M. del C., Schulz, K., Brinkmann, C. K., Ramirez Rojas, A. A., Hochrein, L., & Schindler, D. (2022). L-SCRaMbLE creates large-scale genome rearrangements in synthetic Sc2.0 chromosomes. bioRxiv. https://doi.org/10.1101/2022.12.12.519280


## References
Rozen, S., & Skaletsky, H. (n.d.). Primer3 on the WWW for General Users and for Biologist Programmers. In Bioinformatics Methods and Protocols (pp. 365–386). Humana Press. https://doi.org/10.1385/1-59259-192-2:365

Koressaar, T., & Remm, M. (2007). Enhancements and modifications of primer design program Primer3. In Bioinformatics (Vol. 23, Issue 10, pp. 1289–1291). Oxford University Press (OUP). https://doi.org/10.1093/bioinformatics/btm091

Untergasser, A., Cutcutache, I., Koressaar, T., Ye, J., Faircloth, B. C., Remm, M., & Rozen, S. G. (2012). Primer3—new capabilities and interfaces. In Nucleic Acids Research (Vol. 40, Issue 15, pp. e115–e115). Oxford University Press (OUP). https://doi.org/10.1093/nar/gks596

Langmead, B., Trapnell, C., Pop, M., & Salzberg, S. L. (2009). Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. In Genome Biology (Vol. 10, Issue 3, p. R25). Springer Science and Business Media LLC. https://doi.org/10.1186/gb-2009-10-3-r25

Peter J. A. Cock, Tiago Antao, Jeffrey T. Chang, Brad A. Chapman, Cymon J. Cox, Andrew Dalke, Iddo Friedberg, Thomas Hamelryck, Frank Kauff, Bartek Wilczynski, Michiel J. L. de Hoon, Biopython: freely available Python tools for computational molecular biology and bioinformatics, Bioinformatics, Volume 25, Issue 11, 1 June 2009, Pages 1422–1423, https://doi.org/10.1093/bioinformatics/btp163