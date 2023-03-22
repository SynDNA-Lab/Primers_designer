# qTagGer
## Overview
qTagGer creates a set of diagnostic PCR primers for a user defined motif (i.e. loxPsym) within a user defined input sequence and utilizes [Primer3](https://doi.org/10.1385/1-59259-192-2:365) for primer design. Off-target detection is performed by applying [bowtie](https://doi.org/10.1186/gb-2009-10-3-r25) to detect all sequences with multiple matches in the genome with a defined mismatch value (default value ≤ 4 bp). All off-target and Primer3 config parameters can be customized in the corresponding `config.yaml` and `settings.bak` file. 

## Getting Started

### Dependencies
qTagGer relies on the following software packages and Python libraries:   <br>
<br>
**Programs:**
| Name | Version |
|------|---------|
| [Bowtie](https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.3.1/) | 1.3.1 |
| [Primer3](https://github.com/primer3-org/primer3) | 2.6.1 |


**Python libraries:**

| Name | Version |
|------|---------|
| [pandas](https://pandas.pydata.org/) | ~= 1.4.3 |
| [PyYAML](https://pyyaml.org/) | ~= 6.0 |
| [Biopython](https://biopython.org/) | ~= 1.79 |   

**Tested**   
qTagGer has been tested on Ubuntu 22.04 LTS and Python 3.10.4 with the dependencies mentioned above.

### Docker
Due to the number of dependencies, it is recommended to run qTagGer via the Docker Image we've created. While Docker can be run from the command line interface, we recommend installing the [Docker Desktop](https://www.docker.com/products/docker-desktop/) app.
If you're new to Docker follow the quick start guide that can be found [here](https://docs.docker.com/desktop/get-started/).
The qTagGer Docker image can be downloaded from [**here**](https://example.com/). To load the image please follow the instructions from the official [Docker documentation](https://docs.docker.com/engine/reference/commandline/load/).

**Creating the Docker Container**   
When creating the container, it is important to set the following optional settings:   

**Volumes:**
> **Host path**: File path to an empty folder on your host machine.

> **Container path**: /home/L_SCRaMbLE/qTagGer/input

Providing these settings will allow you to interact with the container by parsing files into the host folder. 

**Running qTagGer in Docker**
1. Start and connect to the container.
2. Switch the shell from `sh` to `bash` by entering and executing the `bash` command.
3. Navigate to the script directory by executing 
```cd home/L_SCRaMbLE/qTagGer/```
4. Change the files as needed (Note: Make sure to change the target path in the config.yml to include "input/" as a prefix)
6. Execute qTagGer by running
```python3 qTagGer.py -c input/config.yml```
7. qTagGer now executes and generates the primers and outputs them to the folder `{jobname}_{datetime}`

### Running qTagGer without Docker
Make sure that you have installed all aforementioned Python libraries as well as primer3 and bowtie correctly. After downloading qTagGer and modifying the config.yml, you can run the program as follows:
> python3 qTagGer.py -c config.yml


## Configuration
| Name | Type | Example | Description |
|------|------|---------|-------------|
| jobname | str | "NeoChrIII" | Jobtitle for the current run|
| top | int | 5 | Number of primer pairs per site |
| target | str | "neochrIII.gb" | Path to target record |
| offtarget_size_cutoff | int | 10000 | Threshold for max considered PCR offtarget product |
| sponge_value | int | 5 | Number of max allowed primer binding sites|
| home | str | "/home/L_SCRaMbLE/qTagGer" | Path to parent folder |
| primer3_path | str | "/home/L_SCRaMbLE/qTagGer/primer3/src/primer3_core" | Path to primer3_core file |
| bowtie_path | str | "/home/L_SCRaMbLE/qTagGer/bowtie/index" | Path to bowtie index files |

## License
This is a repository written under the [CC BY-NC-SA 4.0 license](https://creativecommons.org/licenses/by-nc-sa/4.0/)

## Publication
Lindeboom, T. A., Sanchez Olmos, M. del C., Schulz, K., Brinkmann, C. K., Ramirez Rojas, A. A., Hochrein, L., & Schindler, D. (2022). L-SCRaMbLE creates large-scale genome rearrangements in synthetic Sc2.0 chromosomes. bioRxiv. https://doi.org/10.1101/2022.12.12.519280


## References
Rozen, S., & Skaletsky, H. (n.d.). Primer3 on the WWW for General Users and for Biologist Programmers. In Bioinformatics Methods and Protocols (pp. 365–386). Humana Press. https://doi.org/10.1385/1-59259-192-2:365

Koressaar, T., & Remm, M. (2007). Enhancements and modifications of primer design program Primer3. In Bioinformatics (Vol. 23, Issue 10, pp. 1289–1291). Oxford University Press (OUP). https://doi.org/10.1093/bioinformatics/btm091

Untergasser, A., Cutcutache, I., Koressaar, T., Ye, J., Faircloth, B. C., Remm, M., & Rozen, S. G. (2012). Primer3—new capabilities and interfaces. In Nucleic Acids Research (Vol. 40, Issue 15, pp. e115–e115). Oxford University Press (OUP). https://doi.org/10.1093/nar/gks596

Langmead, B., Trapnell, C., Pop, M., & Salzberg, S. L. (2009). Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. In Genome Biology (Vol. 10, Issue 3, p. R25). Springer Science and Business Media LLC. https://doi.org/10.1186/gb-2009-10-3-r25

Peter J. A. Cock, Tiago Antao, Jeffrey T. Chang, Brad A. Chapman, Cymon J. Cox, Andrew Dalke, Iddo Friedberg, Thomas Hamelryck, Frank Kauff, Bartek Wilczynski, Michiel J. L. de Hoon, Biopython: freely available Python tools for computational molecular biology and bioinformatics, Bioinformatics, Volume 25, Issue 11, 1 June 2009, Pages 1422–1423, https://doi.org/10.1093/bioinformatics/btp163
