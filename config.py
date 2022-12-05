import yaml
import logging


def load(config_path:str) -> None:
    global jobname
    global target
    global chromosome
    global subject

    global pcr_size
    global homepath
    global primer3_path
    global bowtie_path
    global p3_settings
    

    logging.debug("Started config.py load procedure")
    with open(config_path) as file:
        config = yaml.load(file, Loader=yaml.FullLoader)

    jobname = config["jobname"]
    target = config["target"]
    chromosome = config["chromosome"]
    pcr_size = config["pcr_size_target"] # fix this here
    homepath = config["path"]
    primer3_path = config["primer3_path"]
    bowtie_path = config["bowtie_path"]
    p3_settings = config["primer3_settings"] 