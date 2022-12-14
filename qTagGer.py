import argparse
import logging

from config import Config
from target import Target
from selection import RoxP
from primer3 import Primer3Interface
from bowtie import BowtieInterface
from offtarget import OfftargetChecker



logging.basicConfig(filename='report.log', encoding='utf-8', level=logging.DEBUG)


def main(args:argparse.Namespace) -> None:
    config = Config(config_path=args.config)
    target = Target(target_path=config.target_path)
    selection = RoxP(target=target)

    p3interface = Primer3Interface(selection=selection, primer3_path=config.primer3_path)
    p3interface.run()

    bowtie = BowtieInterface(target_path=)
    check_offtargets = OfftargetChecker()

    # Cleanup routine


if __name__ == "__main__":
    # Parser
    parser = argparse.ArgumentParser(description="qTagGer - An Automtic Amplicon Primer Designer")
    parser.add_argument("-c", "--config", required=True, type=str, help="Path to the config file (YAML)")
    args = parser.parse_args()