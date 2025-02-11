#!/usr/bin/env python

import argparse
import codex_pipeline
import os


def initialize_workspace():
    codex_pipeline.dump_parameters()


    os.mkdir("Input_L2")
    os.mkdir("Output_L2")




if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="CODEX Data Pipeline Master")

    parser.add_argument("--init", action="store_true", help="Initialized to workspace")
    
    args = parser.parse_args()


    if( args.init ):    initialize_workspace()





    