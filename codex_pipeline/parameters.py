import json
import os
import codex_pipeline
import sys


dir_cal_pkg = {}
dir_cal_pkg['qe'] = 'QE'
dir_cal_pkg['dem_tens'] = 'DemTensor'
dir_cal_pkg['vign'] = 'Vignetting'

parameter_file = 'parameters.json'

def init_parameters():
    
    
    parameter_dir, _ = os.path.split(__file__)
    parameter_path = os.path.join(parameter_dir, parameter_file)

    # Load CODEX Parameters from JSON file
    with open(parameter_path) as f:
        codex_pipeline.parameters = json.load(f)

    if os.path.exists(parameter_file):
        with open(parameter_path) as f:
            override_parameters = json.load(f)

    codex_pipeline.parameters.update(override_parameters)

    if not os.path.exists(codex_pipeline.parameters['cal_pkg']):
        root_dir, _ = os.path.split(sys.argv[0])
        codex_pipeline.parameters['cal_pkg'] =  os.path.join(root_dir, "CalibrationPackage")
    
    

def dump_parameters():

    with open(parameter_file, "w", encoding="utf-8") as file:
        json.dump(codex_pipeline.parameters, file, indent=4)


def get_param(key):
    return codex_pipeline.parameters[key]

def get_path(key):

    cal_pkg = codex_pipeline.parameters['cal_pkg']

    if key in ["vign"]:
        return os.path.join(cal_pkg, dir_cal_pkg['vign'])
    elif key in ["Darks"]:
        return os.path.join(cal_pkg, 'Darks')
    elif key in ["dem_tens"]:
        return os.path.join(cal_pkg, dir_cal_pkg['dem_tens'])
    elif key in ["qe"]:
        return os.path.join(cal_pkg, dir_cal_pkg['qe'])
    else:
        raise KeyError
    

