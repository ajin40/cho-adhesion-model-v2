import os
import sys
import model
import numpy as np


def parameter_sweep(yaml_file, path="/storage/home/hcoda1/0/ajin40/p-mkemp6-0/ajin/cho_adhesion_model"):
    print(yaml_file + " ...")
    os.chdir(path + '/yaml_parameters')
    model.TestSimulation.start_sweep(path + '/outputs', yaml_file, f"{yaml_file[:-5]}", 0)
    print("Finished")


if __name__ == "__main__":
    path = "/storage/home/hcoda1/0/ajin40/p-mkemp6-0/ajin/cho_adhesion_model"
    process_tag = int(sys.argv[1])
    #replicate_tag = int(sys.argv[2])
    #cells_tag = int(sys.argv[3])
    os.chdir(path)
    yaml_array = os.listdir(path + "/yaml_parameters")
    #parameter_sweep(f'{process_tag}_dox_aba_{replicate_tag}_{cells_tag}_cells.yaml')
    parameter_sweep(f'{process_tag}_sim.yaml')

