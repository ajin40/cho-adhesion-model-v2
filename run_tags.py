import os
import sys
import model
import ast
import numpy as np

def parameter_sweep(model_params, name, path="/storage/scratch1/0/ajin40/ST_cho_model"):
    model.TestSimulation.start_sweep(path + '/outputs', model_params, name)
    print("Finished")

def parameter_sweep_abm(par, directory, RR, YY, RY, dox_ratio, aba_ratio, final_ts=60):
    """ Run model with specified parameters
        :param par: simulation number
        :param directory: Location of model outputs. A folder titled 'outputs' is required
        :param RR: Adhesion value for R-R cell interactions
        :param YY: Adhesion value for Y-Y cell interactions
        :param RY: Adhesion value for R-Y cell interactions
        :param dox_ratio: final ratio of Red cells at simulation end. 1 - (dox_ratio + aba_ratio) = # of remaining uncommitted blue cells.
        :param aba_ratio: final ratio of Yellow cells at simulation end.
        :param final_ts: Final timestep. 60 ts = 96h, 45 ts = 72h
        :type par int
        :type directory: String
        :type RR: float
        :type YY: float
        :type RY: float
        :type dox_ratio: float
        :type aba_ratio: float
        :type final_ts: int
    """
    model_params = {
        "u_bb": 1,
        "u_rb": 1,
        "u_yb": 1,
        "u_rr": RR,
        "u_repulsion": 10000,
        "alpha": 10,
        "gravity": 2,
        "end_step": final_ts,
        "u_yy": YY,
        "u_ry": RY,
        "dox_ratio": dox_ratio,
        "aba_ratio": aba_ratio,
        "PACE": False
    }
    name = f'urr{RR}_uyy{YY}_ury{RY}_dox{dox_ratio}_aba{aba_ratio}'
    sim = model.TestSimulation(model_params)
    sim.start_sweep(directory + '/outputs', model_params, name)
    return par, sim.image_quality, sim.image_quality, 3, final_ts/sim.sub_ts

if __name__ == "__main__":
    ury = [1, 5, 10, 20, 30]
    path = "/storage/scratch1/0/ajin40/ST_cho_model"
    u_rr = int(sys.argv[1])
    u_yy = int(sys.argv[2])
    u_ry = ast.literal_eval(sys.argv[3])
    dox = float(sys.argv[4])
    aba = float(sys.argv[5])
    for value in u_ry:
        model_params = {
            "u_bb": 1,
            "u_rb": 1,
            "u_yb": 1,
            "u_rr": u_rr,
            "u_repulsion": 10000,
            "alpha": 10,
            "gravity": 2,
            "u_yy": u_yy,
            "u_ry": value,
            "dox_ratio": dox,
            "aba_ratio": aba,
            "PACE": True
        }
        os.chdir(path)
        #parameter_sweep(f'{process_tag}_dox_aba_{replicate_tag}_{cells_tag}_cells.yaml')
        parameter_sweep(model_params, f'urr{u_rr}_uyy{u_yy}_ury{value}_dox{dox}_aba{aba}', path)

