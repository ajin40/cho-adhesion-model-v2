import sys
import model
import os
import numpy as np

def dox_aba_matrix(directory, rr, yy, ry, replicate, final_ts=45, sub_ts=960):
    """
    Runs simulations for all dose ratios of a specific experimental replicate (0-2).
    """
    dox_ratio = get_ratios(directory, "dox", replicate)
    aba_ratio = get_ratios(directory, "aba", replicate)
    dox_aba_list = tuple(zip(dox_ratio, aba_ratio))
    for inducer_pair in dox_aba_list:
        model.run_abm(0, path, rr, yy, ry, inducer_pair[0], inducer_pair[1])

def get_ratios(directory, cell_type, replicate):
    os.chdir(directory)
    ratio = np.loadtxt(f"channel_intensities_{cell_type}_{replicate}.csv", delimiter=',', encoding='utf-8-sig')
    return ratio

if __name__ == "__main__":
    path = "[USER SPECIFIES BEFOREHAND]"
    u_rr = int(sys.argv[1])
    u_yy = int(sys.argv[2])
    u_ry = int(sys.argv[3])
    dox = float(sys.argv[4])
    aba = float(sys.argv[5])
    final_ts = int(sys.argv[6])
    sub_ts = int(sys.argv[7])

    model.run_abm(0, path + '/outputs', u_rr, u_yy, u_ry, dox, aba, final_ts=final_ts, sub_ts=sub_ts)
