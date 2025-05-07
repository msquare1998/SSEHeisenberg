"""*************************************************
    SSE for 1D S=1/2 Heisenberg model
    Author: Yi-Ming Ding
    Email: dingyiming@westlake.edu.cn
    Updated: May 07, 2025
*************************************************"""
from hm import Heisenberg
import time
from tqdm import tqdm

def get_system_time_as_seed():
    return int(time.time())

def main():
    # ========================================================
    #  Input the parameters
    # ========================================================
    para_l = 4
    para_beta = 4
    num_thm = 20000
    num_stat = 10000
    num_bins = 5
    para_seed = get_system_time_as_seed()

    # ========================================================
    #  Report the environment
    # ========================================================
    print(f"■ SSE for S=1/2 Heisenberg model, L = {para_l}, beta = {para_beta}")
    print(f"■ num_thm = {num_thm}, num_stat = {num_stat}, num_bins = {num_bins}")

    # ===============================
    #  Prepare for saving data
    # ===============================
    f_energy = open("./data/energy.dat", "w")
    f_zz_corr = open("./data/zz_corr.dat", "w")

    # ===============================================================
    #  Quantum Monte Carlo simulations
    # ===============================================================
    model = Heisenberg(para_l, para_beta, para_seed)

    print("\t---> Thermalization...")
    for _ in tqdm(range(num_thm)):
        model.mc_update()
        model.adjust_m()

    print(f"\t---> Maximum cut-off = {model.m}")
    print("\t---> Measuring...")
    for _ in tqdm(range(num_bins)):
        model.init_measure()
        for __ in range(num_stat):
            model.mc_update()
            model.measure()
        model.statisticize(num_stat)

        f_energy.write(str(model.energy) + "\n")
        f_zz_corr.write("\t".join(str(x) for x in model.zz_correlation) + "\n")

if __name__ == "__main__":
    main()