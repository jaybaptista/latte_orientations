import os.path as path
from jmb_utils import obtainSFR

simulation_list = [
    "../../../data/m12f_sidm1",
    "../../../data/m12i_sidm1",
    "../../../data/m12m_sidm1",
]

rs50_list = [3.8, 4.7, 8.2]  # m12f  # m12i  # m12m

for i, sim in enumerate(simulation_list):
    df = obtainSFR(sim, 60, rmax=5 * rs50_list[i])
    df.to_hdf(f"sidm_sfr_{path.split(sim)[-1]}.hdf", key="w")
