import os.path as path
from jmb_utils import obtainSFR

simulation_list = [
    "../../../data/latte_metaldiff/m12f_res7100",
    "../../../data/latte_metaldiff/m12i_res7100",
    "../../../data/latte_metaldiff/m12m_res7100",
    "../../../data/latte_metaldiff/m12w_res7100",
]

for sim in simulation_list:
    df = obtainSFR(sim, 600)
    df.to_hdf(f"sfr_{path.split(sim)[-1]}.hdf", key="w")
