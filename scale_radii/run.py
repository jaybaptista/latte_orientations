import os.path as path
import jmb_utils as jb

simulation_list = [
    "../../../data/latte_metaldiff/m12f_res7100",
    "../../../data/latte_metaldiff/m12i_res7100",
    "../../../data/latte_metaldiff/m12m_res7100",
    "../../../data/latte_metaldiff/m12w_res7100",
]

for sim in simulation_list:
    df = jb.getScaleRadii(sim, 600)
    df.to_hdf(f"orientations_{path.split(sim)[-1]}.hdf", key="w")
