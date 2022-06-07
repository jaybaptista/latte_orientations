from jmb_utils.galaxy_utils import getHaloOrientation
import numpy as np
import os.path as path

rmax = 500
dr = 2
radii = np.arange(dr, rmax + dr, step=dr)

simulation_list = [
    "../../../data/latte_metaldiff/m12f_res7100",
    "../../../data/latte_metaldiff/m12i_res7100",
    "../../../data/latte_metaldiff/m12m_res7100",
    "../../../data/latte_metaldiff/m12w_res7100",
]

for sim in simulation_list:
    df = getHaloOrientation(sim, 600, radii)
    df.to_hdf(f"radii_{path.split(sim)[-1]}.hdf", key="w")
