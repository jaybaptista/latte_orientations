from astropy.io.misc.hdf5 import read_table_hdf5
import numpy as np
import pandas as pd
import gizmo_analysis as gizmo
from ob import getSymmetryAxes, getMinAngle

gals = ["m12f_res7100", "m12i_res7100", "m12w_res7100", "m12m_res7100"]

for sim in gals:
    sim_dir = "../data/latte_metaldiff/" + sim

    lmc_position = read_table_hdf5("lmc_positions_{}.h5".format(sim))

    z0_lmc = lmc_position[np.argmax(lmc_position["snapshot"])]["position"]

    z0_lmc_hat = z0_lmc / np.linalg.norm(z0_lmc)

    df = {"tensors": [], "angle": []}

    part = gizmo.io.Read.read_snapshots(["dark"], "redshift", 0, sim_dir)

    positions = part["dark"].prop("host.distance")
    dists = part["dark"].prop("host.distance.total")

    for i in np.arange(10, 410, step=10):
        tensor = getSymmetryAxes(positions, dists, radius=i)
        angle = getMinAngle(tensor[0], z0_lmc_hat) * 180 / np.pi
        df["tensors"].append(tensor)
        df["angle"].append(angle)

    df = pd.DataFrame(df)
    df.to_hdf("_data_lmcAlignment_z0_{}.h5".format(sim[:4]), "w")
