import numpy as np
import utilities as ut
import gizmo_analysis as gizmo
from ob import getMinAngle, getAngularMomentum, getSymmetryAxes
import pandas as pd

gals = ["m12f_cdm-only", "m12i_cdm-only", "m12m_cdm-only"]

for sim in gals:

    print(sim)

    df = {"tensors": [], "angle": []}

    sim_dir = "../data/" + sim

    part = gizmo.io.Read.read_snapshots(["dark"], "redshift", 0, sim_dir)

    positions = part["dark"].prop("host.distance")
    dists = part["dark"].prop("host.distance.total")

    t0 = getSymmetryAxes(positions, dists, radius=10)
    host_min_ax = t0[2]

    for i in np.arange(2, 402, step=2):
        tensor = getSymmetryAxes(positions, dists, radius=i)
        angle = getMinAngle(tensor[2], host_min_ax) * 180 / np.pi
        df["tensors"].append(tensor)
        df["angle"].append(angle)

    df = pd.DataFrame(df)
    df.to_hdf("_data_obliquityAngles_z0_DMO_{}.h5".format(sim[:4]), "w")
