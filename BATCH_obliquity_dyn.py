import numpy as np
import utilities as ut
import gizmo_analysis as gizmo
from ob import getMinAngle, getAngularMomentum, getSymmetryAxes
import pandas as pd

gals = ["m12f_res7100", "m12i_res7100", "m12w_res7100", "m12m_res7100"]

for sim in gals:
    sim_dir = "../data/latte_metaldiff/" + sim

    step = 50

    if sim == "m12i_res7100" or sim == "m12m_res7100":
        step = 60

    radii_data = pd.read_hdf("radii_{}.hdf".format(sim))
    r_vir = radii_data["virial"][0]
    rs90 = radii_data["star.radius.90"][0]

    radii = [r_vir, rs90, 5 * rs90]
    snapshot_range = np.arange(0, 600 + step, step=step)

    df = {"virial": [], "disk": [], "5.disk": []}  # disk is star.radius.90

    for k in np.arange(0, 3):

        for snap in snapshot_range:

            part = gizmo.io.Read.read_snapshots(["dark"], "snapshot", snap, sim_dir)
            host_min_ax = part.host["rotation"][0, 2]

            positions = part["dark"].prop("host.distance")
            dists = part["dark"].prop("host.distance.total")

            tensor = getSymmetryAxes(positions, dists, radius=radii[k])
            angle = getMinAngle(tensor[2], host_min_ax) * 180 / np.pi

            df[list(df.keys())[k]].append(angle)

    d = pd.DataFrame(df)
    d["snapshot"] = snapshot_range
    d.to_hdf("data_obliquityAngles_z_dyn_{}.h5".format(sim), "w")
