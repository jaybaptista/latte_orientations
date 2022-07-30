import os.path as path
import orientations as ori
import gizmo_analysis as gizmo
import numpy as np
import pandas as pd
import asdf

simulation_list = [
    "../../../data/latte_metaldiff/m12f_res7100",
    "../../../data/latte_metaldiff/m12i_res7100",
    "../../../data/latte_metaldiff/m12m_res7100",
    "../../../data/latte_metaldiff/m12w_res7100",
]

ds = 10

for sim in simulation_list:

    sim_name = path.split(sim)[-1]

    radii_data = asdf.open("../scale_radii/scale_radii_z_{}.asdf".format(sim_name))

    snapshot_range = np.arange(ds, 600 + ds, step=ds)

    df = {"virial": [], "disk": [], "5.disk": []}  # disk is star.radius.90

    for i, snap in enumerate(snapshot_range):

        r_vir = radii_data["virial"][i]
        rs90 = radii_data["star.radius.90"][i]
        radii = [r_vir, rs90, 5 * rs90]

        part = gizmo.io.Read.read_snapshots(["dark"], "snapshot", snap, sim)
        positions = part["dark"].prop("host.distance")
        host_min_ax = part.host["rotation"][0, 2]
        dists = part["dark"].prop("host.distance.total")

        for k in np.arange(0, 3):

            tensor = ori.getSymmetryAxes(positions, dists, radius=radii[k])
            angle = ori.getMinAngle(tensor[2], host_min_ax) * 180 / np.pi

            df[list(df.keys())[k]].append(angle)

    d = pd.DataFrame(df)
    d["snapshot"] = snapshot_range

    d.to_hdf(f"orientations_over_z_variable_radius_{sim_name}.hdf", key="w")
