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

ds = 50

for sim in simulation_list:

    sim_name = path.split(sim)[-1]

    radii_data = pd.read_hdf("../scale_radii/scale_radii_{}.hdf".format(sim_name))

    r_vir = radii_data["virial"][0]
    rs90 = radii_data["star.radius.90"][0]

    snapshot_range = np.arange(ds, 600 + ds, step=ds)

    part_tmp = gizmo.io.Read.read_snapshots(["dark"], "snapshot", 600, sim)
    host_min_ax = part_tmp.host["rotation"][0, 2]
    del part_tmp

    df = {"virial": [], "disk": [], "5.disk": []}  # disk is star.radius.90

    radii = [r_vir, rs90, 5 * rs90]

    for snap in snapshot_range:

        part = gizmo.io.Read.read_snapshots(["dark"], "snapshot", snap, sim)
        positions = part["dark"].prop("host.distance")
        dists = part["dark"].prop("host.distance.total")
        mass = part["dark"].prop("mass")
        
        for k in np.arange(0, 3):

            shape_df = ori.getShape(positions, radii[k], mass, distances=dists, tolerance=0.001, iters_max=100)
            df[list(df.keys())[k]].append(shape_df["shape"])

    d = asdf.AsdfFile(df)
    d["snapshot"] = snapshot_range

    d.write_to(f"shape_over_z_{sim_name}.asdf")
