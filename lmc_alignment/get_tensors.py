import numpy as np
import gizmo_analysis as gizmo
import pandas as pd
import asdf
import orientations as ori
import os.path as path

config = {
    "m12f_res7100": np.arange(302, 524, step=5),
    "m12i_res7100": np.arange(200, 480, step=5),
    "m12m_res7100": np.arange(153, 564, step=5),
    "m12w_res7100": np.arange(152, 380, step=5),
}

simulation_list = [
    "../../../data/latte_metaldiff/m12f_res7100",
    "../../../data/latte_metaldiff/m12i_res7100",
    "../../../data/latte_metaldiff/m12m_res7100",
    "../../../data/latte_metaldiff/m12w_res7100",
]

for sim in simulation_list:

    sim_name = path.split(sim)[-1]

    print(f"Getting tensor list for {sim_name}")

    radii_data = pd.read_hdf("../scale_radii/scale_radii_{}.hdf".format(sim_name))

    r_vir = radii_data["virial"][0]
    rs90 = radii_data["star.radius.90"][0]

    lmc_af = asdf.open(f"lmc_positions_{sim_name}.asdf")

    radii = [r_vir, rs90, 5 * rs90]

    snapshot_range = config[sim_name]

    df = {"virial": [], "disk": [], "5.disk": []}

    for snap in snapshot_range:
        print(f"{sim_name} @ snapshot {snap}")
        part = gizmo.io.Read.read_snapshots(["dark"], "snapshot", snap, sim)
        positions = part["dark"].prop("host.distance")
        dists = part["dark"].prop("host.distance.total")

        for k in np.arange(0, 3):
            tensor = ori.getSymmetryAxes(positions, dists, radius=radii[k])
            df[list(df.keys())[k]].append(tensor)

    df["snapshot"] = snapshot_range

    print(f"Writing tensor data file for {sim_name}...")

    af = asdf.AsdfFile(df)
    af.write_to(f"tensors_{sim_name}.asdf")
