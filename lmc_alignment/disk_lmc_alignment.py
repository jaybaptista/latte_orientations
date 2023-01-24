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

    print(f"Calculating Disk-LMC alignment for {sim_name}")

    radii_data = pd.read_hdf("../scale_radii/scale_radii_{}.hdf".format(sim_name))

    r_vir = radii_data["virial"][0]
    rs90 = radii_data["star.radius.90"][0]

    lmc_af = asdf.open(f"lmc_positions_{sim_name}.asdf")

    radii = [r_vir, rs90, 5 * rs90]

    snapshot_range = config[sim_name]

    df = {"angle": []}

    for snap in snapshot_range:
        print(f"{sim_name} @ snapshot {snap}")
        part = gizmo.io.Read.read_snapshots(["dark"], "snapshot", snap, sim)
        disk_vector = part.host["rotation"][0, 2]

        to_lmc = lmc_af["position"][np.where(lmc_af["snapshot"] == snap)[0][0]]
        to_lmc_hat = to_lmc / np.linalg.norm(to_lmc)
        
        angle = ori.getMinAngle(disk_vector, to_lmc_hat) * 180 / np.pi

        df["angle"].append(angle)

    df["snapshot"] = snapshot_range

    print(f"Writing alignment data file for {sim_name}...")

    af = asdf.AsdfFile(df)
    af.write_to(f"disk_lmc_alignment_{sim_name}.asdf")
