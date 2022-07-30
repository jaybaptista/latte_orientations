import os.path as path
import orientations as ori
import gizmo_analysis as gizmo
import numpy as np
import pandas as pd

simulation_list = [
    "../../../data/m12f_sidm1",
    "../../../data/m12i_sidm1",
    "../../../data/m12m_sidm1",
]

dr = rmin = 2
rmax = 400 + dr

for sim in simulation_list:

    df = {"tensors": [], "angle": []}

    part = gizmo.io.Read.read_snapshots(
        ["dark", "star"], "redshift", 0, sim, assign_hosts_rotation=True
    )

    positions = part["dark"].prop("host.distance")
    dists = part["dark"].prop("host.distance.total")

    host_min_ax = part.host["rotation"][0, 2]

    for i in np.arange(rmin, rmax, step=dr):
        tensor = ori.getSymmetryAxes(positions, dists, radius=i)
        angle = ori.getMinAngle(tensor[2], host_min_ax) * 180 / np.pi
        df["tensors"].append(tensor)
        df["angle"].append(angle)

    df = pd.DataFrame(df)
    df.to_hdf(f"orientations_present_day_{path.split(sim)[-1]}.hdf", key="w")
