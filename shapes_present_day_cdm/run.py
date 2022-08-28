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

dr = rmin = 2
rmax = 400 + dr
N_pts = 20
radii = np.logspace(np.log10(rmin), np.log10(rmax), num=N_pts)
N = len(radii)


for sim in simulation_list:

    part = gizmo.io.Read.read_snapshots("dark", "snapshot", 600, sim, assign_hosts_rotation=True)

    dist = part["dark"].prop('host.distance.total')
    pos  = part["dark"].prop('host.distance.principal')
    mass = part["dark"]["mass"]
    
    dfs = asdf.AsdfFile({
            "rmax": np.zeros(N),
            "shape": np.zeros(N),
            "iter": np.zeros(N),
            "tensor": np.zeros((N, 3,3)),
            "ortho": np.zeros(N),
            "a": np.zeros(N),
            "b": np.zeros(N),
            "c": np.zeros(N),
            "s": np.zeros(N),
            "p": np.zeros(N),
            "q": np.zeros(N),
        })


    for i, r_i in enumerate(radii):
        print("Running @", r_i)
        df = ori.getShape(pos, r_i, mass, distances=dist, tolerance=0.001, iters_max=100)
        print(df)
        dfs["rmax"][i] = df["rmax"]
        dfs["shape"][i] = df["shape"]
        dfs["iter"][i] = df["iter"]
        dfs["tensor"][i] = df["tensor"]
        dfs["ortho"][i] = df["ortho"]
        dfs["a"][i] = df["a"]
        dfs["b"][i] = df["b"]
        dfs["c"][i] = df["c"]
        dfs["s"][i] = df["s"]
        dfs["p"][i] = df["p"]
        dfs["q"][i] = df["q"]
        print("Radii calculated: ", dfs["rmax"])

    dfs.write_to(f"shapes_present_day_{path.split(sim)[-1]}.asdf")
