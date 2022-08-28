import numpy as np
import gizmo_analysis as gizmo
import utilities as ut
import orientations as ori
import asdf

# Radius of axis ratio calculations (put into logspace) [kpc]
r_ini = 0.1
r_fin = 200

# Number of points calculated for radius in logspace
N_logspace = 10 #1000

# Put radius values into logspace
r_ini_log = np.log10(r_ini)
r_fin_log = np.log10(r_fin)
r = np.logspace(r_ini_log, r_fin_log, N_logspace, base=10)

# Tolerance in axis ratios and maximum iterations (using reduced inertia tensor method)
tolerance = 0.001
iters_max = 1000

part = gizmo.io.Read.read_snapshots("dark", "snapshot", 600, "../../../data/latte_metaldiff/m12f_res7100", assign_hosts_rotation=True)

dist = part["dark"].prop('host.distance.total')
pos  = part["dark"].prop('host.distance.principal')
mass = part["dark"]["mass"]

## VERSION 3
dfs = {
        "rmax": np.zeros(N_logspace),
        "shape": np.zeros(N_logspace),
        "iter": np.zeros(N_logspace),
        "tensor": np.zeros((N_logspace, 3,3)),
        "ortho": np.zeros(N_logspace),
        "a": np.zeros(N_logspace),
        "b": np.zeros(N_logspace),
        "c": np.zeros(N_logspace),
        "s": np.zeros(N_logspace),
        "p": np.zeros(N_logspace),
        "q": np.zeros(N_logspace),
    }

for i, r_i in enumerate(r):
    print("Running @", r_i)
    df = ori.getShapeV3(pos, r_i, mass, distances=dist, tolerance=0.001, iters_max=1000)
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

af = asdf.AsdfFile(dfs)
af.write_to("tmp_3.asdf")