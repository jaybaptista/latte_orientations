import os.path as path
import orientations as ori
import asdf
import numpy as np

simulation_list = [
    "../../../data/latte_metaldiff/m12f_res7100",
    "../../../data/latte_metaldiff/m12i_res7100",
    "../../../data/latte_metaldiff/m12m_res7100",
    "../../../data/latte_metaldiff/m12w_res7100",
]

for sim in simulation_list:

    # present day scales

    df = ori.getScaleRadii(sim, 600)
    df.to_hdf(f"scale_radii_{path.split(sim)[-1]}.hdf", key="w")

    # time-dependent scales

    ds = 10
    snapshot_range = np.arange(ds, 600 + ds, step=ds)

    af = asdf.AsdfFile(
        {
            "snapshot": snapshot_range,
            "virial": np.zeros(snapshot_range.size),
            "star.radius.50": np.zeros(snapshot_range.size),
            "star.radius.90": np.zeros(snapshot_range.size),
        }
    )

    for i, snap in enumerate(snapshot_range):
        df = ori.getScaleRadii(sim, snap)

        af["virial"][i] = df["virial"]
        af["star.radius.50"][i] = df["star.radius.50"]
        af["star.radius.90"][i] = df["star.radius.90"]

    af.write_to(f"scale_radii_z_{path.split(sim)[-1]}.asdf")
