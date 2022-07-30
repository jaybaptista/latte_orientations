import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

m12f_data = pd.read_hdf("orientations_present_day_m12f_res7100.hdf")
m12i_data = pd.read_hdf("orientations_present_day_m12i_res7100.hdf")
m12m_data = pd.read_hdf("orientations_present_day_m12m_res7100.hdf")
m12w_data = pd.read_hdf("orientations_present_day_m12w_res7100.hdf")

dr = rmin = 2
rmax = 400 + dr
radii = np.arange(rmin, rmax, step=dr)

fig, ax = plt.subplots(dpi=200)

ax.plot(radii, m12f_data["angle"], label="m12f", c="blue")
ax.plot(radii, m12i_data["angle"], label="m12i", c="red")
ax.plot(radii, m12m_data["angle"], label="m12m", c="green")
ax.plot(radii, m12w_data["angle"], label="m12w", c="orange")

ax.axhline(20, c="k", ls="--", alpha=0.25)

ax.set_xlabel("Radius [kpc]")
ax.set_ylabel("Orientation w.r.t. to disk [deg]")

ax.legend()

plt.savefig("orientations_present_day_cdm.pdf", bbox_inches="tight")
