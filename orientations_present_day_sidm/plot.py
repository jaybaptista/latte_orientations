import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "serif"})

m12i_data = pd.read_hdf(
    "../orientations_present_day_cdm/orientations_present_day_m12i_res7100.hdf"
)
m12f_data = pd.read_hdf(
    "../orientations_present_day_cdm/orientations_present_day_m12f_res7100.hdf"
)
m12m_data = pd.read_hdf(
    "../orientations_present_day_cdm/orientations_present_day_m12m_res7100.hdf"
)

m12i_data_DMO = pd.read_hdf(
    "../orientations_present_day_dmo/orientations_present_day_m12i_cdm-only.hdf"
)
m12f_data_DMO = pd.read_hdf(
    "../orientations_present_day_dmo/orientations_present_day_m12f_cdm-only.hdf"
)
m12m_data_DMO = pd.read_hdf(
    "../orientations_present_day_dmo/orientations_present_day_m12m_cdm-only.hdf"
)

m12i_data_SIDM = pd.read_hdf("orientations_present_day_m12i_sidm1.hdf")
m12f_data_SIDM = pd.read_hdf("orientations_present_day_m12f_sidm1.hdf")
m12m_data_SIDM = pd.read_hdf("orientations_present_day_m12m_sidm1.hdf")

fig, ax = plt.subplots(3, dpi=200, figsize=(8, 4), sharex=True)

dr = 2
rmax = 400

r = np.arange(dr, rmax + dr, step=dr)

ax_f, ax_i, ax_m = ax

ax_f.plot(r, m12f_data["angle"], ms=2, c="k", alpha=0.5, label="CDM+Stars", lw=3)

ax_f.plot(r, m12f_data_DMO["angle"], ms=2, c="blue", label="DMO", lw=3)

ax_f.plot(r, m12f_data_SIDM["angle"], ms=2, c="dodgerblue", label="SIDM+Stars", lw=3)

ax_i.plot(r, m12i_data["angle"], ms=2, c="k", alpha=0.5, label="m12i", lw=3)

ax_i.plot(r, m12i_data_DMO["angle"], ms=2, c="blue", label="m12i", lw=3)

ax_i.plot(r, m12i_data_SIDM["angle"], ms=2, c="dodgerblue", label="m12i", lw=3)

ax_m.plot(r, m12m_data["angle"], ms=2, c="k", alpha=0.5, label="m12m", lw=3)

ax_m.plot(r, m12m_data_DMO["angle"], ms=2, c="blue", label="m12m", lw=3)

ax_m.plot(r, m12m_data_SIDM["angle"], ms=2, c="dodgerblue", label="m12m", lw=3)

ax_f.legend(loc="upper left", prop={"size": 11}, bbox_to_anchor=(1.0, 1))

ax_f.axhline(20, c="k", alpha=0.25)
ax_i.axhline(20, c="k", alpha=0.25)
ax_m.axhline(20, c="k", alpha=0.25)

ax_f.set_yticklabels([0, 20])
ax_i.set_yticklabels([0, 20])
ax_m.set_yticklabels([0, 20])

for ax_k in ax:
    ax_k.set_ylim(0, 40)
    ax_k.axvline(60, alpha=0.25, c="k", ls="--")
    ax_k.axvline(200, alpha=0.25, c="k", ls="--")
    ax_k.axvline(400, alpha=0.25, c="k", ls="--")


ax_f.text(65, 38, "Gaia", rotation=90, size=11, va="top")
ax_f.text(205, 38, "AS", rotation=90, size=11, va="top")
ax_f.text(405, 38, "LSST", rotation=90, size=11, va="top")

ax_m.set_xlabel("r [kpc]", size=11)

fig.text(
    -0.075,
    0.5,
    "m12f",
    transform=ax_f.transAxes,
    rotation=90,
    ha="center",
    va="center",
    size=11,
    fontweight="bold",
)
fig.text(
    -0.075,
    0.5,
    "m12i",
    transform=ax_i.transAxes,
    rotation=90,
    ha="center",
    va="center",
    size=11,
    fontweight="bold",
)
fig.text(
    -0.125,
    0.5,
    "$\\theta$ [deg]",
    transform=ax_i.transAxes,
    rotation=90,
    ha="center",
    va="center",
    size=11,
)
fig.text(
    -0.075,
    0.5,
    "m12m",
    transform=ax_m.transAxes,
    rotation=90,
    ha="center",
    va="center",
    size=11,
    fontweight="bold",
)

fig.tight_layout()
fig.subplots_adjust(hspace=0)
plt.savefig("present_day_orientations_altdm.png", bbox_inches="tight")
plt.show()
