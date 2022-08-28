import pandas as pd
import asdf
import matplotlib.pyplot as plt
import numpy as np
import utilities as ut
import orientations as ori

plt.rcParams.update({"font.family": "serif"})

m12f_data = pd.read_hdf("orientations_over_z_variable_radius_m12f_res7100.hdf")
m12i_data = pd.read_hdf("orientations_over_z_variable_radius_m12i_res7100.hdf")
m12m_data = pd.read_hdf("orientations_over_z_variable_radius_m12m_res7100.hdf")
m12w_data = pd.read_hdf("orientations_over_z_variable_radius_m12w_res7100.hdf")

m12f_time = ori.getSnapshotData(
    "../../../data/latte_metaldiff/m12i_res7100/", m12f_data["snapshot"].to_numpy()
)[0]
m12i_time = ori.getSnapshotData(
    "../../../data/latte_metaldiff/m12i_res7100/", m12i_data["snapshot"].to_numpy()
)[0]
m12m_time = ori.getSnapshotData(
    "../../../data/latte_metaldiff/m12i_res7100/", m12m_data["snapshot"].to_numpy()
)[0]
m12w_time = ori.getSnapshotData(
    "../../../data/latte_metaldiff/m12i_res7100/", m12w_data["snapshot"].to_numpy()
)[0]

af_md = asdf.open("../lmc_alignment/peri_md.asdf")

fig, ax = plt.subplots(2, 2, dpi=200, figsize=(8, 4))

((ax_f, ax_i), (ax_w, ax_m)) = ax

ax_i.plot(
    m12i_time,
    m12i_data["disk"],
    c="b",
    alpha=0.25,
    label="Stellar Disk ($R_{*,90}$)",
    lw=3,
)
ax_i.plot(
    m12i_time, m12i_data["5.disk"], c="b", alpha=0.5, label="5$ \\times R_{*,90}$", lw=3
)
ax_i.plot(
    m12i_time, m12i_data["virial"], c="b", label="Virial Radius $(R_{200})$", lw=3
)

ax_f.plot(
    m12f_time,
    m12f_data["disk"],
    c="b",
    alpha=0.25,
#     label="Stellar Disk ($R_{*,90}$)",
    lw=3,
)
ax_f.plot(
    m12f_time, m12f_data["5.disk"], c="b", alpha=0.5, lw=3
)
ax_f.plot(
    m12f_time, m12f_data["virial"], c="b", lw=3
)

ax_m.plot(m12m_time, m12m_data["disk"], c="b", alpha=0.25, lw=3)
ax_m.plot(m12m_time, m12m_data["5.disk"], c="b", alpha=0.5, lw=3)
ax_m.plot(m12m_time, m12m_data["virial"], c="b", lw=3)

ax_w.plot(m12w_time, m12w_data["disk"], c="b", alpha=0.25, lw=3)
ax_w.plot(m12w_time, m12w_data["5.disk"], c="b", alpha=0.5, lw=3)
ax_w.plot(m12w_time, m12w_data["virial"], c="b", lw=3)

ax_f.axvline(af_md["m12f_res7100"]["peri.t"], alpha=0.4, c="k")
ax_i.axvline(
    af_md["m12i_res7100"]["peri.t"], alpha=0.4, c="k", label="$t_\mathrm{peri}$"
)
ax_m.axvline(af_md["m12m_res7100"]["peri.t"], alpha=0.4, c="k")
ax_w.axvline(af_md["m12w_res7100"]["peri.t"], alpha=0.4, c="k")

ax_f.axvline(7.6, alpha=0.4, c="r")
ax_i.axvline(8.5, alpha=0.4, c="r", label="$t_\mathrm{bursty}$")
ax_m.axvline(3.5, alpha=0.4, c="r")
ax_w.axvline(8.4, alpha=0.4, c="r")

ax_m.text(0.05, 0.1, "m12m", transform=ax_m.transAxes, bbox=dict(boxstyle="round", fc="#EEE", ec="#DDD", alpha=0.95))
ax_i.text(0.05, 0.1, "m12i", transform=ax_i.transAxes, bbox=dict(boxstyle="round", fc="#EEE", ec="#DDD", alpha=0.95))
ax_f.text(0.05, 0.1, "m12f", transform=ax_f.transAxes, bbox=dict(boxstyle="round", fc="#EEE", ec="#DDD", alpha=0.95))
ax_w.text(0.05, 0.1, "m12w", transform=ax_w.transAxes, bbox=dict(boxstyle="round", fc="#EEE", ec="#DDD", alpha=0.95))

# ax_i.legend(bbox_to_anchor=(1.05, 0.45), fontsize=9)

for a in ax.flatten():
    a.set_ylim(0, 90)

fig.subplots_adjust(hspace=0.2, wspace=0.1)

fig.text(0.05, 0.5, "Orientation [deg]", rotation=90, va="center", ha="center", fontsize=18)
fig.text(0.5, -.025, "Time since beginning [Gyr]", ha="center", fontsize=18)

ax_f.text(
    0.5, 1.06, "Major Mergers", ha="center", transform=ax_f.transAxes, fontweight="bold", fontsize=14
)
ax_i.text(
    0.5, 1.06, "Minor Mergers", ha="center", transform=ax_i.transAxes, fontweight="bold", fontsize=14
)

plt.figlegend(ncol=5, loc="lower center", bbox_to_anchor=(0.5-(.8/2), -0.15, .8, .075))

plt.savefig("orientations_over_z_variable_radius.png", bbox_inches="tight")
plt.show()
