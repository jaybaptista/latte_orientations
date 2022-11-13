import pandas as pd
import asdf
import matplotlib.pyplot as plt
import numpy as np
import utilities as ut
import orientations as ori

plt.rcParams.update({"font.family": "serif"})

m12f_data = asdf.open("shape_over_z_variable_radius_m12f_res7100.asdf")
m12i_data = asdf.open("shape_over_z_variable_radius_m12i_res7100.asdf")
m12m_data = asdf.open("shape_over_z_variable_radius_m12m_res7100.asdf")
m12w_data = asdf.open("shape_over_z_variable_radius_m12w_res7100.asdf")

m12f_time = ori.getSnapshotData(
    "../../../data/latte_metaldiff/m12i_res7100/", np.array(m12f_data["snapshot"])
)[0]
m12i_time = ori.getSnapshotData(
    "../../../data/latte_metaldiff/m12i_res7100/", np.array(m12i_data["snapshot"])
)[0]
m12m_time = ori.getSnapshotData(
    "../../../data/latte_metaldiff/m12i_res7100/", np.array(m12m_data["snapshot"])
)[0]
m12w_time = ori.getSnapshotData(
    "../../../data/latte_metaldiff/m12i_res7100/", np.array(m12w_data["snapshot"])
)[0]

af_md = asdf.open("../lmc_alignment/peri_md.asdf")

############

disk_col = "#fcba04"
disk5_col = "#5c8001"
vir_col = "#00a8e8"

disk_ls = "-"
disk5_ls = "--"
vir_ls = "-."

disk_lw = 1
disk5_lw = 1
vir_lw = 1

dot_size = 5

##########################

fig, ax = plt.subplots(2, 2, dpi=200, figsize=(8, 4))

((ax_f, ax_i), (ax_w, ax_m)) = ax

ax_i.plot(m12i_time, np.array(m12i_data["disk"]), c=disk_col, label="Stellar Disk ($R_{*,90}$)", lw=disk_lw, ls=disk_ls)
ax_i.plot(m12i_time, np.array(m12i_data["5.disk"]), c=disk5_col, label="5$ \\times R_{*,90}$", lw=disk5_lw, ls=disk5_ls)
ax_i.plot(m12i_time, np.array(m12i_data["virial"]), c=vir_col, label="Virial Radius $(R_{200})$", lw=vir_lw, ls=vir_ls)

ax_f.plot(m12f_time, np.array(m12f_data["disk"]), c=disk_col, lw=disk_lw, ls=disk_ls)
ax_f.plot(m12f_time, np.array(m12f_data["5.disk"]), c=disk5_col, lw=disk5_lw, ls=disk5_ls)
ax_f.plot(m12f_time, np.array(m12f_data["virial"]), c=vir_col, lw=vir_lw, ls=vir_ls)

ax_m.plot(m12m_time, np.array(m12m_data["disk"]), c=disk_col, ls=disk_ls, lw=disk_lw)
ax_m.plot(m12m_time, np.array(m12m_data["5.disk"]), c=disk5_col, ls=disk5_ls, lw=disk5_lw)
ax_m.plot(m12m_time, np.array(m12m_data["virial"]), c=vir_col, ls=vir_ls, lw=vir_lw)

ax_w.plot(m12w_time, np.array(m12w_data["disk"]), c=disk_col, ls=disk_ls, lw=disk_lw)
ax_w.plot(m12w_time, np.array(m12w_data["5.disk"]), c=disk5_col, ls=disk5_ls, lw=disk5_lw)
ax_w.plot(m12w_time, np.array(m12w_data["virial"]), c=vir_col, ls=vir_ls, lw=vir_lw)

ax_f.axvline(af_md["m12f_res7100"]["peri.t"], alpha=0.75, c="k")
ax_i.axvline(
    af_md["m12i_res7100"]["peri.t"], alpha=0.75, c="k", label="$t_\mathrm{peri}$"
)
ax_m.axvline(af_md["m12m_res7100"]["peri.t"], alpha=0.75, c="k")
ax_w.axvline(af_md["m12w_res7100"]["peri.t"], alpha=0.75, c="k")

ax_f.axvline(7.6, alpha=0.75, c="r")
ax_i.axvline(8.5, alpha=0.75, c="r", label="$t_\mathrm{bursty}$")
ax_m.axvline(3.5, alpha=0.75, c="r")
ax_w.axvline(8.4, alpha=0.75, c="r")

ax_m.text(0.05, 0.1, "m12m", transform=ax_m.transAxes, bbox=dict(boxstyle="round", fc="#EEE", ec="#DDD", alpha=0.95))
ax_i.text(0.05, 0.1, "m12i", transform=ax_i.transAxes, bbox=dict(boxstyle="round", fc="#EEE", ec="#DDD", alpha=0.95))
ax_f.text(0.05, 0.1, "m12f", transform=ax_f.transAxes, bbox=dict(boxstyle="round", fc="#EEE", ec="#DDD", alpha=0.95))
ax_w.text(0.05, 0.1, "m12w", transform=ax_w.transAxes, bbox=dict(boxstyle="round", fc="#EEE", ec="#DDD", alpha=0.95))

# ax_i.legend(bbox_to_anchor=(1.05, 0.45), fontsize=9)

for a in ax.flatten():
    a.set_ylim(0, 1)

fig.subplots_adjust(hspace=0.2, wspace=0.125)

fig.text(0.05, 0.5, "Shape", rotation=90, va="center", ha="center", fontsize=18)
fig.text(0.5, -.025, "Cosmic time [Gyr]", ha="center", fontsize=18)

ax_f.text(
    0.5, 1.06, "Major Mergers", ha="center", transform=ax_f.transAxes, fontweight="bold", fontsize=14
)
ax_i.text(
    0.5, 1.06, "Minor Mergers", ha="center", transform=ax_i.transAxes, fontweight="bold", fontsize=14
)

plt.figlegend(ncol=5, loc="lower center", bbox_to_anchor=(0.5-(.8/2), -0.15, .8, .075))

plt.savefig("shape_over_z_variable_radius.png", bbox_inches="tight")
plt.show()
