import asdf
import matplotlib.pyplot as plt
import numpy as np
import utilities as ut
import orientations as ori

plt.rcParams.update({"font.family": "serif"})

config = {
    "m12f_res7100": np.arange(302, 524, step=5),
    "m12i_res7100": np.arange(200, 480, step=5),
    "m12m_res7100": np.arange(153, 564, step=5),
    "m12w_res7100": np.arange(152, 380, step=5),
}

############

m12f_data = asdf.open("lmc_alignment_m12f_res7100.asdf")
m12i_data = asdf.open("lmc_alignment_m12i_res7100.asdf")
m12m_data = asdf.open("lmc_alignment_m12m_res7100.asdf")
m12w_data = asdf.open("lmc_alignment_m12w_res7100.asdf")

af_md = asdf.open("peri_md.asdf")

m12f_pos_data = asdf.open("lmc_positions_m12f_res7100.asdf")
m12i_pos_data = asdf.open("lmc_positions_m12i_res7100.asdf")
m12m_pos_data = asdf.open("lmc_positions_m12m_res7100.asdf")
m12w_pos_data = asdf.open("lmc_positions_m12w_res7100.asdf")

m12f_dist = []
m12i_dist = []
m12m_dist = []
m12w_dist = []

for m12f_snap in config["m12f_res7100"]:
    m12f_dist.append(np.sqrt(np.sum(m12f_pos_data["position"][m12f_snap] ** 2)))

for m12i_snap in config["m12i_res7100"]:
    m12i_dist.append(np.sqrt(np.sum(m12i_pos_data["position"][m12i_snap] ** 2)))

for m12m_snap in config["m12m_res7100"]:
    m12m_dist.append(np.sqrt(np.sum(m12m_pos_data["position"][m12m_snap] ** 2)))

for m12w_snap in config["m12w_res7100"]:
    m12w_dist.append(np.sqrt(np.sum(m12w_pos_data["position"][m12w_snap] ** 2)))

m12f_m = np.array(m12f_dist) > 0.0
m12i_m = np.array(m12i_dist) > 0.0
m12m_m = np.array(m12m_dist) > 0.0
m12w_m = np.array(m12w_dist) > 0.0

m12f_dist = np.array(m12f_dist)[m12f_m]
m12i_dist = np.array(m12i_dist)[m12i_m]
m12m_dist = np.array(m12m_dist)[m12m_m]
m12w_dist = np.array(m12w_dist)[m12w_m]

m12f_time = ori.getSnapshotData(
    "../../../data/latte_metaldiff/m12i_res7100/", np.array(m12f_data["snapshot"])
)[0][m12f_m]
m12i_time = ori.getSnapshotData(
    "../../../data/latte_metaldiff/m12i_res7100/", np.array(m12i_data["snapshot"])
)[0][m12i_m]
m12m_time = ori.getSnapshotData(
    "../../../data/latte_metaldiff/m12i_res7100/", np.array(m12m_data["snapshot"])
)[0][m12m_m]
m12w_time = ori.getSnapshotData(
    "../../../data/latte_metaldiff/m12i_res7100/", np.array(m12w_data["snapshot"])
)[0][m12w_m]

############

disk_col = "#0a9396"
disk5_col = "#005f73"
vir_col = "#001219"
dot_size = 5

############ LMC-ALIGNMENT vs. REDSHIFT ##############

fig, ax = plt.subplots(2, 2, dpi=200, figsize=(8, 4))

((ax_f, ax_i), (ax_w, ax_m)) = ax

ax_i.plot(
    m12i_time,
    np.array(m12i_data["disk"])[m12i_m],
    c="b",
    alpha=0.25,
    label="Stellar Disk ($R_{*,90}$)",
    lw=3,
)
ax_i.plot(
    m12i_time,
    np.array(m12i_data["5.disk"])[m12i_m],
    c="b",
    alpha=0.5,
    label="5$ \\times R_{*,90}$",
    lw=3,
)
ax_i.plot(
    m12i_time,
    np.array(m12i_data["virial"])[m12i_m],
    c="b",
    label="Virial Radius $(R_{200})$",
    lw=3,
)

ax_f.plot(
    m12f_time,
    np.array(m12f_data["disk"])[m12f_m],
    c="b",
    alpha=0.25,
    lw=3,
)
ax_f.plot(
    m12f_time,
    np.array(m12f_data["5.disk"])[m12f_m],
    c="b",
    alpha=0.5,
    lw=3,
)
ax_f.plot(
    m12f_time,
    np.array(m12f_data["virial"])[m12f_m],
    c="b",
    lw=3,
)

ax_m.plot(m12m_time, np.array(m12m_data["disk"])[m12m_m], c="b", alpha=0.25, lw=3)
ax_m.plot(m12m_time, np.array(m12m_data["5.disk"])[m12m_m], c="b", alpha=0.5, lw=3)
ax_m.plot(m12m_time, np.array(m12m_data["virial"])[m12m_m], c="b", lw=3)

ax_w.plot(m12w_time, np.array(m12w_data["disk"])[m12w_m], c="b", alpha=0.25, lw=3)
ax_w.plot(m12w_time, np.array(m12w_data["5.disk"])[m12w_m], c="b", alpha=0.5, lw=3)
ax_w.plot(m12w_time, np.array(m12w_data["virial"])[m12w_m], c="b", lw=3)

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

ax_m.text(0.05, 0.1, "m12m", transform=ax_m.transAxes)
ax_i.text(0.05, 0.1, "m12i", transform=ax_i.transAxes)
ax_f.text(0.05, 0.1, "m12f", transform=ax_f.transAxes)
ax_w.text(0.05, 0.1, "m12w", transform=ax_w.transAxes)

# ax_i.legend(bbox_to_anchor=(1.05, 0.45), fontsize=9)

for a in ax.flatten():
    a.set_ylim(0, 90)

fig.subplots_adjust(hspace=0.25, wspace=0.1)

fig.text(0.05, 0.5, "Alignment to Satellite [deg]", rotation=90, va="center", ha="center", fontsize=16)
fig.text(0.5, -.025, "Cosmic time [Gyr]", ha="center", fontsize=18)

ax_f.text(
    0.5, 1.05, "Major Mergers", ha="center", transform=ax_f.transAxes, fontweight="bold", fontsize=14
)
ax_i.text(
    0.5, 1.05, "Minor Mergers", ha="center", transform=ax_i.transAxes, fontweight="bold", fontsize=14
)
plt.figlegend(ncol=5, loc="lower center", bbox_to_anchor=(0.5-(.8/2), -0.15, .8, .075))
plt.savefig("lmc_alignment.png", bbox_inches="tight")

############ LMC-ALIGNMENT vs. SATELLITE DISTANCE ##############

fig, ax = plt.subplots(2, 2, dpi=200, figsize=(8, 4))

((ax_f, ax_i), (ax_w, ax_m)) = ax

ax_i.scatter(
    m12i_dist,
    np.array(m12i_data["disk"])[m12i_m],
    c=disk_col,
    label="Stellar Disk ($R_{*,90}$)",
    ec=disk_col,
    fc="None",
    s=dot_size,
)
ax_i.scatter(
    m12i_dist,
    np.array(m12i_data["5.disk"])[m12i_m],
    c=disk5_col,
    label="5$ \\times R_{*,90}$",
    ec=disk5_col,
    fc="None",
    s=3 * dot_size,
)
ax_i.scatter(
    m12i_dist,
    np.array(m12i_data["virial"])[m12i_m],
    c=vir_col,
    label="Virial Radius $(R_{200})$",
    ec=vir_col,
    fc="None",
    s=6 * dot_size,
)

ax_f.scatter(
    m12f_dist,
    np.array(m12f_data["disk"])[m12f_m],
    c=disk_col,
    ec=disk_col,
    fc="None",
    s=dot_size,
)
ax_f.scatter(
    m12f_dist,
    np.array(m12f_data["5.disk"])[m12f_m],
    c=disk5_col,
    ec=disk5_col,
    fc="None",
    s=3 * dot_size,
)
ax_f.scatter(
    m12f_dist,
    np.array(m12f_data["virial"])[m12f_m],
    c=vir_col,
    ec=vir_col,
    fc="None",
    s=6 * dot_size,
)

ax_m.scatter(
    m12m_dist, np.array(m12m_data["disk"])[m12m_m], ec=disk_col, fc="None", s=dot_size
)
ax_m.scatter(
    m12m_dist,
    np.array(m12m_data["5.disk"])[m12m_m],
    ec=disk5_col,
    fc="None",
    s=3 * dot_size,
)
ax_m.scatter(
    m12m_dist,
    np.array(m12m_data["virial"])[m12m_m],
    ec=vir_col,
    fc="None",
    s=6 * dot_size,
)

ax_w.scatter(
    m12w_dist, np.array(m12w_data["disk"])[m12w_m], ec=disk_col, fc="None", s=dot_size
)
ax_w.scatter(
    m12w_dist,
    np.array(m12w_data["5.disk"])[m12w_m],
    ec=disk5_col,
    fc="None",
    s=3 * dot_size,
)
ax_w.scatter(
    m12w_dist,
    np.array(m12w_data["virial"])[m12w_m],
    ec=vir_col,
    fc="None",
    s=6 * dot_size,
)

ax_m.text(
    0.05,
    0.85,
    "m12m",
    transform=ax_m.transAxes,
    bbox=dict(boxstyle="round", fc="#EEE", ec="#DDD", alpha=0.95),
)
ax_i.text(
    0.05,
    0.85,
    "m12i",
    transform=ax_i.transAxes,
    bbox=dict(boxstyle="round", fc="#EEE", ec="#DDD", alpha=0.95),
)
ax_f.text(
    0.05,
    0.85,
    "m12f",
    transform=ax_f.transAxes,
    bbox=dict(boxstyle="round", fc="#EEE", ec="#DDD", alpha=0.95),
)
ax_w.text(
    0.05,
    0.85,
    "m12w",
    transform=ax_w.transAxes,
    bbox=dict(boxstyle="round", fc="#EEE", ec="#DDD", alpha=0.95),
)

# ax_i.legend(bbox_to_anchor=(1.05, 0.45), fontsize=9)

for a in ax.flatten():
    a.set_ylim(0, 90)

fig.subplots_adjust(hspace=0.25, wspace=0.1)

fig.text(0.05, 0.5, "Alignment to Satellite [deg]", rotation=90, va="center", ha="center", fontsize=16)
fig.text(0.5, -.025, "Satellite Distance [kpc]", ha="center", fontsize=18)

ax_f.text(
    0.5, 1.05, "Major Mergers", ha="center", transform=ax_f.transAxes, fontweight="bold", fontsize=14
)
ax_i.text(
    0.5, 1.05, "Minor Mergers", ha="center", transform=ax_i.transAxes, fontweight="bold", fontsize=14
)
plt.figlegend(ncol=3, loc="lower center", bbox_to_anchor=(0.5-(.5/2), -0.15, .5, .075))
plt.savefig("lmc_alignment_distance.png", bbox_inches="tight")

# diag_af = asdf.AsdfFile({"m12m_dist": m12m_dist})
# diag_af.write_to("diag.asdf")
