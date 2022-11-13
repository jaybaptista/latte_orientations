import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import asdf

plt.rcParams.update({"font.family": "serif"})

m12f_data = pd.read_hdf("orientations_present_day_m12f_res7100.hdf")
m12i_data = pd.read_hdf("orientations_present_day_m12i_res7100.hdf")
m12m_data = pd.read_hdf("orientations_present_day_m12m_res7100.hdf")
m12w_data = pd.read_hdf("orientations_present_day_m12w_res7100.hdf")

m12f_data_shape = asdf.open("../shapes_present_day_cdm/i100 run/shapes_present_day_m12f_res7100.asdf")
m12i_data_shape = asdf.open("../shapes_present_day_cdm/i100 run/shapes_present_day_m12i_res7100.asdf")
m12m_data_shape = asdf.open("../shapes_present_day_cdm/i100 run/shapes_present_day_m12m_res7100.asdf")
m12w_data_shape = asdf.open("../shapes_present_day_cdm/i100 run/shapes_present_day_m12w_res7100.asdf")

dr = rmin = 2
rmax = 400 + dr
radii = np.arange(rmin, rmax, step=dr)

fig, ax = plt.subplots(2, 1, dpi=250)

# Gaia resolvable

gaia_m = radii < 60
as_m = (radii > 60) & (radii < 200)
lsst_m = (radii > 200) & (radii < 400)

ax[0].plot(radii[gaia_m], m12f_data["angle"][gaia_m], label="m12f", c="blue")
ax[0].plot(radii[gaia_m], m12i_data["angle"][gaia_m], label="m12i", c="red")
ax[0].plot(radii[gaia_m], m12m_data["angle"][gaia_m], label="m12m", c="green")
ax[0].plot(radii[gaia_m], m12w_data["angle"][gaia_m], label="m12w", c="orange")

# AS resolvable

ax[0].plot(radii[as_m], m12f_data["angle"][as_m], ls="-.", c="blue")
ax[0].plot(radii[as_m], m12i_data["angle"][as_m], ls="-.", c="red")
ax[0].plot(radii[as_m], m12m_data["angle"][as_m], ls="-.", c="green")
ax[0].plot(radii[as_m], m12w_data["angle"][as_m], ls="-.", c="orange")

# LSST resolvable
ax[0].plot(radii[lsst_m], m12f_data["angle"][lsst_m], ls="dotted", c="blue")
ax[0].plot(radii[lsst_m], m12i_data["angle"][lsst_m], ls="dotted", c="red")
ax[0].plot(radii[lsst_m], m12m_data["angle"][lsst_m], ls="dotted", c="green")
ax[0].plot(radii[lsst_m], m12w_data["angle"][lsst_m], ls="dotted", c="orange")


ax[0].axhline(20, c="k", ls="--", alpha=0.25)

ax[0].set_ylabel("Orientation w.r.t. to disk [deg]")

# ax[0].legend()

gaia_m_t = m12f_data_shape["rmax"] < 60
as_m_t = (m12f_data_shape["rmax"] > 60) & (m12f_data_shape["rmax"] < 200)
lsst_m_t = (m12f_data_shape["rmax"] > 200) & (m12f_data_shape["rmax"] < 400)

# Gaia resolvable
ax[1].plot(0, 0, label="Gaia", c="k")
ax[1].plot(m12f_data_shape["rmax"][gaia_m_t], m12f_data_shape["shape"][gaia_m_t], c="blue")
ax[1].plot(m12i_data_shape["rmax"][gaia_m_t], m12i_data_shape["shape"][gaia_m_t], c="red")
ax[1].plot(m12m_data_shape["rmax"][gaia_m_t], m12m_data_shape["shape"][gaia_m_t], c="green")
ax[1].plot(m12w_data_shape["rmax"][gaia_m_t], m12w_data_shape["shape"][gaia_m_t], c="orange")

# AS resolvable
ax[1].plot(0, 0, label="Asteroseismic", ls="-.", c="k")
ax[1].plot(m12f_data_shape["rmax"][as_m_t], m12f_data_shape["shape"][as_m_t], ls="-.", c="blue")
ax[1].plot(m12i_data_shape["rmax"][as_m_t], m12i_data_shape["shape"][as_m_t], ls="-.", c="red")
ax[1].plot(m12m_data_shape["rmax"][as_m_t], m12m_data_shape["shape"][as_m_t], ls="-.", c="green")
ax[1].plot(m12w_data_shape["rmax"][as_m_t], m12w_data_shape["shape"][as_m_t], ls="-.", c="orange")

# LSST resolvable
ax[1].plot(0, 0, label="LSST", ls="dotted", c="k")
ax[1].plot(m12f_data_shape["rmax"][lsst_m_t], m12f_data_shape["shape"][lsst_m_t], ls="dotted", c="blue")
ax[1].plot(m12i_data_shape["rmax"][lsst_m_t], m12i_data_shape["shape"][lsst_m_t], ls="dotted", c="red")
ax[1].plot(m12m_data_shape["rmax"][lsst_m_t], m12m_data_shape["shape"][lsst_m_t], ls="dotted", c="green")
ax[1].plot(m12w_data_shape["rmax"][lsst_m_t], m12w_data_shape["shape"][lsst_m_t], ls="dotted", c="orange")

ax[1].axhline(1/3, c="k", alpha=0.25)
ax[1].axhline(2/3, c="k", alpha=0.25)

ax[1].text(.025, 1   - .15, "Prolate", transform=ax[1].transAxes, size=12)
ax[1].text(.025, 2/3 - .15, "Triaxial", transform=ax[1].transAxes, size=12)
ax[1].text(.025, 1/3 - .15, "Oblate", transform=ax[1].transAxes, size=12)

ax[1].set_xlabel("Radius [kpc]")
ax[1].set_ylabel("Triaxiality")

plt.figlegend(ncol=4, loc="lower center", bbox_to_anchor=(0.5-(.8/2), -0.15, .8, .075))

fig.subplots_adjust(hspace=0)

plt.savefig("orientations_present_day_cdm_i100.png", bbox_inches="tight")
