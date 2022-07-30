import orientations as ori
import numpy as np
import gizmo_analysis as gizmo
import asdf  # trying this out for the first time...

radii = [10, 30, 50, 100, 150, 200, 250, 300, 350, 400]

sim_dir = "../../../data/latte_metaldiff/"
sims_name = "m12i_res7100"

part = gizmo.io.Read.read_snapshots(
    ["dark", "star"], "redshift", 0, sim_dir + sims_name, assign_hosts_rotation=True
)

pos_dark = part["dark"].prop("host.distance")
pos_star = part["star"].prop("host.distance")

dist_dark = part["dark"].prop("host.distance.total")
dist_star = part["star"].prop("host.distance.total")

vel_dark = part["dark"].prop("host.velocity")
vel_star = part["star"].prop("host.velocity")

####

af = {}

####
# Get the stellar disk angular momentum vector:

print("Calculating angular momentum unit vector...")

L = ori.getAngularMomentum(pos_star, dist_star, vel_star, radius=10)

af["angular.momentum"] = L

print(f"Angular momentum unit vector calculated: {L}")

####
# Get the symmetry axes and eigenvalues

print("Calculating symmetry axes...")

symm_axes = []
eigenvalues = []

for r in radii:

    print(f"... @ {r} kpc")

    axes, eigs = ori.getSymmetryAxes(pos_dark, dist_dark, radius=r, eigenvals=True)

    symm_axes.append(axes)
    eigenvalues.append(eigs)

af["radius"] = radii
af["symmetry.axes"] = symm_axes
af["eigenvalues"] = eigenvalues

####
# Put it all together!

af = asdf.AsdfFile(af)
af.write_to("m12i_symmetry_axes.asdf")
