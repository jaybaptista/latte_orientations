import halo_analysis as halo
import numpy as np
import pandas as pd
import utilities as ut
import gizmo_analysis as gizmo
from astropy.io.ascii import read
from astropy.cosmology import Planck13, z_at_value
import astropy.units as u
import matplotlib.pyplot as plt
from abg_python.smooth_utils import boxcar_average


def getScaleRadii(sim_dir, snapshot):

    hal = halo.io.IO.read_catalogs("snapshot", [snapshot], sim_dir)

    host_idx = hal["host.index"][0]

    r_vir = hal["radius"][host_idx]
    r_star_50 = hal["star.radius.50"][host_idx]
    r_star_90 = hal["star.radius.90"][host_idx]

    df = {
        "virial": [r_vir],
        "star.radius.50": [r_star_50],
        "star.radius.90": [r_star_90],
    }

    return pd.DataFrame(df)


def getAngularMomentum(pos, dist, vel, radius=10):
    # Note—this is normalized angular momentum (this is mass is not included)
    mask = dist < radius

    L_vector = np.median(np.cross(pos[mask], vel[mask]), axis=0)
    L_vector = L_vector / np.linalg.norm(L_vector)

    return L_vector


def getSymmetryAxes(pos, dist, radius, eigenvals=False):

    print(f"Calculating MOI tensor for r = {radius} kpc")
    mask = dist < radius
    pos = pos[mask]

    # Get moment of inertia tensor
    xx = np.sum(pos[:, 0] ** 2 / (pos[:, 0] ** 2 + pos[:, 1] ** 2 + pos[:, 2] ** 2))
    yy = np.sum(pos[:, 1] ** 2 / (pos[:, 0] ** 2 + pos[:, 1] ** 2 + pos[:, 2] ** 2))
    zz = np.sum(pos[:, 2] ** 2 / (pos[:, 0] ** 2 + pos[:, 1] ** 2 + pos[:, 2] ** 2))
    xy = yx = np.sum(
        pos[:, 0] * pos[:, 1] / (pos[:, 0] ** 2 + pos[:, 1] ** 2 + pos[:, 2] ** 2)
    )
    xz = zx = np.sum(
        pos[:, 0] * pos[:, 2] / (pos[:, 0] ** 2 + pos[:, 1] ** 2 + pos[:, 2] ** 2)
    )
    zy = yz = np.sum(
        pos[:, 2] * pos[:, 1] / (pos[:, 0] ** 2 + pos[:, 1] ** 2 + pos[:, 2] ** 2)
    )

    Im = [[xx, xy, xz], [yx, yy, yz], [zx, zy, zz]]

    e_val, e_vec = np.linalg.eig(Im)
    sorted_idx = np.argsort(e_val)[
        ::-1
    ]  # sort from largest to smallest (idx of 2 is minor axis)
    sorted_val = e_val[sorted_idx]
    sorted_vec = np.transpose(e_vec)[sorted_idx]

    if eigenvals:
        return sorted_vec, sorted_val
    else:
        return sorted_vec


def permuteVector(vec, ref):
    if np.dot(ref, vec) < 0.0:
        return -1 * np.array(vec)
    else:
        return vec


def getMinAngle(vec, ref):
    # ref is the reference vector to perform corrections
    # vec is the vector of interest

    vec = vec / np.linalg.norm(vec)
    ref = ref / np.linalg.norm(ref)

    angle = np.min([np.arccos(np.dot(-1 * ref, vec)), np.arccos(np.dot(ref, vec))])

    return angle


def getHaloOrientation(
    simulation_directory, snapshot, radii, reference_vector=None, return_tensor=False
):
    part = gizmo.io.Read.read_snapshots(
        ["dark"], "snapshot", snapshot, simulation_directory
    )
    pos = part["dark"].prop("host.distance")
    dist = part["dark"].prop("host.distance.total")
    host_minor_axis = part.host["rotation"][0, 2]

    if reference_vector == None:
        reference_vector = host_minor_axis

    symmetry_axes = np.array([getSymmetryAxes(pos, dist, r) for r in radii])
    reference_axes = symmetry_axes[:, 2]
    angles = np.array(
        [
            getMinAngle(minor_axis, reference_vector.reshape(1, 3)) * 180 / np.pi
            for minor_axis in reference_axes
        ]
    )

    df = pd.DataFrame({"radius": radii, "angle": angles})

    if return_tensor:
        df["tensor"] = symmetry_axes

    return df


def obtainSFR(
    simulation_directory,
    snapshot,
    rmax=None,
    dt=1 / 1000,  # 1 Myr in Gyr
    cosmology=None,
    tmax=13.7965284579874688,  # this is for Latte galaxies
    tmin_sfr=0.1,  # Gyr, this is here because there is a bug I'm still working out :/
):

    if cosmology == None:
        cosmology = ut.cosmology.CosmologyClass(source="agora")

    if rmax == None:
        hal = halo.io.IO.read_catalogs("snapshot", snapshot, simulation_directory)
        rmax = 5 * hal["star.radius.50"][hal["host.index"][0]]

    part = gizmo.io.Read.read_snapshots(["star"], "redshift", 0, simulation_directory)

    # select stars within a particular radius
    mask = part["star"].prop("host.distance.total") <= rmax
    ages = part["star"].prop("age")[mask]
    formation_time = part["star"].prop("form.time")[mask]
    formation_mass = part["star"].prop("form.mass")[mask]

    time_edges = np.arange(tmax, 0, -dt)[::-1]
    SFRs, time_edges = np.histogram(
        formation_time, weights=formation_mass / (dt * 1e9), bins=time_edges
    )

    t_mask = time_edges[:-1] >= tmin_sfr

    z_edges = cosmology.convert_time("redshift", "time", time_edges[:-1][t_mask])

    df = pd.DataFrame(
        {"time": time_edges[:-1][t_mask], "redshift": z_edges, "sfr": SFRs[t_mask]}
    )

    return df


def calculateRadius(
    positions,
    masses,
    rmax,
    cdf_thresh=0.5,
    bins=5000,
):
    radii = np.sum(positions**2, axis=1) ** 0.5

    m = radii <= rmax

    edges = np.linspace(0, rmax, bins, endpoint=True)

    h, edges = np.histogram(radii[m], bins=edges, weights=masses[m])

    h /= 1.0 * np.sum(h)

    cdf = np.cumsum(h)

    intersection_idx = np.argmin((cdf - cdf_thresh) ** 2)

    r_intersect = edges[1:][intersection_idx]
    cdf_intersect = cdf[intersection_idx]

    return r_intersect, cdf_intersect


def calculateAssemblyTime(part, assembly_radius, threshold=.5, disk_size=30):

    form_z     = part['star'].prop('form.redshift')
    form_radii = part['star'].prop('form.host.distance.total')
    radii      = part['star'].prop('host.distance.total')
    mass       = part['star'].prop('mass')

    redshifts = np.linspace(0.001, 10, 50)

    dmask = (radii < assembly_radius)
    dform_mask = form_radii < disk_size
    
    cdf = []
    
    for z in redshifts:
        zmask = form_z > z

        m1 = dmask & zmask & dform_mask
        m2 = dmask & zmask

        total_z_mass = np.sum(mass[m1])
        total_z0_mass = np.sum(mass[m2])

        f = total_z_mass/total_z0_mass

        cdf.append(f)
     
    intersection_idx = np.argmin((cdf - threshold) ** 2)
    z_intersect = redshifts[1:][intersection_idx]
    cdf_intersect = cdf[intersection_idx]
    
    return z_intersect, cdf_intersect