import numpy as np
import halo_analysis as halo
from astropy.io.ascii import read


def getAngularMomentum(pos, dist, vel, radius=10):

    mask = dist < radius
    pos_sel = pos[mask]

    # L = r x v (direction invariant of mass)
    L_vector = np.median(np.cross(pos, vel), axis=0)
    L_hat = L_vector / np.linalg.norm(L_vector)  # normalize
    return L_vector


def getSymmetryAxes(pos, dist, radius):
    print("Assigning inertia tensor at r = ", radius, "kpc")
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

    # Enforce RHR onto an axis
    #     sorted_vec[2] = np.cross(sorted_vec[0], sorted_vec[1])

    return sorted_vec


def getMinAngle(vec, ref):
    # ref is the reference vector
    # vec is the vector of interest

    vec = vec / np.linalg.norm(vec)
    ref = ref / np.linalg.norm(ref)

    angle = np.min([np.arccos(np.dot(-1 * ref, vec)), np.arccos(np.dot(ref, vec))])

    return angle


# Track a halo's `host.position` given a reference tree index.
def trackHaloPosition(sim_dir, ref_idx):

    df = {
        "fs": [],
        "fc": [],
        "bs": [],
        "bc": [],
    }

    halt = halo.io.IO.read_tree(simulation_directory=sim_dir)

    tidx = ref_idx

    # Forward tracking
    while tidx > 0:
        df["fs"].append(halt["snapshot"][tidx])
        df["fc"].append(halt["host.distance"][tidx])
        tidx = halt["descendant.index"][tidx]

    # Backward tracking
    # Start tracking the first progenitor halo (to prevent duplication of the reference snapshot)
    tidx = halt["progenitor.main.index"][ref_idx]
    while tidx > 0:
        df["bs"].append(halt["snapshot"][tidx])
        df["bc"].append(halt["host.distance"][tidx])
        tidx = halt["progenitor.main.index"][tidx]

    # Reverse order of backward tracked positions/snapshots
    df["bs"] = np.flipud(df["bs"])
    df["bc"] = np.flipud(df["bc"])

    # Combine datasets in chronological order.
    dat = {
        "snapshot": np.append(df["bs"], df["fs"]),
        "position": np.vstack([df["bc"], df["fc"]]),
    }

    return dat


# Returns a simulation time in Gyr from given redshift or snapshot
def cvtRedshift(z, sim="../../data/latte_metaldiff/m12i_res7100"):
    times = sim + "/snapshot_times.txt"
    snapshot_data = read(times, format="commented_header", header_start=2)
    return snapshot_data[np.argmin(abs(snapshot_data["redshift"] - z))]["time[Gyr]"]


def cvtSnapshot(snapshot, sim="../../data/latte_metaldiff/m12i_res7100"):
    times = sim + "/snapshot_times.txt"
    snapshot_data = read(times, format="commented_header", header_start=2)
    return snapshot_data[snapshot]["time[Gyr]"]
