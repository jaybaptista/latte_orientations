import halo_analysis as halo
import gizmo_analysis as gizmo
import numpy as np
import asdf
import os.path as path


# Track a halo's `host.position` given a reference tree index.
def trackHaloPosition(sim_dir, ref_idx, write=False):

    sim_name = path.split(sim_dir)[-1]

    df = {
        "fs": [],  # "forward snap"
        "fc": [],  # "forward coordinates"
        "bs": [],  # "backward snap"
        "bc": [],  # "backward coordinates"
    }

    halt = halo.io.IO.read_tree(simulation_directory=sim_dir)

    tidx = ref_idx

    # Forward tracking
    print("Forward tracking halo from tree index...")
    while tidx > 0:
        df["fs"].append(halt["snapshot"][tidx])
        df["fc"].append(halt["host.distance"][tidx])
        tidx = halt["descendant.index"][tidx]

    # Backward tracking
    # Start tracking the first progenitor halo (to prevent duplication of the reference snapshot)
    print("Backward tracking halo from tree index...")
    tidx = halt["progenitor.main.index"][ref_idx]
    while tidx > 0:
        df["bs"].append(halt["snapshot"][tidx])
        df["bc"].append(halt["host.distance"][tidx])
        tidx = halt["progenitor.main.index"][tidx]

    # Reverse order of backward tracked positions/snapshots
    df["bs"] = np.flipud(df["bs"])
    df["bc"] = np.flipud(df["bc"])

    # Combine datasets in chronological order.
    af = asdf.AsdfFile(
        {
            "snapshot": np.append(df["bs"], df["fs"]),
            "position": np.vstack([df["bc"], df["fc"]]),
        }
    )

    if write:
        print(f"Writing LMC positions for {sim_name}...")
        af.write_to(f"lmc_positions_{sim_name}.asdf")

    return af


def assemble_test(sim_dir, assembly_radius):
    
    redshifts = [
        0, 0.1, 0.2, 0.3, 0.4,
        0.5, 0.6, 0.7, 0.8, 0.9,
        1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7,
        1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4,
        2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4, 4.2, 4.4, 4.6, 4.8,
        5, 5.2, 5.4, 5.6, 5.8, 6, 7, 8, 9,
    ]
    
    part = gizmo.io.Read.read_snapshots('star', 'redshift', 0, simulation_directory=sim_dir, assign_host_rotation=True)
    hal = halo.io.IO.read_catalogs('redshift', redshifts, file_kind='hdf5', simulation_directory=sim_dir, all_snapshot_list=False)
    
    pass
