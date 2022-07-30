import halo_analysis as halo
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
