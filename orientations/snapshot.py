import utilities as ut
from astropy.io.ascii import read
import os.path as path


def getSnapshotData(sim_dir, snapshot):
    # For Latte simulations
    cosmology = ut.cosmology.CosmologyClass(source="agora")

    snapshot_path = path.join(sim_dir, "snapshot_times.txt")

    snapshot_data = read(snapshot_path, format="commented_header", header_start=2)

    snapshot_time = snapshot_data[snapshot]["time[Gyr]"]

    snapshot_redshift = snapshot_data[snapshot]["redshift"]

    return snapshot_time, snapshot_redshift
