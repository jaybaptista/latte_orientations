import utilities as ut
from astropy.io.ascii import read
from astropy.io.ascii import NoHeader
import os.path as path


def getSnapshotData(sim_dir, snapshot):
    # For Latte simulations
    cosmology = ut.cosmology.CosmologyClass(source="agora")

    snapshot_path = path.join(sim_dir, "snapshot_times.txt")

    snapshot_data = read(snapshot_path, format="commented_header", header_start=2)

    snapshot_time = snapshot_data[snapshot]["time[Gyr]"]

    snapshot_redshift = snapshot_data[snapshot]["redshift"]

    return snapshot_time, snapshot_redshift

def getSnapshotData_SIDM(sim_dir, snapshot):
    # For Latte simulations
    cosmology = ut.cosmology.CosmologyClass(source="agora")
    
    snapshot_path = path.join(sim_dir, "snapshot_scale-factors.txt")
    
    snapshot_data = read(snapshot_path, Reader=NoHeader)
    
    snapshot_scale_factor = snapshot_data[snapshot]
    
    snapshot_time = cosmology.convert_time('time', 'scalefactor', np.array(snapshot_scale_factor).astype('float'))
    
    snapshot_redshift = cosmology.convert_time('redshift', 'time', snapshot_time)
    
    return snapshot_time, snapshot_redshift