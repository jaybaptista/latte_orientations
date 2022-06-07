import gizmo_analysis as gizmo
import numpy as np
from astropy.io.ascii import read
import jmb_utils as jb
import pandas as pd

sim_dir = '../../../data/latte_metaldiff/m12i_res7100/'
part = gizmo.io.Read.read_snapshots(['star'], 'redshift', 0, sim_dir, assign_hosts_rotation=True, assign_formation_coordinates=True,)

ratios = jb.calculateAssemblyTime(
    part,
    assembly_radius=15
)

df = pd.DataFrame({
    'ratios': ratios
})

df.to_hdf('m12i.hdf', key='w')