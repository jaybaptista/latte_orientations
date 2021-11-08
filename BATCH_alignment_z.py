from astropy.io.misc.hdf5 import read_table_hdf5
import numpy as np
import pandas as pd
import gizmo_analysis as gizmo
from ob import getSymmetryAxes, getMinAngle

gals = [
    'm12f_res7100',
    'm12i_res7100',
    'm12w_res7100',
    'm12m_res7100'
]

config = {
    'm12f': np.arange(302, 524, step=6),
    'm12i': np.arange(200, 480, step=6),
    'm12m': np.arange(153, 564, step=10),
    'm12w': np.arange(152, 380, step=6)
}

for sim in gals:
    sim_dir = '../data/latte_metaldiff/' + sim

    radii_data = pd.read_hdf('radii_{}.hdf'.format(sim))
    r_vir = radii_data['virial'][0]
    rs90 = radii_data['star.radius.90'][0]
    
    lmc_position = read_table_hdf5('lmc_positions_{}.h5'.format(sim))
    
#     max_snap = np.max(lmc_position['snapshot'])
#     min_snap = np.min(lmc_position['snapshot'])
    
#     step = 50
    
#     if (sim == 'm12i_res7100' or  sim == 'm12m_res7100'):
#         step = 60
    
#     if sim == 'm12i_res7100':
#         min_snap = min_snap + 1

    radii = [r_vir, rs90, 5*rs90]
#     snapshot_range = np.arange(min_snap, max_snap, step = step)

    snapshot_range = config[sim[:4]]

    df = {
        'virial': [],
        'disk': [], # disk is star.radius.90
        '5.disk': []
    }

    for k in np.arange(0,3):

        for snap in snapshot_range:

            part = gizmo.io.Read.read_snapshots(['dark'], 'snapshot', snap, sim_dir)

            positions = part['dark'].prop('host.distance')
            dists = part['dark'].prop('host.distance.total')
            
            z_lmc = lmc_position[np.where(lmc_position['snapshot'] == snap)[0][0]]['position']
            z_lmc_hat = z_lmc / np.linalg.norm(z_lmc)

            tensor = getSymmetryAxes(positions, dists, radius = radii[k])
            angle = getMinAngle(tensor[0], z_lmc_hat) * 180/np.pi

            df[list(df.keys())[k]].append(angle)

    d = pd.DataFrame(df)
    d['snapshot'] = snapshot_range
    d.to_hdf('data_lmcAlignment_z_{}.h5'.format(sim), 'w')