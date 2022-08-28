import asdf
import gizmo_analysis as gizmo
import os.path as path
import orientations as ori
import numpy as np

simulation_list = [
    "../../../data/latte_metaldiff/m12f_res7100",
    "../../../data/latte_metaldiff/m12i_res7100",
    "../../../data/latte_metaldiff/m12m_res7100",
    "../../../data/latte_metaldiff/m12w_res7100",
]

radii_list = [
    "../scale_radii/scale_radii_z_m12f_res7100.asdf",
    "../scale_radii/scale_radii_z_m12i_res7100.asdf",
    "../scale_radii/scale_radii_z_m12m_res7100.asdf",
    "../scale_radii/scale_radii_z_m12w_res7100.asdf",
]

for k, sim in enumerate(simulation_list):
    sim_name = path.split(sim)[-1]
    radii = asdf.open(radii_list[k])
    snapshot_range = radii["snapshot"]
        
    af = asdf.AsdfFile({
        "snapshot": snapshot_range,
        "v.sigma.ratio": np.zeros(len(snapshot_range)),
    })
    
    for i, snap in enumerate(snapshot_range):

        print(f"Calculating ratio for {sim_name} at snapshot {snap}...")
        
        virial_radius = radii["virial"][i]
        
        part = gizmo.io.Read.read_snapshots("gas", "snapshot", snap, sim)
        
        v_phi = part["gas"].prop('host.velocity.principal.cylindrical')[:, 1]
        mass  = part["gas"].prop("mass.hydrogen.neutral")
        dist  = part["gas"].prop("host.distance.principal.total")

        
        ratio = ori.vsigma_ratio(v_phi, mass, dist, virial_radius)
        
        af["v.sigma.ratio"][i] = ratio
        
        print(f"Calculated {sim_name} @ {snap}: v.sigma.ratio = {ratio}")
    
    af.write_to(f"v_sigma_ratios_{sim_name}.asdf")