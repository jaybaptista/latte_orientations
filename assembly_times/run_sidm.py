import os.path as path
import orientations as ori
import gizmo_analysis as gizmo
import numpy as np
import asdf

simulation_list = [
    "../../../data/m12f_sidm1",
    "../../../data/m12i_sidm1",
    "../../../data/m12m_sidm1",
]

af = asdf.AsdfFile({})

for sim in simulation_list:

    sim_name = path.split(sim)[-1]

    part = gizmo.io.Read.read_snapshots(
        ["star"],
        "redshift",
        0,
        sim,
        assign_hosts_rotation=True,
        assign_formation_coordinates=True,
    )

    assembly_z15, assembly_f15 = ori.calculateAssemblyTime(part, assembly_radius=15)
    assembly_z2, assembly_f2 = ori.calculateAssemblyTime(part, assembly_radius=2)

    af[sim_name] = {
        "z.15": assembly_z15,
        "frac.15": assembly_f15,
        "z.2": assembly_z2,
        "frac.2": assembly_f2,
    }

af.write_to("assembly_times_sidm.asdf")
