import asdf
import utilities as ut
from astropy.cosmology import Planck13, z_at_value
import os.path as path
import halo_analysis as halo
import numpy as np

cosmology = ut.cosmology.CosmologyClass(source="agora")

pericenter_data = {
    "m12f_res7100": {
        "snapshot": 463,
        "tree.idx": 221716,
        "peri.z": 0.26,
        "peri.t": cosmology.convert_time("time", "redshift", 0.26),
    },
    "m12i_res7100": {
        "snapshot": 356,
        "tree.idx": 422845,
        "peri.z": 0.6,
        "peri.t": cosmology.convert_time("time", "redshift", 0.6),
    },
    "m12m_res7100": {
        "snapshot": 290,
        "tree.idx": 1821315,
        "peri.z": 0.92,
        "peri.t": cosmology.convert_time("time", "redshift", 0.92),
    },
    "m12w_res7100": {
        "snapshot": 359,
        "tree.idx": 403617,
        "peri.z": 0.59,
        "peri.t": cosmology.convert_time("time", "redshift", 0.59),
    },
}

af_md = asdf.AsdfFile(pericenter_data)

af_md.write_to("peri_md.asdf")


### Get LMC pericenter characteristics (stellar mass, infall mass)

simulation_list = [
    "../../../data/latte_metaldiff/m12f_res7100",
    "../../../data/latte_metaldiff/m12i_res7100",
    "../../../data/latte_metaldiff/m12m_res7100",
    "../../../data/latte_metaldiff/m12w_res7100",
]

char_af = asdf.AsdfFile({})

for sim in simulation_list:
    sim_name = path.split(sim)[-1]
    hal = halo.io.IO.read_catalogs("snapshot", af_md[sim_name]["snapshot"], sim)

    lmc_idx = np.where(hal["tree.index"] == af_md[sim_name]["tree.idx"])[0][0]

    tmp_dict = {
        "star.mass": np.format_float_scientific(hal["star.mass"][lmc_idx]),
        "infall.mass": np.format_float_scientific(hal["infall.mass"][lmc_idx]),
        "mass.peak": np.format_float_scientific(hal["mass.peak"][lmc_idx]),
        "mass": np.format_float_scientific(hal["mass"][lmc_idx]),
        "mass.bound": np.format_float_scientific(hal["mass.bound"][lmc_idx]),
        "mass.vir": np.format_float_scientific(hal["mass.vir"][lmc_idx]),
        "mass.200c": np.format_float_scientific(hal["mass.200c"][lmc_idx]),
        "mass.200m": np.format_float_scientific(hal["mass.200m"][lmc_idx]),
    }

    char_af[sim_name] = tmp_dict

char_af.write_to("peri_lmc_details.asdf")
