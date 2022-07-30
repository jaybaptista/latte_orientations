import orientations as ori
import os.path as path
import asdf

simulation_list = [
    "../../../data/latte_metaldiff/m12f_res7100",
    "../../../data/latte_metaldiff/m12i_res7100",
    "../../../data/latte_metaldiff/m12m_res7100",
    "../../../data/latte_metaldiff/m12w_res7100",
]

peri_md = asdf.open("peri_md.asdf")

for sim in simulation_list:

    sim_name = path.split(sim)[-1]

    tidx = peri_md[sim_name]["tree.idx"]

    ori.trackHaloPosition(sim, tidx, write=True)
