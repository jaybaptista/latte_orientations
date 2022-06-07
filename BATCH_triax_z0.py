import gizmo_analysis as gizmo
import utilities as ut
import numpy as np
import matplotlib.pyplot as plt
import importlib
from matplotlib import colors
from astropy.table import Table, Column
from astropy.io import ascii
import pandas as pd

sim_dir = "../data/latte_metaldiff/"
sims = ["m12f_res7100", "m12i_res7100", "m12w_res7100", "m12m_res7100"]
rs = np.arange(10, 410, step=10)


def axis_ratio_calculation(tolerance, iters_max, r_max):
    x = pos_dark[:, 0]
    y = pos_dark[:, 1]
    z = pos_dark[:, 2]
    # Initialize s = c/a, p = c/b, and q = b/a
    s = 1
    p = 1
    q = 1
    # Select all particles with total distances less that 2*r_max
    n_sel = d_tot < 2 * r_max
    x = x[n_sel]
    y = y[n_sel]
    z = z[n_sel]
    d = (x**2 + y**2 / q**2 + z**2 / s**2) ** (1 / 2)
    # Initialize difference values
    ds = 2 * tolerance
    dp = 2 * tolerance
    dq = 2 * tolerance
    # Loop-counter
    loop_counter = 0
    # While-loop to check if
    while not ((ds < tolerance) and (dp < tolerance) and (dq < tolerance)):
        # If-statement for loop-counter
        loop_counter = loop_counter + 1
        if loop_counter == iters_max + 1:
            print("- Iterations maxed out.")
            loop_counter = loop_counter - 1
            break
        # Define new as old values
        s_old = s
        p_old = p
        q_old = q
        d_old = d
        # Select particles less than r_max
        n_select = d_old < r_max
        x_select = x[n_select]
        y_select = y[n_select]
        z_select = z[n_select]
        m_select = 1
        # Calculate new distance measure for this iteration
        d_select = (
            x_select**2 + y_select**2 / q_old**2 + z_select**2 / s_old**2
        ) ** (1 / 2)
        # Calculate parts and total reduced inertia tensor
        Ixx_terms = m_select * x_select * x_select / d_select**2
        Ixy_terms = m_select * x_select * y_select / d_select**2
        Ixz_terms = m_select * x_select * z_select / d_select**2
        Iyy_terms = m_select * y_select * y_select / d_select**2
        Iyz_terms = m_select * y_select * z_select / d_select**2
        Izz_terms = m_select * z_select * z_select / d_select**2
        Ixx = Ixx_terms.sum()
        Ixy = Ixy_terms.sum()
        Ixz = Ixz_terms.sum()
        Iyy = Iyy_terms.sum()
        Iyz = Iyz_terms.sum()
        Izz = Izz_terms.sum()
        I_m = [[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]]
        #         Eigen values and eigen vectors of reduced inertia tensor matrix
        eigen_values, eigen_vectors = np.linalg.eigh(I_m)

        #         order = np.argsort(eigen_values)
        #         eigen_values = eigen_values[order]

        # Sort eigen vectors to find new x, y, z (where x<->a, y<->b, z<->c)
        # teigen_vectors = np.transpose(eigen_vectors)
        # sorted_vectors = teigen_vectors[indices]
        XX = eigen_vectors[:, 2]
        YY = eigen_vectors[:, 1]
        ZZ = eigen_vectors[:, 0]
        sorted_vectors = [XX, YY, ZZ]
        #         sorted_vectors = np.transpose(sorted_vectors)[order]
        sorted_vectors = np.transpose(sorted_vectors)
        # If-statement to skip non-orthogonal eigen vectors
        dot1 = np.around(np.dot(XX, YY), decimals=5)
        dot2 = np.around(np.dot(YY, ZZ), decimals=5)
        dot3 = np.around(np.dot(ZZ, XX), decimals=5)

        # FORCE RHR COMPLIANCE
        #         sorted_vectors[2] = np.cross(sorted_vectors[0], sorted_vectors[1])

        if (dot1 == 0) and (dot2 == 0) and (dot3 == 0):
            X = ut.coordinate.get_coordinates_rotated(
                np.transpose([x, y, z]), sorted_vectors
            )
            nonorthogonal_counter = 0
        else:
            nonorthogonal_counter = 1
            break
        # Sort eigen values to find new a, b, c (where a>=b>=c)
        # indices = np.unravel_index(np.argsort(eigen_values), eigen_values.shape)
        # sorted_values = eigen_values[indices]
        sorted_values = eigen_values ** (1 / 2)
        a = sorted_values[2]
        b = sorted_values[1]
        c = sorted_values[0]
        # X = ut.coordinate.get_coordinates_rotated(np.transpose([x,y,z]), sorted_vectors)
        x_new = X[:, 0]
        y_new = X[:, 1]
        z_new = X[:, 2]
        # New axis-ratios
        s = c / a
        p = c / b
        q = b / a
        d = (x_new**2 + y_new**2 / q**2 + z_new**2 / s**2) ** (1 / 2)
        # Differences in axis-ratios
        ds = np.abs(s - s_old)
        dp = np.abs(p - p_old)
        dq = np.abs(q - q_old)
    # Calculate the triaxiality
    T = (a**2 - b**2) / (a**2 - c**2)
    # Finally set axis ratios, triaxiality, loop_counter of comp_type, non-orthogonal eigen_counter
    axis_ratios_comp = np.array(
        [s, p, q, a, b, c, T, loop_counter, nonorthogonal_counter]
    )
    print(axis_ratios_comp)
    # Print for confirmation of Pooling completion
    print(
        "FINISHED axis ratio calculations for: "
        + sims_name
        + " @ r = {} kpc".format(r_max)
    )
    return sorted_vectors, axis_ratios_comp, T


for sims_name in sims:

    df = {"t": [], "vec": [], "meta": []}

    part = gizmo.io.Read.read_snapshots(
        ["dark"], "redshift", 0, sim_dir + sims_name, assign_hosts_rotation=True
    )
    pos_dark = part["dark"].prop("host.distance")
    vel_dark = part["dark"].prop("host.velocity.principal")
    d_tot = part["dark"].prop("host.distance.total")

    for r_i in rs:
        print("Calculating triaxiality at: ", r_i)
        vec, ratio, t_i = axis_ratio_calculation(0.001, 1, r_i)

        df["t"].append(t_i)
        df["vec"].append(vec)
        df["meta"].append(ratio)

    df = pd.DataFrame(df)
    df.to_hdf("_data_triax_z0_{}_i1.h5".format(sims_name[:4]), "w")
