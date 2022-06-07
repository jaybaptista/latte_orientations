def axis_ratio_calculation(comp, tolerance, iters_max, r_max):

    # If-statement to determine x, y, z, m, and d_tot
    if comp == "star":
        x = x_star
        y = y_star
        z = z_star
        m = m_star
        d_tot = d_tot_star
    elif comp == "gass":
        x = x_gass
        y = y_gass
        z = z_gass
        m = m_gass
        d_tot = d_tot_gass
    elif comp == "dark":
        x = x_dark
        y = y_dark
        z = z_dark
        m = m_dark
        d_tot = d_tot_dark
    elif comp == "totl":
        x = np.hstack((x_star, x_gass, x_dark))
        y = np.hstack((y_star, y_gass, y_dark))
        z = np.hstack((z_star, z_gass, z_dark))
        m = np.hstack((m_star, m_gass, m_dark))
        d_tot = np.hstack((d_tot_star, d_tot_gass, d_tot_dark))

    # Initialize s = c/a, p = c/b, and q = b/a
    s = 1
    p = 1
    q = 1

    # Select all particles with total distances less that 2*r_max
    n_sel = d_tot < 2 * r_max
    x = x[n_sel]
    y = y[n_sel]
    z = z[n_sel]
    m = m[n_sel]
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
        m_select = m[n_select]

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

        # Eigen values and eigen vectors of reduced inertia tensor matrix
        eigen_values, eigen_vectors = np.linalg.eigh(I_m)

        # Sort eigen vectors to find new x, y, z (where x<->a, y<->b, z<->c)
        # teigen_vectors = np.transpose(eigen_vectors)
        # sorted_vectors = teigen_vectors[indices]
        XX = eigen_vectors[:, 2]
        YY = eigen_vectors[:, 1]
        ZZ = eigen_vectors[:, 0]
        sorted_vectors = [XX, YY, ZZ]
        sorted_vectors = np.transpose(sorted_vectors)

        # If-statement to skip non-orthogonal eigen vectors
        dot1 = np.around(np.dot(XX, YY), decimals=5)
        dot2 = np.around(np.dot(YY, ZZ), decimals=5)
        dot3 = np.around(np.dot(ZZ, XX), decimals=5)
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
        [
            s,
            p,
            q,
            a,
            b,
            c,
            T,
            loop_counter,
            nonorthogonal_counter,
            sorted_vectors[0, 0],
            sorted_vectors[0, 1],
            sorted_vectors[0, 2],
            sorted_vectors[1, 0],
            sorted_vectors[1, 1],
            sorted_vectors[1, 2],
            sorted_vectors[2, 0],
            sorted_vectors[2, 1],
            sorted_vectors[2, 2],
        ]
    )
    print(axis_ratios_comp)

    # Print for confirmation of Pooling completion
    print("FINISHED axis ratio calculations for: " + sims_name + " - " + comp)

    return axis_ratios_comp
