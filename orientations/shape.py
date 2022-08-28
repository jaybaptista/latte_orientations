import numpy as np
import utilities as ut


def getShape(positions, rmax, mass, distances=None, tolerance=0.001, iters_max=1000, rfin=200):
    x = positions[:, 0]
    y = positions[:, 1]
    z = positions[:, 2]

    if distances is None:
        distances = np.sqrt(np.sum(positions**2, axis=1))

    s = p = q = 1

    m_select = distances < 2 * rfin

    x = x[m_select]
    y = y[m_select]
    z = z[m_select]
    mass = mass[m_select]

    d = np.sqrt((x**2) + ((y/q) ** 2) + ((z/s) ** 2))

    ds = dp = dq = 2 * tolerance

    nonorthogonal_counter = 2

    i = 0

    while not ((ds < tolerance) and (dp < tolerance) and (dq < tolerance)):
        i += 1

        if (i % 10) == 0:
            print(f"Shape iteration = {i}")

        if i == iters_max + 1:
            print("Reached maximum iteration")
            i += -1
            break

        s_old = s
        p_old = p
        q_old = q
        d_old = d

        m_sel = d_old < rmax
        
        if m_sel.sum() == 0:
            print("Iteration region too small, zero particles selected. Ending iteration...")
            break
        
        x_sel = x[m_sel]
        y_sel = y[m_sel]
        z_sel = z[m_sel]
        mass_sel = mass[m_sel]

        d_sel = ((x_sel**2) + ((y_sel/q_old)**2) + ((z_sel / s_old)**2))**(1/2)

        IXX = (mass_sel*(x_sel**2 / d_sel**2)).sum()
        IYY = (mass_sel*(y_sel**2 / d_sel**2)).sum()
        IZZ = (mass_sel*(z_sel**2 / d_sel**2)).sum()

        IXY = (mass_sel*x_sel * y_sel / d_sel**2).sum()
        IXZ = (mass_sel*x_sel * z_sel / d_sel**2).sum()
        IYZ = (mass_sel*z_sel * y_sel / d_sel**2).sum()

        I_M = np.array([[IXX, IXY, IXZ], [IXY, IYY, IYZ], [IXZ, IYZ, IZZ]]) / mass_sel.sum()

        evals, evecs = np.linalg.eigh(I_M)

        XX = evecs[:, 2]
        YY = evecs[:, 1]
        ZZ = evecs[:, 0]
        sorted_vectors = np.transpose([XX, YY, ZZ])

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

        sorted_values = evals ** (1 / 2)
        a = sorted_values[2]
        b = sorted_values[1]
        c = sorted_values[0]

        x_new = X[:, 0]
        y_new = X[:, 1]
        z_new = X[:, 2]

        s = c / a
        p = c / b
        q = b / a
        d = ((x_new**2) + (y_new**2 / q**2) + (z_new**2 / s**2)) ** (1 / 2)

        ds = abs(s - s_old)
        dp = abs(p - p_old)
        dq = abs(q - q_old)

    T = (a**2 - b**2) / (a**2 - c**2)

    df = {
        "rmax": rmax,
        "shape": T,
        "iter": i,
        "tensor": sorted_vectors,
        "ortho": nonorthogonal_counter,
        "a": a,
        "b": b,
        "c": c,
        "s": s,
        "p": p,
        "q": q,
    }

    return df

# reorient everything into new frame
def getShapeV2(positions, rmax, mass, distances=None, tolerance=0.001, iters_max=1000, rfin=200):
    x = positions[:, 0]
    y = positions[:, 1]
    z = positions[:, 2]

    if distances is None:
        distances = np.sqrt(np.sum(positions**2, axis=1))

    s = p = q = 1

    m_select = distances < 2 * rfin

    x = x[m_select]
    y = y[m_select]
    z = z[m_select]
    mass = mass[m_select]

    d = np.sqrt((x**2) + ((y/q) ** 2) + ((z/s) ** 2))

    ds = dp = dq = 2 * tolerance

    nonorthogonal_counter = 2

    i = 0

    while not ((ds < tolerance) and (dp < tolerance) and (dq < tolerance)):
        i += 1

        if (i % 10) == 0:
            print(f"Shape iteration = {i}")

        if i == iters_max + 1:
            print("Reached maximum iteration")
            i += -1
            break

        s_old = s
        p_old = p
        q_old = q
        d_old = d

        m_sel = d_old < rmax
        
        if m_sel.sum() == 0:
            print("Iteration region too small, zero particles selected. Ending iteration...")
            break
        
        x_sel = x[m_sel]
        y_sel = y[m_sel]
        z_sel = z[m_sel]
        mass_sel = mass[m_sel]

        d_sel = ((x_sel**2) + ((y_sel/q_old)**2) + ((z_sel / s_old)**2))**(1/2)

        IXX = (mass_sel*(x_sel**2 / d_sel**2)).sum()
        IYY = (mass_sel*(y_sel**2 / d_sel**2)).sum()
        IZZ = (mass_sel*(z_sel**2 / d_sel**2)).sum()

        IXY = (mass_sel*x_sel * y_sel / d_sel**2).sum()
        IXZ = (mass_sel*x_sel * z_sel / d_sel**2).sum()
        IYZ = (mass_sel*z_sel * y_sel / d_sel**2).sum()

        I_M = np.array([[IXX, IXY, IXZ], [IXY, IYY, IYZ], [IXZ, IYZ, IZZ]]) / mass_sel.sum()

        evals, evecs = np.linalg.eigh(I_M)

        XX = evecs[:, 2]
        YY = evecs[:, 1]
        ZZ = evecs[:, 0]
        sorted_vectors = np.transpose([XX, YY, ZZ])

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

        sorted_values = evals ** (1 / 2)
        a = sorted_values[2]
        b = sorted_values[1]
        c = sorted_values[0]

        x_new = X[:, 0]
        y_new = X[:, 1]
        z_new = X[:, 2]

        s = c / a
        p = c / b
        q = b / a
        d = ((x_new**2) + (y_new**2 / q**2) + (z_new**2 / s**2)) ** (1 / 2)

        ds = abs(s - s_old)
        dp = abs(p - p_old)
        dq = abs(q - q_old)
        
        # UPDATE COORDINATE SYSTEM TO NEW ROTATED FRAME
        x = x_new
        y = y_new
        z = z_new

    T = (a**2 - b**2) / (a**2 - c**2)

    df = {
        "rmax": rmax,
        "shape": T,
        "iter": i,
        "tensor": sorted_vectors,
        "ortho": nonorthogonal_counter,
        "a": a,
        "b": b,
        "c": c,
        "s": s,
        "p": p,
        "q": q,
    }

    return df
