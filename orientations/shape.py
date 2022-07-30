import numpy as np
import utilities as ut


def getShape(positions, rmax, distances=None, tolerance=0.001, iters_max=1000):
    x = positions[:, 0]
    y = positions[:, 1]
    z = positions[:, 2]

    if distances is None:
        distances = np.sqrt(np.sum(positions**2, axis=1))

    s = p = q = 1

    m_sel = distances < 2 * rmax

    x_sel = x[m_sel]
    y_sel = y[m_sel]
    z_sel = z[m_sel]

    d = np.sqrt(x_sel**2 + (y_sel / q) ** 2 + (z_sel / s) ** 2)

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
        x_sel = x_sel[m_sel]
        y_sel = y_sel[m_sel]
        z_sel = z_sel[m_sel]

        d_sel = np.sqrt(x_sel**2 + (y_sel / q_old) ** 2 + (z_sel / s_old) ** 2)

        IXX = np.sum((x_sel**2 / d_sel**2))
        IYY = np.sum((y_sel**2 / d_sel**2))
        IZZ = np.sum((z_sel**2 / d_sel**2))

        IXY = np.sum(-x_sel * y_sel / d_sel**2)
        IXZ = np.sum(-x_sel * z_sel / d_sel**2)
        IYZ = np.sum(-z_sel * y_sel / d_sel**2)

        I_M = np.array([[IXX, IXY, IXZ], [IXY, IYY, IYZ], [IXZ, IYZ, IZZ]])

        evals, evecs = np.linalg.eigh(I_M)
        sort_idx = np.argsort(evecs)[::-1]

        ####### red tape begins here #######
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
        ####### red tape ends here ########

        X = ut.coordinate.get_coordinates_rotated(
            np.transpose([x_sel, y_sel, z_sel]), sorted_vectors
        )
        x_sel = X[:, 0]
        y_sel = X[:, 1]
        z_sel = X[:, 2]

        s = c / a
        p = c / b
        q = b / a
        d = (x_sel**2 + y_sel**2 / q**2 + z_sel**2 / s**2) ** (1 / 2)

        ds = np.abs(s - s_old)
        dp = np.abs(p - p_old)
        dq = np.abs(q - q_old)

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
