import ntpath
import os
from tkinter import Tk
from tkinter.filedialog import askopenfilename

import numpy as np
from scipy import interpolate

np.seterr("raise")


def calcular_IJ(XC, YC, XB, YB, phi, S):
    num_pan = len(XC)
    I = np.zeros([num_pan, num_pan])
    J = np.zeros([num_pan, num_pan])

    for i in range(num_pan):
        for j in range(num_pan):
            if j != i:
                A = -(XC[i] - XB[j]) * np.cos(phi[j]) - (
                    YC[i] - YB[j]
                ) * np.sin(phi[j])
                B = (XC[i] - XB[j]) ** 2 + (YC[i] - YB[j]) ** 2
                Cn = np.sin(phi[i] - phi[j])
                Dn = -(XC[i] - XB[j]) * np.sin(phi[i]) + (
                    YC[i] - YB[j]
                ) * np.cos(phi[i])
                Ct = -np.cos(phi[i] - phi[j])
                Dt = (XC[i] - XB[j]) * np.cos(phi[i]) + (
                    YC[i] - YB[j]
                ) * np.sin(phi[i])
                E = np.sqrt(B - A**2)
                if E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E):
                    I[i, j] = 0
                    J[i, j] = 0
                else:
                    term1 = (
                        0.5 * Cn * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B)
                    )
                    term2 = ((Dn - A * Cn) / E) * (
                        np.arctan2((S[j] + A), E) - np.arctan2(A, E)
                    )
                    I[i, j] = term1 + term2
                    term1 = (
                        0.5 * Ct * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B)
                    )
                    term2 = ((Dt - A * Ct) / E) * (
                        np.arctan2((S[j] + A), E) - np.arctan2(A, E)
                    )
                    J[i, j] = term1 + term2
            if np.iscomplex(I[i, j]) or np.isnan(I[i, j]) or np.isinf(I[i, j]):
                I[i, j] = 0
            if np.iscomplex(J[i, j]) or np.isnan(J[i, j]) or np.isinf(J[i, j]):
                J[i, j] = 0
    return I, J


def calcular_streamline(XP, YP, XB, YB, phi, S):
    num_pan = len(XB) - 1
    Mx = np.zeros(num_pan)
    My = np.zeros(num_pan)

    for j in range(num_pan):
        A = -(XP - XB[j]) * np.cos(phi[j]) - (YP - YB[j]) * np.sin(phi[j])
        B = (XP - XB[j]) ** 2 + (YP - YB[j]) ** 2
        Cx = -np.cos(phi[j])
        Dx = XP - XB[j]
        Cy = -np.sin(phi[j])
        Dy = YP - YB[j]
        E = np.sqrt(B - A**2)
        if E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E):
            Mx[j] = 0
            My[j] = 0
        else:
            term1 = 0.5 * Cx * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B)
            term2 = ((Dx - A * Cx) / E) * (
                np.arctan2((S[j] + A), E) - np.arctan2(A, E)
            )
            Mx[j] = term1 + term2
            term1 = 0.5 * Cy * np.log((S[j] ** 2 + 2 * A * S[j] + B) / B)
            term2 = ((Dy - A * Cy) / E) * (
                np.arctan2((S[j] + A), E) - np.arctan2(A, E)
            )
            My[j] = term1 + term2

        if np.iscomplex(Mx[j]) or np.isnan(Mx[j]) or np.isinf(Mx[j]):
            Mx[j] = 0
        if np.iscomplex(My[j]) or np.isnan(My[j]) or np.isinf(My[j]):
            My[j] = 0

    return Mx, My


def xfoil(NACA, PPAR, AoA, flagAirfoil):
    xFoilResults = list(range(9))

    if flagAirfoil[0] == 1:
        airfoilName = NACA
        xFoilResults[0] = airfoilName
    elif flagAirfoil[1] == 1:
        root = Tk()
        ftypes = [("dat file", "*.dat")]
        ttl = "Select Airfoil File"
        dir1 = "/Airfoil_DAT_Selig/"
        root.withdraw()
        root.update()
        root.fileName = askopenfilename(
            filetypes=ftypes, initialdir=dir1, title=ttl
        )
        root.destroy()

        head, tail = ntpath.split(root.fileName)
        airfoilName = tail[0 : len(tail) - 4]
        xFoilResults[0] = airfoilName

    saveFlnm = "Save_" + airfoilName + ".txt"
    saveFlnmCp = "Save_" + airfoilName + "_Cp.txt"
    saveFlnmPol = "Save_" + airfoilName + "_Pol.txt"

    if os.path.exists(saveFlnm):
        os.remove(saveFlnm)
    if os.path.exists(saveFlnmCp):
        os.remove(saveFlnmCp)
    if os.path.exists(saveFlnmPol):
        os.remove(saveFlnmPol)

    fid = open("xfoil_input.inp", "w")
    if flagAirfoil[0] == 1:
        fid.write("NACA " + NACA + "\n")
    elif flagAirfoil[1] == 1:
        fid.write("LOAD " + "./Airfoil_DAT_Selig/" + tail + "\n")

    fid.write("PPAR\n")
    fid.write("N " + PPAR[0] + "\n")
    fid.write("P " + PPAR[1] + "\n")
    fid.write("T " + PPAR[2] + "\n")
    fid.write("R " + PPAR[3] + "\n")
    fid.write("XT " + PPAR[4] + "\n")
    fid.write("XB " + PPAR[5] + "\n")
    fid.write("\n")
    fid.write("\n")

    fid.write("PSAV " + saveFlnm + "\n")

    fid.write("OPER\n")
    fid.write("Pacc 1 \n")
    fid.write("\n\n")
    fid.write("Alfa " + str(AoA) + "\n")
    fid.write("CPWR " + saveFlnmCp + "\n")
    fid.write("PWRT\n")
    fid.write(saveFlnmPol + "\n")
    if os.path.exists(saveFlnmPol):
        fid.write("y \n")

    fid.close()

    os.system("xfoil.exe < xfoil_input.inp")

    if os.path.exists("xfoil_input.inp"):
        os.remove("xfoil_input.inp")

    dataBufferCp = np.loadtxt(saveFlnmCp, skiprows=3)
    xFoilResults[1] = dataBufferCp[:, 0]
    xFoilResults[2] = dataBufferCp[:, 1]
    xFoilResults[3] = dataBufferCp[:, 2]

    if os.path.exists(saveFlnmCp):
        os.remove(saveFlnmCp)

    dataBuffer = np.loadtxt(saveFlnm, skiprows=0)
    xFoilResults[4] = dataBuffer[:, 0]
    xFoilResults[5] = dataBuffer[:, 1]

    if os.path.exists(saveFlnm):
        os.remove(saveFlnm)

    dataBufferPol = np.loadtxt(saveFlnmPol, skiprows=12)
    xFoilResults[6] = dataBufferPol[1]
    xFoilResults[7] = dataBufferPol[2]
    xFoilResults[8] = dataBufferPol[4]

    if os.path.exists(saveFlnmPol):
        os.remove(saveFlnmPol)

    return xFoilResults


def calcular_circ(a, b, x0, y0, numT, Vx, Vy, X, Y):
    t = np.linspace(0, 2 * np.pi, numT)
    xC = a * np.cos(t) + x0
    yC = b * np.sin(t) + y0
    fx = interpolate.RectBivariateSpline(Y, X, Vx)
    fy = interpolate.RectBivariateSpline(Y, X, Vy)
    VxC = fx.ev(yC, xC)
    VyC = fy.ev(yC, xC)
    Gamma = -(np.trapz(VxC, xC) + np.trapz(VyC, yC))
    return Gamma, xC, yC, VxC, VyC
