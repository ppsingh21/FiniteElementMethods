import numpy as np
import matplotlib.pyplot as plt

def shape_functions(h):
    """Compute coefficients of shape functions."""
    co=[[0.5, -0.75, 0, 0.25],[-0.125*h, 0.125*h, 0.125*h, -0.125*h],[0.5, 0.75, 0, -0.25],[0.125*h, 0.125*h, -0.125*h, -0.125*h]]
    cod=[[-0.75, 0, 0.75], [0.125*h, 0.25*h, -0.375*h],[0.75, 0, -0.75],[0.125*h, -0.25*h, -0.375*h]]
    codd=[[0, 1.5], [0.25*h, -0.75*h], [0, -1.5], [-0.25*h, -0.75*h]]
    coddd=[[1.5],[-0.75*h],[-1.5],[-0.75*h]]
    return co, cod, codd, coddd

def integration_approximation(re):
    """Approximate integration coefficients."""
    coi=np.array([[0.23862,0.46791],[-0.23862,0.46791],[0.66121,0.36076],[-0.66121,0.36076],[0.93247,0.17132],[-0.93247,0.17132]])
    s=0
    for i in range(0, 6):
        su = 0
        for j in range(0, int(re.shape[0])):
            su += re[j] * coi[i][0] ** j
        s += su * coi[i][1]
    return s

def elemental_matrix(f, co, h):
    """Compute elemental matrix."""
    fe = np.zeros((4, 1))
    for i in range(0, 4):
        fe[i][0] = (h / 2) * integration_approximation(np.polynomial.polynomial.polymul(co[i], f[0]))
    return fe

def get_value(co, x, p, d):
    """Map element position to global matrix position."""
    s = 0
    for i in range(0, p - d +1):
        s += co[i] * x ** i
    return s

def main():
    n = int(input("Enter number of elements: "))
    h = 1 / n
    co, cod, codd, coddd = shape_functions(h)
    bc1 = int(input("Type of BC at x=0 (1 or 2): "))
    bc2 = int(input("Type of BC at x=1 (1 or 2): "))
    if bc1 == 2:
        force1 = float(input("Enter force: "))
        moment1 = float(input("Enter moment: "))
    else:
        dis1 = float(input("Enter displacement: "))
        slope1 = float(input("Enter slope: "))
    if bc2 == 2:
        force2 = float(input("Enter force: "))
        moment2 = float(input("Enter moment: "))
    else:
        dis2 = float(input("Enter displacement: "))
        slope2 = float(input("Enter slope: "))
    fp = int(input("Enter order of t: "))
    fco = np.zeros((1, fp+1))
    for i in range(0, fp + 1):
        fco[0][i] = int(input("Enter coeff: "))
    if n % 2 == 0:
        forcemid = int(input("Enter force at midpoint: "))
        momentmid = int(input("Enter moment at midpoint: "))
    nodloc = np.zeros((n, 2))
    for i in range(0, n):
        for j in range(0, 2):
            if i == j == 0:
                continue
            elif j % 2 == 0:
                nodloc[i][j] = nodloc[i - 1][j + 1]
            else:
                nodloc[i][j] = nodloc[i][j - 1] + h
    K = np.zeros((2 * n + 2, 2 * n + 2))
    F = np.zeros((2 * n + 2, 1))
    Q = np.zeros((2 * n + 2, 1))
    ke = (1000 / (3 * h ** 3)) * np.array([[6, -3 * h, -6, -3 * h],
                                           [-3 * h, 2 * h ** 2, 3 * h, h ** 2],
                                           [-6, 3 * h, 6, 3 * h],
                                           [-3 * h, h ** 2, 3 * h, 2 * h ** 2]])
    re = np.zeros((1, 5))
    for i in range(0, n):
        fe = np.zeros((4, 1))
        sunod = (nodloc[i][0] + nodloc[i][1]) / 2
        fcof = np.array([[fco[0][0] + fco[0][1] * sunod + fco[0][2] * sunod * 2, fco[0][1] * (h / 2) + fco[0][2] * h * sunod,
                          fco[0][2] * (h / 4)]])
        fe = elemental_matrix(fcof, co, h)
        if i == 0:
            K[:4, :4] += ke
            F[:4, 0:1] += fe
        else:
            K[2 * i:2 * (i + 2), 2 * i:2 * (i + 2)] += ke
            F[2 * i:2 * (i + 2), 0:1] += fe
    KF = K
    if bc1 == 1:
        for i in range(1, 2 * n + 2):
            F[i] = F[i] - dis1 * KF[i][0]
        for i in range(0, 2 * n + 2):
            for j in range(0, n * 2 + 2):
                if i == 0 or j == 0:
                    KF[i][j] = 0
        KF[0][0] = 1
        F[0][0] = dis1
        Q[0][0] = 0
        for i in range(2, 2 * n + 2):
            F[i] = F[i] - slope1 * KF[i][1]
        for i in range(0, 2 * n + 2):
            for j in range(0, 2 * n + 2):
                if i == 1 or j == 1:
                    KF[i][j] = 0
        KF[1][1] = 1
        F[1][0] = slope1
        Q[1][0] = 0
    if bc2 == 1:
        for i in range(0, 2 * n + 1):
            F[i] = F[i] - slope2 * KF[i][2 * n + 1]
        for i in range(0, n * 2 + 2):
            for j in range(0, 2 * n + 2):
                if i == 2 * n + 1 or j == n * 2 + 1:
                    KF[i][j] = 0
        KF[n * 2 + 1][n * 2 + 1] = 1
        F[n * 2 + 1][0] = slope2
        Q[n * 2 + 1][0] = 0
        for i in range(0, 2 * n):
            F[i] = F[i] - dis2 * KF[i][2 * n]
        for i in range(0, n * 2 + 2):
            for j in range(0, 2 * n + 2):
                if i == 2 * n or j == n * 2:
                    KF[i][j] = 0
        KF[n * 2][n * 2] = 1
        F[n * 2][0] = dis2
        Q[n * 2][0] = 0
    if bc1 == 2:
        Q[0][0] = force1
        Q[1][0] = -moment1
    if bc2 == 2:
        Q[2 * n][0] = force2
        Q[2 * n + 1][0] = -moment2
    if n % 2 == 0:
        Q[n][0] = forcemid
        Q[n + 1][0] = -momentmid
    U = np.linalg.inv(KF) @ (F + Q)

    # Plotting
    x = np.linspace(0, 1, 1000)
    yh, yhslope, moment, shear, stress = [], [], [], [], []

    for i in range(n):
        for j in range(1000):
            if nodloc[i][0] <= x[j] <= nodloc[i][1]:  # Check if x[j] is within the element range
                aux = (2 * x[j] - (nodloc[i][0] + nodloc[i][1])) / h
                w, wslope, mom, sh, st, b = 0, 0, 0, 0, 0, 0
                for m in range(i * 2, i * 2 + 4):
                    w=w+(U[m][0]*get_value(co[b],aux,3,0))
                    wslope=wslope+(U[m][0]*get_value(cod[b],aux,3,1))*(2/h)
                    mom=mom+(500/3)*(U[m][0]*get_value(codd[b],aux,3,2))*(2/h)**2
                    sh=sh+(-500/3)*(U[m][0]*get_value(coddd[b],aux,3,3))*(2/h)**3
                    st=st+(-1000)*(U[m][0]*get_value(codd[b],aux,3,2))*(2/h)**2
                    b=b+1
                yh.append(w)
                yhslope.append(wslope)
                moment.append(mom)
                shear.append(sh)
                stress.append(st)

    # Ensure yh and other lists have the same length as x
    if len(yh) < 1000:
        yh.extend([yh[-1]] * (1000 - len(yh)))
        yhslope.extend([yhslope[-1]] * (1000 - len(yhslope)))
        moment.extend([moment[-1]] * (1000 - len(moment)))
        shear.extend([shear[-1]] * (1000 - len(shear)))
        stress.extend([stress[-1]] * (1000 - len(stress)))
    elif len(yh) > 1000:
        yh = yh[:1000]
        yhslope = yhslope[:1000]
        moment = moment[:1000]
        shear = shear[:1000]
        stress = stress[:1000]
        
    ye, yed = [], []
    i = 0
    while i <= 1:
        ye.append((3 / 500) * ((i ** 4) / 24 - (i ** 3) / 3 + 5 * i * i / 4))
        yed.append(((i ** 3) / 6 - i * i + 5 * i / 2) * (3 / 500))
        i += 0.001
    er = 0
    for i in range(0, 1000):
        er += ((ye[i] - yh[i]) ** 2)
    e = np.sqrt(er / 1000) * 100
    print("RMSE error%:", e, "%")

    # Plotting results
    line1, = plt.plot(x, yh, label="FEM")
    line2, = plt.plot(x, ye, label="Exact")
    plt.legend(handles=[line1, line2])
    plt.title("Plot of deflection with x")
    plt.ylabel("Deflection")
    plt.xlabel("x")
    plt.show()

    line1, = plt.plot(x, yhslope, label="FEM")
    line2, = plt.plot(x, yed, label="Exact")
    plt.legend(handles=[line1, line2])
    plt.title("Plot of slope with x")
    plt.ylabel("Slope")
    plt.xlabel("x")
    plt.show()

    plt.plot(x, moment)
    plt.title("Plot of moment with x")
    plt.ylabel("Moment")
    plt.xlabel("x")
    plt.show()

    plt.plot(x, shear)
    plt.title("Plot of shear force with x")
    plt.ylabel("Shear Force")
    plt.xlabel("x")
    plt.show()

    plt.plot(x, stress)
    plt.title("Plot of stress with x")
    plt.ylabel("Stress")
    plt.xlabel("x")
    plt.show()

if __name__ == "__main__":
    main()

