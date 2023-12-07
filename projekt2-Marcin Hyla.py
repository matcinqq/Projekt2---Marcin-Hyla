import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as sp

data = np.array([[1.0, 3.0], [2.0, 1.0], [3.5, 4.0], [5.0, 0.0], [6.0, 0.5], [9.0, -2.0], [9.5, -3.0]])
x = np.linspace(0, 10, 100)
N = data.shape[0]
hi = np.zeros([N - 1])
bi = np.zeros([N - 1])
u = np.zeros([N - 1])
v = np.zeros([N - 1])
z = np.zeros([N])

for i in range(len(hi)):  #hi
    hi[i] = data[i + 1, 0] - data[i, 0]

for i in range(len(bi)):  #bi
    bi[i] = 6 / hi[i] * (data[i + 1][1] - data[i][1])

for i in range(len(u)): #u
    if i == 0:
        u[i] = 2 * (hi[0] + hi[1])
    else:
        u[i] = 2 * (hi[i - 1] + hi[i]) - ((hi[i - 1]) ** 2) / (u[i - 1])

for i in range(len(v)): #v
    if i == 0:
        v[i] = bi[1] - bi[0]
    else:
        v[i] = bi[i] - bi[i - 1] - (hi[i - 1] * (v[i - 1] / u[i - 1]))

for i in range(len(z) - 2, -1, -1): #z
    z[i] = (v[i] - hi[i] * z[i + 1]) / u[i]


#Funkcje obliczające współczynniki Ai, Bi i Ci
def Ai(a, b):
    A = np.zeros([len(a)])
    for i in range(len(a)):
        A[i] = (1/(6*a[i]))*(b[i+1]-b[i])
    return A


def Bi(b):
    B = np.zeros([len(b)])
    for i in range(len(b)):
        B[i] = (b[i])/2
    return B


def Ci(a, b, c):
    C = np.zeros([len(a)])
    for i in range(len(a)):
        C[i] = (-a[i]/6)*(b[i+1]+2*b[i])+1/a[i]*(c[i+1][1]-c[i][1])
    return C


#Funkcja licząca Si
def Si(hi, z, data, x):
    S = np.zeros(x.shape)
    for i in range(len(x)):
        for j in range(len(hi)):
            if data[j, 0] <= x[i] <= data[j + 1, 0]:
                xi = x[i]
                ti = data[j,0]
                yi = data[j, 1]
                S[i] = yi + (xi - ti) * (Ci(hi, z, data)[j] + (xi - ti) * (Bi(z)[j] + (xi - ti) * Ai(hi, z)[j]))

    return S


#Dla komfortu, zamiana na odzdzielną zmienną
S = Si(hi,z,data,x)
print(S)


#Wielomiany Lagrange
def l(x, data):
    poly = np.ones((data.shape[0], x.shape[0]))
    for i in range(data.shape[0]):
        for j in range(data.shape[0]):
            if i != j:
                poly[i, :] *= (x - data[j, 0]) / (data[i, 0] - data[j, 0])
    return poly


#Funkcją znajdująca wielomian interpolujący
def poly_interpolation(x, data):
    p = np.zeros(x.shape[0])
    basis = l(x, data)
    for n in range(data.shape[0]):
        p += basis[n, :] * data[n, 1]
    return p

#Potrzebne do sklejenia 1. stopnia
def Si_1(x, data):
    S = np.zeros(x.shape)
    for i in range(data.shape[0] - 1):
        for m in range(len(x)):
            try:
                if data[i, 0] <= x[m] <= data[i + 1, 0]:
                    S[m] = ((data[i + 1, 1] - data[i, 1]) / (data[i + 1, 0] - data[i, 0])) * (x[m] - data[i, 0]) + data[i, 1]
            except:
                continue
    return S

for i in range(len(x) - 1):
    xdd = np.zeros(x.shape)
    try:
        xdd[i] += Si(x, data)[i] * (x[i] <= data[i + 1, 0]) * (x[i] >= data[i, 0])
    except:
        xdd[i] = 0
        continue

#zad 3
x_data = data[:, 0]
y_data = data[:, 1]
cs = sp.CubicSpline(x_data, y_data)
wartosci_sklejenia = cs(x)

#zad.1
fig = plt.figure()
axes = fig.add_subplot(1, 1, 1)
axes.plot(data[:, 0], data[:, 1], 'ro', label="węzły")
plt.plot(x,S, color = "c", linewidth = 1.5, label = "sklejenie st.3")
plt.legend()
plt.show()

#Zad.2
fig = plt.figure()
axes = fig.add_subplot(1, 1, 1)
axes.plot(x, poly_interpolation(x, data), 'b', label="Lagrange")
axes.plot(data[:, 0], data[:, 1], 'ro', label="węzły")
plt.plot(x,S, color = "c", linewidth = 1.5, label = "sklejenie st.3")
plt.plot(x,Si_1(x,data), color = "yellow", label = "sklejenie st.1")
plt.legend()
plt.show()

#Zad.3
fig = plt.figure()
axes = fig.add_subplot(1, 1, 1)
plt.plot(x,S, color = "c", linewidth = 1.5, label = "sklejenie st.3")
plt.plot(x_data, y_data, "ro", label='wiezy')
plt.plot(x, wartosci_sklejenia, color = "blue", label='CubicSpline')
plt.legend()
plt.show()

'''
Funkcja narysowana przez CubicSpline różni się trochę od wyznaczonej samodzielnie
'''
