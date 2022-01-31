INFILE_LOC = "fullCBM"
MOD = 150
CUT = False

kx, ky, kz, e = [], [], [], []

#Read
with open(INFILE_LOC, 'r') as infile:
    for n, lin in enumerate(infile):
        if(n%MOD != 0):
            continue

        line = lin.split()

        kx.append(float(line[0]))
        ky.append(float(line[1]))
        kz.append(float(line[2]))
        e.append(float(line[3]))

    infile.close()

#Remove some unwanted stuff for better visualization
if(CUT):
    kx_, ky_, kz_, e_ = [], [], [], []
    for i in range(0, len(e)):
        if(kz[i] > -0.001):
            continue
        kx_.append(kx[i])
        ky_.append(ky[i])
        kz_.append(kz[i])
        e_.append(e[i])
    kx, ky, kz, e = kx_, ky_, kz_, e_

#Plot
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

plt.title(INFILE_LOC)
img = ax.scatter(kx, ky, kz, c=e, cmap=plt.jet(), s=3.5)
fig.colorbar(img)

ax.set_xlabel(r"$\vec{k_\mathrm{x}}$")
ax.set_ylabel(r"$\vec{k_\mathrm{y}}$")
ax.set_zlabel(r"$\vec{k_\mathrm{z}}$")

plt.show()
