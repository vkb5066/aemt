
INFILE_LOC = "vz"

ks, nrgs, nrgsS = [], [], []

with open(INFILE_LOC, 'r') as infile:
    for lin in infile:
        line = lin.split()

        ks.append(float(line[0]))
        nrgs.append(float(line[1]))
        nrgsS.append(float(line[2]))

    infile.close()

ks, nrgs, nrgsS = (list(t) for t in zip(*sorted(zip(ks, nrgs, nrgsS))))

from matplotlib import pyplot as plt
plt.title(INFILE_LOC)
#plt.plot(ks, nrgs, 'k.')
plt.plot(ks, nrgsS, 'k.')

#plt.ylim(0, 10)
#plt.axvline(0)
print(len(ks))

plt.show()
