import headerPoscar as hp
import headerEigenval as he
import headerKTransform as hk

PATH_TO_POSCAR = "CONTCAR"
PATH_TO_EIGENVAL = "EIGENVAL"
PATH_TO_OUTFILE = "aemEig"

#Get reciprocal lattice
poscar = hp.Poscar(PATH_TO_POSCAR)
A = [poscar.superCellVecA,
     poscar.superCellVecB,
     poscar.superCellVecC]
B = hk.RealToRecip(A)

#Transform k from direct to cartesian
origKps = he.ParseEigenval("EIGENVAL", occupanciesListed=True)[:400]
cartKps = hk.TransformKps(origKps=origKps, B=B, isym=0, cubeLen=0.471009)

#Write new formatting
with open(PATH_TO_OUTFILE, 'w') as outfile:
     outfile.write(str(len(cartKps)) + ' ' + str(len(cartKps[0].bands)) + "\n")
     for i in range(0, len(cartKps)):
          outfile.write(str(cartKps[i].a) + ' ' + str(cartKps[i].b) + ' ' + str(cartKps[i].c) + "\n")
          for j in range(0, len(cartKps[i].bands)):
               outfile.write(str(cartKps[i].bands[j][1]) + "\n")
     outfile.write("STOP")
     outfile.close()
