#Header for writing the data needed for computing the effective mass tensor
from copy import deepcopy as Deepcopy
import headerEigenval as he
from scipy.spatial.kdtree import KDTree
pi = 3.141592654

def KpsDirToCart(B, kpLis):
    for i in range(0, len(kpLis)):
        tmp = [kpLis[i].a, kpLis[i].b, kpLis[i].c]
        tmp = FracToCart(frac=tmp, B=B)
        kpLis[i].a = tmp[0]
        kpLis[i].b = tmp[1]
        kpLis[i].c = tmp[2]
        kpLis[i].abc = [tmp[0], tmp[1], tmp[2]]
    return kpLis

#From a set of original kps, returns a set of new kps spanning the entire BZ + a bit more in cartesian
#(recip) coords.
def TransformKps(origKps, isym, B, cubeLen):
    #Deal with symmetry: get k points into the full BZ
    additionalList = None
    if(isym == -1): ##no symmetry at all
        additionalList = [-1, 0, +1]
    elif(isym == 0): ##only inversion symmetry
        additionalList = [-1, 0, +1]
        newKps_ = []
        for k in origKps:
            newK = Deepcopy(k)
            newK.a = -newK.a
            newK.b = -newK.b
            newK.c = -newK.c
            newKps_.append(newK)
        origKps += newKps_
    else:
        print("I don't know how to support this yet! Use ISYM = 0/-1")
        exit(1)

    #Deal with periodic images: get images into list 'additionalKs'
    additionalKs = []
    for kp in origKps:
        for i in additionalList:
            for j in additionalList:
                for k in additionalList:
                    tstK = Deepcopy(kp)
                    tstK.a += i
                    tstK.b += j
                    tstK.c += k

                    #We only want to include new k points that are within a resonably sized cube surrounding
                    #the first BZ
                    tstKCart = Deepcopy(tstK)
                    tmp = [tstKCart.a, tstKCart.b, tstKCart.c]
                    tmp = FracToCart(frac=tmp, B=B)
                    tstKCart.a = tmp[0]
                    tstKCart.b = tmp[1]
                    tstKCart.c = tmp[2]

                    if (abs(tstKCart.a) < cubeLen and abs(tstKCart.b) < cubeLen and abs(tstKCart.c) < cubeLen):
                        additionalKs.append(tstKCart)

    #Get the original k points into cartesian coordinates to compare vs the additional list
    origKps = KpsDirToCart(B=B, kpLis=origKps)

    #Now, dispose of any duplicate k-points
    # !!! NOTE !!! This might be unnecessary with RBFs.  We could save a lot of time by being smart about
    #this instead
    #I do this with KD trees so that it dosn't take 20 years to finish for dense meshes
    potentialCopiedKLis = origKps + additionalKs
    lis = [[k.a, k.b, k.c] for k in potentialCopiedKLis]
    node = KDTree(lis, copy_data=True)

    ret = []
    indsToIgnore = []
    for i in range(0, len(lis)):
        if(i in indsToIgnore):
            continue
        inds = node.query_ball_point(x=lis[i], r=1E-6, return_sorted=True)
        ret.append(potentialCopiedKLis[inds[0]])
        if(len(inds) > 1):
            for j in range(1, len(inds)):
                indsToIgnore.append(inds[j])

    return ret





#Returns a fractional vector frac to cartesian coordinates, given B = [[b1x, b1y, b1z], [b2x, b2y, b2z],
#[b3x, b3y, b3z]]
def FracToCart(frac, B):
    cart = [None, None, None]
    cart[0] = B[0][0]*frac[0] + B[1][0]*frac[1] + B[2][0]*frac[2]
    cart[1] = B[0][1]*frac[0] + B[1][1]*frac[1] + B[2][1]*frac[2]
    cart[2] = B[0][2]*frac[0] + B[1][2]*frac[1] + B[2][2]*frac[2]
    return cart

def dot(u, v):
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
def cross(u, v):
    return [u[1]*v[2] - u[2]*v[1], -(u[0]*v[2] - u[2]*v[0]), u[0]*v[1] - u[1]*v[0]]
def RealToRecip(A):
    a1, a2, a3 = A[0], A[1], A[2]
    V = dot(a1, cross(a2, a3))

    return  [[2*pi/V * b for b in cross(a2, a3)],
             [2*pi/V * b for b in cross(a3, a1)],
             [2*pi/V * b for b in cross(a1, a2)]]
