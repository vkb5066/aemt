import re
from copy import deepcopy

#Holds info on individual atoms--------------------------------------------------------------------
class Atom:
    atomType = None
    idNum = None
    a, b, c = None, None, None
    flagA, flagB, flagC = None, None, None
    equivPositions = None
    bondedTo = None
    misc = None

    #I can't have overloaded init functions in python so this is what we're stuck with :(
    def __init__(self):
        self.atomType = None
        self.idNum = None
        self.a, self.b, self.c = None, None, None
        self.flagA, self.flagB, self.flagC = None, None, None
        self.equivPositions = []
        self.bondedTo = None
        self.misc = None

    def SetAtomType(self, atomType):
        self.atomType = atomType
    def SetAtomId(self, idNum):
        self.idNum = idNum
    def SetAtomPositions(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c
    def SetAtomPositionsByList(self, lis):
        self.a = lis[0]
        self.b = lis[1]
        self.c = lis[2]
    def SetAtomFlags(self, a, b, c):
        self.flagA = a
        self.flagB = b
        self.flagC = c
    def SetAtomFlagsByList(self, lis):
        self.flagA = lis[0]
        self.flagB = lis[1]
        self.flagC = lis[2]

    def __eq__(self, atom):
        return (Dist(self, atom) < 0.01)
    def __ne__(self, atom):
        return (Dist(self, atom) >= 0.01)


    def GetDeepCopy(self):
        return deepcopy(self)

    def Print(self, string = False):
        if(not string):
            print(self.idNum, self.atomType, self.a, self.b, self.c, self.flagA, self.flagB,
                  self.flagC)
        else:
            return("%.10f"%(self.a) + "   " + "%.10f"%(self.b) + "   " + "%.10f"%(self.c) + "   " +
                   str(self.flagA) + "   " + str(self.flagB) + "   " + str(self.flagC))

def Dist(atomA, atomB):
    return ((atomB.a - atomA.a)**2 + (atomB.b - atomA.b)**2 + (atomB.c - atomA.c)**2)**(1/2)


#Holds the info on an entire POSCAR (or CONTCAR) file---------------------------------------------
class Poscar:
    #Misc. Stuff
    filePath = None

    #Things explicitly written in the file
    comment = None
    univScaleFactor = None
    superCellVecA = None
    superCellVecB = None
    superCellVecC = None
    atomTypes = None
    atomTypeNums = None
    atomTypesAndNums = None
    selectiveDynamicsTag = None
    directTag, cartesianTag = None, None
    atoms =[]
    velocities = None

    #Possibly useful calculations from file info
    volume = None
    elemRanges = None

    #Initialization Functions
    def FetchComment(self):
        self.comment = ReadLines(self.filePath, 0, isNum = False)[0]
    def FetchUnivScaleFactor(self):
        self.univScaleFactor = float(ReadLines(self.filePath, 1)[0])
    def FetchSuperCellVecs(self):
        self.superCellVecA = [float(a) for a in ReadLines(self.filePath, 2)[0].split()]
        self.superCellVecB = [float(b) for b in ReadLines(self.filePath, 3)[0].split()]
        self.superCellVecC = [float(c) for c in ReadLines(self.filePath, 4)[0].split()]
    def FetchAtomTypes(self):
        self.atomTypes = ReadLines(self.filePath, 5, isNum = False, stripSpecial = False)[0].split()
    def FetchAtomTypeNums(self):
        self.atomTypeNums = [int(a) for a in ReadLines(self.filePath, 6)[0].split()]
    def FetchAtomTypesAndNums(self):
        self.atomTypesAndNums = dict(zip(self.atomTypes, self.atomTypeNums))
    def FetchSelTag(self):
        string = ReadLines(self.filePath, 7, isNum = False)[0]
        if(string[0] == 's' or string[0] == 'S'):
            self.selectiveDynamicsTag = True
        else:
            self.selectiveDynamicsTag = False
    def FetchDirectCartesianFlag(self):
        lineNum = 7
        if(self.selectiveDynamicsTag):
            lineNum = lineNum + 1
        string = ReadLines(self.filePath, lineNum, isNum = False)[0]
        if(string[0] == 'd' or string[0] == 'D'):
            self.cartesianTag = False
            self.directTag = True
        else:
            self.directTag = False
            self.cartesianTag = True
    def FetchAtoms(self):
        lineNum = 8
        if(self.selectiveDynamicsTag):
            lineNum = lineNum + 1
        i, idN = 0, 1
        while(i < len(self.atomTypes)):
            j = self.atomTypeNums[i]
            while(j > 0): #this is a LOT of file accesses.  TODO:  make a better way to do this
                line = ReadLines(self.filePath, lineNum, isNum = False,
                                 stripSpecial = False)[0].split()
                thisAtom = Atom()
                thisAtom.SetAtomId(int(idN))
                thisAtom.SetAtomType(str(self.atomTypes[i]))
                thisAtom.SetAtomPositionsByList([float(v) for v in line[0:3]])
                try:
                    thisAtom.SetAtomFlagsByList([str(f) for f in line[3:6]])
                except(IndexError):
                    thisAtom.SetAtomFlags('', '', '')
                self.atoms.append(thisAtom)
                j = j - 1
                lineNum = lineNum + 1
                idN = idN + 1
            i = i + 1

    def FetchElementRanges(self):
        tot = sum(self.atomTypeNums)
        self.elemRanges = []
        self.elemRanges.append([self.atomTypes[0], [1, self.atomTypeNums[0]]])
        tot -= self.atomTypeNums[0]
        i = 1
        while(tot > 0):
            self.elemRanges.append([self.atomTypes[i], [self.elemRanges[i-1][1][1] + 1,
                                                       self.elemRanges[i-1][1][1] +
                                                       self.atomTypeNums[i]]])
            tot = tot - self.atomTypeNums[i]
            i = i + 1

    def __init__(self, path):
        #Make absolutly positivly 100% sure that new instances of this class do not have old items
        #in it.  This is very annoying to have to worry about
        self.filePath = None
        self.univScaleFactor = None
        self.superCellVecA = None
        self.superCellVecB = None
        self.superCellVecC = None
        self.atomTypes = None
        self.atomTypeNums = None
        self.atomTypesAndNums = None
        self.selectiveDynamicsTag = None
        self.directTag, self.cartesianTag = None, None
        self.atoms = []
        self.velocities = []
        self.volume = None
        self.elemRanges = None

        #Now I can re-set (or initially set) the items
        self.filePath = path
        self.FetchComment()
        self.FetchUnivScaleFactor()
        self.FetchSuperCellVecs()
        self.FetchAtomTypes()
        self.FetchAtomTypeNums()
        self.FetchAtomTypesAndNums()
        self.FetchElementRanges()
        self.FetchSelTag()
        self.FetchDirectCartesianFlag()
        self.FetchAtoms()

    #Other Functions
    def GetDeepCopy(self):
        return deepcopy(self)

    def RecalculateAtomTypesAndNums(self):
        self.atomTypesAndNums = {}
        #Initial pass
        for atom in self.atoms:
            self.atomTypesAndNums[str(atom.atomType)] = 0
        #Counting pass
        for atom in self.atoms:
            self.atomTypesAndNums[str(atom.atomType)] += 1

        self.atomTypes, self.atomTypeNums = [], []
        for k, v in zip(self.atomTypesAndNums.keys(), self.atomTypesAndNums.values()):
            self.atomTypes.append(k)
            self.atomTypeNums.append(v)

    def RecalculateAtomOrder(self, resetIDs = True):
        newAtoms = []
        count = 1
        for elem in self.atomTypes:
            for i in range(0, len(self.atoms)):
                if (self.atoms[i].atomType == elem):
                    if(resetIDs):
                        self.atoms[i].idNum = count
                    newAtoms.append(self.atoms[i].GetDeepCopy())
                    count += 1
        self.atoms = newAtoms

    def Refresh(self, resetIDs = True):
        self.RecalculateAtomTypesAndNums()
        self.RecalculateAtomOrder(resetIDs)
        self.FetchElementRanges()

    def ChangeAtomOrder(self, order, resetIDs = True):
        if(order == self.atomTypes):
            return
        if(set(self.atomTypes[:]) != set(order[:])):
            print("WARNING:  Original elements do not match new elements.  Original elements are: ")
            for elem in self.atomTypes:
                print(elem)
            print("And new elements are: ")
            for elem in order:
                print(elem)
            print("If you know that the new atom list reflects this new element order, you may"
                  + " ignore this warning.")
        self.atomTypes = order
        self.RecalculateAtomOrder(resetIDs)
        self.RecalculateAtomTypesAndNums()
        self.FetchElementRanges()

    def GetAtomEquivPositions(self):
        self.ConvertToDirect()
        for atom in self.atoms:
            for addA in [-1, 0, 1]:
                for addB in [-1, 0, 1]:
                    for addC in [-1, 0, 1]:
                        new = Atom()
                        new.SetAtomType(atom.atomType)
                        new.SetAtomId(atom.idNum)
                        pos = [atom.a + addA, atom.b + addB, atom.c + addC]
                        new.SetAtomPositionsByList(pos)
                        new.SetAtomFlags(atom.flagA, atom.flagB, atom.flagC)
                        atom.equivPositions.append(new)


    def ConvertToDirect(self):
        if(self.directTag):
            return

        #Apply univ scale factor to supercell vectors
        for i in range(0, 3):
            self.superCellVecA[i] = self.univScaleFactor * self.superCellVecA[i]
            self.superCellVecB[i] = self.univScaleFactor * self.superCellVecB[i]
            self.superCellVecC[i] = self.univScaleFactor * self.superCellVecC[i]
        #Apply univ scale factor to atom coords
        for atom in self.atoms:
            atom.a = self.univScaleFactor * atom.a
            atom.b = self.univScaleFactor * atom.b
            atom.c = self.univScaleFactor * atom.c
        self.univScaleFactor = 1.0

        #Find the determinant of the supercell matrix
        def det2x2(tl, tr, bl, br):
            return tl*br - tr*bl
        det3x3 =   self.superCellVecA[0]*det2x2(self.superCellVecB[1], self.superCellVecB[2],
                                                self.superCellVecC[1], self.superCellVecC[2]) \
                 - self.superCellVecA[1]*det2x2(self.superCellVecB[0], self.superCellVecB[2],
                                                self.superCellVecC[0], self.superCellVecC[2]) \
                 + self.superCellVecA[2]*det2x2(self.superCellVecB[0], self.superCellVecB[1],
                                                self.superCellVecC[0], self.superCellVecC[1])
        if(det3x3 == 0.0):
            print("ConvertToDirect(): Found a zero determinant.")
            return

        #Get the inverse matrix by dividing the co-factor matrix by the determinent
        inv = [0*i for i in range(0, 9)]
        inv[0] = det2x2(self.superCellVecB[1], self.superCellVecB[2], self.superCellVecC[1],
                        self.superCellVecC[2])/det3x3
        inv[1] = -det2x2(self.superCellVecB[0], self.superCellVecB[2], self.superCellVecC[0],
                         self.superCellVecC[2])/det3x3
        inv[2] = det2x2(self.superCellVecB[0], self.superCellVecB[1], self.superCellVecC[0],
                        self.superCellVecC[1])/det3x3
        inv[3] = -det2x2(self.superCellVecA[1], self.superCellVecA[2], self.superCellVecC[1],
                        self.superCellVecC[2])/det3x3
        inv[4] = det2x2(self.superCellVecA[0], self.superCellVecA[2], self.superCellVecC[0],
                         self.superCellVecC[2])/det3x3
        inv[5] = -det2x2(self.superCellVecA[0], self.superCellVecA[1], self.superCellVecC[0],
                        self.superCellVecC[1])/det3x3
        inv[6] = det2x2(self.superCellVecA[1], self.superCellVecA[2], self.superCellVecB[1],
                        self.superCellVecB[2])/det3x3
        inv[7] = -det2x2(self.superCellVecA[0], self.superCellVecA[2], self.superCellVecB[0],
                         self.superCellVecB[2])/det3x3
        inv[8] = det2x2(self.superCellVecA[0], self.superCellVecA[1], self.superCellVecB[0],
                        self.superCellVecB[1])/det3x3

        #Update the atomic coordinates
        for atom in self.atoms:
            aOrig, bOrig, cOrig = atom.a, atom.b, atom.c
            atom.a = inv[0]*aOrig + inv[1]*bOrig + inv[2]*cOrig
            atom.b = inv[3]*aOrig + inv[4]*bOrig + inv[5]*cOrig
            atom.c = inv[6]*aOrig + inv[7]*bOrig + inv[8]*cOrig
            for equiv in atom.equivPositions:
                aOrig, bOrig, cOrig = equiv.a, equiv.b, equiv.c
                equiv.a = inv[0]*aOrig + inv[1]*bOrig + inv[2]*cOrig
                equiv.b = inv[3]*aOrig + inv[4]*bOrig + inv[5]*cOrig
                equiv.c = inv[6]*aOrig + inv[7]*bOrig + inv[8]*cOrig


        self.cartesianTag = False
        self.directTag = True

    def ConvertToCartesian(self):
        if(self.cartesianTag):
            return

        #Apply univ scale factor to atom coords and do the conversion all at once
        for atom in self.atoms:
            aOrig, bOrig, cOrig = atom.a, atom.b, atom.c
            atom.a = self.univScaleFactor * ((self.superCellVecA[0]*aOrig) +
                                             (self.superCellVecB[0]*bOrig) +
                                             (self.superCellVecC[0]*cOrig))
            atom.b = self.univScaleFactor * ((self.superCellVecA[1]*aOrig) +
                                             (self.superCellVecB[1]*bOrig) +
                                             (self.superCellVecC[1]*cOrig))
            atom.c = self.univScaleFactor * ((self.superCellVecA[2]*aOrig) +
                                             (self.superCellVecB[2]*bOrig) +
                                             (self.superCellVecC[2]*cOrig))
            for equiv in atom.equivPositions:
                aOrig, bOrig, cOrig = equiv.a, equiv.b, equiv.c
                equiv.a = self.univScaleFactor * ((self.superCellVecA[0]*aOrig) +
                                                 (self.superCellVecB[0]*bOrig) +
                                                 (self.superCellVecC[0]*cOrig))
                equiv.b = self.univScaleFactor * ((self.superCellVecA[1]*aOrig) +
                                                 (self.superCellVecB[1]*bOrig) +
                                                 (self.superCellVecC[1]*cOrig))
                equiv.c = self.univScaleFactor * ((self.superCellVecA[2]*aOrig) +
                                                 (self.superCellVecB[2]*bOrig) +
                                                 (self.superCellVecC[2]*cOrig))
        #Apply univ scale factor to supercell vectors
        for i in range(0, 3):
            self.superCellVecA[i] = self.univScaleFactor * self.superCellVecA[i]
            self.superCellVecB[i] = self.univScaleFactor * self.superCellVecB[i]
            self.superCellVecC[i] = self.univScaleFactor * self.superCellVecC[i]

        self.univScaleFactor = 1.0

        self.directTag = False
        self.cartesianTag = True

    #Moves all atoms to inside unit cell all direct coords between 0 and 1
    def MoveAtomsToUnitCell(self):
        toCart = False
        if(self.cartesianTag == True):
            toCart = True

        self.ConvertToDirect()
        for i in range(0, len(self.atoms)):
            while(self.atoms[i].a < 0.):
                self.atoms[i].a = self.atoms[i].a + 1.0
            while(self.atoms[i].a > 1.):
                self.atoms[i].a = self.atoms[i].a - 1.0
            while(self.atoms[i].b < 0.):
                self.atoms[i].b = self.atoms[i].b + 1.0
            while(self.atoms[i].b > 1.):
                self.atoms[i].b = self.atoms[i].b - 1.0
            while(self.atoms[i].c < 0.):
                self.atoms[i].c = self.atoms[i].c + 1.0
            while(self.atoms[i].c > 1.):
                self.atoms[i].c = self.atoms[i].c - 1.0

        if(toCart):
            self.ConvertToCartesian()

    def Write(self, path):
        outfile = open(path, 'w')
        outfile.write(str(self.comment) + "\n")
        outfile.write("%.5f"%(self.univScaleFactor) + "\n")
        outfile.write("%.10f"%(self.superCellVecA[0]) + "   " + "%.10f"%(self.superCellVecA[1]) +
                      "   " + "%.10f"%(self.superCellVecA[2]) + "\n")
        outfile.write("%.10f"%(self.superCellVecB[0]) + "   " + "%.10f"%(self.superCellVecB[1]) +
                      "   " + "%.10f"%(self.superCellVecB[2]) + "\n")
        outfile.write("%.10f"%(self.superCellVecC[0]) + "   " + "%.10f"%(self.superCellVecC[1]) +
                      "   " + "%.10f"%(self.superCellVecC[2]) + "\n")
        for a in self.atomTypes:
            outfile.write(str(a) + "   ")
        outfile.write("\n")
        for a in self.atomTypeNums:
            outfile.write(str(a) + "   ")
        outfile.write("\n")
        if(self.selectiveDynamicsTag):
            outfile.write("Sel\n")
        if(self.directTag):
            outfile.write("Direct\n")
        if(self.cartesianTag):
            outfile.write("Cartesian\n")
        j = 0
        while(j < len(self.atomTypes)):
            i = 0
            while(i < len(self.atoms)):
                if(self.atoms[i].atomType == self.atomTypes[j]):
                    outfile.write(self.atoms[i].Print(string = True) + "\n")
                i = i + 1
            j = j + 1
        outfile.close()


#NON_CLASS FUNCTIONS--------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def ReadLines(fileLoc, nums, isNum = True, stripSpecial = True):
    if(not isinstance(nums, list)):
        nums = [nums]

    ret = []
    infile = open(fileLoc, 'r')
    for lineNum, line in enumerate(infile):
        for n in nums:
            if(lineNum == n):
                ret.append(line)
    infile.close()

    if(isNum):
        ret = [re.sub('[^0-9,. ]','', r) for r in ret]
    elif(stripSpecial):
        ret = [re.sub('[^a-zA-Z ]', '', r) for r in ret]
    return ret
