import numpy as np

TOL = 0.000001
VBM_MIN_OCC = 0.99 #if VBM_MIN_OCC <= occupancy <= 1.0, occupied by an electron
CBM_MAX_OCC = 0.01 #if 0.0 <= occupancy <= CBM_MAX_OCC, occupied by a hole


def Approx(a, b):
    return abs(a - b) < TOL

def AlmostEqual(u, v):
    return ((Approx(u[0], v[0])) and (Approx(u[1], v[1])) and (Approx(u[2], v[2])))

class kPoint:
    idChar = None
    a, b, c = None, None, None
    weight = None
    bands = None
    isMag = None

    def __init__(self, header, bandInfo, isMag=False):
        self.idChar = None #set this later - no way to obtain from EIGENVAL file alone
        self.a, self.b, self.c = None, None, None
        self.weight = None
        self.bands = []
        self.isMag = None

        self.a, self.b, self.c = float(header.split()[0]), float(header.split()[1]), \
                                 float(header.split()[2])
        self.weight = float(header.split()[3])

        if(not isMag):
            for i in range(0, len(bandInfo)):
                self.bands.append([int(bandInfo[i].split()[0]),     #Band ID
                                   float(bandInfo[i].split()[1]),   #Band Energy
                                   float(bandInfo[i].split()[2])])  #Occupancy (not int bc of part. occ)
        else:
            for i in range(0, len(bandInfo)):
                self.bands.append([int(bandInfo[i].split()[0]),     #Band ID
                                   float(bandInfo[i].split()[1]),   #Band Energy up
                                   float(bandInfo[i].split()[2]),   #Band Energy dwn
                                   float(bandInfo[i].split()[3]),   #Occupancy up
                                   float(bandInfo[i].split()[4])])  #Occupancy dwn

        self.isMag = isMag

    def __eq__(self, other):
        return ((self.a - other.a)*(self.a - other.a) + \
                (self.b - other.b)*(self.b - other.b) + \
                (self.c - other.c)*(self.c - other.c) < 1E-4)


    def SetId(self, char):
        self.idChar = char

    def RemoveLowOcc(self, rmBelow=0.5):
        newBands = []
        for b in self.bands:
            if(b[2] > rmBelow):
                newBands.append(b)
        self.bands = newBands[:]

    #As of now, only shifts the c.b. up by specified amount
    def ApplyScissor(self, scissor = 0, tol = 0.5):
        for b in self.bands:
            if(b[2] < tol): #if occupancy less than tolerance
                b[1] = b[1] + scissor #increase energy by scissor in eV


#Returns a list of k-point instances from the EIGENVAL file
def ParseEigenval(infileLoc, occupanciesListed=True, eFermi=None, isMag=False, autoDetermineFormat=True):
    if((not occupanciesListed) and (eFermi == None)):
        print("ParseEigenval: ERR: if the occupancies are not explicitly written to EIGENVAL, you need to")
        print("provide the fermi energy to this function!")
        return [None]

    #Length of lines in eigenval  Initial non-mag defaults:
    occLisBandLineLen = 3
    nonOccLisBandLineLen = 2
    if(isMag): ##new length if magnetized
        occLisBandLineLen = 5
        nonOccLisBandLineLen = 3

    kpts = []
    with open(infileLoc, 'r') as infile:
        readingBands = False
        head, info = None, []
        for num, line in enumerate(infile):
            #Remove EOL characters from line
            line.replace("\n", '')
            line.replace("\r", '')

            #Obtain important params before reading in the k-points
            if(num == 5):
                ##nElectrons = int(line.split()[0])
                nKpts = int(line.split()[1])
                nBands = int(line.split()[2])
                bandCounter = int(line.split()[2])

            #Check for reading k-point heads
            if(num > 5 and not readingBands and len(line.split()) == 4):
                head = line
                readingBands = True

            #Check for reading in that k-point's band info
            if(occupanciesListed):
                if(num > 5 and readingBands and len(line.split()) == occLisBandLineLen):
                    info.append(line)
                    bandCounter = bandCounter - 1
            else:
                if(num > 5 and readingBands and len(line.split()) ==  nonOccLisBandLineLen):
                    if(isMag):
                        ##spin-up occ
                        if (float(line.split()[1]) <= eFermi):  ##if e < eFermi, assume fully occupied
                            line += " 1.00000"
                        else: ##otherwise no occ
                            line += " 0.00000"
                        ##spin-dn occ
                        if (float(line.split()[2]) <= eFermi):  ##if e < eFermi, assume fully occupied
                            line += " 1.00000"
                        else: ##otherwise no occ
                            line += " 0.00000"
                    else:
                        ##general occ
                        if (float(line.split()[1]) <= eFermi):
                            line += " 1.00000"
                        else:
                            line += " 0.00000"

                    info.append(line)
                    bandCounter = bandCounter - 1

            #Check for ending that k-point's info
            if(num > 5 and readingBands and (len(line.split()) == 0 or bandCounter == 0)):
                kpts.append(kPoint(head, info, isMag))
                bandCounter = nBands
                head, info = None, []
                readingBands = False

        infile.close()

        if(len(kpts) == nKpts):
            return kpts
        else:
            print("I couldn't read in all the k-points.  SAD!!!")
        return []

#Feed me a specific k-point.
#if spin=None, returns the highest of spin up and down
#if spin=up, returns only the highest spin up.  Same for spin=down
def GetHighestOccupiedBand(thisKpt, spin=None):
    if(thisKpt.isMag):
        #Spin polarized
        if(spin==None):
            for i in range(0, len(thisKpt.bands)):
                if ((VBM_MIN_OCC <= thisKpt.bands[i][3] or VBM_MIN_OCC <= thisKpt.bands[i][4]) and \
                (thisKpt.bands[i + 1][3] < VBM_MIN_OCC or thisKpt.bands[i + 1][4] < VBM_MIN_OCC)):
                    return thisKpt.bands[i]
        #Spin polarized, only want spin up valence band
        elif(spin[0] == 'u' or spin[0] == 'U'):
            for i in range(0, len(thisKpt.bands)):
                if ((VBM_MIN_OCC <= thisKpt.bands[i][3]) and (thisKpts.bands[i + 1][3] < VBM_MIN_OCC)):
                    return thisKpt.bands[i]
        #Spin polarized, only want spin down valence band
        elif (spin[0] == 'd' or spin[0] == 'D'):
            for i in range(0, len(thisKpt.bands)):
                if ((VBM_MIN_OCC <= thisKpt.bands[i][4]) and (thisKpts.bands[i + 1][4] < VBM_MIN_OCC)):
                    return thisKpt.bands[i]
    #non spin-pol
    else:
        for i in range(0, len(thisKpt.bands)):
            if((VBM_MIN_OCC <= thisKpt.bands[i][2]) and (thisKpt.bands[i+1][2] < VBM_MIN_OCC)):
                return thisKpt.bands[i]

    return None ##uh oh

#Feed me a specific k-point.
#if spin=None, returns the highest of spin up and down
#if spin=up, returns only the highest spin up.  Same for spin=down
def GetLowestUnoccupiedBand(thisKpt, spin=None):
    if (thisKpt.isMag):
        # Spin polarized
        if (spin == None):
            for i in range(0, len(thisKpt.bands)):
                if ((thisKpt.bands[i][3] <= CBM_MAX_OCC or thisKpt.bands[i][4] <= CBM_MAX_OCC) and \
                        (CBM_MAX_OCC < thisKpt.bands[i-1][3] or CBM_MAX_OCC < thisKpt.bands[i-1][4])):
                    return thisKpt.bands[i]
        # Spin polarized, only want spin up valence band
        elif (spin[0] == 'u' or spin[0] == 'U'):
            for i in range(0, len(thisKpt.bands)):
                if ((thisKpt.bands[i][3] <= CBM_MAX_OCC) and (CBM_MAX_OCC < thisKpt.bands[i - 1][3])):
                    return thisKpt.bands[i]
        # Spin polarized, only want spin down valence band
        elif (spin[0] == 'd' or spin[0] == 'D'):
            for i in range(0, len(thisKpt.bands)):
                if ((thisKpt.bands[i][4] <= CBM_MAX_OCC) and (CBM_MAX_OCC < thisKpt.bands[i - 1][4])):
                    return thisKpt.bands[i]
    # non spin-pol
    else:
        for i in range(0, len(thisKpt.bands)):
            if ((thisKpt.bands[i][2] <= CBM_MAX_OCC) and (CBM_MAX_OCC < thisKpt.bands[i - 1][2])):
                return thisKpt.bands[i]

    return None  ##uh oh


#Get the k-point closest to the top of the valence band
def GetVBMKpoint(kpLis, spin=None):
    bestKp, bandMaxEnergy = kpLis[0], kpLis[0].bands[0][1] ##initialize maxEnergy to the bottom of the vb

    if(bestKp.isMag):
        #spin polarized
        if(spin == None):
            for kp in kpLis:
                hob = GetHighestOccupiedBand(kp)
                if (hob[1] > bandMaxEnergy or hob[2] > bandMaxEnergy):
                    bestKp = kp
                    bandMaxEnergy = max(hob[1], hob[2])
                return bestKp
        #spin polarized, only consider spin up
        elif(spin[0] == 'u' or spin[0] == 'U'):
            for kp in kpLis:
                hob = GetHighestOccupiedBand(kp)
                if (hob[1] > bandMaxEnergy):
                    bestKp = kp
                    bandMaxEnergy = hob[1]
                return bestKp
        #spin polarized, only consider spin down
        elif(spin[0] == 'd' or spin[0] == 'D'):
            for kp in kpLis:
                hob = GetHighestOccupiedBand(kp)
                if (hob[2] > bandMaxEnergy):
                    bestKp = kp
                    bandMaxEnergy = hob[2]
                return bestKp
    #non spin-pol
    else:
        for kp in kpLis:
            if(GetHighestOccupiedBand(kp)[1] > bandMaxEnergy):
                bestKp = kp
                bandMaxEnergy = GetHighestOccupiedBand(kp)[1]
        return bestKp

    return None ##super uh oh

#Get the k-point closest to the bottom of the conduction band
def GetCBMKpoint(kpLis, spin=None):
    bestKp, bandMinEnergy = kpLis[-1], kpLis[-1].bands[-1][1] ##initialize minEnergy to the top of the cb

    if(bestKp.isMag):
        #spin polarized
        if(spin == None):
            for kp in kpLis:
                hob = GetLowestUnoccupiedBand(kp)
                if (hob[1] < bandMinEnergy or hob[2] < bandMinEnergy):
                    bestKp = kp
                    bandMinEnergy = min(hob[1], hob[2])
                return bestKp
        #spin polarized, only consider spin up
        elif(spin[0] == 'u' or spin[0] == 'U'):
            for kp in kpLis:
                hob = GetLowestUnoccupiedBand(kp)
                if (hob[1] < bandMinEnergy):
                    bestKp = kp
                    bandMinEnergy = hob[1]
                return bestKp
        #spin polarized, only consider spin down
        elif(spin[0] == 'd' or spin[0] == 'D'):
            for kp in kpLis:
                hob = GetLowestUnoccupiedBand(kp)
                if (hob[2] < bandMinEnergy):
                    bestKp = kp
                    bandMinEnergy = hob[2]
                return bestKp
    #non spin-pol
    else:
        for kp in kpLis:
            if(GetLowestUnoccupiedBand(kp)[1] < bandMinEnergy):
                bestKp = kp
                bandMinEnergy = GetLowestUnoccupiedBand(kp)[1]
        return bestKp

    return None ##super uh oh

#Returns whether the current list of k-points signifies a direct bandgap semiconductor
def IsDirectBg(kpLis):
    vbm = GetVBMKpoint(kpLis)
    cbm = GetCBMKpoint(kpLis)
    return AlmostEqual([vbm.a, vbm.b, vbm.c][:], [cbm.a, cbm.b, cbm.c][:])

#Determines whether the material has a gap about eFermi or not (is a conductor)
def IsConductor(kpLis):
    mag = kpLis[0].isMag
    vbm = GetVBMKpoint(kpLis)
    cbm = GetCBMKpoint(kpLis)
    vbmBand = GetHighestOccupiedBand(vbm, spin=None) ##don't care about spin.  any vbm > cbm
    cbmBand = GetLowestUnoccupiedBand(cbm, spin=None) #indicates a conductor

    #if vbm is higher than cbm, then the bands cross and the material is a conductor
    if(mag):
        return (min(cbmBand[1], cbmBand[2]) - max(vbmBand[1], vbmBand[2]) <= 0.)
    return (cbmBand[1] - vbmBand[1] <= 0.)

#Gets the direct/indirect bandgap for a given kpoint.
def GetBandgap(instructions, dirPoint=None, vbmPoint=None, cbmPoint=None, spin=None):
    mag = -1 ##this should be initialized in instructions check.  If left as -1, will throw error
    ##for direct
    if(instructions[0] == 'd' or instructions[0] == 'D'):
        if(dirPoint == None):
            print("GetBandgap: Supply a k-point to use for the bandgap calculation\n")
            return None
        mag = int(dirPoint.isMag)
        vbmBand = GetHighestOccupiedBand(dirPoint, spin)
        cbmBand = GetLowestUnoccupiedBand(dirPoint, spin)
    ##for indirect
    if(instructions[0] == 'i' or instructions[0] == 'I'):
        if(vbmPoint == None or cbmPoint == None):
            print("GetBandgap: Supply a vbm and cbm k-point to use for the bandgap calculation\n")
            return None
        mag = int(vbmPoint.isMag)
        vbmBand = GetHighestOccupiedBand(vbmPoint, spin)
        cbmBand = GetLowestUnoccupiedBand(cbmPoint, spin)

    if(mag == 1):
        if(spin == None):
            return (min(cbmBand[1], cbmBand[2]) - max(vbmBand[1], vbmBand[2]))
        elif(spin[0] == 'u' or spin[0] == 'U'):
            return (cbmBand[1] - vbmBand[1])
        elif(spin[0] == 'd' or spin[0] == 'D'):
            return (cbmBand[2] - vbmBand[2])
    if(mag == 0):
        return (cbmBand[1] - vbmBand[1])

    print("GetBandgap: Invalid arguments")
    return None
