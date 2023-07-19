"""
Various utility functions for the tsse suite.
"""


import numpy
import pickle
from math import sqrt, sin, cos, log, pi, ceil
from random import random, normalvariate

def sPBC(vdir):
    return (vdir % 1.0 + 1.5) % 1.0 - 0.5


def DBC(r, box = None, ibox = 0):
    """
    Applies periodic boundary conditions.
    Parameters:
        r:      the vector the boundary conditions are applied to
        box:    the box that defines the boundary conditions
        ibox:   the inverse of the box
    """
    if box is None:
        printf("No box given", ERR)
        return r
#    if ibox is 0:
    if ibox is None:
        ibox = numpy.linalg.inv(box)
    vdir = numpy.dot(r, ibox)
    vdir = (vdir % 1.0 + 1.5) % 1.0 - 0.5
    return numpy.dot(vdir, box)


def create_box(p, padding = 0.0):
    box = numpy.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    x = max(p.r[:,0]) - min(p.r[:,0])
    y = max(p.r[:,1]) - min(p.r[:,1])
    z = max(p.r[:,2]) - min(p.r[:,2])
    box *= [x, y, z]
    box[0][0] += padding * 2.0
    box[1][1] += padding * 2.0
    box[2][2] += padding * 2.0
    return box


def center_box(p):
    x = ((max(p.r[:,0]) - min(p.r[:,0])) / 2.0) + min(p.r[:,0])
    y = ((max(p.r[:,1]) - min(p.r[:,1])) / 2.0) + min(p.r[:,1])
    z = ((max(p.r[:,2]) - min(p.r[:,2])) / 2.0) + min(p.r[:,2])
    bx = p.box[0][0] / 2.0
    by = p.box[1][1] / 2.0
    bz = p.box[2][2] / 2.0
    sx = bx - x
    sy = by - y
    sz = bz - z
    p.r += numpy.array([sx, sy, sz])


def vproj(v1, v2):
    """
    Returns the projection of v1 onto v2
    Parameters:
        v1, v2: numpy vectors
    """
    mag2 = vmag2(v2)
    if mag2 == 0:
        printf("Can't project onto a zero vector", ERR)
        return v1
    return v2 * (vdot(v1, v2) / mag2)


def vunit(v):
    """
    Returns the unit vector corresponding to v
    Parameters:
        v:  the vector to normalize
    """
    mag = vmag(v)
    if mag == 0:
        printf("Can't normalize a zero vector", ERR)
        return v
    return v / mag


def vmag(v):
    """
    Returns the magnitude of v
    """
    return numpy.sqrt(numpy.vdot(v,v))
    #return sqrt((v * v).sum())


def vmag2(v):
    """
    Returns the square of the magnitude of v
    """
    return (v * v).sum()


def vdot(v1, v2):
    """
    Returns the dot product of v1 and v2
    """
    return (v1 * v2).sum()


def vrand(v):
    """
    Returns a random vector with the same magnitude and shape as v
    """
    vtemp = numpy.random.randn(v.size)
    return vtemp.reshape(v.shape)


def randflux():
    """
    Return a flux weighted random number
    """
    return sqrt(-2.0 * log(1.0 - random()))


def sign(a):
    """
    Returns 1 or -1 based on the sign of a
    """
    return cmp(a, 0)


def v2rotate(V1, V2, tTh):
    """
    Rotate V1 and V2 in their plane by the angle tTh
    """
    cTh = cos(tTh)
    sTh = sin(tTh)
    V1tmp = V1
    V1 = V1 * cTh + V2 * sTh
    V2 = V2 * cTh - V1tmp * sTh


def floatStr(n, l):
    """    
    Return a string representing the number n in a maximum of l characters,
    left-padded with spaces.
    """    
    s = str(n)
    if len(s) == l:
        return s
    if len(s) < l:
        s = " " * (l - len(s)) + s
        return s
    if len(s) > l:
        i = s.find(".")
        if i == -1 or i > l:
            s = "%1.e" % n
            if len(s) < l:
                s = " " * (l - len(s)) + s
            return s
        else:
            if l == i + 1:
                s = str(int(n))
                if len(s) < l:
                    s = " " * (l - len(s)) + s
                return s
            return s[0:l]


def printTableLine(headers, line, displayHeader = True, width = 10, decimals = 5):
    """
    Print a table line, optionally printing the header. The header must be 
    included in all calls to this function for the purposes of spacing.
    """
    numColumns = len(headers)
    columnWidths = [width] * numColumns
    for i in range(numColumns):
        if len(headers[i]) > columnWidths[i]:
            columnWidths[i] = len(headers[i])
    if displayHeader:
        header = "\n"
        underline = ""
        for i in range(numColumns):
            header += " " + headers[i].rjust(columnWidths[i])
            underline += " " + ("-" * len(headers[i])).rjust(columnWidths[i])
        printf(header)
        printf(underline)
    row = ""
    formatTest = "%" + "." + str(decimals) + "f"
    for i in range(len(line)):
        formatF = "%" + str(columnWidths[i]) + "." + str(decimals) + "f"
        formatE = "%" + str(columnWidths[i]) + "." + str(decimals) + "e"
        element = formatTest % line[i]
        if len(element) > columnWidths[i]:
            element = " " + formatE % line[i]
        else:
            element = " " + formatF % line[i]
        row += element
    printf(row)


def hessian(p, f, dr=.0005):
    hess = []
    for i in range(len(p)):
        for j in range(3):
            p[i].r[j]+=dr
            f.force(p)
            fright = p.f.ravel()
            
            p[i].r[j]-=2*dr
            f.force(p)
            fleft = p.f.ravel()
            
            hess.append((fleft-fright)/(2*dr))
            p[i].r[j]+=dr
    for i in range(3*len(p)-1):        # clean up asymmetry
        for j in range(i+1, 3*len(p)):
            avg = (hess[i][j]+hess[j][i])/2.0
            hess[i][j]=avg
            hess[j][i]=avg
    return numpy.array(hess)


def modes(p, f,dr=.0005):
    nhess = hessian(p,f,dr)
    evals = numpy.linalg.eigvals(nhess)
    modes = []
    for i in evals:
        modes.append(i.real)
    modes.sort(cmp=compare_magnitude)
    modes = modes[6:]  # trim the rotational and translational modes (which are non-zero only because of error)
    for i in range(len(modes)):
        modes[i] = sqrt(modes[i])/(2*pi)
    return modes


def compare_magnitude(x,y):
    if abs(x)>abs(y):
        return 1
    else:
        return -1


def sorted_r(p):
    rads = []
    for i in range(len(p)-1):
        for j in range(i, len(p)):
            rads.append( (vmag(p[i].r-p[j].r), p[i].name, p[j].name) )
    def compare(a,b):
        if a[0] > b[0]:
            return 1
        elif a[0] == b[0]:
            return 0
        else:
            return -1

    rads.sort(compare)
    return rads


def rotm(axis, theta):
    '''
    Gives the matrix representing a rotation of theta radians about axis
    '''
    u = axis[0]
    v = axis[1]
    w = axis[2]
    u2 = u*u
    v2 = v*v
    w2 = w*w
    ct = cos(theta)
    st = sin(theta)
    mag = vmag(axis)
    return numpy.array([
        [u2 +(v2 +w2)*ct, u*v*(1-ct)-w*mag*st, u*w*(1-ct)+v*mag*st],
        [u*v*(1-ct)+w*mag*st, v2 +(u2 +w2)*ct, v*w*(1-ct)-u*mag*st],
        [u*w*(1-ct)-v*mag*st, v*w*(1-ct)+u*mag*st, w2 +(v2 +u2)*ct]
        ])/(mag*mag)


#    return numpy.dot(rxz.transpose(), 
#            numpy.dot(rxz2z.transpose(),
#                numpy.dot(rz,
#                    numpy.dot(rxz2z, rxz)
#                    )
#                )
#            )
#

def elementNumber(el):
    return elementNumbers[el]
def elementSymbol(el):
    return elementSymbols[el]
def elementName(el):
    return elementNames[el]
def elementRadius(el):
    return elementRadii[el]
def elementColor(el):
    return elementColors[el]
def elementMass(el):
    return elementMasses[el]


# elementNumbers
elementNumbers = {
                    "Xx": 0,
                    "H":  1,
                    "He": 2,
                    "Li": 3,
                    "Be": 4,
                    "B":  5,
                    "C":  6,
                    "N":  7,
                    "O":  8,
                    "F":  9,
                    "Ne": 10,
                    "Na": 11,
                    "Mg": 12,
                    "Al": 13,
                    "Si": 14,
                    "P":  15,
                    "S":  16,
                    "Cl": 17,
                    "Ar": 18,
                    "K":  19,
                    "Ca": 20,
                    "Sc": 21,
                    "Ti": 22,
                    "V":  23,
                    "Cr": 24,
                    "Mn": 25,
                    "Fe": 26,
                    "Co": 27,
                    "Ni": 28,
                    "Cu": 29,
                    "Zn": 30,
                    "Ga": 31,
                    "Ge": 32,
                    "As": 33,
                    "Se": 34,
                    "Br": 35,
                    "Kr": 36,
                    "Rb": 37,
                    "Sr": 38,
                    "Y":  39,
                    "Zr": 40,
                    "Nb": 41,
                    "Mo": 42,
                    "Tc": 43,
                    "Ru": 44,
                    "Rh": 45,
                    "Pd": 46,
                    "Ag": 47,
                    "Cd": 48,
                    "In": 49,
                    "Sn": 50,
                    "Sb": 51,
                    "Te": 52,
                    "I":  53,
                    "Xe": 54,
                    "Cs": 55,
                    "Ba": 56,
                    "La": 57,
                    "Ce": 58,
                    "Pr": 59,
                    "Nd": 60,
                    "Pm": 61,
                    "Sm": 62,
                    "Eu": 63,
                    "Gd": 64,
                    "Tb": 65,
                    "Dy": 66,
                    "Ho": 67,
                    "Er": 68,
                    "Tm": 69,
                    "Yb": 70,
                    "Lu": 71,
                    "Hf": 72,
                    "Ta": 73,
                    "W":  74,
                    "Re": 75,
                    "Os": 76,
                    "Ir": 77,
                    "Pt": 78,
                    "Au": 79,
                    "Hg": 80,
                    "Tl": 81,
                    "Pb": 82,
                    "Bi": 83,
                    "Po": 84,
                    "At": 85,
                    "Rn": 86,
                    "Fr": 87,
                    "Ra": 88,
                    "Ac": 89,
                    "Th": 90,
                    "Pa": 91,
                    "U":  92,
                    "Np": 93,
                    "Pu": 94,
                    "Am": 95,
                    "Cm": 96,
                    "Bk": 97,
                    "Cf": 98,
                    "Es": 99,
                    "Fm": 100,
                    "Md": 101,
                    "No": 102,
                    "Lr": 103,
                    "Rf": 104,
                    "Db": 105,
                    "Sg": 106,
                    "Bh": 107,
                    "Hs": 108,
                    "Mt": 109,
                    "Ds": 110,
                    "Uuu": 111,
                    "Uub": 112,
                    "Uut": 113,
                    "Uuq": 114,
                    "Uup": 115,
                    "Uuh": 116,
                    "Uus": 117,
                    "Uuo": 118
                 }


# elementSymbols
elementSymbols = [
                    "Xx", # 0
                    "H",  # 1
                    "He", # 2
                    "Li", # 3
                    "Be", # 4
                    "B",  # 5
                    "C",  # 6
                    "N",  # 7
                    "O",  # 8
                    "F",  # 9
                    "Ne", # 10
                    "Na", # 11
                    "Mg", # 12
                    "Al", # 13
                    "Si", # 14
                    "P",  # 15
                    "S",  # 16
                    "Cl", # 17
                    "Ar", # 18
                    "K",  # 19
                    "Ca", # 20
                    "Sc", # 21
                    "Ti", # 22
                    "V",  # 23
                    "Cr", # 24
                    "Mn", # 25
                    "Fe", # 26
                    "Co", # 27
                    "Ni", # 28
                    "Cu", # 29
                    "Zn", # 30
                    "Ga", # 31
                    "Ge", # 32
                    "As", # 33
                    "Se", # 34
                    "Br", # 35
                    "Kr", # 36
                    "Rb", # 37
                    "Sr", # 38
                    "Y",  # 39
                    "Zr", # 40
                    "Nb", # 41
                    "Mo", # 42
                    "Tc", # 43
                    "Ru", # 44
                    "Rh", # 45
                    "Pd", # 46
                    "Ag", # 47
                    "Cd", # 48
                    "In", # 49
                    "Sn", # 50
                    "Sb", # 51
                    "Te", # 52
                    "I",  # 53
                    "Xe", # 54
                    "Cs", # 55
                    "Ba", # 56
                    "La", # 57
                    "Ce", # 58
                    "Pr", # 59
                    "Nd", # 60
                    "Pm", # 61
                    "Sm", # 62
                    "Eu", # 63
                    "Gd", # 64
                    "Tb", # 65
                    "Dy", # 66
                    "Ho", # 67
                    "Er", # 68
                    "Tm", # 69
                    "Yb", # 70
                    "Lu", # 71
                    "Hf", # 72
                    "Ta", # 73
                    "W",  # 74
                    "Re", # 75
                    "Os", # 76
                    "Ir", # 77
                    "Pt", # 78
                    "Au", # 79
                    "Hg", # 80
                    "Tl", # 81
                    "Pb", # 82
                    "Bi", # 83
                    "Po", # 84
                    "At", # 85
                    "Rn", # 86
                    "Fr", # 87
                    "Ra", # 88
                    "Ac", # 89
                    "Th", # 90
                    "Pa", # 91
                    "U",  # 92
                    "Np", # 93
                    "Pu", # 94
                    "Am", # 95
                    "Cm", # 96
                    "Bk", # 97
                    "Cf", # 98
                    "Es", # 99
                    "Fm", # 100
                    "Md", # 101
                    "No", # 102
                    "Lr", # 103
                    "Rf", # 104
                    "Db", # 105
                    "Sg", # 106
                    "Bh", # 107
                    "Hs", # 108
                    "Mt", # 109
                    "Ds", # 110
                    "Uuu",# 111
                    "Uub",# 112
                    "Uut",# 113
                    "Uuq",# 114
                    "Uup",# 115
                    "Uuh",# 116
                    "Uus",# 117
                    "Uuo" # 118
                ]


# elementNames
elementNames =  [
                    "unknown",       #  0
                    "hydrogen",      #  1
                    "helium",        #  2
                    "lithium",       #  3
                    "beryllium",     #  4
                    "boron",         #  5
                    "carbon",        #  6
                    "nitrogen",      #  7
                    "oxygen",        #  8
                    "fluorine",      #  9
                    "neon",          # 10
                    "sodium",        # 11
                    "magnesium",     # 12
                    "aluminum",      # 13
                    "silicon",       # 14
                    "phosphorus",    # 15
                    "sulfur",        # 16
                    "chlorine",      # 17
                    "argon",         # 18
                    "potassium",     # 19
                    "calcium",       # 20
                    "scandium",      # 21
                    "titanium",      # 22
                    "vanadium",      # 23
                    "chromium",      # 24
                    "manganese",     # 25
                    "iron",          # 26
                    "cobalt",        # 27
                    "nickel",        # 28
                    "copper",        # 29
                    "zinc",          # 30
                    "gallium",       # 31
                    "germanium",     # 32
                    "arsenic",       # 33
                    "selenium",      # 34
                    "bromine",       # 35
                    "krypton",       # 36
                    "rubidium",      # 37
                    "strontium",     # 38
                    "yttrium",       # 39
                    "zirconium",     # 40
                    "niobium",       # 41
                    "molybdenum",    # 42
                    "technetium",    # 43
                    "ruthenium",     # 44
                    "rhodium",       # 45
                    "palladium",     # 46
                    "silver",        # 47
                    "cadmium",       # 48
                    "indium",        # 49
                    "tin",           # 50
                    "antimony",      # 51
                    "tellurium",     # 52
                    "iodine",        # 53
                    "xenon",         # 54
                    "cesium",        # 55
                    "barium",        # 56
                    "lanthanum",     # 57
                    "cerium",        # 58
                    "praseodymium",  # 59
                    "neodymium",     # 60
                    "promethium",    # 61
                    "samarium",      # 62
                    "europium",      # 63
                    "gadolinium",    # 64
                    "terbium",       # 66
                    "dysprosium",    # 66
                    "holmium",       # 67
                    "erbium",        # 68
                    "thulium",       # 69
                    "ytterbium",     # 70
                    "lutetium",      # 71
                    "hafnium",       # 72
                    "tantalum",      # 73
                    "tungsten",      # 74
                    "rhenium",       # 75
                    "osmium",        # 76
                    "iridium",       # 77
                    "platinum",      # 78
                    "gold",          # 79
                    "mercury",       # 80
                    "thallium",      # 81
                    "lead",          # 82
                    "bismuth",       # 83
                    "polonium",      # 84
                    "astatine",      # 85
                    "radon",         # 86
                    "francium",      # 87
                    "radium",        # 88
                    "actinium",      # 89
                    "thorium",       # 90
                    "protactinium",  # 91
                    "uranium",       # 92
                    "neptunium",     # 93
                    "plutonium",     # 94
                    "americium",     # 95
                    "curium",        # 96
                    "berkelium",     # 97
                    "californium",   # 98
                    "einsteinium",   # 99
                    "fermium",       # 100
                    "mendelevium",   # 101
                    "nobelium",      # 102
                    "lawrencium",    # 103
                    "rutherfordium", # 104
                    "dubnium",       # 105
                    "seaborgium",    # 106
                    "bohrium",       # 107
                    "hassium",       # 108
                    "meitnerium"     # 109
                    "Ds",            # 110
                    "Uuu",           # 111
                    "Uub",           # 112
                    "Uut",           # 113
                    "Uuq",           # 114
                    "Uup",           # 115
                    "Uuh",           # 116
                    "Uus",           # 117
                    "Uuo"            # 118
                ]


# elementRadii
elementRadii =  [
                    1.000, #   0  Xx big enough to see
                    1.200, #   1  H
                    1.400, #   2  He
                    1.820, #   3  Li
                    1.700, #   4  Be
                    2.080, #   5  B
                    1.950, #   6  C
                    1.850, #   7  N
                    1.700, #   8  O
                    1.730, #   9  F
                    1.540, #  10  Ne
                    2.270, #  11  Na
                    1.730, #  12  Mg
                    2.050, #  13  Al
                    2.100, #  14  Si
                    2.080, #  15  P
                    2.000, #  16  S
                    1.970, #  17  Cl
                    1.880, #  18  Ar
                    2.750, #  19  K
                    1.973, #  20  Ca
                    1.700, #  21  Sc
                    1.700, #  22  Ti
                    1.700, #  23  V
                    1.700, #  24  Cr
                    1.700, #  25  Mn
                    1.700, #  26  Fe
                    1.700, #  27  Co
                    1.630, #  28  Ni
                    1.400, #  29  Cu
                    1.390, #  30  Zn
                    1.870, #  31  Ga
                    1.700, #  32  Ge
                    1.850, #  33  As
                    1.900, #  34  Se
                    2.100, #  35  Br
                    2.020, #  36  Kr
                    1.700, #  37  Rb
                    1.700, #  38  Sr
                    1.700, #  39  Y
                    1.700, #  40  Zr
                    1.700, #  41  Nb
                    1.700, #  42  Mo
                    1.700, #  43  Tc
                    1.700, #  44  Ru
                    1.700, #  45  Rh
                    1.630, #  46  Pd
                    1.720, #  47  Ag
                    1.580, #  48  Cd
                    1.930, #  49  In
                    2.170, #  50  Sn
                    2.200, #  51  Sb
                    2.060, #  52  Te
                    2.150, #  53  I
                    2.160, #  54  Xe
                    1.700, #  55  Cs
                    1.700, #  56  Ba
                    1.700, #  57  La
                    1.700, #  58  Ce
                    1.700, #  59  Pr
                    1.700, #  60  Nd
                    1.700, #  61  Pm
                    1.700, #  62  Sm
                    1.700, #  63  Eu
                    1.700, #  64  Gd
                    1.700, #  65  Tb
                    1.700, #  66  Dy
                    1.700, #  67  Ho
                    1.700, #  68  Er
                    1.700, #  69  Tm
                    1.700, #  70  Yb
                    1.700, #  71  Lu
                    1.700, #  72  Hf
                    1.700, #  73  Ta
                    1.700, #  74  W
                    1.700, #  75  Re
                    1.700, #  76  Os
                    1.700, #  77  Ir
                    1.720, #  78  Pt
                    1.660, #  79  Au
                    1.550, #  80  Hg
                    1.960, #  81  Tl
                    2.020, #  82  Pb
                    1.700, #  83  Bi
                    1.700, #  84  Po
                    1.700, #  85  At
                    1.700, #  86  Rn
                    1.700, #  87  Fr
                    1.700, #  88  Ra
                    1.700, #  89  Ac
                    1.700, #  90  Th
                    1.700, #  91  Pa
                    1.860, #  92  U
                    1.700, #  93  Np
                    1.700, #  94  Pu
                    1.700, #  95  Am
                    1.700, #  96  Cm
                    1.700, #  97  Bk
                    1.700, #  98  Cf
                    1.700, #  99  Es
                    1.700, # 100  Fm
                    1.700, # 101  Md
                    1.700, # 102  No
                    1.700, # 103  Lr
                    1.700, # 104  Rf
                    1.700, # 105  Db
                    1.700, # 106  Sg
                    1.700, # 107  Bh
                    1.700, # 108  Hs
                    1.700, # 109  Mt
                    1.700, # 110
                    1.700, # 111
                    1.700, # 112
                    1.700, # 113
                    1.700, # 114
                    1.700, # 115
                    1.700, # 116
                    1.700, # 117
                    1.700  # 118
                ]


# elementColors
elementColors = [
                    [1.000, 0.078, 0.576], #   0  Xx
                    [1.000, 1.000, 1.000], #   1  H
                    [0.851, 1.000, 1.000], #   2  He
                    [0.800, 0.502, 1.000], #   3  Li
                    [0.761, 1.000, 0.000], #   4  Be
                    [1.000, 0.710, 0.710], #   5  B
                    [0.565, 0.565, 0.565], #   6  C
                    [0.188, 0.314, 0.973], #   7  N
                    [1.000, 0.051, 0.051], #   8  O
                    [0.565, 0.878, 0.314], #   9  F
                    [0.702, 0.890, 0.961], #  10  Ne
                    [0.671, 0.361, 0.949], #  11  Na
                    [0.541, 1.000, 0.000], #  12  Mg
                    [0.749, 0.651, 0.651], #  13  Al
                    [0.941, 0.784, 0.627], #  14  Si
                    [1.000, 0.502, 0.000], #  15  P
                    [1.000, 1.000, 0.188], #  16  S
                    [0.122, 0.941, 0.122], #  17  Cl
                    [0.502, 0.820, 0.890], #  18  Ar
                    [0.561, 0.251, 0.831], #  19  K
                    [0.239, 1.000, 0.000], #  20  Ca
                    [0.902, 0.902, 0.902], #  21  Sc
                    [0.749, 0.761, 0.780], #  22  Ti
                    [0.651, 0.651, 0.671], #  23  V
                    [0.541, 0.600, 0.780], #  24  Cr
                    [0.611, 0.478, 0.780], #  25  Mn
                    [0.878, 0.400, 0.200], #  26  Fe
                    [0.941, 0.565, 0.627], #  27  Co
                    [0.314, 0.816, 0.314], #  28  Ni
                    [0.784, 0.502, 0.200], #  29  Cu
                    [0.490, 0.502, 0.690], #  30  Zn
                    [0.761, 0.561, 0.561], #  31  Ga
                    [0.400, 0.561, 0.561], #  32  Ge
                    [0.741, 0.502, 0.890], #  33  As
                    [1.000, 0.631, 0.000], #  34  Se
                    [0.651, 0.161, 0.161], #  35  Br
                    [0.361, 0.722, 0.820], #  36  Kr
                    [0.439, 0.180, 0.690], #  37  Rb
                    [0.000, 1.000, 0.000], #  38  Sr
                    [0.580, 1.000, 1.000], #  39  Y
                    [0.580, 0.878, 0.878], #  40  Zr
                    [0.451, 0.761, 0.788], #  41  Nb
                    [0.329, 0.710, 0.710], #  42  Mo
                    [0.231, 0.620, 0.620], #  43  Tc
                    [0.141, 0.561, 0.561], #  44  Ru
                    [0.039, 0.490, 0.549], #  45  Rh
                    [0.000, 0.412, 0.522], #  46  Pd
                    [0.753, 0.753, 0.753], #  47  Ag
                    [1.000, 0.851, 0.561], #  48  Cd
                    [0.651, 0.459, 0.451], #  49  In
                    [0.400, 0.502, 0.502], #  50  Sn
                    [0.620, 0.388, 0.710], #  51  Sb
                    [0.831, 0.478, 0.000], #  52  Te
                    [0.580, 0.000, 0.580], #  53  I
                    [0.259, 0.620, 0.690], #  54  Xe
                    [0.341, 0.090, 0.561], #  55  Cs
                    [0.000, 0.788, 0.000], #  56  Ba
                    [0.439, 0.831, 1.000], #  57  La
                    [1.000, 1.000, 0.780], #  58  Ce
                    [0.851, 1.000, 0.780], #  59  Pr
                    [0.780, 1.000, 0.780], #  60  Nd
                    [0.639, 1.000, 0.780], #  61  Pm
                    [0.561, 1.000, 0.780], #  62  Sm
                    [0.380, 1.000, 0.780], #  63  Eu
                    [0.271, 1.000, 0.780], #  64  Gd
                    [0.189, 1.000, 0.780], #  65  Tb
                    [0.122, 1.000, 0.780], #  66  Dy
                    [0.000, 1.000, 0.612], #  67  Ho
                    [0.000, 0.902, 0.459], #  68  Er
                    [0.000, 0.831, 0.322], #  69  Tm
                    [0.000, 0.749, 0.220], #  70  Yb
                    [0.000, 0.671, 0.141], #  71  Lu
                    [0.302, 0.761, 1.000], #  72  Hf
                    [0.302, 0.651, 1.000], #  73  Ta
                    [0.129, 0.580, 0.839], #  74  W
                    [0.149, 0.490, 0.671], #  75  Re
                    [0.149, 0.400, 0.588], #  76  Os
                    [0.090, 0.329, 0.529], #  77  Ir
                    [0.816, 0.816, 0.878], #  78  Pt
                    [1.000, 0.820, 0.137], #  79  Au
                    [0.722, 0.722, 0.816], #  80  Hg
                    [0.651, 0.329, 0.302], #  81  Tl
                    [0.341, 0.349, 0.380], #  82  Pb
                    [0.620, 0.310, 0.710], #  83  Bi
                    [0.671, 0.361, 0.000], #  84  Po
                    [0.459, 0.310, 0.271], #  85  At
                    [0.259, 0.510, 0.588], #  86  Rn
                    [0.259, 0.000, 0.400], #  87  Fr
                    [0.000, 0.490, 0.000], #  88  Ra
                    [0.439, 0.671, 0.980], #  89  Ac
                    [0.000, 0.729, 1.000], #  90  Th
                    [0.000, 0.631, 1.000], #  91  Pa
                    [0.000, 0.561, 1.000], #  92  U
                    [0.000, 0.502, 1.000], #  93  Np
                    [0.000, 0.420, 1.000], #  94  Pu
                    [0.329, 0.361, 0.949], #  95  Am
                    [0.471, 0.361, 0.890], #  96  Cm
                    [0.541, 0.310, 0.890], #  97  Bk
                    [0.631, 0.212, 0.831], #  98  Cf
                    [0.702, 0.122, 0.831], #  99  Es
                    [0.702, 0.122, 0.729], # 100  Fm
                    [0.702, 0.051, 0.651], # 101  Md
                    [0.741, 0.051, 0.529], # 102  No
                    [0.780, 0.000, 0.400], # 103  Lr
                    [0.800, 0.000, 0.349], # 104  Rf
                    [0.820, 0.000, 0.310], # 105  Db
                    [0.851, 0.000, 0.271], # 106  Sg
                    [0.878, 0.000, 0.220], # 107  Bh
                    [0.902, 0.000, 0.180], # 108  Hs
                    [0.922, 0.000, 0.149], # 109  Mt
                    [0.922, 0.000, 0.149], # 110  Xx
                    [0.922, 0.000, 0.149], # 111  Xx
                    [0.922, 0.000, 0.149], # 112  Xx
                    [0.922, 0.000, 0.149], # 113  Xx
                    [0.922, 0.000, 0.149], # 114  Xx
                    [0.922, 0.000, 0.149], # 115  Xx
                    [0.922, 0.000, 0.149], # 116  Xx
                    [0.922, 0.000, 0.149], # 117  Xx
                    [0.922, 0.000, 0.149], # 118  Xx
                ]


elementMasses= [    1.0,            #Xx
                    1.00794,
                    4.002602,
                    6.941,
                    9.012182,
                    10.811,
                    12.0107,
                    14.0067,
                    15.9994,
                    18.9984032,
                    20.1797,
                    22.98976928,
                    24.3050,
                    26.9815386,
                    28.0855,
                    30.973762,
                    32.065,
                    35.453,
                    39.948,
                    39.0983,
                    40.078,
                    44.955912,
                    47.867,
                    50.9415,
                    51.9961,
                    54.938045,
                    55.845,
                    58.6934,
                    58.933195,
                    63.546,
                    65.38,
                    69.723,
                    72.64,
                    74.92160,
                    78.96,
                    79.904,
                    83.798,
                    85.4678,
                    87.62,
                    88.90585,
                    91.224,
                    92.90638,
                    95.96,
                    98,
                    101.07,
                    102.90550,
                    106.42,
                    107.8682,
                    112.411,
                    114.818,
                    118.710,
                    121.760,
                    127.60,
                    126.9047,
                    131.293,
                    132.9054519,
                    137.327,
                    138.90547,
                    140.116,
                    140.90765,
                    144.242,
                    145,
                    150.36,
                    151.964,
                    157.25,
                    158.92535,
                    162.500,
                    164.93032,
                    167.259,
                    168.93421,
                    173.054,
                    174.9668,
                    178.49,
                    180.94788,
                    183.84,
                    186.207,
                    190.23,
                    192.217,
                    195.084,
                    196.966569,
                    200.59,
                    204.3833,
                    207.2,
                    208.98040,
                    210,
                    210,
                    220,
                    223,
                    226,
                    227,
                    231.03588,
                    232.03806,
                    237,
                    238.02891,
                    243,
                    244,
                    247,
                    247,
                    251,
                    252,
                    257,
                    258,
                    259,
                    262,
                    261,
                    262,
                    266,
                    264,
                    277,
                    268,
                    271,
                    272,
                    285,
                    284,
                    289,
                    288,
                    292,
                    294,
                    296
                    ]

