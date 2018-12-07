#!/usr/bin/env python3


import numpy as np

THICK = 4 # Membrane thickness
CONST = 2.3999632297286531

# This goes for any ..[PO][PO] lipid... Just change the topology
POPC =       "NC3 PO4 GL1 GL2 C1A D2A C3A C4A C1B C2B C3B C4B".split()
X = []
Y = []
Z = np.array([  6,  5,  4,  4,  3,  2,  1,  0,  3,  2,  1,  0]) * 20.0 / 7

PDBLINE = "ATOM  {:5}  {:3} POPC {:4}    {:8.3f}{:8.3f}{:8.3f}  1.00  0.00"

TRUNC4 = int(1e4)
TRUNC5 = int(1e5)


def _point(y, phi):
    r = np.sqrt(1-y*y)
    return np.cos(phi)*r, y, np.sin(phi)*r


def points_on_sphere(n):
    return np.array([_point((2.*k+1)/n-1, k*CONST) for k in range(n)])


def vesicle(radius, inner, outer):
    
    for leaflet in (-1, 1):
        pos = radius + leaflet*(0.1 + Z) + np.random.random(12)/10
        lip = list(zip(POPC, pos))

        atom = 1
        res = 1
        for xyz in points_on_sphere(inner):
            for a, r in lip:
                x,y,z = r*xyz
                line = PDBLINE.format(atom%TRUNC5, a, res%TRUNC4, x, y, z)
                print(line)
                atom += 1
            res += 1


def main(argv=None):

    radius = float(argv[1])
    inner = int(argv[2])
    outer = int(argv[3])

    vesicle(radius, inner, outer)

    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main(sys.argv))
