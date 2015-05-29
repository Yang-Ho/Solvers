import random
import sys

def rank(coord = [4, 1, 2, 3]):
    thisProc = "P.coord.rank"
    ABOUT = """
This function takes a permutation coordinate as a list such as [4, 1, 2, 3] and
returns its inversion number, which can also be interpreted as the distance 
from [1, 2, 3, 4] or also as 'the rank' of the coordinate in the Hasse graph 
with respect to its coordinate in natural order, here [1, 2, 3, 4].
"""
    if coord == "??":
        print ABOUT
        return
    elif coord == "?":
        print "Valid query is '{} ??'".format(thisProc)
        return
        #raise Exception("Valid query is '{} ??'".format(thisProc))
        #sys.stderr.write("Valid query is '{}' ??'\n".format(thisProc))
        #sys.exit(1)

    # Standard definition for permutation inversion: given a permutation 
    # b_1, b_2, b_3,..., b_n of the n integers 1, 2, 3, ..., n, 
    # an inversion is a pair (b_i, b_j) where i < j and b_i > b_j.
    size = len(coord)
    inversion = 0
    for i in xrange(size):
        b_i = coord[i]
        for j in xrange(i+1, size):
            b_j = coord[j]
            if (i < j and b_i > b_j):
                inversion += 1
    
    return inversion

def IS(coordP = [3, 1, 2, 4]):
    thisProc = "P.coord.IS"
    ABOUT = """
This proc takes a permutation coordinate as a list such as **4 1 2 3** and
returns 1 if coordinate is a valid permutation and returns 0 otherwise.
"""
    if coordP == "??":
        print ABOUT
        return
    elif coordP == "?":
        print "Valid query is '{} ??'".format(thisProc)
        return
        #raise Exception("Valid query is '{} ??'".format(thisProc))
        #sys.stderr.write("Valid query is '{}' ??'\n".format(thisProc))
        #sys.exit(1)

    elif coordP == "NA":
        return False 
    size = len(coordP)
    coordP = set(coordP)
    if (size == len(coordP)):
        # coordP is a valid permutation coordinate
        return True
    else:
        return False

def rand(nDim = 21):
    thisProc = "P.coord.rand"
    ABOUT = """
This proc take the length of the permutation coordinate (nDim) and returns
a random permutation of coordinate (coordType = P) of specified length.
"""
    if nDim == "??":
        print ABOUT
        return
    elif nDim == "?":
        print "Valid query is '{} ??'".format(thisProc)
        return
        #raise Exception("Valid query is '{} ??'".format(thisProc))
        #sys.stderr.write("Valid query is '{}' ??'\n".format(thisProc))
        #sys.exit(1)

    coord = [i+1 for i in range(nDim)]
    random.shuffle(coord)
    return coord
    
def distance(coord0 = [3, 1, 2, 3], coord1 = [4, 1, 2, 3]):
    thisProc = "P.coord.distance"
    ABOUT = """This procedure takes two permutations such as **3 1 2 4** and **4 1 2 3**
and returns the value of permutation distance, which is defined by the 
crossing number of edges in a two layer graph formed by the permutation at 
each layer. This formulation is a special case of formulas introduced in 
the  article 
JOHN N. WARFIELD
Crossing Theory and Hierarchy Mapping
IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS, 
VOL. SMC-6, NO. 7, JuLY 1977, pp. 505--524"""

    if coord0 == "??":
        print ABOUT
        return
    elif coord0 == "?":
        print "Valid query is '{} ??'".format(thisProc)
        return
        #raise Exception("Valid query is '{} ??'".format(thisProc))
        #sys.stderr.write("Valid query is '{}' ??'\n".format(thisProc))
        #sys.exit(1)
    distance = 0
    # create a position indices of elements in one of the coordinates 
    size = len(coord1)
    # parray pos ;# positions array 
    pos = [None] * (size + 1)

    for i in xrange(size):
        pos[coord0[i]] = i

    for i in xrange(size):
        a_i = pos[coord0[i]]
        b_j = pos[coord1[i]]
        #puts i,j...a_i,b_j=$i,$j...$a_i,$b_j

        for j in xrange(size):
            a_i1 = pos[coord0[j]]
            b_j1 = pos[coord1[j]]
            #puts i,j...a_i,b_j=$i,$j...$a_i1,$b_j1

            if a_i1 > a_i and b_j1 < b_j or a_i1 < a_i and b_j1 > b_j:
                distance += 1

    return distance / 2
    # the computed distance represents total number of edge crossings and
    # must be divided by 2

def neighborhood(coordP = ["a", "b", "c", "d", "e"]):
    thisProc = "P.coord.neighborhood"
    ABOUT = """
This proc takes a permutation coordinate such as **a b c d e** (here, of size 
L = 5) and returns the complete set of adjacent coordinates (i.e. coordinates
with the permutation distance of 1 from the reference permutation. The size 
of this set of adjacent (or neighborhood coordinates) is L-1."""
    if coordP == "??":
        print ABOUT
        return
    elif coordP == "?":
        print "Valid query is '{} ??'".format(thisProc)
        return
        #raise Exception("Valid query is '{} ??'".format(thisProc))
        #sys.stderr.write("Valid query is '{}' ??'\n".format(thisProc))
        #sys.exit(1)

    coordAdj = []
    L = len(coordP)
    Lm1 = L - 1
    elm_i = coordP[0]

    for i in xrange(Lm1):
        ip1 = i + 1
        swapL = list(coordP)
        # swap elements at i & ip1
        elm_ip1 = coordP[ip1]
        swapL[i] = elm_ip1

        if ip1 <= Lm1:
            swapL[ip1] = elm_i
            coordAdj.append(list(swapL))
            elm_i = coordP[ip1]

    print "coordAdj"
    for coord in coordAdj:
        print coord
    print

    return coordAdj

if __name__ == '__main__':
    inversion([8,2,5,3,1,7,4,6])
"""
    distance()
    neighborhood()
    print (inversion([4, 3, 2, 1]))
    print (inversion())
    print (inversion([3, 1, 2, 4]))
    print (inversion([4, 1, 2, 3]))
    print (is())
    print (is([10, 5, 7, 21, 17, 1, 19, 15, 13, 3, 11, 20, 4, 2, 12, 9, 8, 18, 16, 6, 13]))
    print (rand())
    print (rand(4))
"""
