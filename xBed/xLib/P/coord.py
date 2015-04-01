import random

def inversion(coord = [4, 1, 2, 3]):
    """
    This proc takes a permutation coordinate as a list such as **4 1 2 3**
    and returns its inversion number, which can also be interpreted as the 
    distance from **1 2 3 4** or also as the rank of the coordinate in the
    Hasse graph with respect to its bottom coordinate, here **1 2 3 4**.
    Standard definition for permutation inversion:
    Given a permutation b_1, b_2, b_3,..., b_n of the n integers 1, 2, 3, ..., n, 
    an inversion is a pair (b_i, b_j) where i < j and b_i > b_j.
    """
    size = len(coord)
    inversion = 0
    for i in xrange(size):
        b_i = coord[i]
        for j in xrange(i+1, size):
            b_j = coord[j]
            if (i < j and b_i > b_j):
                inversion += 1
    
    return inversion

def IS(coordP = [10, 5, 7, 21, 17, 1, 19, 15, 13, 3, 11, 20, 4, 2, 12, 9, 8, 18, 16, 6, 14]):
    if coordP == "NA":
        return 0
    size = len(coordP)
    coordP = set(coordP)
    if (size == len(coordP)):
        return 1
    else:
        return 0

def rand(nDim = 41):
    """
    This proc take the dimension number (nDim) and returns
    a random permutation coordinate (coordType = P) of length nDim.
    """
    coord = [i+1 for i in range(nDim)]
    random.shuffle(coord)
    return coord
    
def distance(coord0 = [2, 1, 4, 3, 6, 5, 8, 7, 9], coord1 = [9, 7, 8, 5, 6, 3, 4, 1, 2]):
    '''
    This procedure takes two permutations such as **3 1 2 4** and
    **4 1 2 3** and returns the value of permutation distance, which is
    defined by the crossing number of edges in a two layer graph formed
    by the permutation at each layer. For details and a citation, see the
    example under the comments at the tail of this proc.
    '''
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

    print(distance / 2)
    return distance / 2
    # the computed distance represents total number of edge crossings and
    # must be divided by 2

def neighborhood(coordP = ["a", "b", "c", "d", "e"]):
    '''
    This proc takes a permutation coordinate such as **a b c d e**
    (here of size L = 5) and returns a set of permutation_distance=1 neighborhood
    coordinates. The size of this set is L-1. See the example below.
    '''
    print("coordAdj")
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

        print(swapL)
    
    print "\n", coordAdj
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
