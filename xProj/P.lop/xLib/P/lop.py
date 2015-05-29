from itertools import imap
import P.coord
import time
import sys
import os
import pwd
import platform
import random
import util
import core

aLits = None
aAdjacent = None
aPI = {}
aV = {}
aValueBest = {}
aWalkBest = {}
aWalk = {}
aWalkProbed = {}
aCoordHash0 = {}
reqAry = {}
optAry = {}
optAryNames = []
optAryNamesBool = []
optAryList = []

def saw_pivot_bee():
    return

def saw_pivot_ant():
    return

def saw_pivot_simple(coordPiv=[5,3,2,1,4], valuePiv=-46):
    thisProc = "P.lop.saw.pivot.simple"
    ABOUT = (
        "This procedure takes a pivot coordinate/value, probes the distance=1"
        " neighborhood of a 'lop' (a linear ordering problem), subject to the "
        "constraints of a SAW (self-avoiding walk) -- i.e. the best "
        "coord/value it returns has not been yet been selected as the pivot "
        "for the next step. Neighborhood size of 0 signifies that the next "
        "step of a SAW is blocked.\n\tThis implementation is 'simple', i.e. "
        "for each pivot coordinate of length L, there are up to L-1 explicit "
        "probes of the function P.lop.f"
        )
    if coordPiv == "??":
        print ABOUT
        return
    elif coordPiv == "?":
        print "Valid query is '{} ??'".format(thisProc)
        return

    #info global variables
    global all_info
    global all_valu
    global aV
    #instance global variables
    global aStruc
    #solver global variables
    global aCoordHash0
    """
    global aHashTmp
    global aHashNeighb
    global aWalk 
    global aHashWalk 
    global aWalkProbed 
    global aWalkBest 
    global aValueBest 
    global aAdjacent
    """

    try:
        aV["writeVar"] = int(aV["writeVar"])
    except:
        print "writeVar is not an int"

    coordBest = None
    coordBestList = []
    valueProbedList = []
    coordProbedList = []
    valueBest = 2147483641
    neighbSize = 0

    L = aV["nDim"]
    Lm1 = L -1
    swapList = []
    coordAdj = {}

    elm_i = coordPiv[0]
    for i in range(Lm1):
        ip1 = i + 1
        swapL = coordPiv[:]
        elm_ip1 = coordPiv[ip1]
        swapL[i] = elm_ip1
        if ip1 <= Lm1:
            swapL[ip1] = elm_i
            coordAdj[ip1] = swapL
            elm_i = coordPiv[ip1]
    
    if aV["writeVar"] == 3:
        rank = P.coord.distance(coordPiv, aV["coordRef"])
        rowLines = ("\nFROM: {}"
                "\nProbing ALL distance=1 neighbors of the\ninitial pivot"
                " coordinate = {}"
                "\npair\tcoord\tf\tfBest\trank\tcntNeighb\tcntProbe"
                "\n--\t{}\t{}\t{}\t{}\t0\t{}\n".format(thisProc, coordPiv, 
                    coordPiv, valuePiv, valuePiv, rank, aV["cntProbe"]))

    for i in range(Lm1):
        ip1 = i + 1
        coordA = coordAdj[ip1]
        if tuple(coordA) not in aCoordHash0:
            neighbSize += 1
            valueA = f(coordA)
            aV["cntProbe"] += 1
            if aV["isWalkTables"]:
                coordProbedList.append(coordA)
                valueProbedList.append(valueA)
            if valueA <= valueBest:
                if valueA < valueBest:
                    coordBestList = []
                valueBest = valueA
                coordBest = coordA
                #aAdjacent[valueA] = coordA
                coordBestList.append(coordBest)
            if aV['writeVar'] == 3:
                iP = coordA[i]
                iP1 = coordA[ip1]
                pair = (iP, iP1)
                rank = P.coord.distance(coordA, aV["coordRef"])
                rowLines += "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(pair,
                        coordA, valueA, valueBest, rank, neighbSize, aV["cntProbe"])

    if aV["writeVar"] == 3:
        print rowLines

    idx = int(len(coordBestList)*random.random())
    if len(coordBestList) > 0:
        coordBest = coordBestList[idx]
    else:
        coordBest = None
    return (coordBest, valueBest, neighbSize, coordProbedList, valueProbedList)

def saw_pivot(coordPiv=[5,3,2,1,4], valuePiv=-46):
    thisProc = "P.lop.saw.pivot"
    ABOUT = ( "This procedure takes a pivot coordinate/value, probes the "
        "distance=1 neighborhood of a 'lop' (a linear ordering probelm), "
        "subject to the constraints of a SAW (self-avoiding walk) -- i.e. "
        "the best coord/value it returns has not been yet been selected as "
        "the pivot for the next step. Neighborhood size of 0 signifies that "
        "the next step of a SAW is blocked.\n This implementation is 'FAST', "
        "i.e. for each pivot coordinate of length L, there are up to L-1 FAST "
        "tableau-based probes of each pivot coordinate.")

    if coordPiv == "??":
        print ABOUT
        return
    elif coordPiv == "?":
        print "Valid query is '{} ??'".format(thisProc)
        return
        #raise Exception("Valid query is '{} ??'".format(thisProc))
        #sys.stderr.write("Valid query is '{}' ??'\n".format(thisProc))
        #sys.exit(1)

    # info global variables
    global all_info
    global all_valu
    global aV
    # instance global variables
    global aStruc
    # Solver global variables
    global aCoordHash0
    global aWalkProbed

    coordBest = None
    valueBest = 2147483641 #@TODO: change to max int
    valueProbedList = []
    coordBestList = []
    neighbSize = 0
    coordProbedList = []

    # PASS 1: given coordPiv, get sumP(i), valuePiv, coordAdj
    L = aV["nDim"]
    Lm1 = L - 1
    swapList = []
    coordAdj = {}

    # needed for tableau probing
    valuePiv = 0
    # first elm_i to be swapped 
    elm_i = coordPiv[0]

    sumP = {}
    for i in range(L):
        #compute valPiv (pivot value)
        ip1 = i + 1
        iP = coordPiv[i]
        sumP[iP] = 0 # needed for PASS 2
        for j in range(ip1, L):
            jP = coordPiv[j]
            sumP[iP] -= aStruc[jP-1][iP-1]

        valuePiv += sumP[iP]

        # initialize for coordAdj for phase 2
        swapL = coordPiv[:]
        elm_ip1 = coordPiv[ip1-1]
        swapL[i] = elm_ip1
        if ip1 <= Lm1:
            swapL[ip1] = elm_i
            coordAdj[ip1] = swapL
            elm_i = coordPiv[ip1]
    aV["cntProbe"] += 1
    
    if aV["writeVar"] == 3:
        rank = P.coord.distance(coordPiv, aV["coordRef"])
        rowLines = ("\nFROM: {}"
            "\nProbing ALL distance=1 neighbors of the"
            "\ninitial pivot coordinate = {}"
            "\npair\tcoord\tf\tfBest\trank\tcntNeighb\tcntprobe"
            "\n--\t{}\t{}\t{}\t{}\0\t{}\n".format(thisProc, coordPiv, coordPiv,
                valuePiv,valuePiv,rank,aV["cntProbe"]))

    # PASS 2: get valueA for ea. admissible adjacent coordinate
    # Find adjacent coordinates (all coordinates in the walk are excluded from
    #       the neighborhood of the current pivot
    for i in range(Lm1):
        ip1 = i + 1
        coordA = coordAdj[ip1]
        if tuple(coordA) not in aCoordHash0:
            neighbSize += 1
            iP = coordA[i]
            dif_ij = -sumP[iP]
            iP1 = coordA[ip1]
            dif_ji = -sumP[iP1]

            for j in range(ip1, L):
                ij = coordAdj[ip1][j-1]
                dif_ij -= aStruc[iP-1][ij-1]
                ji = coordAdj[ip1][j-1]
                dif_ji -= aStruc[iP1-1][ji-1]

            valueA = valuePiv + dif_ij + dif_ji
            aV["cntProbe"] += 1

            # aggregate coordBestList for random selection
            if valueA <= valueBest:
                if valueA < valueBest:
                    coordBestList = []
                valueBest = valueA
                coordBest = coordA
                coordBestList.append(coordBest)
            if aV["writeVar"] == 3:
                pair = "{},{}".format(iP,iP1)
                print "coordA={}".format(coordA)
                rank = P.coord.distance(coordA, aV["coordRef"])
                rowLines += "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(pair,coordA,
                        valueA,valueBest,rank,neighbSize,aV["cntProbe"])

    if aV["writeVar"] == 3:
        print rowLines

    idx = int(len(coordBestList)*random.random())
    if len(coordBestList) > 0:
        coordBest = coordBestList[idx]
    else:
        coordBest = None
    return (coordBest, valueBest, neighbSize, coordProbedList, valueProbedList)

def saw(arg=""):
    thisProc = "P.lop.saw"
    ABOUT = """
Procedure $thisProc takes global array values initialized under P.lop.init
and constructs a segment of a self-avoiding walk (SAW). Either P.lop.saw.pivot.simple 
or the significantly more efficient procedure P.lop.saw.pivot.ant is invoked.
More to come ....
"""
    if arg == "??":
        print ABOUT
        return
    elif arg == "?":
        print "Valid query is '{} ??'".format(thisProc)
        return
        #raise Exception("Valid query is '{} ??'".format(thisProc))
        #sys.stderr.write("Valid query is '{}' ??'\n".format(thisProc))
        #sys.exit(1)

    # info global variables
    global all_info
    global all_valu
    global aV
    # instance global variables
    global aStruc
    # solver global variables
    global aCoordHash0
    global aWalkProbed

    # primary input variables
    functionID = aV["functionID"]
    runtimeLmt = aV["runtimeLmt"] 
    cntProbeLmt = aV["cntProbeLmt"]
    #walkRepeatsLmt = aV["walkRepeatsLmt"]
    #walkIntervalLmt = aV["walkIntervalLmt"]
    walkSegmLmt = aV["walkSegmLmt"]
    valueTarget = aV["valueTarget"]

    if aV["isSimple"]:
        procPivotNext = saw_pivot_simple
    else:
        procPivotNext = saw_pivot
    print "# FROM: {}, searching for pivotBest via {}".format(thisProc, 
            procPivotNext.__name__)
    # auxiliary variables
    aV["coordPivot"] = aV["coordInit"][:]
    aV["valuePivot"] = aV["valueInit"]
    aV["coordBest"] = aV["coordPivot"][:]
    aV["valueBest"] = aV["valuePivot"]
    step = 0

    while True:
        # Timing
        microSecs = time.time()

        # UPDATE walkLength
        step += 1
        aV["walkLength"] += 1
        # PROBE neighborhood of current pivot
        bestNeighb = procPivotNext(aV["coordPivot"], aV["valuePivot"])

        # SELECT next pivot
        aV["coordPivot"] = bestNeighb[0]
        aV["valuePivot"] = bestNeighb[1]
        aV["neighbSize"] = bestNeighb[2]
        if aV["isWalkTables"]:
            aV["coordProbedList"] = bestNeighb[3]
            aV["valueProbedList"] = bestNeighb[4]

            walkLengthM1 = aV["walkLength"] - 1
            neighbSize = 0
            isPivot = 0
            for coord, value in aV["coordProbedList"], aV["valueProbedList"]:
                neighbSize += 1
                rank = P.coord.rank(coord)
                aWalkProbed[(walkLengthM1,neighbSize)] = ( walkLengthM1,
                        aV["cntRestart"], coord, value, rank, isPivot, 
                        neighbSize, None)
            isPivot = 1
            aV["rankPivot"] = P.coord.rank(aV["coordPivot"])
            aWalkProbed[(aV["walkLength"],0)] = (aV["walkLength"], 
                    aV["cntRestart"], aV["coordPivot"], aV["valuePivot"],
                    aV["rankPivot"], isPivot, aV["neighbSize"], aV["cntProbe"])

        if aV["coordPivot"] is not None:
            aCoordHash0[tuple(aV["coordPivot"])] = []

        # UPDATE valueBest, aValueBest, aWalkBest
        if aV["valuePivot"] <= aV["valueBest"]:
            aV["valueBest"] = aV["valuePivot"]
            aV["coordBest"] = aV["coordPivot"]
        # CHECK the nighboroodSize
        if aV["neighbSize"] == 0:
            aV["isBlocked"] = 1
            aV["speedProbe"] = int(aV["cntProbe"]/aV["runtime"])
            print ("WARNING from {}: isBlocked=1, aV['neighbSize']={} ..."
                    "no available neighborhood coordinates ...".format(thisProc,
                        aV["neighbSize"]))
            return
        #end timing
        microSecs = time.time() - microSecs

        # Record runtime for step
        aV["runtime"] += microSecs
        aV["speedProbe"] = int(aV["cntProbe"]/aV["runtime"])

        if aV["valueBest"] <= valueTarget:
            break

        if aV["cntProbe"] > cntProbeLmt:
            aV["isCensored"] = 1
            aV["speedProbe"] = int(aV["cntProbe"]/aV["runtime"])
            print ("WARNING from {}: isCensored=1, cntProbe={} > cntProbeLmt"
                    "={}\n".format(thisProc, aV["cntProbe"], aV["cntProbeLmt"]))
            return
        if aV["runtime"] > runtimeLmt:
            aV["isCensored"] = 1
            aV["speedProbe"] = int(aV["cntProbe"]/aV["runtime"])
            print ("WARNING from {}: isCensored=1, runtime={} > runtimeLmt"
                    "={}\n".format(thisProc, aV["runtime"], aV["runtimeLmt"]))
            return

    if aV["valueBest"] == aV["valueTarget"]:
        aV["targetReached"] = 1
    elif aV["valueBest"] < aV["valueTarget"]:
        aV["targetReached"] = 2
    else:
        aV["targetReached"] = 0

    return aV["targetReached"]

def pFile_read(fileName):
    thisProc = "pFile_read"
    ABOUT = """
This proc takes the path of an instance file associated with the 
sandbox P.lop, reads the file contents and returns parameter values and data 
structures expected by the combinatorial solver under this sandbox.
"""
    if fileName == "??":
        print ABOUT
        return
    elif fileName == "?":
        print "Valid query is '{} ??'".format(thisProc)
        return
        #raise Exception("Valid query is '{} ??'".format(thisProc))
        #sys.stderr.write("Valid query is '{}' ??'\n".format(thisProc))
        #sys.exit(1)
    rList = core.file_read(fileName).split("\n")

    nDim = int(rList[0])
    nItems = 0
    Ary = [[0 for i in range(nDim)] for j in range(nDim)]
    i = 0
    for line in rList[1:]:
        line = line.split()
        if len(line) > 0:
            for j in range(nDim):
                val = int(line[j])
                Ary[i][j] = val
                if val != 0:
                    nItems += 1
            i += 1

    varList = []
    for i in range(nDim):
        varList.append(i+1)

    density = float(nItems)/float(nDim * nDim - nDim)
    return (nDim, Ary, varList, density)

def pFile_write(coordPerm = [5, 3, 4, 2, 1], fileName = "../xBenchm/lop/tiny/i-5-book.lop"): 
    thisProc = "pFile_write"
    ABOUT = """
This proc takes a permutation coordinate and the path of an instance file 
associated with the sandbox P.lop, reads the file contents and returns 
an isomorph instance file created by permuting rows and columns of the
original instance file.
"""
    if coordPerm == "??":
        print ABOUT
        return
    elif coordPerm == "?":
        print "Valid query is '{} ??'".format(thisProc)
        return

    rList = core.file_read(fileName).split("\n")

    nDim = int(rList[0])
    Ary = [[0 for i in range(nDim + 1)] for j in range(nDim + 1)]
    i = 0
    for line in rList[1:]:
        if line != "":
            i += 1
            for j in range(nDim):
                Ary[i][j + 1] = line.split()[j]
    
    pAry = [[0 for i in range(nDim + 1)] for j in range(nDim + 1)]
    
    for i in range(1, nDim + 1):
        iP = coordPerm[i - 1]
        for j in range(1, nDim + 1):
            jP = coordPerm[j - 1]
            pAry[i][j] = Ary[iP][jP]
                    
    rowLines = str(nDim) + "\n"
    for i in range(1, nDim + 1):
        row = []
        for j in range(1, nDim + 1):
            jP = coordPerm[j - 1]
            row.append(str(pAry[i][j]))
            pAry[i][j] = Ary[iP][jP]
        rowLines += " " + str(row).translate(None, "[,\']") + "\n"
    filePerm = fileName.split(".lop")[0] + "-" + str(coordPerm).translate(None, "[ ]") + ".lop"
    file_write(filePerm, rowLines)
    print ".. created file " + filePerm

def f(coord = (1, 2, 3, 4, 5)):
    thisProc = "P.lop.f"
    ABOUT = """
This procedure takes a permutation coordinate (passed as an argument) and  
the parameters from the instance file associated with the sandbox P.lop 
(passed in a global array aStruc). It computes and returns the instance 
function value, given this coordinate.
"""
    if coord == "??":
        print ABOUT
        return
    elif coord == "?":
        print "Valid query is '{} ??'".format(thisProc)
        return

    #instance global variables
    global aStruc

    L = len(coord) 
    Lm1 = L - 1
    sumTot = 0
    for i in range(Lm1):
        iP = coord[i]
        ip1 = i + 1
        sumRow = 0
        for j in range(ip1, L):
            jP = coord[j]
            sumRow -= aStruc[iP - 1][jP - 1]
        sumTot += sumRow
    return sumTot

def exhA(instanceFile = "../xBenchm/lop/tiny/i-4-test1.lop"):
    thisProc = "exhA"
    ABOUT = """
This proc reads an instance file associated with the sandbox P.lop and a file 
of all coordinates (here, permutations) induced by the dimension defined by 
this instance. It then performs an exhaustive evaluation for function values 
defined by the instance and returns the list of coordinate:value pairs as
solutions with the function value minima. A rank value (in context of the 
underlying Hasse graph) is also associated with each coordinate. For coordinate
lengths of size <= 5, the procedure returns the exhaustive solution set
and a data structure that can be passed on to a follow-up tcl procedure which
creates a file of vertices and a file of edges annotated with x-y coordinates
for plotting of Hasse graphs under R.
"""
    if instanceFile == "??":
        print ABOUT
        return
    elif instanceFile == "?":
        print "Valid query is '{} ??'".format(thisProc)
        return

    global all_info
    global all_valu
    global aV
    global aStruc

    hasseAry = {}
    print core.file_read(instanceFile)
    aInstance = pFile_read(instanceFile)
    aPI = {"nDim": aInstance[0], "varList": aInstance[2], "density": aInstance[3]}
    aStruc = aInstance[1]
    nDim = aPI["nDim"]
    coordRef = [i for i in range(nDim + 1)]
    permFile = "../../../xBed/xBenchm/perm/perm.0" + str(aPI["nDim"]) + ".txt"
    coordList = core.file_read(permFile)
    print ".. " + str(len(coordList.split("\n")) - 1) + " permutations read from file " + permFile + "\n"
    valueMin = int(1e30)
    list = coordList.split("\n")
    bestAry = ""
    bestRank = 0
    for coord in list[:]:
        if len(coord) > 0:
            coord = map(int, coord.split(","))
            value = f(coord)
            rank = P.coord.rank(coord)
            if not hasseAry.has_key(rank):
                hasseAry[rank] = []
            hasseAry[rank].append(str(coord) + ":" + str(value))
            if value < valueMin:
                valueMin = value
                bestAry = str(coord) +  ":" + str(value)
                bestRank = rank
    rankMax = len(hasseAry) - 1
    print "rank\tsize\tcoord_value_pairs"
    widthMax = 0
    for i in range(rankMax):
        hasseAry[i].sort()
        print str(i) + "\t" + str(len(hasseAry[i])) + "\t" + str(hasseAry[i])
    print "\nvalueBest=-" + str(valueMin)
    if len(bestAry) > 0:
        print "bestAry(" + str(bestRank) + ")\t= " + str(bestAry)

    return [aPI["nDim"], rankMax, widthMax, coordList, bestAry, hasseAry]

def exhB( instanceFile = "../xBenchm/lop/tiny/i-4-test1.lop" ):
    global aStruc
    thisProc = "exhB"
    ABOUT = """
This proc takes an instance file (instanceDef) under the sanbox P.lop, 
defined on the space of permutation coordinates of length L. The procedure 
iteratively generates (L!) permutation coordinates to perform an exhaustive
evaluation of the function specified by the instance. The principle behind 
this coordinate generation is repeated re-use of the associative array 
aCoordHash0. Given this array, the generation proceeds from all coordinates 
at rank k to all coordinates at rank k+1. The value of k is in the range 
[0, L*(L-1)/2]. The exhaustive evaluation includes comprehensive  
instrumentation to measure the computational cost and the efficiency of 
the procedure.

For a stdout query, use one of these these commands:
          P.lop.exhB  ??  (sourced under tclsh)
  ../xBin/P.lop.exhBT     (executable under bash) 

"""
    if instanceFile == "??":
        print ABOUT
        return
    elif instanceFile == "?":
        print "Valid query is '{} ??'".format(thisProc)
        return

    # Read the instance
    print core.file_read(instanceFile)
    aInstance = pFile_read( instanceFile )
    aPI = { "nDim" : aInstance[0], "varList" : aInstance[2], "density" : aInstance[3] }
    aStruc = aInstance[1]

    # Setup for permutation coordinates
    L = aPI["nDim"]
    Lm1 = L - 1
    rankMax = L * ( L - 1 ) / 2

    # Define coordinate:value pair with rank=0
    coordRef = aPI["varList"]
    coordList0 = [coordRef]
    value = f(coordRef)
    valueBest = 1e30
    bestAry = {value: ["000_"+",".join(imap(str,coordRef))+":"+str(value)]}
    coordDistrib = {"000":1}
    sizeTot = 1

    # For each rank, create aCoordHash and aggregate coordList1, then probe
    # P.lop.f for function value
    aCoordHash = {}
    runtimeCoord = 0.0
    runtimeProbe = 0.0
    if ( L <= 5 ):
        hasseVerticies = {(0,1): [",".join(imap(str,coordRef)) + ":" + str(value)] }

    sizeRank = 0

    Lm1Range = xrange(Lm1)
    for rank in range(1, rankMax+1):
        coordList1 = []
        # Timing
        microSecs = time.time() 
        for coord in coordList0:
            elm_i = coord[0]
            for i in Lm1Range:
                ip1 = i + 1
                swapL = coord[:]
                elm_ip1 = coord[ip1]
                swapL[i] = elm_ip1
                if ip1 <= Lm1:
                    swapL[ip1] = elm_i
                    coordAdj = swapL[:]
                    elm_i = coord[ip1]

                inversion = P.coord.rank( coordAdj )
                if inversion == rank and not (",".join(imap(str,coordAdj)) in aCoordHash):
                    aCoordHash[",".join(imap(str,coordAdj))] = []
                    coordList1.append(coordAdj)
                    sizeTot += 1
                    sizeRank += 1
        #end timing
        microSecs = time.time() - microSecs
        runtimeCoord += microSecs 

        # Timing
        microSecs = time.time() 
        for coord in coordList1:
            value = f( coord )
            solutionString = ",".join(imap(str,coord))+":"+str(value)
            if L <= 5:
                if (rank,sizeRank) in hasseVerticies:
                    hasseVerticies[(rank,sizeRank)].append(solutionString)
                else:
                    hasseVerticies[(rank,sizeRank)] = [solutionString]
            if value < valueBest:
                coordString = "{:03}_".format(rank)+solutionString
                if bestAry.has_key(value):
                    bestAry[value].append(coordString)
                else:
                    bestAry[value] = [coordString]
        # end timing
        microSecs = time.time() - microSecs
        runtimeProbe += microSecs

        coordDistrib["{:03}".format(rank)] = sizeRank
        coordList0 = coordList1[:]
        aCoordHash = {}
        coordList1 = []
        sizeRank = 0

    print   
    if L <= 5:
        for key in hasseVerticies:
            print "hasseVerticies" + str(key) + " =", hasseVerticies[key]
    
    valueMin = min(bestAry)
    print "\nThere are {} valueBest = {} solutions:\nrank\tsolution".format(len(bestAry[valueMin]), valueMin)

    for item in bestAry[valueMin]:
        item = item.split("_")
        rank = item[0]
        solution = item[1]
        print rank+"\t"+solution

    # Out put file stuff
    print "\n".join([
        "\n",
        "instanceFile = {}".format(instanceFile),
        "There are {} valueBest = {} solutions".format(len(bestAry[valueMin]),
            valueMin),
        "rankMax = {}".format(rankMax),
        "coordLength = {}".format(L),
        "coordTotal = {}".format(sizeTot),
        "runtimeCoord = {:6.5}".format(runtimeCoord),
        "runtimeProbe = {:6.5}".format(runtimeProbe),
        "runtimeRatio = {:6.5}".format(runtimeCoord/runtimeProbe),
        "hostID = {}@{}-{}-{}".format(pwd.getpwuid(os.getuid())[0],
            os.uname()[1], platform.system(), os.uname()[2]),
        "dateLine = {}".format(time.strftime("%a %b %H:%M:%S %Z %Y")),
        "thisProc = {}\n".format(thisProc)])

    for key in sorted(coordDistrib):
        print "coordDistrib({})".format(key)+" =",coordDistrib[key]

def init( instanceDef, args = [] ):
    thisProc = "P.lop.init"
    ABOUT = """
Proc {} takes a variable number of arguments: 'instanceDef' as the
required argument and a number of optional arguments (in any order).
Proc $thisProc initializes all global data structures expected by the solver and 
also explicitly returns values of these three variables:

  targetReached coordInit valueInit  
       
Proc {} is invoked by P.lop.main, for command-line examples, querry

  P.lop.main ?? and/or P.lop.info 1
""".format(thisProc,thisProc)

    if instanceDef == "??":
        print ABOUT
        return
    elif instanceDef == "?":
        print "Valid query is '{} ??'".format(thisProc)
        return

    # info global variables
    global all_info
    global all_valu
    global aV
    # instance global variables 
    global aStruc
    # solver global variables
    global aCoordHash0
    global aWalkProbed
 
    argsOptions = args
    # TEMPORARY UNTIL I FIGURE OUT GLOBAL ISSUES
    tempInfoGlobals()

    # (0A) Phase 0A: initialize global variables
    aV = {}
    aStruc = {}
    aCoordHash0 = {}
    aWalkProbed = {}

    print "\n".join([
        "# ** from: {}:".format(thisProc),
        "# instanceDef={}".format(instanceDef),
        "# argsOptions={}".format(argsOptions)
        ])

    # (0B) Phase 0B: extract variable groups from all_valu
    namesRequired = []
    namesInternal = []
    namesOptionalBool = []
    namesOptional = []
    for name in all_valu.keys():
    
        val = all_valu[name]
    
        if val == "required":
            namesRequired.append(name)
        elif val == "internal":
            namesInternal.append(name)
            aV[name] = None
        elif val == "FALSE":
            namesOptionalBool.append(name)
            aV[name] = False 
        elif val == "TRUE":
            namesOptionalBool.append(name)
            aV[name] = True
        else:
            namesOptional.append(name)
            try:
                aV[name] = float(val)
            except:
                aV[name] = val
  
    # (1) Phase 1: initialize required commandline variables
    aV["instanceDef"] = instanceDef
    
    # Timing
    microSecs = time.time() 
    rList = pFile_read(aV["instanceDef"])
    aV["nDim"] = rList[0]
    aStruc = rList[1]
    aV["varList"] = rList[2]
    aV["density"] = rList[3]
    aV["coordRef"] = aV["varList"]
    aV["instanceID"] = os.path.basename(aV["instanceDef"])[:-4] 

    infoSolutionsDir = os.path.dirname(os.path.realpath(aV["instanceDef"]))
    infoSolutionsFile = all_info["sandboxName"] + ".info_solutions.txt"
    infoSolutionsFile = "/".join([all_info["sandboxPath"], "xBenchm", "lop",
        infoSolutionsFile])
    aV["infoSolutionsFile"] = infoSolutionsFile

    if not os.path.isfile(infoSolutionsFile):
        print "\nERROR from {}:\nfile {} is missing!\n".format(thisProc, infoSolutionsFile)
        return
    rList = core.file_read(infoSolutionsFile).split('\n')
    rList.pop()
    rListTmp = []
    for line in rList:
        firstChar = line[0]
        if firstChar != "#" and len(line) > 0:
            rListTmp.append(line)
    isFound = False 
    for line in rListTmp:
        line = line.split()
        varName = line[0]
        if varName == aV["instanceID"]:
            aV["valueTarget"] = int(line[1])
            try:
                aV["isProven"] = int(line[2].strip('-'))
            except:
                aV["isProven"] = False
            isFound = True
        if isFound:
            break
    if not isFound:
        print ("\nERROR from {}:"
                "\n .. instance {} was not found in file"
                "\n     {}\n".format(thisProc, aV["instanceID"], infoSolutionsFile))
        return
    #end timing
    microSecs = time.time() - microSecs

    aV["runtimeRead"] = microSecs
    aV["infoSolutionsFile"] = infoSolutionsFile
    aV["commandName"] = all_info["sandboxName"] + ".main"
    aV["commandLine"] = "{} {} {}".format(aV["commandName"] , instanceDef , argsOptions)
    aV["valueTarget"] = aV["valueTarget"] * (1 - aV["valueTol"])

    # (2A) Phase 2A: initialize optional command line variables
    if len(argsOptions) > 0:
        tmpList = argsOptions
        while len(tmpList) > 0:
            name = tmpList[0].strip("-")
            if name in namesOptional:
                aV[name] = tmpList[1]
                tmpList = tmpList[2:]
            elif name in namesOptionalBool:
                aV[name] = True
                tmpList = tmpList[1:]
            elif not name:
                print ("\nERROR from {}:"
                        "\n.. option name {} is not either of two lists below:"
                        "\n{}"
                        "\n\nor\n"
                        "\n{}".format(thisProc, name, namesOptional, namesOptionalBool))
                return
            else:
                aV[name] = int(tmpList[1])
                tmpList = tmpList[2:]

    # (2B) Phase 2B: continue to initialize the optional command line variables
    # aV["seedInit"] needs to be initialized first
    if aV["seedInit"] == "NA":
        # initialize RNG with random seed
        aV["seedInit" ] = 1e9 * random.random()
        random.seed(aV["seedInit"])
    else:
        # initialize RNG with user provided seed
        try:
            aV["seedInit"] = int(aV["seedInit"])
            random.seed(aV["seedInit"])
        except:
            print ("ERROR from {}:"
                    ".. only -seedInit NA or -ssedInit <int> are valid assignments,"
                    "not -seedInit {}\n".format(thisProc, aV['seedInit']))
            return
   
    # initialize permutation coordinate
    if aV["coordInit"] == "NA":
        # generate a random permutation coordinate
        aV['coordInit'] = P.coord.rand(aV["nDim"])
        aV['rankInit'] = P.coord.rank(aV["coordInit"])
    else:
        # check if user provided coordInit is the valid length
        aV["coordInit"] = [int(c) for c in aV["coordInit"].split(",")]
        if len(aV["coordInit"]) != aV["nDim"]:
            print ("\nERROR from {}:"
                    "\nThe permutation coordinate is of length {},"
                    "not the expected length {}\n".format(thisProc, aV["coordinit"], aV["nDim"]))
            return
        aV["rankInit"] = P.coord.rank(aV["coordInit"])
        """
    if aV["walkIntervalLmt"] == "NA" and aV["walkIntervalCoef"] != "NA":
        try:
            walkIntervalCoef = float(aV["walkIntervalCoef"])
            if walkIntevalCoef > 0.:
                aV["walkIntervalLmt"] = int(walkIntervalCoef * aV["nDim"])
                aV["walkIntervalCoef"] = walkIntervalCoef
        except:
            print ("\nERROR from {}:"
                    "The walkIntervalCoef can only be assigned a value of NA"
                    "or a positive number, not {} \n".format(thisProc, aV["walkIntervalCoef"]))
            return
    elif aV["walkIntervalLmt"] != "NA" and aV["walkIntervalCoef"] == "NA":
        try:
            walkIntervalLmt = int(aV["walkIntervalLmt"])
            if walkIntervalLmt > 0:
                aV["walkIntervalLmt"] = walkIntervalLmt
        except:
            print ("\nERROR from {}:"
                    "The walkIntervalLmt can only be assigned a value of NA"
                    "or a positive number, not {} \n".format(thisProc, aV["walkIntervalLmt"]))
            return
    elif aV["walkIntervalLmt"] != "NA" and aV["walkIntervalCoef"] != "NA":
        print ("ERROR from {}:"
                "The walkIntervalLmt and walkIntervalCoef can only be assigned"
                "pairwise values of\n(NA NA) (default) (NA double) or (integer NA)"
                "not ({} {})\n".format(thisProc, aV["walkIntervalLmt"], aV["walkIntervalCoef"]))
        return """

    if aV["walkSegmLmt"] == "NA" and aV["walkSegmCoef"] != "NA":
        try:
            walkSegmCoef = float(aV["walkSegmCoef"])
            if walkIntevalCoef > 0.:
                aV["walkSegmLmt"] = int(walkSegmCoef * aV["nDim"])
                aV["walkSegmCoef"] = walkSegmCoef
        except:
            print ("\nERROR from {}:"
                    "The walkSegmCoef can only be assigned a value of NA"
                    "or a positive number, not {} \n".format(thisProc, aV["walkSegmCoef"]))
            return
    elif aV["walkSegmLmt"] != "NA" and aV["walkSegmCoef"] == "NA":
        try:
            walkSegmLmt = int(aV["walkSegmLmt"])
            if walkSegmLmt > 0:
                aV["walkSegmLmt"] = walkSegmLmt
        except:
            print ("\nERROR from {}:"
                    "The walkSegmLmt can only be assigned a value of NA"
                    "or a positive number, not {} \n".format(thisProc, aV["walkSegmLmt"]))
            return
    elif aV["walkSegmLmt"] != "NA" and aV["walkSegmCoef"] != "NA":
        print ("ERROR from {}:"
                "The walkSegmLmt and walkSegmCoef can only be assigned"
                "pairwise values of\n(NA NA) (default) (NA double) or (integer NA)"
                "not ({} {})\n".format(thisProc, aV["walkSegmLmt"], aV["walkSegmCoef"]))
        return

    # (3A) Phase 3A: directly initialize the remaining internal variables
    aV["coordType"] = "P"
    aV["functionDomain"] = "P.lop"
    aV["functionID"] = "lop"

    # define solverID
    if aV["solverMethod"] == "ant" and aV["isSimple"]:
        aV["solverID"] = "ant_saw_simple"
    elif aV["solverMethod"] == "ant":
        aV["solverID"] = "ant_saw"
    elif aV["solverMethod"] == "bee":
        aV["solverID"] = "bee_saw"
    elif aV["solverMethod"] == "rw" and not aV["notSAW"]:
        aV["solverID"] = "rw_saw"
    elif aV["solverMethod"] == "rw" and aV["notSAW"]:
        aV["solverID"] = "rw"
    elif aV["solverMethod"] == "sa" and not aV["notSAW"]:
        aV["solverID"] = "sa_saw"
    elif aV["solverMethod"] == "sa" and aV["notSAW"]:
        aV["solverID"] = "sa"
    else:
        aV["solverID"] = saw

    # Time stamp (14-digit GMT formatting)
    aV["solverVersion"] = time.strftime("%Y %m %d %H:%M:%S")
    aV["timeStamp"] =  time.strftime("%Y%m%d%H%M%S")
    aV["dateLine"] = time.strftime("%a %b %d %H:%M:%S %Z %Y") 
    #@TODO HostID stuff (lines 583-486)
    aV["hostID"] = "{}@{}-{}-{}".format(pwd.getpwuid(os.getuid())[0],
            os.uname()[1], platform.system(), os.uname()[2])
    aV["compiler"] = "python-"+".".join(imap(str,sys.version_info[:3]))

    # find aV["valueInit"] by doing the first probe for function value
    # Timing
    microSecs = time.time() 
    aV["valueInit"] = f(aV["coordInit"])
    # end timing
    microSecs = time.time() - microSecs

    # initialize associated variables for initial probe
    aV["runtime"] = microSecs
    aV["cntProbe"] = 1
    aV["cntStep"] = 0
    aV["coordPivot"] = aV["coordInit"]
    aV["coordBest"] = aV["coordInit"]
    aV["valuePivot"] = aV["valueInit"]
    aV["valueBest"] = aV["valueInit"]
    aCoordHash0= {}
    aCoordHash0[tuple(aV["coordInit"])] = []

    # (4) Phase 4: check if valueTarget has been reached, return to main if > 0
    if aV["valueInit"] == aV["valueTarget"]:
        aV["targetReached"] = 1
    elif aV["valueInit"] < aV["valueTarget"]:
        aV["targetReached"] = 2
    else:
        aV["targetReached"] = 0

    # (5) Phase 5: complete initialization of aV before the first step
    #aHashWalk[tuple(aV["coordInit"])] = []
    aV["isCensored"] = 0
    aV["cntRestart"] = 0
    aV["walkLength"] = aV["cntStep"]
    """
    if aV["neighbDist"] == 1:
        aV["neighbSize"] = aV["nDim"] - 1
    else:
        aV["neighbSize"] = "dynamic"
    """

    # (6) Phase 6: initialize special arrays that can be selected w/ arguments
    #       from command line
    if aV["writeVar"] >= 4:
        aV["isSimple"] = 1
    aValueBest[aV["valueInit"]] = [0,0,aV["coordInit"]]
    aWalkBest[aV["valueInit"]] = [0,0,aV["coordInit"],0,0]

    aWalk[aV["cntStep"]] = "{} {} {} {} {} {}".format(aV["cntStep"], aV["cntRestart"],
            aV["coordPivot"], aV["valuePivot"], aV["neighbSize"], aV["cntProbe"])

    isPivot = 1
    aV["rankPivot"] = P.coord.rank(aV["coordPivot"])
    aWalkProbed[(aV["walkLength"],0)] = (aV["walkLength"], aV["cntRestart"], 
            aV["coordPivot"], aV["valuePivot"], aV["rankPivot"], isPivot, 
            aV["neighbSize"], aV["cntProbe"])

    # (7) Phase 7: verify that all variables under aV have been initialized in
    #       solverDomain table all_info["infoVariablesFile"]
    errorItems = []
    errorLines = []
    for name in aV.keys():
        if name not in all_valu:
            errorLines.append( ("{} -- this variable is missing from the solver"
                " domain table in the file {}\n".format(name, all_info["infoVariablesFile"])))
            errorItems.append(name)
    if len(errorItems) > 0:
        print "\nWarning from {}\n{}".format(thisProc, errorLines)
        print "Missing variables\n{}\n".format(errorItems)

    # (8) Phase 8: check whether coordInit caused targetReached to be > 1
    if aV["targetReached"] > 0:
        print "# BINGO: valueTarget has been reached or exceeded with coordInit"
        stdout(withWarning=1)
        return

    # (9) Phase 9: write to stdout based on writeVar
    if aV["writeVar"] == 1:
        print "\n** Final values of initialized variables (array aV) **"
        print aV
        print "\n** Values associated with instance array aStruc **"
        print aStruc
        print "\n** as reported on {}, returning".format("@TODO: TIME STAMP")
        print "targetReached\tvalueInit\tcoordInit"

    #result = "{} {} {}".format(aV["targetReached"], aV["valueInit"], aV["coordInit"])
    result = (aV["targetReached"], aV["valueInit"], aV["coordInit"])
    return result

def main( instanceDef, args=[] ):
    thisProc = "P.lop.main"
    ABOUT = ("Procedure {} takes a variable number of arguments: "
            "instanceDef as required argument and args (a reserved-name variable)."
            "It outputs either (1) a response to a command line query "
            "(see below) or a summary of initialized variables only, and (2) "
            "completes variable initialization and invokes the solver under "
            "the default or user-specified argument values.\n"
            "For a stdout of command line details, query with command shown "
            "below:\n\t\tP.lop.main ?\t(under python interpreter)\n"
            "\tor\n"
            "\t\t../xBin/P.lopP (under bash)".format(thisProc))

    if instanceDef == "??":
        print ABOUT
        return
    elif instanceDef == "?":
        print "Valid query is '{} ??'".format(thisProc)
        return
    
    #info global variables
    global all_info
    global all_valu
    global aV
    #instance global variables
    global aStruc
    #solver global variables
    global aHashTmp
    global aHashNeighb
    global aHashWalk
    global aWalkBest
    global aWalk
    global aWalkProbed

    # TEMPORARY UNTIL I FIGURE OUT GLOBAL ISSUES
    """
    thisDir = os.getcwd()
    sandboxPath = os.path.dirname(thisDir) 
    sandboxName = os.path.basename(sandboxPath)
    infoVariablesFile = sandboxName + ".info_variables.txt"
    infoVariablesFile = "/".join([sandboxPath,"xLib",infoVariablesFile])

    all_info = {}
    all_info["infoVariablesFile"] = infoVariablesFile
    all_info["sandboxName"] = sandboxName 
    all_info["sandboxPath"] = sandboxPath 
    """
    tempInfoGlobals()

    # (1) Phase 1: query about the commandLine **OR** 
    #       to initialize all variables
    if instanceDef == "?":
        # read solver domain table and return query about commandline
        info(1,all_info["infoVariablesFile"])
        return
    else:
        # variable initialization
        rList = init(instanceDef,args)

    if rList is None:
        return
    elif aV["isInitOnly"]:
        print ("\n{}"
            "\n.. completed initialization of all variables,"
            "\n   exiting the solver since option -{} has been asserted."
            "\n{}\n".format("-"*78,aV["isInitOnly"],"-"*78))
        print "targetReached coordInit valueInit = {}\n".format(rList)
        for key, value in aV.items():
            print "aV({}) = {}".format(key,value)
        print
        for i in range(len(aStruc)):
            print "aStruc[{}] = {}".format(i, aStruc[i])
        print
        for key, value in aCoordHash0.items():
            print "aCoordHash0({}) = {}".format(key,value)
        print
        for key, value in aWalkProbed.items():
            print "aWalkProbed({}) = {}".format(key,value)
        return
    else:
        print ("\n{}"
            "\n.. completed initialization of all variables,"
            "\n   proceeding with the search under solverID = {}"
            "\n{}\n".format("-"*78,aV["solverID"],"-"*78))
    
    # (2) Phase 2: proceed with the combinatorial search
    if aV["solverMethod"] == "saw":
        aV["solverID"] = saw
    else:
        print "\nERROR from {}:\nsolverMethod = {} is not implemented\n".format(thisProc, aV["solverMethod"])
        return
    print "#    Proceeding with the search under solverID = {}".format(aV["solverID"])
    print "#"*78
    aV["solverID"]()

    stdout()
    if aV["isWalkTables"]:
        print "TODO: walk.tables method"
    return
    
def stdout( withWarning = 1 ):
    ABOUT = ("This procedure outputs results afer a successful completion of "
            "a combinatorial solver. The output is directed to 'stdout' and "
            "includes a solution (a coordinate-value pair) and the observed "
            "performance values. The format consists of a few comment lines, "
            "followed by a tabbed name-value pairs. The first piar is always\n"
            "instanceDef <value>\n"
            "This procedure is universal under any function coordType=P!")
    thisProc = "P.lop.stdout"
    
    #info global vairables
    global all_info
    global all_value
    global aV

    print ("# \n# FROM {}: A SUMMARY OF NAME-VALUE PAIRS"
            "\n# commandLine = {}"
            "\n#    dateLine = {}"
            "\n#   timeStamp = {}"
            "\n#".format(thisProc, aV["commandLine"], aV["dateLine"], aV["timeStamp"]))

    stdoutNames = ("instanceDef", "solverID", "coordInit", "coordBest", "nDim",
            "walkLengthLmt", "walkLength", "cntRestartLmt", "cntRestart", 
            "cntProbeLmt", "cntProbe", "runtimeLmt", "runtime", "runtimeRead", 
            "speedProbe", "hostID", "isSimple", "solverMethod", "compiler",
            "walkSegmLmt", "walkSegmCoef", "seedInit", "valueInit", 
            "valueBest", "valueTarget", "valueTol", "targetReached", 
            "isCensored")

    for name in stdoutNames:
        if name in aV:
            print "{}\t\t{}".format(name, aV[name])
        else:
            if withWarning:
                print "# WARNING: no value exist for {}".format(name)

def info( isQuery=0, infoVariablesFile="../xLib/P.lop.info_variables.txt"):
    thisProc = "P.lop.info"
    ABOUT = """
This proc takes a variable 'isQuery' and the hard-wired path to file
infoVariablesFile *info_variables.txt.  

  if isQuery == 0    then $thisProc ONLY reads infoVariablesFile and 
                     initializes global arrays 'all_info' and 'all_valu'
              
  if isQuery == 1    then $thisProc initializes the global arrays
                     'all_info' and 'all_valu' and then outputs to stdout  
                     the complete information about the command line for 
                     P.lop.main. The information about the command line 
                     is auto-generated within $thisProc from the
                     tabulated data which is read from infoVariablesFile 
                     defined above.
                   
   if isQuery == ??  then $thisProc responds to stdout with information 
                     about all three case of valid arguments. 

            P.lop.main ??   (under tclsh)
       or
            ../xBin/P.lopT  (under bash)
"""
    if isQuery == "??":
        print ABOUT
        return
    elif isQuery == "?":
        print "Valid query is '{} ??'".format(thisProc)
        return

    #info global variables
    global all_info
    global all_valu
    global aV

    # read *.info_variables.txt for this solver domain
    rList = util.table_info_variables(infoVariablesFile) 
    all_info = rList[0]
    all_valu = rList[1]

    if not isQuery:
        return (all_info, all_valu)
    elif isQuery != 1:
        print ("\nERROR from {}: incorrect arguement value!"
                "\nFor more information, use the command P.lop.info ??\n")

    # Preferred order of optional commandLine arguements
    optInfoList = [ "runtimeLmt", "cntProbeLmt", "cntRestartLmt", 
            "walkLengthLmt", "seedInit", "coordInit", "valueInit", 
            "valueTarget", "valueTol", "walkSegmLmt",
            "walkSegmCoef", "isSimple", "writeVar" ]

    print "\n".join([
        "USAGE:\n",
        "under TkCon shell (which has sourced ../xLib/all_python.py",
        "\tP.lop.main instanceFile [optional_arguments]",
        "",
        "under bash, invoking the 'tcl executable P.lopT' which sources" + 
            " libraries directly",
        "\t../xBin/P.lopT instanceFile [optional_arguments]",
        "",
        "under bash, invoking the 'python executable P.lopP' which sources" + 
            " libraries directly",
        "\t../xBin/P.lopP instanceFile [optional_arguments]",
        "",
        "under bash, invoking the 'compiled C++ code as a binary P.lopX",
        "\t../xBin/P.lopX instanceFile [optional_arguments]",
        "",
        "EXAMPLES:",
        "\tP.lop.main ../xBenchm/lop/tiny/i-5-book2.lop -initOnly",
        "\tP.lop.main ../xBenchm/lop/tiny/i-5-book2.lop -valueTarget -51",
        "\t../xBin/P.lopT ../xBenchm/lop/tiny/i-5-book2.lop -valueTarget" + 
            " -51 -runtimeLmt 5 -seedInit 1789",
        "\t../xBin/P.lopP ../xBenchm/lop/tiny/i-5-book2.lop -valueTarget" +
            " -51 -coordInit 1,2,3,4,5 -runtimeLmt 5",
        "\t../xBin/P.lopX ../xBenchm/lop/tiny/i-5-book2.lop -valueTarget" +
            " -51 -coordInit 5,3,2,1,4 -seedinit 1914",
        "DESCRIPTION:",
        "P.lop.main, P.lopT, P.lopP, or P.lopX take one REQUIRED argument",
        "",
        "\tinstanceDef  (Here a filePath with an extension .lop " + 
            "or NO extension)",
        "",
        "and a number of OPTIONAL arguments in any order. The most" + 
            " significant parameter,",
        "extracted from the instanceDef is",
        "\tnDim ... coordinate size,",
        "\t\t\ti.e. the number of variables (columns in the square matrix)",
        "Here is a complete list of pairs 'name defaultValue', with short",
        "in-line descriptions:"
        ])

    # create nameList and valuelist
    for name in optInfoList:
        value = all_valu[name]
        nameList = []
        valueList = []
        nameList.append(name)
        valueList.append(value)
        # pad with blank spaces
        for i in range(len(all_info[name])):
            nameList.append(" "*len(name))
            valueList.append(" "*len(value))

        optArgData = [(nameList[i], valueList[i], all_info[name][i]) for i in range(len(all_info[name]))]
        for (nam, val, descr) in optArgData:
            cnt = 17 - len(nam)
            space1 = " "*cnt
            cnt = 12 - len(val)
            space2 = " "*cnt
            if len(nam) > 0 and nam.strip():
                # prefix with -
                nam1 = "-"+nam
            else:
                nam1 = " "+nam

            print "\t{}{}{}{}\t{}".format(nam1,space1,val,space2,descr)

    print "\n".join([
        "\nDETAILS:",
        "This solver reads an instance of the 'linear ordering problem in a " +
            "matrix format",
        "and returns a column/row ordering that minimizes the negative sum " +
            "of matrix",
        "elements above the diagonal. The example below shows an instance " +
            "of such a matrix",
        "with sum = -8 under its 'natural order', and an instance under an " +
            "optimal",
        "permutation of 3,1,4,2 with a sum of -13. For this matrix, there " +
            "are two more",
        "such optimal permutations: 2,3,1,4 and 4,2,3,1.",
        "",
        "natural order   under permutation",
        "  1,2,3,4          3,1,4,2",
        "  sum = -8         sum = -13",
        "------------    ------------",
        "4               4           ",
        "  0 0 0 5         0 4 1 1",
        "  1 0 2 0         0 0 5 0",
        "  4 1 0 1         1 3 0 2",
        "  3 2 1 0         2 1 0 0"
        ])

    return

def tempInfoGlobals():
    global all_info
    global all_valu
    # TEMPORARY UNTIL I FIGURE OUT GLOBAL ISSUES
    thisDir = os.getcwd()
    sandboxPath = os.path.dirname(thisDir) 
    sandboxName = os.path.basename(sandboxPath)
    infoVariablesFile = sandboxName + ".info_variables.txt"
    infoVariablesFile = "/".join([sandboxPath,"xLib",infoVariablesFile])

    all_info = {}
    all_info["infoVariablesFile"] = infoVariablesFile
    rList = info(0, all_info["infoVariablesFile"])

    all_info = rList[0]
    all_valu = rList[1]
    all_info["sandboxName"] = sandboxName 
    all_info["sandboxPath"] = sandboxPath 
    all_info["infoVariablesFile"] = infoVariablesFile

if __name__ == "__main__":
    #pFile_write()
   
    #print f()
    #print "========\n"
    info()
    #main("../xBenchm/lop/tiny/i-5-book.lop","-writeVar 3")
    main()
    saw_pivot_simple()
    #exhB( "../xBenchm/lop/tiny/i-7-test1.lop")
    #print "========\n"
    #exhB( "../xBenchm/lop/tiny/i-8-test1.lop")
    #print "========\n"
    #exhB( "../xBenchm/lop/tiny/i-9-test1.lop")
    #print "========\n"
    #exhB( "../xBenchm/lop/tiny/i-10-test1.lop")
    #main(sys.argv[1],sys.argv[2:])
