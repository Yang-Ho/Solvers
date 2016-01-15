#! /usr/bin/env python2
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
    if coordPiv == "?":
        #Error
        print "Valid query is '{} ??'".format(thisProc)
        return

    #info global variables
    global all_info
    global all_valu
    global aV
    #instance global variables
    global aStruc
    #solver globla variables
    global aHashTmp
    global aHashNeighb
    global aWalk 
    global aHashWalk 
    global aWalkProbed 
    global aWalkBest 
    global aValueBest 
    global aAdjacent

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
        coordA = tuple(coordAdj[ip1])
        if coordA not in aHashTmp:
            neighbSize += 1
            valueA = f(coordA)
            aV["cntProbe"] += 1
            if aV["writeVar"] == 6:
                coordProbedList.append(coordA)
                valueProbedList.append(valueA)
            if valueA <= valueBest:
                if valueA < valueBest:
                    coordBestList = []
                valueBest = valueA
                coordBest = coordA
                aAdjacent[valueA] = coordA
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
    coordBest = coordBestList[idx]
    return (coordBest, valueBest, neighbSize, coordProbedList, valueProbedList)

def saw():
    ABOUT = (
        "This procedure takes values initialized under 'global arrays:' by the"
        " proedure P.lop.init and constructs a segment of a self-avoiding walk"
        " (SAW). This procedure invokes either P.lop.saw.pivot.simple or the "
        "significantly more efficient procedure P.lop.saw.pivot.ant" 
        )
    return

def pFile_read(fileName):
    thisProc = "pFile_read"
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
    rList = file_read(fileName).split("\n")

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

def f(coord = [5, 3, 4, 2, 1]):
    global aStruc
    thisProc = "f"
    nDim = len(coord)
    L = nDim
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
    hasseAry = {}
    print file_read(instanceFile)
    aInstance = pFile_read(instanceFile)
    aPI = {"nDim": aInstance[0], "varList": aInstance[2], "density": aInstance[3]}
    aStruc = aInstance[1]
    nDim = aPI["nDim"]
    coordRef = [i for i in range(nDim + 1)]
    permFile = "..//xBenchm/perm/perm.0" + str(aPI["nDim"]) + ".txt"
    coordList = file_read(permFile)
    print ".. " + str(len(coordList.split("\n")) - 1) + " permutations read from file " + permFile + "\n"
    valueMin = int(1e30)
    list = coordList.split("\n")
    bestAry = ""
    bestRank = 0
    for coord in list[:]:
        if len(coord) > 0:
            coord = map(int, coord.split(","))
            value = f(coord, aStruc)
            rank = P_coord.inversion(coord)
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

    # Read the instance
    print file_read(instanceFile)
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
    # P_lop.f for function value
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

                inversion = P_coord.inversion( coordAdj )
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

def initOld( instanceDef, argsOptions ):
    thisProc = "p_lop.{}".format("init")

    global aPI
    global reqAry
    global optAry
    global optAryNames
    global optAryNamesBool
    global optAryList
   
    aPI = {} 

    #if len(argsoptions) == 1:
        #argsoptions = argsoptions[0]
    
    aPI["coordType"] = "P"
    aPI["functionBase"] = "lop"
    aPI["commandName"] = aPI["coordType"]+"_"+aPI["functionBase"]+".main"
    aPI["solverId"] = aPI["commandName"]

    numAllConstants = 4

    aPI["instanceDef"] = instanceDef
    
    for item in optAry:
        name = item[0]
        valDefault = item[1]
        aPI[name] = valDefault

    print argsOptions

    if len(argsOptions) > 0:
        tmpList = argsOptions
        #print optAryNames
        while len(tmpList) > 0:
            name = tmpList[0].strip("-")
            if name not in optAryNames:
                #error
                print "\n".join([
                "\nERROR from {}:".format(thisProc),
                ".. option name {} not in this list".format(name),
                "\n.. {}".format(optAryNames.sort())])
            if name not in optAryNamesBool:
                aPI[name] = True
                del tmpList[0]
            else:
                aPI[name] = tmpList[1]
                del tmpList[0]
                del tmpList[0]


    nCheck = numAllConstants + len(reqAry) + len(optAry)
    if len(aPI) != nCheck:
        #Error
        print "\n".join([
        "\nERROR from {}:".format(thisProc),
        "\nsize of reqAry+optAry ({}) is not matched to".format(nCheck),
        "\nsize of aPI ({})".format(len(aPI))])

    aPI["commandLine"] = "{} {} {}".format(aPI["commandName"], instanceDef, argsOptions)
    aPI["commandLineOptionNames"] = optAryList
    if abs(int(aPI["writeVar"])) == 1:
        print "\n".join([
        "\n** TRACE FROM {} **".format(thisProc),
        "** ALL commandLine variable names ** with either default or \
        user-assigned values"])
        print aPI

    microSecs = time.time()
    aPI["instanceFile"] = aPI["instanceDef"]

    aInstance = pFile_read(aPI["instanceFile"])
    aPI["nDim"] = aInstance[0]
    aStruc = aInstance[1]
    aPI["varList"] = aInstance[2]
    aPI["density"] = aInstance[3]
    aPI["coordRef"] = aPI["varList"]
    microSecs = time.time() - microSecs  
    aV["runtimeRead"] = microSecs

    if aPI["seedInit"] == "NA":
        seedInit = 1e9 * random.random()
        random.seed(seedInit)
        init = random.random()
        aPI["seedInit"] = seedInit
    elif isinstance(aPI["seedInit"], (long, int)):
        random.seed(seedInit)
        init = random.random()
    else:
        #Error
        print "\n".join([
            "ERROR from {}:".format(thisProc),
            ".. only -seedInit NA or -seedInit <int> are valid assignments, not \
            -seedInit {}\ ".format(aPI["seedInit"])])

    #ABOUT Coord Init
    if aPI["coordInit"] == "NA":
        aPI["coordInit"] = P_coord.rand(aPI["nDim"])
        aPI["rankInit"] = P_coord.inversion(aPI["coordInit"])
    else:
        coordTmp = aPI["coordInit"]
        coordTmp.sort()
        if len(coordTmp) != aPI["nDim"]:
            #Error
            print "\n".join([
                "\nERROR from {}:".format(thisProc),
                "The permutation coordinate is of length {},".format(len(coordTmp)),
                "not the expected length {}\n".format(aPI["nDim"])])
        aPI["rankInit"] = P_coord.inversion(aPI["coordInit"])

    # ABOUT walkSegmLmt and walkSegmCoef
    aPIwalkSegmLmt = aPI["walkSegmLmt"]
    aPIwalkSegmCoef = aPI["walkSegmLmt"]
    if aPIwalkSegmLmt == "NA" and aPIwalkSegmLmt == "NA":
        pass
    elif aPIwalkSegmLmt == "NA":
        try:
            aPIwalkSegmCoef = float(aPIwalkSegmCoef)
            if aPIwalkSegmCoef > 0:
                aPI["walkSegmLmt"] = int(aPIwalkSegmCoef * aPI["nDim"])
            else:
                #Error
                print "\n".join([
                    "\nERROR from {}:".format(thisProc),
                    "The walkSegmCoef can only be assigned a value of NA or a\
                    positive number, not {}".format(aPI["walkSegmCeof"])])
        except ValueError:
            #Error
            print "\n".join([
                "\nERROR from {}:".format(thisProc),
                "The walkSegmCoef can only be assigned a value of NA or a\
                positive number, not {}".format(aPI["walkSegmCeof"])])
    elif aPIwalkSegmCoef == "NA": 
        try:
            aPIwalkSegmLmt = float(aPIWalkSegmLmt)
            if aPIwalkSegmLmt > 0:
                aPI["walkSegmLmt"] = int(aPIWalkSegmLmt)
            else:
                #Error
                print "\n".join([
                    "\nERROR from {}:".format(thisProc),
                    "The walkSegmLmt can only be assigned a value of NA or a\
                    positive number, not {}".format(aPI["walkSegmLmt"])])
        except ValueError:
            #Error
            print "\n".join([
                "\nERROR from {}:".format(thisProc),
                "The walkSegmLmt can only be assigned a value of NA or a\
                positive number, not {}".format(aPI["walkSegmLmt"])])
    else:
        #Error
        print "\n".join([
            "\nERROR from {}:".format(thisProc),
            "The walkSegmLmt and walkSegmCoef can only be assigned pairwise \
            values of",
            "(NA NA) (default), (NA double), or (integer, NA)"])

    # ABOUT walkIntervalLmt and walkIntervalCoef
    aPIwalkIntervalLmt = aPI["walkIntervalLmt"]
    aPIwalkIntervalCoef = aPI["walkIntervalCoef"]
    if aPIwalkIntervalLmt == "NA" and aPIwalkIntervalCoef == "NA":
        pass
    elif aPIwalkIntervalLmt == "NA": 
        try:
            aPIwalkIntervalCoef = float(aPIWalkIntervalCoef)
            if aPIwalkIntervalCoef > 0:
                aPI["walkIntervalLmt"] = int(aPI["walkIntervalCoef"] * aPI["nDim"])
            else:
                #Error
                print "\n".join([
                    "\nERROR from {}:".format(thisProc),
                    "The walkIntervalCoef can only be assigned a value of NA or a\
                    positive number > 0, not {}".format(aPIwalkIntervalCoef)])
        except ValueError:
            #Error
            print "\n".join([
                "\nERROR from {}:".format(thisProc),
                "The walkIntervalCoef can only be assigned a value of NA or a\
                positive number > 0, not {}".format(aPIwalkIntervalCoef)])
    elif aPIwalkSegmCoef == "NA":
        try:
            aPIwalkIntervalLmt = float(aPIWalkIntervalLmt)
            if aPIwalkIntervalLmt > 0:
                aPI["walkIntervalLmt"] = int(aPIWalkIntervalLmt)
            else:
                #Error
                print "\n".join([
                    "\nERROR from {}:".format(thisProc),
                    "The walkIntervalLmt can only be assigned a value of NA or a\
                    positive number > 0, not {}".format(aPIwalkIntervalLmt)])
        except ValueError:
            #Error
            print "\n".join([
                "\nERROR from {}:".format(thisProc),
                "The walkIntervalLmt can only be assigned a value of NA or a\
                positive number > 0, not {}".format(aPIwalkIntervalLmt)])
    else:
        #Error
        print "\n".join([
            "ERROR from {}:".format(thisProc),
            "The walkIntervalLmt and walkIntervalCoef can only be assigned pairwise values of",
            "(NA NA) (default), (NA double), or (integer, NA)"])

    aPI["valueMin"] = aPI["valueTarget"]
    aPI["valueTarget"] = aPI["valueMin"] * ( 1 + aPI["valueTol"] )
    aV["valueTarget"] = aPI["valueTarget"]
    aPI["coordType"] = "P"
    aPI["functionBase"] = "lop"
    aPI["dateLine"] = time.strftime("%a %b %H:%M:%S %Z %Y")
    aPI["timeStamp"] = time.strftime("%Y %m %d %H %M %S")
    aPI["hostID"] = "{}@{}-{}-{}".format(pwd.getpwuid(os.getuid())[0], os.uname()[1], platform.system(), os.uname()[2])

    if aPI["neighbDist"] == 1:
        aPI["walkMethod"] = "ant"
    else:
        aPI["walkMethod"] = "bee"

    microSecs = time.time()
    rList = [f(aPI["coordInit"],aStruc)]
    microSecs = time.time() - microSecs 
    
    aV["runtime"] = microSecs
    aV["coordInit"] = aPI["coordInit"]
    aV["valueInit"] = rList[0]
    aV["coordPivot"] = aV["coordInit"]
    aV["valuePivot"] = aV["valueInit"]
    
    if aV["valueInit"] == aV["valueTarget"]:
        aV["targetReached"] = 1
    elif aV["valueInit"] < aV["valueTarget"]:
        aV["targetReached"] = 2
    else:
        aV["targetReached"] = 0
    
    if aV["targetReached"] > 0:
        aV["coordInit"] = aV["coordInit"]
        print "\n".join([
            "\# BINGO, targetReached = {} for seed = {},".format(aV["targetReached"], aPI["seedInit"]),
            "\# coordInit = {}, valueInit = {} ".format(aV["coordinit"], aV["valueInit"])])
        aV["cntProbe"] = 1
        aV["walkLength"] = 0
        aV["isCensored"] = 0
        aV["speed"] = "NA"
        aV["restart"] = 0
        aV["restartUniq"] = 0
        return "{} {} {}".format(aV["coordInit"], aV["valueInit"], aV["targetReached"])
    
    aV["cntProbe"] = 1
    aCoordHash0 = {",".join(imap(str,aV["coordInit"])):[]}
    aCoordHash1 = {",".join(imap(str,aV["coordInit"])):[]}
    aV["isCensored"] = 0
    aV["isBlocked"] = 0
    aV["walkLength"] = 0
    aV["cntRestart"] = 0
    aV["walkRepeats"] = 0
    aV["walkInterval"] = 0
    aValueBest[aV["valueInit"]] = (0,0,aV["coordInit"])
    if aPI["neighbDist"] == 1:
        aV["neighbSize"] = aPI["nDim"] - 1
    else:
        aV["neighbSize"] = "dynamic"
    
    aWalkBest[aV["valueInit"]] = (0,0,aV["coordInit"],0,0)
    
    aWalk[aV["walkLength"]] = "{} {} {} {} {} {} {}".format(aV["walkLength"], aV["cntRestart"],aV["coordPivot"], aV["valuePivot"], aV["neighbSize"], aV["cntProbe"], aV["targetReached"])

    aWalkProbed[(aV["walkLength"],0)] = "{} {} {} {} {} {} {}".format(aV["walkLength"], aV["cntRestart"],aV["coordPivot"], aV["valuePivot"], aV["neighbSize"], aV["cntProbe"], 1)
    
    if aPI["writeVar"] == 1 or aPI["writeVar"] == -1:
        print "\n".join([
            "** Final values of initialized primary input variables (array aPI) **",
            "{}".format(aPI),
            "** Final values of initialized auxiliary variables (array aV) **",
            "{}".format(aV),
            "** as reported on {}, returning".format(aPI["dateLine"]),
            "coordInit\tvalueInit\ttargetReached"])
    
    return (aV["coordInit"], aV["valueInit"], aV["targetReached"])

def init( instanceDef, args ):
    thisProc = "P.lop.init"

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
    global aWalk
    global aWalkProbed
    global aWalkBest
    global aValueBest
    global aAdjacent
 
    argsOptions = args

    # TEMPORARY UNTIL I FIGURE OUT GLOBAL ISSUES
    thisDir = os.getcwd()
    sandboxPath = os.path.dirname(thisDir) 
    sandboxName = os.path.basename(sandboxPath)
    infoVariablesFile = sandboxName + ".info_variables.txt"
    infoVariablesFile = "/".join([sandboxPath,"xLib",infoVariablesFile])

    all_info = {}
    all_info["infoVariablesFile"] = infoVariablesFile
    rList = info(all_info["infoVariablesFile"], 0)

    all_info = rList[0]
    all_valu = rList[1]
    all_info["sandboxName"] = sandboxName 
    all_info["sandboxPath"] = sandboxPath 
    all_info["infoVariablesFile"] = infoVariablesFile

    #@TODO dynamically set all these dictionaries (see line 147 - 149 in tcl)
    aV = {}
    aStruc = {}
    aHashTmp = {}
    aHashNeighb = {}
    aHashWalk = {}
    aWalk = {}
    aWalkProbed = {}
    aWalkBest = {}
    aValueBest = {}
    aAdjacent = {}

    print "\n".join([
        "# ** from: {}:".format(thisProc),
        "# instanceDef={}".format(instanceDef),
        "# argsOptions={}".format(argsOptions)
        ])
    
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
   
    aV["instanceDef"] = instanceDef
    
    # Timing
    microSecs = time.time() 
    rList = pFile_read(aV["instanceDef"])
    aV["nDim"] = rList[0]
    aStruc = rList[1]
    aV["varList"] = rList[2]
    aV["density"] = rList[3]
    aV["coordRef"] = aV["varList"]
    aV["instanceID"] = os.path.basename(aV["instanceDef"]).strip(".lop") 

    infoSolutionsDir = os.path.dirname(os.path.realpath(aV["instanceDef"]))
    infoSolutionsFile = all_info["sandboxName"] + ".info_solutions.txt"
    infoSolutionsFile = "/".join([all_info["sandboxPath"], "xBenchm", "lop",
        infoSolutionsFile])
    aV["infoSolutionsFile"] = infoSolutionsFile

    if not os.path.isfile(infoSolutionsFile):
        #Error
        print "\nERROR from {}:\nfile {} is missing!\n".format(thisProc, infoSolutionsFile)
    rList = core.file_read(infoSolutionsFile).split('\n')
    rList.pop()
    rListTmp = []
    for line in rList:
        firstChar = line[0]
        if firstChar != "#" and len(line) > 0:
            rListTmp.append(line)
    isFound = 0
    for line in rListTmp:
        line = line.split()
        varName = line[0]
        if varName == aV["instanceID"]:
            aV["valueTarget"] = int(line[1].strip('-'))
            aV["isProven"] = line[2].strip('-')
            isFound = 1
        if isFound:
            break
    if not isFound:
        #Error
        print ("\nERROR from {}:"
                "\n .. instance {} was not found in file"
                "\n     {}\n".format(thisProc, aV["instanceID"], infoSolutionsFile))
    #end timing
    microSecs = time.time() - microSecs

    aV["runtimeRead"] = microSecs
    aV["infoSolutionsFile"] = infoSolutionsFile
    aV["commandName"] = all_info["sandboxName"] + ".main"
    aV["commandLine"] = "{} {} {}".format(aV["commandName"] , instanceDef , argsOptions)
    aV["valueTarget"] = aV["valueTarget"] * (1 - aV["valueTol"])

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
                #Error
                print ("\nERROR from {}:"
                        "\n.. option name {} is not either of two lists below:"
                        "\n{}"
                        "\n\nor\n"
                        "\n{}".format(thisProc, name, namesOptional, namesOptionalBool))
                break

    if aV["seedInit"] == "NA":
        aV["seedInit" ] = 1e9 * random.random()
        random.seed(aV["seedInit"])
    elif isinstance(aV["seedInit"], (int,long)):
        random.seed(aV["seedInit"])
    else:
        #Error
        print ("ERROR from {}:"
                ".. only -seedInit NA or -ssedInit <int> are valid assignments,"
                "not -seedInit {}\n".format(thisProc, aV['seedInit']))
    
    if aV["coordInit"] == "NA":
        aV['coordInit'] = P.coord.rand(aV["nDim"])
        aV['rankInit'] = P.coord.rank(aV["coordInit"])
    else:
        aV["coordInit"] = [int(c) for c in aV["coordInit"].split(",")]
        if len(aV["coordInit"]) != aV["nDim"]:
            #Error
            print ("\nERROR from {}:"
                    "\nThe permutation coordinate is of length {},"
                    "not the expected length {}\n".format(thisProc, aV["coordinit"], aV["nDim"]))
        aV["rankInit"] = P.coord.rank(aV["coordInit"])

    if aV["walkIntervalLmt"] == "NA" and aV["walkIntervalCoef"] != "NA":
        try:
            walkIntervalCoef = float(aV["walkIntervalCoef"])
            if walkIntevalCoef > 0.:
                aV["walkIntervalLmt"] = int(walkIntervalCoef * aV["nDim"])
                aV["walkIntervalCoef"] = walkIntervalCoef
        except:
            #Error
            print ("\nERROR from {}:"
                    "The walkIntervalCoef can only be assigned a value of NA"
                    "or a positive number, not {} \n".format(thisProc, aV["walkIntervalCoef"]))
    elif aV["walkIntervalLmt"] != "NA" and aV["walkIntervalCoef"] == "NA":
        try:
            walkIntervalLmt = int(aV["walkIntervalLmt"])
            if walkIntervalLmt > 0:
                aV["walkIntervalLmt"] = walkIntervalLmt
        except:
            #Error
            print ("\nERROR from {}:"
                    "The walkIntervalLmt can only be assigned a value of NA"
                    "or a positive number, not {} \n".format(thisProc, aV["walkIntervalLmt"]))
    elif aV["walkIntervalLmt"] != "NA" and aV["walkIntervalCoef"] != "NA":
        #Error
        print ("ERROR from {}:"
                "The walkIntervalLmt and walkIntervalCoef can only be assigned"
                "pairwise values of\n(NA NA) (default) (NA double) or (integer NA)"
                "not ({} {})\n".format(thisProc, aV["walkIntervalLmt"], aV["walkIntervalCoef"]))

    if aV["walkSegmLmt"] == "NA" and aV["walkSegmCoef"] != "NA":
        try:
            walkSegmCoef = float(aV["walkSegmCoef"])
            if walkIntevalCoef > 0.:
                aV["walkSegmLmt"] = int(walkSegmCoef * aV["nDim"])
                aV["walkSegmCoef"] = walkSegmCoef
        except:
            #Error
            print ("\nERROR from {}:"
                    "The walkSegmCoef can only be assigned a value of NA"
                    "or a positive number, not {} \n".format(thisProc, aV["walkSegmCoef"]))
    elif aV["walkSegmLmt"] != "NA" and aV["walkSegmCoef"] == "NA":
        try:
            walkSegmLmt = int(aV["walkSegmLmt"])
            if walkSegmLmt > 0:
                aV["walkSegmLmt"] = walkSegmLmt
        except:
            #Error
            print ("\nERROR from {}:"
                    "The walkSegmLmt can only be assigned a value of NA"
                    "or a positive number, not {} \n".format(thisProc, aV["walkSegmLmt"]))
    elif aV["walkSegmLmt"] != "NA" and aV["walkSegmCoef"] != "NA":
        #Error
        print ("ERROR from {}:"
                "The walkSegmLmt and walkSegmCoef can only be assigned"
                "pairwise values of\n(NA NA) (default) (NA double) or (integer NA)"
                "not ({} {})\n".format(thisProc, aV["walkSegmLmt"], aV["walkSegmCoef"]))

    aV["coordType"] = "P"
    aV["functionDomain"] = "P.lop"
    aV["functionID"] = "lop"

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

    aV["solverVersion"] = time.strftime("%Y %m %d %H %M %S")
    aV["timeStamp"] =  time.strftime("%Y %m %d %H %M %S")
    aV["dateLine"] = "@TODO: DATE LINE (line 582)"
    #@TODO HostID stuff (lines 583-486)

    # Timing
    microSecs = time.time() 
    aV["valueInit"] = f(aV["coordInit"])
    #end timing
    microSecs = time.time() - microSecs

    aV["runtime"] = microSecs
    aV["cntProbe"] = 1
    aV["cntStep"] = 0
    aV["coordPivot"] = aV["coordInit"]
    aV["coordBest"] = aV["coordInit"]
    aV["valuePivot"] = aV["valueInit"]
    aV["valueBest"] = aV["valueInit"]
    aHashTmp = {}
    aHashTmp[tuple(aV["coordInit"])] = []

    if aV["valueInit"] == aV["valueTarget"]:
        aV["targetReached"] = 1
    elif aV["valueInit"] < aV["valueTarget"]:
        aV["targetReached"] = 2
    else:
        aV["targetReached"] = 0

    aHashWalk[tuple(aV["coordInit"])] = []
    aV["isCensored"] = 0
    aV["cntRestart"] = 0
    aV["walkLength"] = aV["cntStep"]
    if aV["neighbDist"] == 1:
        aV["neighbSize"] = aV["nDim"] - 1
    else:
        aV["neighbSize"] = "dynamic"

    if aV["writeVar"] >= 4:
        aV["isSimple"] = 1
    aValueBest[aV["valueInit"]] = [0,0,aV["coordInit"]]
    aWalkBest[aV["valueInit"]] = [0,0,aV["coordInit"],0,0]

    aWalk[aV["cntStep"]] = "{} {} {} {} {} {}".format(aV["cntStep"], aV["cntRestart"],
            aV["coordPivot"], aV["valuePivot"], aV["neighbSize"], aV["cntProbe"])

    isPivot = 1
    aV["rankPivot"] = P.coord.rank(aV["coordPivot"])
    aWalkProbed[(aV["walkLength"],0)] = "{} {} {} {} {} {} {} {}".format(
            aV["walkLength"], aV["cntRestart"], aV["coordPivot"],
            aV["valuePivot"], aV["rankPivot"], isPivot, aV["neighbSize"], aV["cntProbe"])

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

   
    if aV["targetReached"] > 0:
        print "# BINGO: valueTarget has been reached or exceeded with coordInit"
        stdout(withWarning=1)
        return

    if aV["writeVar"] == 1:
        print "@TODO: Lines 675-678"

    return "{} {} {}".format(aV["targetReached"], aV["valueInit"], aV["coordInit"])

def main( instanceDef, args ):
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

    # (1) Phase 1: to query about the commandLine **OR** 
    #   to initialize all variables
    if instanceDef == "?":
        info(all_info["infoVariablesFile"],1)
        return
    else:
        rList = init(instanceDef,args)

    if len(rList) == 0:
        return
    elif aV["isInitOnly"]:
        print ("\n{}"
            "\n.. completed initialization of all variables,"
            "\n   exiting the solver since option -{} has been asserted."
            "\n{}\n".format("-"*78,aV["isInitOnly"],"-"*78))
    else:
        print ("\n{}"
            "\n.. completed initialization of all variables,"
            "\n   proceeding with the search under solverID = {}"
            "\n{}\n".format("-"*78,aV["solverID"],"-"*78))
    
    # (2) Phase 2: proceed with the combinatorial search
    cntRestart = 0
    coordInitList = aV["coordInit"]
    aV["cntRestartUniq"] = 0
    aV["cntRestart"] = cntRestart
    print ( "# FROM {}: initialized for restart={}"
        "\n# coordInit={}, valueInit={}".format(thisProc, cntRestart,
            aV["coordInit"], aV["valueInit"]))

    while cntRestart < aV["cntRestartLmt"]:
        # invoke the self-avoiding walk and return aV["targetReached"] (global)
        #@TODO Dynamic function calling? (P.lop.saw())
        saw(cntRestart)
        if aV["targetReached"] > 0:
            break
        else:
            cntRestart += 1
            aV["cntRestart"] = cntRestart
            aV["coordInit"] = P.coord.rand(aV["nDim"])

            #@TODO Dynamic function calling? (P.lop.f())
            rList = f(aV["coordInit"])
            aV["valueInit"] = rList[0]
            coordInitList.append(aV["coordInit"])

        aV["cntRestartUniq"] = len(coordInitList) - 1
        if aV["valueInit"] == aV["valueTarget"]:
            aV["targetReached"] = 1
        elif aV["valueInit"] < aV["valueTarget"]:
            aV["targetReached"] = 2
        else:
            aV["targetReached"] = 0

        if aV["targetReached"] > 0:
            break
        elif aV["writeVar"] == 2:
            print ("# FROM {}: initialized for restart={}"
                "\n# coordInit={}, valueInit={}".format(thisProc, cntRestart,
                    aV["coordInit"], aV["valueInit"]))

        if aV["runtime"] > aV["runtimelmt"]:
            break
    
    stdout(withWarning=1)
    return
    
def stdout( withWarning = 1 ):
    return

def info( infoVariablesFile="../xLib/P.lop.info_variables.txt", isQuery=0 ):
    thisProc = "P.lop.info"
    ABOUT = (
        "This proc reads contents of *info_variables.txt and generates a "
        "stdout response to user's query such as\n"
        "\tP.lop.main ? (under tclsh)\n"
        "or\n"
        "\t../xBin/P.lopT (under bash)"
    )
    if isQuery == "??":
        print ABOUT
        return

    #info global variables
    global all_info
    global all_valu
    global aV

    rList = util.table_info_variables(infoVariablesFile) 
    all_info = rList[0]
    all_valu = rList[1]

    if not isQuery:
        return (all_info, all_valu)

    optInfoList = [ "runtimeLmt", "cntProbeLmt", "cntRestartLmt", 
            "walkLengthLmt", "seedInit", "coordInit", "valueInit", 
            "valueTarget", "valueTol", "neighbDist", "walkSegmLmt",
            "walkSegmCoef", "walkIntervalLmt", "walkIntervalCoef",
            "walkRepeatsLmt", "isSimple", "writeVar" ]

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

    for name in optInfoList:
        value = all_valu[name]
        nameList = []
        valueList = []
        nameList.append(name)
        valueList.append(value)
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
