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
        sys.stderr.write("\n".join([
        "\nERROR from {}:".format(thisProc),
        "\nsize of reqAry+optAry ({}) is not matched to".format(nCheck),
        "\nsize of aPI ({})".format(len(aPI))]))
        sys.exit(1)

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
        sys.stderr.write("\n".join([
            "ERROR from {}:".format(thisProc),
            ".. only -seedInit NA or -seedInit <int> are valid assignments, not \
            -seedInit {}\ ".format(aPI["seedInit"])]))
        sys.exit(1)

    #ABOUT Coord Init
    if aPI["coordInit"] == "NA":
        aPI["coordInit"] = P_coord.rand(aPI["nDim"])
        aPI["rankInit"] = P_coord.inversion(aPI["coordInit"])
    else:
        coordTmp = aPI["coordInit"]
        coordTmp.sort()
        if len(coordTmp) != aPI["nDim"]:
            sys.stderr.write("\n".join([
                "\nERROR from {}:".format(thisProc),
                "The permutation coordinate is of length {},".format(len(coordTmp)),
                "not the expected length {}\n".format(aPI["nDim"])]))
            sys.exit(1)
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
                sys.stderr.write("\n".join([
                    "\nERROR from {}:".format(thisProc),
                    "The walkSegmCoef can only be assigned a value of NA or a\
                    positive number, not {}".format(aPI["walkSegmCeof"])]))
                sys.exit(1)
        except ValueError:
            sys.stderr.write("\n".join([
                "\nERROR from {}:".format(thisProc),
                "The walkSegmCoef can only be assigned a value of NA or a\
                positive number, not {}".format(aPI["walkSegmCeof"])]))
            sys.exit(1)
    elif aPIwalkSegmCoef == "NA": 
        try:
            aPIwalkSegmLmt = float(aPIWalkSegmLmt)
            if aPIwalkSegmLmt > 0:
                aPI["walkSegmLmt"] = int(aPIWalkSegmLmt)
            else:
                sys.stderr.write("\n".join([
                    "\nERROR from {}:".format(thisProc),
                    "The walkSegmLmt can only be assigned a value of NA or a\
                    positive number, not {}".format(aPI["walkSegmLmt"])]))
                sys.exit(1)
        except ValueError:
            sys.stderr.write("\n".join([
                "\nERROR from {}:".format(thisProc),
                "The walkSegmLmt can only be assigned a value of NA or a\
                positive number, not {}".format(aPI["walkSegmLmt"])]))
            sys.exit(1)
    else:
        sys.stderr.write("\n".join([
            "\nERROR from {}:".format(thisProc),
            "The walkSegmLmt and walkSegmCoef can only be assigned pairwise \
            values of",
            "(NA NA) (default), (NA double), or (integer, NA)"]))
        sys.exit(1)

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
                sys.stderr.write("\n".join([
                    "\nERROR from {}:".format(thisProc),
                    "The walkIntervalCoef can only be assigned a value of NA or a\
                    positive number > 0, not {}".format(aPIwalkIntervalCoef)]))
                sys.exit(1)
        except ValueError:
            sys.stderr.write("\n".join([
                "\nERROR from {}:".format(thisProc),
                "The walkIntervalCoef can only be assigned a value of NA or a\
                positive number > 0, not {}".format(aPIwalkIntervalCoef)]))
            sys.exit(1)
    elif aPIwalkSegmCoef == "NA":
        try:
            aPIwalkIntervalLmt = float(aPIWalkIntervalLmt)
            if aPIwalkIntervalLmt > 0:
                aPI["walkIntervalLmt"] = int(aPIWalkIntervalLmt)
            else:
                sys.stderr.write("\n".join([
                    "\nERROR from {}:".format(thisProc),
                    "The walkIntervalLmt can only be assigned a value of NA or a\
                    positive number > 0, not {}".format(aPIwalkIntervalLmt)]))
                sys.exit(1)
        except ValueError:
            sys.stderr.write("\n".join([
                "\nERROR from {}:".format(thisProc),
                "The walkIntervalLmt can only be assigned a value of NA or a\
                positive number > 0, not {}".format(aPIwalkIntervalLmt)]))
            sys.exit(1)
    else:
        sys.stderr.write("\n".join([
            "ERROR from {}:".format(thisProc),
            "The walkIntervalLmt and walkIntervalCoef can only be assigned pairwise values of",
            "(NA NA) (default), (NA double), or (integer, NA)"]))
        sys.exit(1)

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
