# xPeriment configuration file:  xPeriment-config-P.lopT-test1.tcl
# ------------------------------------------------------------------
# these variables can be entered in any order
set xP(solverDir)      ../xBin/

set xP(solverOptions)  {
    -runtimeLmt 7200 -seedInit 1903 
}
set xP(sampleSize)   10

set xP(sampleColumns) {
    instanceID instanceDim hostID sandboxID solverName 
    valueTarget valueBest targetReached isCensored runtimeLmt runtime 
    cntProbe walkLength cntRestart speedProbe coordBest
}

set xP(stdoutNames) {
    instanceDef solverID coordInit coordBest nDim  
    walkLengthLmt walkLength  cntRestartLmt cntRestart   
    cntProbeLmt cntProbe runtimeLmt runtime runtimeRead speed hostID 
    isSimple neighbDist solverMethod 
    walkSegmLmt walkSegmCoef walkIntervalLmt walkIntervalCoef 
    walkRepeatsLmt seedInit  valueInit 
    valueBest  valueTarget valueTol targetReached isCensored
} 


set xP(instanceDefList) {
    ../xBenchm/lop/pal/pal-11-35.lop
    ../xBenchm/lop/pal/pal-13-55.lop
    ../xBenchm/lop/pal/pal-19-107.lop
}


set xP(instanceDefListFull) {
    ../xBenchm/lop/pal/pal-11-35.lop
    ../xBenchm/lop/pal/pal-13-55.lop
    ../xBenchm/lop/pal/pal-19-107.lop
    ../xBenchm/lop/pal/pal-23-161.lop
    ../xBenchm/lop/pal/pal-27-252.lop
    ../xBenchm/lop/pal/pal-31-297.lop
    ../xBenchm/lop/pal/pal-31-300.lop
    ../xBenchm/lop/pal/pal-43-558.lop
    ../xBenchm/lop/pal/pal-43-597.lop
    ../xBenchm/lop/pal/pal-55-1054.lop
    ../xBenchm/lop/pal/pal-55-1084.lop
}
