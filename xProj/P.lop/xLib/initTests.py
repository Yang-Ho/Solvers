execfile('all_python.py')
#P.lop.init('../xBenchm/lop/tiny/i-5-book2.lop',['-isInitOnly'])
P.lop.init('../xBenchm/lop/tiny/i-5-book2.lop',['-isInitOnly','-writeVar','3','-coordInit','5,3,2,1,4'])
P.lop.saw_pivot_simple([5,3,2,1,4],-46)
