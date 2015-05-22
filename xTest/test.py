from Globals import *

def f():
    global one
    print numbers.one
    one = numbers.one
    numbers.one = 2
    print one, numbers.one

if __name__ == "__main__":
    f()
