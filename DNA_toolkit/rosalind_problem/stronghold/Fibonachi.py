def fibonachi_loop(num):
    old = 1
    new = 1
    for _ in range(num -1):
        temp = new
        new = old
        old = old + temp
    return new

#print(fibonachi_loop(6))

def fibonachi_loopP(num):
    old, new = 1,1
    for _ in range(num-1):
        new, old = old, old +new 
    return new

print(fibonachi_loopP(24))