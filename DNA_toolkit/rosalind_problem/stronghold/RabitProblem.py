def fibonachi_loopP(months, offsprings):
    parent, child = 1,1
    for _ in range(months-1):
        child, parent = parent, parent +(child * offsprings) 
    return child

print(fibonachi_loopP(34 ,3))