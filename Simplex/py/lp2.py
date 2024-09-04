from pulp import *
m = LpProblem(sense=LpMinimize);
x0 = LpVariable("x0", lowBound=0)
x1 = LpVariable("x1", lowBound=0)
x2 = LpVariable("x2", lowBound=0)
x3 = LpVariable("x3", lowBound=0)
x4 = LpVariable("x4", lowBound=0)
x5 = LpVariable("x5", lowBound=0)
x6 = LpVariable("x6", lowBound=0)
x7 = LpVariable("x7", lowBound=0)
x8 = LpVariable("x8", lowBound=0)
x9 = LpVariable("x9", lowBound=0)
x10 = LpVariable("x10", lowBound=0)
x11 = LpVariable("x11", lowBound=0)
x12 = LpVariable("x12", lowBound=0)
x13 = LpVariable("x13", lowBound=0)
x14 = LpVariable("x14", lowBound=0)
x15 = LpVariable("x15", lowBound=0)
x16 = LpVariable("x16", lowBound=0)
x17 = LpVariable("x17", lowBound=0)
x18 = LpVariable("x18", lowBound=0)
x19 = LpVariable("x19", lowBound=0)
x20 = LpVariable("x20", lowBound=0)
x21 = LpVariable("x21", lowBound=0)
x22 = LpVariable("x22", lowBound=0)
x23 = LpVariable("x23", lowBound=0)
x24 = LpVariable("x24", lowBound=0)
x25 = LpVariable("x25", lowBound=0)
x26 = LpVariable("x26", lowBound=0)
x27 = LpVariable("x27", lowBound=0)
x28 = LpVariable("x28", lowBound=0)
x29 = LpVariable("x29", lowBound=0)
x30 = LpVariable("x30", lowBound=0)

m += 3.66 -1.02 *x0 + 0.07*x1 -1.02*x4 + 0.01*x7 + 0.07*x8 -0.21*x11 -0.11*x15 + 0.01*x17  -0.21*x18  -0.11*x19  -0.58*x20 + 0.65*x23 -1.25*x24  -1.02*x25 + 0.06*x27 + 0.31*x29
m += 0 <= 1 + 1*x0 + 1*x4  -1*x11 -1*x14  -1*x15 -1*x18  -1*x19 + 0.5*x20  -0.5*x23 + -1*x24  -0.5*x25 + 0.5*x27 + 0.5*x29
m += 0 <= 0  -1*x20 + 1*x24 + 1*x25
m += 0 <= 0 -1*x23 + 1*x24 + 1*x25
m += 0<= 2  -1*x7  -1*x11  -1*x15  -1*x16  -1*x17  -1*x18  -1*x19  -1*x24  -1*x25 + 1*x29
m += 0<= 1 + 1*x1  -1*x7 + 1*x8  -1*x13  -1*x15  -1*x17  -1*x19 + 0.5*x20 + 0.5*x23  -1*x24  -1.5*x25 + 0.5*x27 + 0.5*x29
m += 0<= 0  -1*x0  -1*x1  -1*x4 + 1*x7  -1*x8 + 1*x11  -1*x12 + 1*x15 + 1*x17 + 1*x18 + 1*x19  -1*x20 + 1*x24 + 1*x25  -1*x29
m += 0 <= 0 + 1*x24 + 1*x25 + -1*x27
m += 0<= 1 -1*x0 -1*x1 -1*x4 -1*x8 -1*x9 + 1*x15 + 1*x19  -0.5*x20 + 0.5*x23 + 0.5*x25  -0.5*x27  -0.5*x29
m += 0 <= 0 + 1*x24 + 1*x25 -1*x29
m  += 0 <= 0 + 1*x24 + 1*x25
      
m.solve()
print(value(x29))
