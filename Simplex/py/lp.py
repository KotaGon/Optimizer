from pulp import *
import math
m = LpProblem(sense=LpMinimize);
N = int(math.sqrt(sum([1 for _ in open('../out/dist.dat')])))
print(N)
with open("../out/dist.dat", "r") as f:
  dist = [[0.0 for j in range(N)] for i in range(N)]
  index = 0
  for line in f:
    dist[int(index/N)][int(index%N)] = float(line[:-1])
    index += 1

  x = [[0 for j in range(N)] for i in range(N)]
  for i in range(N):
    for j in range(i+ 1, N):
      x[i][j] = LpVariable("x" + str(i) + "_" + str(j), lowBound=0, upBound=1)

  f = 0
  for i in range(N):
    for j in range(i +1, N):
      f += dist[i][j] * x[i][j]

  m += f

  for i in range(N):
    f = 0
#for j in range(i + 1, N):
#     m += x[i][j] <= 1.0

    for j in range(0, i):
      f += x[j][i]
    for j in range(i+1,N):
      f += x[i][j]
#m += f <= 2
    m += f == 2

#solver = PULP_CBC_CMD(presolve=0, msg = 1) 
result = m.solve(pulp.PULP_CBC_CMD( msg=1, threads=1, presolve=False,timeLimit=10000, mip=False,fracGap = 0.0, gapAbs = 0.0))
m.writeLP("test.lp")
