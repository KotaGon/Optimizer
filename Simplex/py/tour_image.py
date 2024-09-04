import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import math

pos = []
from_list = [] 
to_list = []

with open("../out/pos.dat", "r") as fp:
  for line in fp:
    pos.append(line[:-1].split("\t")[0:2])

N = len(pos)

flow = [0] * N

opt_sol = 0.0
with open("../out/sol.dat", "r") as fp:
  for line in fp:
   elm = line[:-1].split("\t")
   _from = int(elm[0])
   _to = int(elm[1])
   sol = float(elm[2])
   dist = math.sqrt( (float(pos[_from][0]) - float(pos[_to][0]))**2 + (float(pos[_from][1]) - float(pos[_to][1]))**2)
#opt_sol += sol * float(elm[3])
   opt_sol += sol * dist

   if(sol > 1):
    print("sol > 1")

   flow[_from] += sol
   flow[_to] += sol

   if(sol > 1.0e-6): 
     from_list.append(_from)
     to_list.append(_to)

for val in flow:
  if val != 2:
    print("wa")
print("opt sol = {0}".format(opt_sol))
    
start_pos = []
end_pos = []

N = len(from_list)

for i in range(N):
  start_pos.append(pos[from_list[i]])
for i in range(N):
  end_pos.append(pos[to_list[i]])

lines = [ [sp,ep] for sp,ep in zip(start_pos, end_pos)]

lc = LineCollection(lines)

fig, ax = plt.subplots()
ax.add_collection(lc)
ax.autoscale()
plt.show()
