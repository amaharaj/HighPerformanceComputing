import matplotlib.pyplot as plt
import numpy as np
import math
import pylab

f1 = open('numThreadsData.txt','r')
lines = f1.readlines()

nthreads = []
walltime = []

fig = plt.figure()

for line in lines:
   p = line.split()
   walltime.append(float(p[1]))
   nthreads.append(float(p[0]))

plt.plot(walltime,nthreads)

plt.xlabel("Wall Time")
plt.ylabel("Number of Threads")
plt.title("Wall Time w.r.t. Number of Threads")

plt.show()
