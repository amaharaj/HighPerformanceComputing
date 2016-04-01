import matplotlib.pyplot as plt
import numpy as np
import math
import pylab

f1 = open('numThreadsData2-2.txt','r')
f1_lines = f1.readlines()
f2 = open('numThreadsData1-2.txt','r')
f2_lines = f2.readlines()
f3 = open('numThreadsData3-2.txt','r')
f3_lines = f3.readlines()
f4 = open('numThreadsData4-2.txt','r')
f4_lines = f4.readlines()


nthreads = []
walltime1 = []
walltime2 = []
walltime3 = []
walltime4 = []

fig = plt.figure()

for l1 in f1_lines:
   p = l1.split()
   walltime1.append(float(p[1]))
   nthreads.append(float(p[0]))

for l2 in f2_lines:
   q = l2.split()
   walltime2.append(float(q[1]))

for l3 in f3_lines:
   q = l3.split()
   walltime3.append(float(q[1]))

for l4 in f4_lines:
   q = l4.split()
   walltime4.append(float(q[1]))

a, = plt.plot(nthreads,walltime1,label='n=10000, t_max=100')
b, = plt.plot(nthreads,walltime2,label='n=10000, t_max=1000')
c, = plt.plot(nthreads,walltime3,label='n=100000, t_max=1000')
d, = plt.plot(nthreads,walltime4,label='n=100000, t_max=100')

plt.legend(handles=[a,b,c,d])

plt.xlabel("Number of Threads")
plt.ylabel("Wall Time")
plt.title("Wall Time w.r.t. Number of Threads")

plt.show()
