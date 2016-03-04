#!/usr/bin/python

from ellipticBeam import ellipticBeam
from matplotlib import pyplot as plt

t = 0.4
c = 0.009
beta = 0.7

myBunchGenerator = ellipticBeam(t, c, beta)
bunch = myBunchGenerator.generateBunch(1.e-5, 10000)
xArray = []
yArray = []
for idx in range(10000):
    xArray.append(bunch[idx][0])
    yArray.append(bunch[idx][2])

plt.scatter(xArray, yArray)
plt.show()
plt.clf()
