#! /usr/bin/env python3.6
from pylab import *

data = genfromtxt("output")
rank = data.shape[1]
time = data.shape[0]
print(rank, time)

nrows = int(sqrt(rank))
print(nrows)
result = zeros((nrows, nrows))
xx = linspace(0, 1, nrows + 2, endpoint = True)

for i in range(0, time):
    result = data[i, :].reshape((nrows, nrows))
    #contourf(result.T)
    plot(xx[1:-1], result[:,nrows//2])
    plot(xx[1:-1:2], (exp(10 * xx[1:-1:2]) - 1) / (exp(10) - 1), 'ro')
    show()
