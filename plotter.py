# Рисует получившиеся многочлены
# Нужно переместить в build

import sys
import matplotlib.pyplot as plt
import numpy as np

# lagrange/spline
xrange = []
fxrange = []
with open('functionInfo.txt') as funcinfo:
    for eg in funcinfo:
        dta = eg.split(" ")
        xrange.append(float(dta[0]))
        fxrange.append(float(dta[1]))

plt.plot(xrange, fxrange)

tablex, tabley = [], []

with open('tableInfo.txt') as tableinfo:
    for eg in tableinfo:
        if len(eg) > 1:
            dta = eg.split(" ")
            tablex.append(float(dta[1]))
            tabley.append(float(dta[3]))

plt.scatter(tablex, tabley, c='r')

if sys.argv[1] == 'lagrange':
    def builder(degs, point):
        ansv = 0
        for (num, coef) in enumerate(degs):
            ansv += coef * point ** num
        return ansv


    # with open('LagrangeInterpolationInfo.txt') as Linfo:
    #     degs = [float(i) for i in list(Linfo)[0].split(" ")]
    # degs.reverse()
    # lxrange = [builder(degs, pt) for pt in xrange]
    # plt.plot(xrange, lxrange)
    lrange = []
    lxrange = []
    with open('FakeLagrangeInfo.txt') as Linfo:
        for eg in Linfo:
            dta = eg.split(" ")
            lrange.append(float(dta[0]))
            lxrange.append(float(dta[1]))
    plt.plot(lrange, lxrange)


if sys.argv[1] == 'spline':

    def builder(degs, point, num):
        return degs[0] + degs[1] * (point - tablex[num]) + degs[2] * (point - tablex[num]) ** 2 + degs[3] * (
                    point - tablex[num]) ** 3


    list_degs = []
    with open('SplineInterpolationInfo.txt') as Sinfo:
        for eg in Sinfo:
            if len(eg) > 1:
                list_degs.append([float(i) for i in eg.split(" ")[1:]])
    for (num, degs) in enumerate(list_degs):
        degs.reverse()
        xtemp = np.arange(tablex[num], tablex[num + 1], (tablex[num + 1] - tablex[num]) / 100)
        fxtemp = [builder(degs, pt, num) for pt in xtemp]
        plt.plot(xtemp, fxtemp, c='orange')

plt.grid()
plt.show()
