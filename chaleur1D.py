#!/usr/bin/python3.4
#-*-coding:utf-8-*

import numpy as np
from numpy import linalg
from valeurs import*
import matplotlib.pyplot as plt

#définition de la solution
U=[U0]

#la matrice du schéma à inverser
A=I(nx)+nu*B(nx)

#on implémente le schéma d'euler
for n in range(1,nt):
	Un=linalg.solve(A,U[n-1])
	Un+=Ubord
	U.append(Un)

#on trace la solution
x=range(0,nx)
for n in range(0,naff-1):
	plt.plot(x,U[n*nt/naff][x])
	plt.show()
