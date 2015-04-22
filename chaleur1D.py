#!/usr/bin/python3.4
#-*-coding:utf-8-*

import numpy as np
from numpy import linalg
from valeurs import*

#définition de la solution
U=[U0]

#la matrice du schéma à inverser
A=I(nx)+nu*B(nx)

#nombre d'itérations
N=10
for n in range(1,N):
	U.append(linalg.solve(A,U[n-1]))

