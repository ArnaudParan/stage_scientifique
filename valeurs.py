#!/usr/bin/python3.4
#-*-coding:utf-8-*

import numpy as np

def B(n):
	B=np.diag(np.ones(n-1),1)
	B+=B.transpose()
	B+=(-2)*np.identity(n)
	return (-1)*B

def I(n):
	return np.identity(n)

#coeff diffusion
D=1E1

#définition du maillage et du temps
lx=10
dx=0.1
tmax=0.3
dt=0.01

#nombre d'affichages
naff=5

#grandeurs relatives à l'équation
nu=D*dt/(dx*dx)
nx=int(lx/dx)
nt=int(tmax/dt)

#condition initiale intérieure
U0=np.array([0.]*nx)
U0[nx/4]=10.
U0[3*nx/4]=10.

#conditions de bord
Ubord=np.array([0.]*nx)
Ubord[0]=10.
Ubord[nx-1]=10.


U0+=Ubord
