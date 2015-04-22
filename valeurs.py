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
D=0.1

#définition du maillage
lx=10
dx=0.1
tmax=10
dt=0.1

#nombre d'affichages
naff=5

nx=int(lx/dx)
nt=int(tmax/dt)

#condition initiale intérieure
U0=np.array([0.]*nx)
U0[nx/4]=10.
U0[3*nx/4]=10.

#conditions de bord ici, on fixe u à 0. sur le bord
Ubord=np.array([0.]*nx)
Ubord[0]=0.
Ubord[nx-1]=0.


U0+=Ubord



#grandeurs relatives à l'équation
nu=D*dt/(dx*dx)
