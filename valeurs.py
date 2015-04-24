#!/usr/bin/python3.4
#-*-coding:utf-8-*

import numpy as np
from scipy.sparse import csr_matrix

#matrice pour la diffusion thermique
def B(n):
	B=np.diag(np.ones(n-1),1)
	B+=B.transpose()
	B+=(-2)*np.identity(n)
	B[0][n-1]=1.
	B[n-1][0]=1.
	return csr_matrix((-1)*B)

#matrice pour l'advection
def C(n):
	C=np.diag(-np.ones(n-1),-1)
	C+=np.diag(np.ones(n),0)
	C[0][n-1]=-1
	return csr_matrix(C)

def I(n):
	return csr_matrix(np.identity(n))

#si on intègre advection ou diffusion ou les deux
diffusion=False
advection=True

#coeff diffusion
D=1E1
a=1E1

#définition du maillage et du temps
lx=1E1
tmax=5*1e-1

#nombre d'affichages
naff=10

#température maximale (pour l'affichage des graphes)
Tmax=10.

#condition initiale intérieure
def u0(x):
	epsilon=lx/20.
	if(x<=lx/4.+epsilon and x>=lx/4.-epsilon):
		return 10.
	elif(x<=3.*lx/4.+epsilon and x>=3.*lx/4.-epsilon):
		return 10.
	return 0.

#conditions de bord
def ubord(x):
	"""if(x==0.):
		return 10.
	elif(x==lx):
		return 10."""
	return 0.


dx=lx/100.
dt=tmax/100.
