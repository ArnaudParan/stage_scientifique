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
lx=10.
tmax=0.3

#nombre d'affichages
naff=5



#condition initiale intérieure
def u0(x):
	epsilon=lx/20.
	if(x<=lx/4.+epsilon and x>=lx/4.-epsilon):
		return 10.
	elif(x<=3.*lx/4.+epsilon and x>=3.*lx/4.-epsilon):
		return 10.
	else:
		return 0.

#conditions de bord
def ubord(x):
	if(x==0.):
		return 10.
	elif(x==lx):
		return 10.
	else:
		return 0.


dx=lx/100.
dt=tmax/100.
