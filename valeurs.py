#!/usr/bin/python3.4
#-*-coding:utf-8-*

import numpy as np
from scipy.sparse import csr_matrix

#si on intègre advection ou diffusion ou les deux
diffusion=True
advection=False

#coeff diffusion
D=1E1
a=1E1

#définition du maillage et du temps
lx=1E1
tmax=1*1e-1

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
