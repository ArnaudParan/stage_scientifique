#!/usr/bin/python3.4
#-*-coding:utf-8-*

import numpy as np
from scipy.sparse import csr_matrix

#si on intègre advection ou diffusion ou les deux
diffusion=False
advection=True

#coeff diffusion
D=1E1
a=-1E-1

#définition du maillage et du temps
lx=1E1
ly=1E1
tmax=1*1e-1

#nombre d'affichages
naff=10

#température maximale (pour l'affichage des graphes)
Tmax=10.

#condition initiale intérieure
def u0_2d(x,y):
	epsilon=lx/5.
	if(x<=lx/2.+epsilon and x>=lx/2.-epsilon and y<=ly/2.+epsilon and y>=ly/2.-epsilon):
		return 10.
	return 0.

#conditions de bord
def ubord_2d(x,y):
	if(x==0.):
		return 10.
	elif(x==lx):
		return 10.
	return 0.


dx=lx/100.
dt=tmax/100.
