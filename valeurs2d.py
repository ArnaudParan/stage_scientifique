#!/usr/bin/python3.4
#-*-coding:utf-8-*

import numpy as np
from scipy.sparse import csr_matrix

#si on intègre advection ou diffusion ou les deux
diffusion=False
advection=False
dif_autre=True
#up down ou center
schema_adv='up'
schema_advx=schema_adv
schema_advy=schema_adv

#coeff diffusion
D=1E1
a=1E0
ax=-a
ay=a
s1=1e0
s2=-1e0

#définition du maillage et du temps
lx=1E1
ly=1E1
tmax=1*1e0

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
	"""if(x==0.):
		return 10.
	elif(x==lx):
		return 10.
	if(y==0.):
		return 10.
	elif(y==ly):
		return 10."""
	return 0.	


dx=lx/20.
dy=ly/20.
dt=tmax/100.
options=[diffusion,advection,schema_advx,schema_advx]
