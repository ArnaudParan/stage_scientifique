#!/usr/bin/python3.4
#-*-coding:utf-8-*

import numpy as np
import pylab as pl
from numpy import linalg
from pylab import*
from valeurs import*
from gradientConjugue import*
import matplotlib.pyplot as plt

#fonction annexe multipliant tous les termes d'une liste par un scalaire
def mult(liste, scalaire):
	retour=[]
	for i in range(0,len(liste)):
		retour.append(scalaire*liste[i])
	return retour


#grandeurs relatives à l'équation
nu=D*dt/(dx*dx)
nx=1+int(lx/dx)
nt=1+int(tmax/dt)

#condition initiale intérieure
U0=np.array([0.]*nx)
for i in range(0,nx-1):
	U0[i]=u0(i*dx)

#condition de bord
Ubord=np.array([0.]*nx)
Ubord[0]=ubord(0.)
Ubord[nx-1]=ubord(lx)

U0+=Ubord

#définition de la solution
U=[U0]

#la matrice du schéma à inverser
A=I(nx)+nu*B(nx)

#on implémente le schéma d'euler
for n in range(1,nt):
	Un=conjugate_gradient(A,U[n-1],U[n-1])
	#Un=linalg.solve(A,U[n-1])
	Un[0]=Ubord[0]
	Un[nx-1]=Ubord[nx-1]
	U.append(Un)

#on trace la solution
x=range(0,nx)
"""for n in range(0,naff-1):
	plt.plot(mult(x,dx),U[n*nt/naff][x])
	plt.ylabel('T')
	plt.xlabel('x, t='+str(dt*n*nt/naff)+', D='+str(D))
	plt.show()"""

#autre méthode, on trace sur la même fenêtre
subplot(2,2,1)
plt.plot(mult(x,dx),U[0][x])
plt.ylabel('T')
plt.xlabel('x, t=0, D='+str(D))

subplot(2,2,2)
plt.plot(mult(x,dx),U[nt/3][x])
plt.xlabel('t='+str(dt*nt/3.))

subplot(2,2,3)
plt.plot(mult(x,dx),U[2*nt/3][x])
plt.xlabel('t='+str(dt*2.*nt/3.))

subplot(2,2,4)
plt.plot(x,U[nt-1][x])
plt.xlabel('t='+str(nt*dt))

show()
