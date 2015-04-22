#!/usr/bin/python3.4
#-*-coding:utf-8-*

import numpy as np
import pylab as pl
from numpy import linalg
from pylab import*
from valeurs import*
import matplotlib.pyplot as plt

#fonction annexe multipliant tous les termes d'une liste par un scalaire
def mult(liste, scalaire):
	retour=[]
	for i in range(0,len(liste)):
		retour.append(scalaire*liste[i])
	return retour


#définition de la solution
U=[U0]

#la matrice du schéma à inverser
A=I(nx)+nu*B(nx)

#on implémente le schéma d'euler
for n in range(1,nt):
	Un=linalg.solve(A,U[n-1])
	Un[0]=Ubord[0]
	Un[nx-1]=Ubord[nx-1]
	U.append(Un)

#on trace la solution
x=range(0,nx)
"""jfor n in range(0,naff-1):
	plt.plot(mult(x,dx),U[n*nt/naff][x])
	plt.ylabel('T')
	plt.xlabel('x, t='+str(dt*n*nt/naff)+', D='+str(D))
	plt.show()i"""

#autre méthode, on trace sur la même fenêtre
plt.title('blabla')
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

# savefig('../figures/subplot-grid.png', dpi=64)
show()
