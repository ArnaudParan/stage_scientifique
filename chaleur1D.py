#!/usr/bin/python3.4
#-*-coding:utf-8-*

import numpy as np
import pylab as pl
from pylab import*
from valeurs import*
import matplotlib.pyplot as plt
import scipy as scipy
from scipy.sparse import csr_matrix
from scipy.sparse import*
from scipy.sparse.linalg.isolve.iterative import cg

#fonction annexe multipliant tous les termes d'une liste par un scalaire
def mult(liste, scalaire):
	retour=[]
	for i in range(0,len(liste)):
		retour.append(scalaire*liste[i])
	return retour

def euler1D(lx,tmax,D,dx,dt):
	#vérification des paramètres
	if dt>dx/(a*1.):
		print('Système mal paramétré, l\'algorithme procède à un reparamétrage')
		dt=dx/(2.*a)

	#grandeurs relatives à l'équation
	nuD=D*dt/(dx*dx)
	nuA=a*dt/dx
	nx=1+int(lx/dx)
	nt=1+int(tmax/dt)
	print(dx)
	print(dt)

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

	if diffusion and advection:
		if dt>dx/(a*1.):
			print('Système mal paramétré, l\'algorithme procède à un reparamétrage')
			dt=dx/(2.*a)
		A=I(nx)+nuD*B(nx)+nuA*C(nx)
	#la matrice du schéma de diffusion pure
	elif diffusion:
		A=I(nx)+nuD*B(nx)
	elif advection:
		if dt>dx/(a*1.):
			print('Système mal paramétré, l\'algorithme procède à un reparamétrage')
			dt=dx/(2.*a)
		A=I(nx)+nuA*C(nx)

	#on implémente le schéma d'euler
	for n in range(1,nt):
		Un=scipy.sparse.linalg.isolve.iterative.bicg(A,U[n-1])[0]
		#Un=linalg.solve(A,U[n-1])
		Un[0]=Ubord[0]
		Un[nx-1]=Ubord[nx-1]
		#Un[0]=30.*sin(n*pi/10.)
		#Un[nx-1]=0.
		U.append(Un)
	
	return U

def tracer(U,lx,tmax,dx,dt,n=naff, T=Tmax):
	#on trace la solution
	nx=1+int(lx/dx)
	nt=1+int(tmax/dt)
	x=range(0,nx)
	#première méthode, affiche des pages
	"""for n in range(0,n):
		plt.plot(mult(x,dx),U[n*nt/n][x])
		plt.axis([0,int(nx*dx),0,T])
		plt.ylabel('T')
		plt.xlabel('x, t='+str(dt*n*nt/n)+', D='+str(D))
		plt.show()"""

	#autre méthode, on trace sur la même fenêtre
	subplot(2,2,1)
	plt.plot(mult(x,dx),U[0][x])
	plt.axis([0,int(nx*dx),0,T])
	plt.ylabel('T')
	plt.xlabel('x, t=0, D='+str(D))

	subplot(2,2,2)
	plt.plot(mult(x,dx),U[nt/3][x])
	plt.axis([0,int(nx*dx),0,T])
	plt.xlabel('t='+str(dt*nt/3.))

	subplot(2,2,3)
	plt.plot(mult(x,dx),U[2*nt/3][x])
	plt.axis([0,int(nx*dx),0,T])
	plt.xlabel('t='+str(dt*2.*nt/3.))
	
	subplot(2,2,4)
	plt.plot(mult(x,dx),U[nt-1][x])
	plt.axis([0,int(nx*dx),0,T])
	plt.xlabel('t='+str(nt*dt))
	
	show()
