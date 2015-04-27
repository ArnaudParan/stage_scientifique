#!/usr/bin/python3.4
#-*-coding:utf-8-*

from matplotlib import cm
import numpy as np
import pylab as pl
from pylab import*
from valeurs import*
from valeurs2d import*
import matplotlib.pyplot as plt
import scipy as scipy
from scipy.sparse import csr_matrix
from scipy.sparse import*
from scipy.sparse.linalg.isolve.iterative import cg
from pylab import *
from mpl_toolkits.mplot3d import Axes3D


#fait la transition de la représentation machine en 1D à la représentation humaine en 2D
def to1d(i,j,n):
	return i+n*j

def to2d(a,n):
	return [a%n,int(a/n)]

def to2dArray(u,n,m):
	retour=[]
	for i in range(n):
		retour.append([0]*m)
	for i in range(n*m):
		retour[(i%n)][(i/n)]=u[i]
	return retour

#matrice pour la diffusion thermique
def Bx(n,m):
	B=(2.)*np.identity(n*m)
	for y in range(0,m):
		for x in range(1,n-1):
			B[to1d(x+1,y,n)][to1d(x,y,n)]=-1.
			B[to1d(x-1,y,n)][to1d(x,y,n)]=-1.
		B[to1d(0,y,n)][to1d(n-1,y,n)]=-1.
		B[to1d(n-1,y,n)][to1d(0,y,n)]=-1.
		B[to1d(n-2,y,n)][to1d(n-1,y,n)]=-1.
		B[to1d(1,y,n)][to1d(0,y,n)]=-1.
	return csr_matrix(B)
	
def By(n,m):
	B=(2.)*np.identity(n*m)
	for x in range(0,n):
		for y in range(1,m-1):
			B[to1d(x,y+1,n)][to1d(x,y,n)]=-1.
			B[to1d(x,y-1,n)][to1d(x,y,n)]=-1.
		B[to1d(x,0,n)][to1d(x,m-1,n)]=-1.
		B[to1d(x,m-1,n)][to1d(x,0,n)]=-1.
		B[to1d(x,m-2,n)][to1d(x,m-1,n)]=-1.
		B[to1d(x,1,n)][to1d(x,0,n)]=-1.
	return csr_matrix(B)

#matrice pour l'advection"
def Cx(n,m):
	for y in range(0,m):
		if up:
			Cx=(1.)*np.identity(n*m)
			for x in range(1,n-1):
				Cx[to1d(x-1,y,n)][to1d(x,y,n)]=-1.
			Cx[to1d(n-1,y,n)][to1d(0,y,n)]=-1.
			Cx[to1d(n-2,y,n)][to1d(n-1,y,n)]=-1.
		else:
			Cx=(-1.)*np.identity(n*m)
			for x in range(1,n-1):
				Cx[to1d(x+1,y,n)][to1d(x,y,n)]=1.
			Cx[to1d(0,y,n)][to1d(n-1,y,n)]=1.
			Cx[to1d(1,y,n)][to1d(0,y,n)]=1.
	return csr_matrix(Cx)
	
def Cy(n,m):
	for x in range(0,n):
		if up:
			B=(1.)*np.identity(n*m)
			for y in range(1,m-1):
				B[to1d(x,y-1,n)][to1d(x,y,n)]=-1.
			B[to1d(x,m-1,n)][to1d(x,0,n)]=-1.
			B[to1d(x,m-2,n)][to1d(x,n-1,n)]=-1.
		else:
			B=(-1.)*np.identity(n*m)
			for y in range(1,m-1):
				B[to1d(x,y+1,n)][to1d(x,y,n)]=1.
			B[to1d(x,0,n)][to1d(x,n-1,n)]=1.
			B[to1d(x,1,n)][to1d(x,0,n)]=1.
	return csr_matrix(B)

def I(n,m):
	return csr_matrix(np.identity(m*n))

#fonction annexe multipliant tous les termes d'une liste par un scalaire
def mult(liste, scalaire):
	retour=[]
	for i in range(0,len(liste)):
		retour.append(scalaire*liste[i])
	return retour

#réinitialise le bord dans le code
def initBord(U, Ubord, nx, ny):
	for x in range(0,nx):
		U[to1d(x,0,nx)]=Ubord[to1d(x,0,nx)]
		U[to1d(x,ny-1,nx)]=Ubord[to1d(x,ny-1,nx)]
	for y in range(0,ny):
		U[to1d(0,y,nx)]=Ubord[to1d(0,y,nx)]
		U[to1d(nx-1,y,nx)]=Ubord[to1d(nx-1,y,nx)]

def euler2D(lx,ly,tmax,D,ax,ay,dx,dy,dt):
	#vérification des paramètres
	if advection and (dt>dx/abs(a*10.) or dt>dy/abs(a*10.)):
		print('Système mal paramétré, l\'algorithme procède à un reparamétrage')
		dt=min([dx,dy])/abs(10.*a)
	#vérification des paramètres
	if (dt>dx/abs(a*10.) or dt>dy/abs(a*10.)):
		print('Système mal paramétré, l\'algorithme procède à un reparamétrage')
		dt=min([dx,dy])/abs(10.*a)

	#grandeurs relatives à l'équation
	nuDx=D*dt/(dx*dx)
	nuDy=D*dt/(dy*dy)
	nuAx=ax*dt/dx
	nuAy=ay*dt/dy
	nx=1+int(lx/dx)
	ny=1+int(ly/dy)
	nt=1+int(tmax/dt)

	#condition initiale intérieure
	U0=np.array([0.]*nx*ny)
	for i in range(0,nx*ny):
		U0[i]=u0_2d(to2d(i,nx)[0]*dx,to2d(i,nx)[1]*dy)
	
	#condition de bord
	Ubord=np.array([0.]*nx*ny)
	for i in range(0,nx*ny):
		Ubord[i]=ubord_2d(to2d(i,nx)[0]*dx,to2d(i,nx)[1]*dy)

	initBord(U0,Ubord,nx,ny)

	#définition de la solution
	U=[U0]

	#définition du schéma
	if diffusion and advection:
		A=I(nx,ny)+nuDx*Bx(nx,ny)+nuDy*By(nx,ny)+nuAx*Cx(nx,ny)+nuAy*Cy(nx,ny)
	elif diffusion:
		A=I(nx,ny)+nuDx*Bx(nx,ny)+nuDy*By(nx,ny)
	elif advection:
		A=I(nx,ny)+nuAx*Cx(nx,ny)+nuAy*Cy(nx,ny)

	#on implémente le schéma d'euler
	for n in range(1,nt):
		Un=scipy.sparse.linalg.isolve.iterative.bicg(A,U[n-1])[0]
		#Un=linalg.solve(A,U[n-1])
		initBord(Un,Ubord,nx,ny)
		U.append(Un)
	
	return U

def tracer(U,lx,ly,tmax,dx,dy,dt,n=naff, T=Tmax):
	#on trace la solution
	nx=1+int(lx/dx)
	ny=1+int(ly/dy)
	nt=1+int(tmax/dt)
	x=mult(np.arange(0,nx),dx)
	y=mult(np.arange(0,ny),dy)
	x,y=np.meshgrid(x,y)
	#première méthode, affiche des pages
	"""for n in range(0,n):
		fig=figure()
		ax = Axes3D(fig)
		ax.plot_surface(x,y,to2dArray(U[n*nt/naff],nx,ny), rstride=1, cstride=1, cmap='hot')
		ax.set_zlim(0,10)
		show()"""

	#autre méthode, on trace sur la même fenêtre
	fig1=figure(1)
	ax = Axes3D(fig1)
	ax.plot_surface(x,y,to2dArray(U[0],nx,ny), rstride=1, cstride=1, cmap=cm.jet)
	ax.set_zlim(0,10)
	show()

	fig2=figure(2)
	ax = Axes3D(fig2)
	ax.plot_surface(x,y,to2dArray(U[nt/3],nx,ny), rstride=1, cstride=1, cmap=cm.jet)
	ax.set_zlim(0,10)
	show()

	fig3=figure(3)
	ax = Axes3D(fig3)
	ax.plot_surface(x,y,to2dArray(U[2*nt/3],nx,ny), rstride=1, cstride=1, cmap=cm.jet)
	ax.set_zlim(0,10)
	show()
	
	fig4=figure(4)
	ax = Axes3D(fig4)
	ax.plot_surface(x,y,to2dArray(U[nt-1],nx,ny), rstride=1, cstride=1, cmap=cm.jet)
	ax.set_zlim(0,10)

	show()
