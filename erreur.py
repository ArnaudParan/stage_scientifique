#!/usr/bin/python3.4
#-*-coding:utf-8-*

from chaleur1D import*

#la solution finale qui sera la référence
dxFin=lx/100.
dtFin=tmax/100.
UFin=euler1D(lx,tmax,D,dxFin,dtFin)

#des fonctions d'interpolation 1D pour ne pas avoir de mauvaise surprise
def interpolx(u,x,t,dx):
	return u[t][int(x)]+(x-int(x))*(u[t][int(x)+1]-u[t][int(x)])

def interpolt(u,x,t,dt):
	return u[int(t)][x]+(t-int(t))*(u[int(t)+1][x]-u[int(t)][x])

#calcule l'erreur en espace quand les temps sont identiques
def erreurx(Uapp,Uf,dx,dt,dxf,dtf):
	err=0.
	for i in range(0,len(Uapp)-1):
		for j in range(0,len(Uapp[i])-1):
			if(err< abs(interpolx(Uf,j*dx/dxf,i,dxf)-Uapp[i][j])):
				err= abs(interpolx(Uf,j*dx/dxf,i,dxf)-Uapp[i][j])
	return err

#calcule l'erreur en temps quand les maillages sont identiques
def erreurt(Uapp,Uf,dx,dt,dxf,dtf):
	err=0.
	for i in range(0,len(Uapp)-1):
		for j in range(0,len(Uapp[i])-1):
			if(err< abs(interpolt(Uf,j,i*dt/dtf,dxf)-Uapp[i][j])):
				err= abs(interpolt(Uf,j,i*dt/dtf,dxf)-Uapp[i][j])
	return err

#listes dans lesquelles on va stocker les erreurs
errx=[]
errt=[]

#représente des delta qu'on va utiliser
fact=range(10,110,10)

#calcule l'erreur en x
for i in fact:
	dx=lx/i
	dt=dtFin
	#tracer(euler1D(lx,tmax,D,dx,dt),lx,tmax,dx,dt,naff)
	errx.append(erreurx(euler1D(lx,tmax,D,dx,dt),UFin,dx,dt,dxFin,dtFin))

#calcule l'erreur en t
for i in fact:
	dx=dxFin
	dt=tmax/i
	#tracer(euler1D(lx,tmax,D,dx,dt),lx,tmax,dx,dt,naff)
	errt.append(erreurt(euler1D(lx,tmax,D,dx,dt),UFin,dx,dt,dxFin,dtFin))

#graphe solution finale
tracer(UFin,lx,tmax,dx,dt,naff)
#graphe erreur en x
plt.plot(fact,errx)
plt.title("Erreur en x")
#graphe erreur en t
plt.show()
plt.plot(fact,errt)
plt.title("Erreur en t")
plt.show()
