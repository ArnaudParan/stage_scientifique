#!/usr/bin/python3.4
#-*-coding:utf-8-*

from chaleur1D import*

dxFin=lx/100.
dtFin=tmax/100.
UFin=euler1D(lx,tmax,D,dxFin,dtFin)

def interpolx(u,x,t,dx):
	return u[t][int(x)]+(x-int(x))*(u[t][int(x)+1]-u[t][int(x)])

def interpolt(u,x,t,dt):
	return u[int(t)][x]+(t-int(t))*(u[int(t)+1][x]-u[int(t)][x])

def erreurx(Uapp,Uf,dx,dt,dxf,dtf):
	err=0.
	for i in range(0,len(Uapp)-1):
		for j in range(0,len(Uapp[i])-1):
			if(err< abs(interpolx(Uf,j*dx/dxf,i,dxf)-Uapp[i][j])):
				err= abs(interpolx(Uf,j*dx/dxf,i,dxf)-Uapp[i][j])
	return err

def erreurt(Uapp,Uf,dx,dt,dxf,dtf):
	err=0.
	for i in range(0,len(Uapp)-1):
		for j in range(0,len(Uapp[i])-1):
			if(err< abs(interpolt(Uf,j,i*dt/dtf,dxf)-Uapp[i][j])):
				err= abs(interpolt(Uf,j,i*dt/dtf,dxf)-Uapp[i][j])
	return err

errx=[]
errt=[]

fact=range(10,110,10)

for i in fact:
	dx=lx/i
	dt=dtFin
	#tracer(euler1D(lx,tmax,D,dx,dt),lx,tmax,dx,dt,naff)
	errx.append(erreurx(euler1D(lx,tmax,D,dx,dt),UFin,dx,dt,dxFin,dtFin))

for i in fact:
	dx=dxFin
	dt=tmax/i
	#tracer(euler1D(lx,tmax,D,dx,dt),lx,tmax,dx,dt,naff)
	errt.append(erreurt(euler1D(lx,tmax,D,dx,dt),UFin,dx,dt,dxFin,dtFin))

tracer(UFin,lx,tmax,dx,dt,naff)
plt.plot(fact,errx)
plt.title("Erreur en x")
plt.show()
plt.plot(fact,errt)
plt.title("Erreur en t")
plt.show()
