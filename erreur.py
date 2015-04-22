#!/usr/bin/python3.4
#-*-coding:utf-8-*

from chaleur1D import*

dxFin=lx/100.
dtFin=tmax/100.
UFin=euler1D(lx,tmax,D,dxFin,dtFin)


def erreur(Uapp,Uf,dx,dt,dxf,dtf):
	err=0.
	for i in range(0,len(Uapp)-1):
		for j in range(0,len(Uapp[i])-1):
			if(err< abs(Uf[int(i*dt/dtf)][int(j*dx/dxf)]-Uapp[i][j])):
				err= abs(Uf[int(i*dt/dtf)][int(j*dx/dxf)]-Uapp[i][j])
	return err

err=[]

fact=range(1,100,10)

for i in fact:
	dx=lx/i
	dt=tmax/i
	tracer(euler1D(lx,tmax,D,dx,dt),lx,tmax,dx,dt,naff)
	err.append(erreur(euler1D(lx,tmax,D,dx,dt),UFin,dx,dt,dxFin,dtFin))
	print(err[i/10-1])

print(err[1])
plt.plot(fact,err)
plt.show()
