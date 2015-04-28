#!/usr/bin/python3.4
#-*-coding:utf-8-*

from matplotlib import cm
import numpy as np
import pylab as pl
from pylab import*
from valeurs import*
from matrices import*
from valeurs2d import*
import matplotlib.pyplot as plt
import scipy as scipy
from scipy.sparse import csr_matrix
from scipy.sparse import*
from scipy.sparse.linalg.isolve.iterative import cg
from pylab import *
from mpl_toolkits.mplot3d import Axes3D



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

def euler2D(lx,ly,tmax,D,ax,ay,dx,dy,dt,options):
    diffusion=options[0]
    advection=options[1]
    schema_advx=options[2]
    schema_advy=options[3]
    #vérification des paramètres
    if advection and (dt>dx/abs(ax*10.) or dt>dy/abs(ay*10.)):
        print('Système mal paramétré, l\'algorithme procède à un reparamétrage')
        dt=min([dx/abs(10.*ax),dy/abs(10.*ay)])

    if (dt>dx/abs(ax*10.) or dt>dy/abs(ay*10.)):
        print('Système mal paramétré, l\'algorithme procède à un reparamétrage')
        dt=min([dx/abs(10.*ax),dy/abs(10.*ay)])

    #vérification de l'instabilité
    if advection and ax<0 and schema_advx=='down':
        print('Attention, votre schéma étant trop instable, nous en avons choisi un autre')
        schema_advx='up'
    if advection and ay<0 and schema_advy=='down':
        print('Attention, votre schéma étant trop instable, nous en avons choisi un autre')
        schema_advy='up'
    if advection and ax>0 and schema_advx=='up':
        print('Attention, votre schéma étant trop instable, nous en avons choisi un autre')
        schema_advx='down'
    if advection and ay>0 and schema_advy=='up':
        print('Attention, votre schéma étant trop instable, nous en avons choisi un autre')
        schema_advy='down'

    #grandeurs relatives à l'équation
    nuDx=D*dt/(dx*dx)
    nuDy=D*dt/(dy*dy)
    nuAx=ax*dt/dx
    nuAy=ay*dt/dy
    nx=1+int(lx/dx)
    ny=1+int(ly/dy)
    nt=1+int(tmax/dt)
    print('dx :'+str(dx)+' ; dy :'+str(dy)+' ; dt :'+str(dt))
    print('erreur :'+str(dx+dy+dt))
    if advection:
        print('facteur d\'instabilité a :('+str(dt*abs(ax)/dx)+','+str(dt*abs(ay)/dy)+')')
    if diffusion:
        print('facteur d\'instabilité D :('+str(dt*abs(2.*D)/(dx*dx))+','+str(dt*abs(2.*D)/(dy*dy))+')')

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
        A=I(nx,ny)+nuDx*Bx(nx,ny)+nuDy*By(nx,ny)+nuAx*Cx(nx,ny,schema_advx)+nuAy*Cy(nx,ny,schema_advy)
    elif diffusion:
        A=I(nx,ny)+nuDx*Bx(nx,ny)+nuDy*By(nx,ny)
    elif advection:
        A=I(nx,ny)+nuAx*Cx(nx,ny,schema_advx)+nuAy*Cy(nx,ny,schema_advy)
    elif test:
        A=I(nx,ny)+s1*nuAx/dx*H11(nx,ny)+s2*nuAy/dy*H22(nx,ny)+(2*nuAx/dx+s1*nuAy/dy)*H12(nx,ny)

    #on implémente le schéma d'euler
    for n in range(1,nt):
        Un=scipy.sparse.linalg.isolve.iterative.bicg(A,U[n-1])[0]
        #Un=linalg.solve(A,U[n-1])
        initBord(Un,Ubord,nx,ny)
        U.append(Un)
    
    return [U,lx,ly,tmax,dx,dy,dt]

def tracer(sol,n=naff, T=Tmax):
    U=sol[0]
    lx=sol[1]
    ly=sol[2]
    tmax=sol[3]
    dx=sol[4]
    dy=sol[5]
    dt=sol[6]
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
