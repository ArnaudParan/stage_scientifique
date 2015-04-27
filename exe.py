#!/usr/bin/python3.4
#-*-coding:utf-8-*

from chaleur2D import*

tracer(euler2D(lx,ly,tmax,D,ax,ay,dx,dy,dt),lx,ly,tmax,dx,dy,dt,naff,Tmax)
