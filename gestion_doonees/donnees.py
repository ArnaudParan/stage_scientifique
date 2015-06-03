#!/usr/bin/python3.4
#-*-coding:utf-8-*

##
# @file
# @brief la prise des données de data_lanzoni.dat
# @author Arnaud Paran
# @version 0.1
# @date 7 mai 2015

import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from matplotlib import collections  as mc

from os import chdir

##
# @brief lit une ligne de points du fichier de données
# @param fichier le fichier dont on lit les données
# @return un tableau de float représentant le point
def lit_ligne_pts(fichier):
        ligne=fichier.readline().split()
        retour=[]
        for donnee in ligne:
                retour.append(float(donnee))
        return retour

##
# @brief lit une ligne de la matrice de connectivité dans le .dat
# @param fichier le fichier dont on lit les données
# @return un tableau d'int représentant les indices des points du simplex
def lit_ligne_connect(fichier):
        ligne=fichier.readline().split()
        retour=[]
        for donnee in ligne:
                retour.append(int(float(donnee))-1)
        return retour

chdir('/home/arnaud/Documents/stage_scientifique/gestion_doonees/')

#Cet objet contien le fichier à ouvrir
Fdata=open('data_lanzoni.dat','r')

#positionne le curseur au début des points
for i in range(14):
        Fdata.readline()

#lit les points
points=[]
connectivite=[]
segments=[]
corresp=[]

#début de la zone qu'on considère
L=60

point=0
for i in range(7500):
	ligne=lit_ligne_pts(Fdata)
	if ligne[0]>=L:
		points.append(ligne)
		corresp.append(point)
		point=point+1
	else:
		corresp.append(-1)


#lit la connectivité
for i in range(13884):
	ligne=lit_ligne_connect(Fdata)
	ligne[0]=corresp[ligne[0]]
	ligne[1]=corresp[ligne[1]]
	ligne[2]=corresp[ligne[2]]
	if ligne[0]!=-1 and ligne[1]!=-1 and ligne[2]!=-1:
		connectivite.append(ligne)
		segments.append([(points[ligne[0]][0],points[ligne[0]][1]),(points[ligne[1]][0],points[ligne[1]][1])])
		segments.append([(points[ligne[1]][0],points[ligne[1]][1]),(points[ligne[2]][0],points[ligne[2]][1])])
		segments.append([(points[ligne[2]][0],points[ligne[2]][1]),(points[ligne[0]][0],points[ligne[0]][1])])

Fdata.close()


lc = mc.LineCollection(segments, linewidths=1)
fig, ax = pl.subplots()
ax.add_collection(lc)
ax.margins(0.1)
ax.set_autoscaley_on(False)
ax.set_ylim([-60,60])
plt.show()
