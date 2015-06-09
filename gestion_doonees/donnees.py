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
# @brief lit une ligne de matrice_points du fichier de données
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
# @return un tableau d'int représentant les indices des matrice_points du simplex
def lit_ligne_connect(fichier):
        ligne=fichier.readline().split()
        retour=[]
        for donnee in ligne:
                retour.append(int(float(donnee))-1)
        return retour

chdir('/home/arnaud/Documents/stage_scientifique/gestion_doonees/')

#Cet objet contien le fichier à ouvrir
Fdata=open('data_lanzoni.dat','r')

#positionne le curseur au début des matrice_points
for i in range(14):
        Fdata.readline()

#lit les matrice_points
matrice_points=[]
matrice_connectivite=[]
segments=[]
corresp=[]

#début de la zone qu'on considère
L=60

point=0
for i in range(7500):
	ligne=lit_ligne_pts(Fdata)
	if ligne[0]>=L:
		matrice_points.append(ligne)
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
		matrice_connectivite.append(ligne)
		segments.append([(matrice_points[ligne[0]][0],matrice_points[ligne[0]][1]),(matrice_points[ligne[1]][0],matrice_points[ligne[1]][1])])
		segments.append([(matrice_points[ligne[1]][0],matrice_points[ligne[1]][1]),(matrice_points[ligne[2]][0],matrice_points[ligne[2]][1])])
		segments.append([(matrice_points[ligne[2]][0],matrice_points[ligne[2]][1]),(matrice_points[ligne[0]][0],matrice_points[ligne[0]][1])])

#Fdata.close()


#lc = mc.LineCollection(segments, linewidths=1)
#fig, ax = pl.subplots()
#ax.add_collection(lc)
#ax.margins(0.1)
#ax.set_autoscaley_on(False)
#ax.set_ylim([-60,60])
#plt.show()

def ligne_point(simplex, point) :
        return matrice_connectivite[simplex][point]

def x (simplex, point) :
        return matrice_points [matrice_connectivite[simplex][point]] [0]

def y (simplex, point) :
        return matrice_points [matrice_connectivite[simplex][point]] [1]

def u (simplex, point) :
        return matrice_points [matrice_connectivite[simplex][point]] [2]

def v (simplex, point) :
        return matrice_points [matrice_connectivite[simplex][point]] [3]

def h (simplex, point) :
        return matrice_points [matrice_connectivite[simplex][point]] [4]

def z (simplex, point) :
        return matrice_points [matrice_connectivite[simplex][point]] [5]
