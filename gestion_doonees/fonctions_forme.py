#!/usr/bin/python3.4
#-*-coding:utf-8-*

"""
 @file fonctions_forme.py
 @brief tout ce qui concerne la manipulation des fonctions forme
 @author Arnaud Paran
 @version 0.1
 @date 6 mai 2015
"""

#ici, on s'intéresse au calcul des fonctions formes, de leur gradient, etc

"""
 @fn mult(x,l)
 @brief multiplie un vecteur par un scalaire
 @param x le vecteur
 @param l le scalaire
 @return le produit
"""
def mult(l,x):
	retour=[]
	for i in range(len(x)):
		retour.append(x[i]*l)
	return retour

"""
 @fn fforme_eval(gradient,x0,y0,x,y)
 @brief evalue une fonction forme en un point
 @param gradient le gradient de la fonction
 @param x0 abcisse du point de valeur 1
 @param y0 ordonnée du point de valeur 1
 @param x abscisse du point auquel on évalue
 @param y ordonnée du point auquel on évalue
 @return la valeur de la fonction au point x,y
"""
def fforme_eval(g,x0,y0,x,y):
	return 1.+g[0]*(x-x0)+g[1]*(y-y0)

"""
 @fn fforme_gradient(points, connectivite, i)
 @brief Renvoie les gradients des fonctions forme des points du simplex
 @param points matrice des points qui vient de donnees.py
 @param connectivite matrice de connectivite qui vient de donnees.py
 @param i nous dit la ligne de la matrice de connectivité à lire
 @return une matrice contenant les gradients
"""
def fforme_gradient(points, connectivite, i):
	#les abscisses et coordonnées des points
	x1=points[connectivite[i][0]][0]
	x2=points[connectivite[i][1]][0]
	x3=points[connectivite[i][2]][0]
	y1=points[connectivite[i][0]][1]
	y2=points[connectivite[i][1]][1]
	y3=points[connectivite[i][2]][1]
	#les normales non normalisées
	u1=[y3-y2,x2-x3]
	u2=[y3-y1,x1-x3]
	u3=[y1-y2,x2-x1]
	#les coefficients
	mu1=1./((x2-x3)*(y1-y3)-(x1-x3)*(y2-y3))
	mu2=1./((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))
	mu3=1./((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))
	#les gradients
	g1=mult(mu1,u1)
	g2=mult(mu2,u2)
	g3=mult(mu3,u3)

	return [g1,g2,g3]

