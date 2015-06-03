#!/usr/bin/python3.4
#-*-coding:utf-8-*

##
# @file
# @brief tout ce qui concerne la manipulation des fonctions forme
# @author Arnaud Paran
# @version 0.1
# @date 6 mai 2015
#ici, on s'intéresse au calcul des fonctions formes, de leur gradient, etc
import mpl_toolkits.mplot3d as a3
import matplotlib.colors as colors
import pylab as pl
import scipy as sp
import math as Math
from reg_lin import *

g=9.81
K=1.

##
# @brief crée la moyenne d'un tableau
# @param tab le tableau qu'on moyenne
# @return la moyenne
def moy(tab):
	l=len(tab)
	moyenne=0.
	for elem in tab:
		moyenne+=elem/l
	return moyenne

##	
# @brief multiplie un vecteur par un scalaire
# @param x le vecteur
# @param l le scalaire
# @return le produit
def mult(l,x):
	produit=[]
	for i in range(len(x)):
		produit.append(x[i]*l)
	return produit

def somme (x1,x2):
        somme=x1
        for cpt in range(len(somme)):
                somme[cpt]+=x2[cpt]
        return somme

##	
# @brief additionne trois vecteurs
# @param x1 le premier vecteur
# @param x2 le deuxième vecteur
# @param x3 le troisième vecteur
# @return la somme
def add(x1,x2,x3):
	somme=[]
	#envoie une exception si les vecteurs n'ont pas la meme taille
	if len(x1)!=len(x2) or len(x2)!=len(x3) or len(x3)!=len(x1):
		raise "trois vecteurs de taille différente additionnés"
	for i in range(len(x1)):
		somme.append(x1[i]+x2[i]+x3[i])
	return somme

##	
# @brief la norme infinie d'un vecteur
# @param x le vecteur
# @return la norme
def norminf(x):
	sup=max(x)
	inf=min(x)
	return max([abs(sup),abs(inf)])

##
# @brief evalue une fonction forme en un point
# @param gradient le gradient de la fonction
# @param x0 abcisse du point de valeur 1
# @param y0 ordonnée du point de valeur 1
# @param x abscisse du point auquel on évalue
# @param y ordonnée du point auquel on évalue
# @return la valeur de la fonction au point x,y
def fforme_eval(g,x0,y0,x,y):
	return (1.+g[0]*(x-x0)+g[1]*(y-y0))

##
# @brief Renvoie les gradients des fonctions forme des points du simplex
# @param points matrice des points qui vient de donnees.py
# @param connectivite matrice de connectivite qui vient de donnees.py
# @param i nous dit la ligne de la matrice de connectivité à lire soit l'indice du simplex
# @return une matrice contenant les gradients
def fforme_gradient(points,connectivite,i):
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

##
# @brief vérifie les gradients calculés
# @param points la matrice de points définie dans donnees.py
# @param connectivite la matrice de connectivité définie dans donneer.py
# @return l'erreur du calcul des gradients
# on calcule la somme des gradients sur chaque simplex et on en prend la norme infinie, sachant que cette somme doit être nulle
def verif_gradient(points, connectivite):
	erreur=0.
	for i in range(len(connectivite)):
		gradients=fforme_gradient(points,connectivite,i)
		g1=gradients[0]
		g2=gradients[1]
		g3=gradients[2]
		erreur = max([erreur,norminf(add(g1,g2,g3))])
	return erreur

##
# @brief vérifie l'équation de la conservation de la masse sur un simplex
# @param points la matrice de points définie dans donnees.py
# @param connectivite la matrice de connectivité définie dans donneer.py
# @param i l'indice du simplex considéré
# @return les coordonnées de la divergence sur la base des fonctions forme
def div_equa_masse(points, connectivite,i):
	#les paramètres des points du simplex
	#h
	h=[points[connectivite[i][0]][4]]
	h.append(points[connectivite[i][1]][4])
	h.append(points[connectivite[i][2]][4])
	H=(h[0]+h[1]+h[2])/3.
	#z
	z=[points[connectivite[i][0]][5]]
	z.append(points[connectivite[i][1]][5])
	z.append(points[connectivite[i][2]][5])
	#bottom
	b=[points[connectivite[i][0]][5]]
	b.append(points[connectivite[i][1]][5])
	b.append(points[connectivite[i][2]][5])
	#vitesse u
	u=[points[connectivite[i][0]][2]]
	u.append(points[connectivite[i][1]][2])
	u.append(points[connectivite[i][2]][2])
	#vitesse v
	v=[points[connectivite[i][0]][3]]
	v.append(points[connectivite[i][1]][3])
	v.append(points[connectivite[i][2]][3])
	#gradients
	g=fforme_gradient(points, connectivite, i)
	#les coordonnées de la divergence dans la base des fonctions forme
	div=[0.,0.,0.]
	#la coordonnée selon la iè fonction forme
	for l in range(3):
		for k in range(3):
			#dérivée par rapport à x
			div[l]+=h[k]*u[l]*g[k][0]+h[l]*u[k]*g[k][0]
			#dérivée par rapport à y
			div[l]+=h[k]*v[l]*g[k][1]+h[l]*v[k]*g[k][1]
			#on divise par la moyenne des h
			div[l]*=1./H
        return div;

##
# @brief trace un graphe montrant si la solution est stationnaire
# @param points la matrice de points définie dans donnees.py
# @param connectivite la matrice de connectivité définie dans donneer.py
def trace_stationnaire(points, connectivite):
        div=[]
        for point in range(len(points)):
                div.append([])
        ax=a3.Axes3D(pl.figure())
        ax.set_xlim3d(0,120)
        ax.set_ylim3d(-5,5)
        ax.set_zlim3d(0.3,10000)
	#ajout des termes
	for simplex in range(len(connectivite)):
                div_elem=div_equa_masse(points,connectivite,simplex)
		div[connectivite[simplex][0]].append(div_elem[0])
		div[connectivite[simplex][1]].append(div_elem[1])
		div[connectivite[simplex][2]].append(div_elem[2])
	#moyennage
        for point in range(len(points)):
		div[point]=moy(div[point])
	#affichage
	for simplex in range(len(connectivite)):
                sommets=[]
                sommets.append([points[connectivite[simplex][0]][0],points[connectivite[simplex][0]][1],div[connectivite[simplex][0]]])
                sommets.append([points[connectivite[simplex][1]][0],points[connectivite[simplex][1]][1],div[connectivite[simplex][1]]])
                sommets.append([points[connectivite[simplex][2]][0],points[connectivite[simplex][2]][1],div[connectivite[simplex][2]]])
                tri=a3.art3d.Poly3DCollection([sommets])
                tri.set_color(matplotlib.colors.rgb2hex([1.,0.,0.]))
                tri.set_edgecolor('k')
                ax.add_collection3d(tri)
        pl.show()

##
# @brief évalue la véracité de la conservation de la masse
# @param points la matrice de points définie dans donnees.py
# @param connectivite la matrice de connectivité définie dans donneer.py
# @return renvoie l'erreur
def equa_masse(points,connectivite):
	erreur=0.
	for i in range(len(connectivite)):
		erreur=max([erreur,norminf(div_equa_masse(points,connectivite,i))])
	return erreur


def lieu_stat(points,connectivite):
	m=1.
	M=0.
	vals=[20,40,60,80,100]
	Mm=[[M,m],[M,m],[M,m],[M,m],[M,m]]
	for j in range(len(vals)):
		for i in range(len(connectivite)):
			if points[connectivite[i][0]][0]<vals[j]:
				m=min([m]+div_equa_masse(points,connectivite,i))
			if points[connectivite[i][0]][0]>vals[len(vals)-j-1]:
				M=max([M]+div_equa_masse(points,connectivite,i))
		Mm[j][1]=m
		Mm[len(vals)-j-1][0]=M
	return Mm

def moments(points, connectivite,i):
	#les paramètres des points du simplex
	#h
	h=[points[connectivite[i][0]][4]]
	h.append(points[connectivite[i][1]][4])
	h.append(points[connectivite[i][2]][4])
	H=(h[0]+h[1]+h[2])/3.
	#z
	z=[points[connectivite[i][0]][5]]
	z.append(points[connectivite[i][1]][5])
	z.append(points[connectivite[i][2]][5])
	#bottom
	b=[points[connectivite[i][0]][5]]
	b.append(points[connectivite[i][1]][5])
	b.append(points[connectivite[i][2]][5])
	#vitesse u
	u=[points[connectivite[i][0]][2]]
	u.append(points[connectivite[i][1]][2])
	u.append(points[connectivite[i][2]][2])
	#vitesse v
	v=[points[connectivite[i][0]][3]]
	v.append(points[connectivite[i][1]][3])
	v.append(points[connectivite[i][2]][3])
	#gradients
	g=fforme_gradient(points, connectivite, i)
	#les coordonnées de la divergence dans la base des fonctions forme
	div=[[],[],[]]
	#la coordonnée selon la iè fonction forme
	for l in range(3):
                S0=[0.,0.]
                Sf=[u[l]*Math.sqrt((u[l]**2)+(v[l]**2))/(h[l]**(4./3.)),v[l]*Math.sqrt((u[l]**2)+(v[l]**2))/(h[l]**(4./3.))]
		for k in range(3):
                        S0[0]+=z[k]*g[k][0]
                        S0[1]+=z[k]*g[k][1]
                        #S0[0]+=(z[k]+h[k])*g[k][0]
                        #S0[1]+=(z[k]+h[k])*g[k][1]
                div[l]=[S0,Sf]
        return div;

def reg_moments(points, connectivite):
        x1=[]
        y1=[]
        x2=[]
        y2=[]
        for point in range(len(points)):
                x1.append([])
                y1.append([])
                x2.append([])
                y2.append([])
	#ajout des termes
        C=calc_C(points,connectivite)
        vgv=calcule_v_grad_v(points,connectivite)
	for simplex in range(len(connectivite)):
                div_elem=moments(points,connectivite,simplex)
		y1[connectivite[simplex][0]].append(-div_elem[0][0][0]+C[simplex][0][0]+vgv[simplex][0][0])
		y1[connectivite[simplex][1]].append(-div_elem[1][0][0]+C[simplex][1][0]+vgv[simplex][1][0])
		y1[connectivite[simplex][2]].append(-div_elem[2][0][0]+C[simplex][2][0]+vgv[simplex][2][0])
		x1[connectivite[simplex][0]].append(div_elem[0][1][0])
		x1[connectivite[simplex][1]].append(div_elem[1][1][0])
		x1[connectivite[simplex][2]].append(div_elem[2][1][0])

		y2[connectivite[simplex][0]].append(-div_elem[0][0][1]+C[simplex][0][1]+vgv[simplex][0][1])
		y2[connectivite[simplex][1]].append(-div_elem[1][0][1]+C[simplex][1][1]+vgv[simplex][1][1])
		y2[connectivite[simplex][2]].append(-div_elem[2][0][1]+C[simplex][2][1]+vgv[simplex][2][1])
		x2[connectivite[simplex][0]].append(div_elem[0][1][1])
		x2[connectivite[simplex][1]].append(div_elem[1][1][1])
		x2[connectivite[simplex][2]].append(div_elem[2][1][1])
	#moyennage
        for point in range(len(points)):
		x1[point]=moy(x1[point])
		y1[point]=moy(y1[point])
		x2[point]=moy(x2[point])
		y2[point]=moy(y2[point])
        print("selon x")
        regx=reg1D(x1,y1)[1]
        if regx < 0 :
                print("Erreur de signe")
        else :
                print("K="+str(Math.sqrt(1./regx)))
        print("selon y")
        regy=reg1D(x2,y2)[1]
        if regy < 0 :
                print("Erreur de signe")
        else :
                print("K="+str(Math.sqrt(1./regy)))
        print("selon x et y")
        regxpy=reg1D(x1+x2,y1+y2)[1]
        if regxpy<0:
                print("Erreur de signe")
        else:
                print("K="+str(Math.sqrt(1./regxpy)))

##
# @brief renvoie la matrice gradient de la vitesse du point pt du simplex
# @param points la matrice de points définie dans donnees.py
# @param connectivite la matrice de connectivité définie dans donneer.py
# @param simplex le simplex considéré
# @param pt le point du simplex (entre 0 et 2)
# @return renvoie la matrice gradient de la vitesse
def grad_vitesse(points, connectivite, simplex):
        gradient=[[0.,0.],[0.,0.]]
        grad_fforme=fforme_gradient(points,connectivite,simplex)
	#vitesse u
	u=[points[connectivite[simplex][0]][2]]
	u.append(points[connectivite[simplex][1]][2])
	u.append(points[connectivite[simplex][2]][2])
	#vitesse v
	v=[points[connectivite[simplex][0]][3]]
	v.append(points[connectivite[simplex][1]][3])
	v.append(points[connectivite[simplex][2]][3])
        for cpt_point in range(3):
                #du/dx
                gradient[0][0]+=u[cpt_point]*grad_fforme[cpt_point][0]
                #du/dy
                gradient[0][1]+=u[cpt_point]*grad_fforme[cpt_point][1]
                #dv/dx
                gradient[1][0]+=v[cpt_point]*grad_fforme[cpt_point][0]
                #dv/dy
                gradient[1][1]+=v[cpt_point]*grad_fforme[cpt_point][1]
        return csr_matrix(gradient)

def calcule_v_grad_v(points,connectivite):
        v_grad_v=[]
        for simplex in range(len(connectivite)):
                G=grad_vitesse(points,connectivite,simplex)
                V1=[points[connectivite[simplex][0]][2],points[connectivite[simplex][0]][3]]
                V2=[points[connectivite[simplex][1]][2],points[connectivite[simplex][1]][3]]
                V3=[points[connectivite[simplex][2]][2],points[connectivite[simplex][2]][3]]
                v_grad_v.append([G.dot(V1),G.dot(V2),G.dot(V3)])
        return v_grad_v

def calcule_h_v_grad_v(points,connectivite):
        hvgv=calcule_v_grad_v(points,connectivite)
        for simplex in range(len(connectivite)):
                hvgv[simplex][0]*=points[connectivite[simplex][0]][4]
                hvgv[simplex][1]*=points[connectivite[simplex][1]][4]
                hvgv[simplex][2]*=points[connectivite[simplex][2]][4]
        return hvgv

def gradh(points,connectivite,simplex):
        gradient=[0.,0.]
        grad_fforme=fforme_gradient(points,connectivite,simplex)
	#h
	h=[points[connectivite[simplex][0]][4]]
	h.append(points[connectivite[simplex][1]][4])
	h.append(points[connectivite[simplex][2]][4])
        for cpt_point in range(3):
                gradient[0]+=h[cpt_point]*grad_fforme[cpt_point][0]
                gradient[1]+=h[cpt_point]*grad_fforme[cpt_point][1]
        return gradient

def calc_Ci(points,connectivite,simplex):
        G=gradh(points,connectivite,simplex)
        h=[points[connectivite[simplex][0]][4]]
        h.append(points[connectivite[simplex][1]][4])
        h.append(points[connectivite[simplex][2]][4])
        return [mult(h[0]*g,G),mult(h[1]*g,G),mult(h[2]*g,G)]

def calc_C(points,connectivite):
        C=[]
        for simplex in range(len(connectivite)):
                C.append(calc_Ci(points,connectivite,simplex))
        return C

def calc_Di(points,connectivite,simplex):
        moment=moments(points,connectivite,simplex)
        S0=[moment[0][0],moment[1][0],moment[2][0]]
        Sf=[mult(1./K**2,moment[0][1]),mult(1./K**2,moment[1][1]),mult(1./K**2,moment[2][1])]
        S0Sf= [somme(S0[0],Sf[0]),somme(S0[1],Sf[1]),somme(S0[2],Sf[2])]
        for i in range(len(S0Sf)):
                mult(g*points[connectivite[simplex][i]][4],S0Sf[i])
        return S0Sf

def calcule_h_v_grad_vi(points,connectivite,simplex):
        G=grad_vitesse(points,connectivite,simplex)
        V1=[points[connectivite[simplex][0]][2],points[connectivite[simplex][0]][3]]
        V2=[points[connectivite[simplex][1]][2],points[connectivite[simplex][1]][3]]
        V3=[points[connectivite[simplex][2]][2],points[connectivite[simplex][2]][3]]
        hvgv = [G.dot(V1),G.dot(V2),G.dot(V3)]
        hvgv[0]*=points[connectivite[simplex][0]][4]
        hvgv[1]*=points[connectivite[simplex][1]][4]
        hvgv[2]*=points[connectivite[simplex][2]][4]
        return hvgv

def exe(points,connectivite):
        A=[]
        M=0.
        for i in range(len(connectivite)):
                c=calc_Ci(points,connectivite,i)
                d=calc_Di(points,connectivite,i)
                hvgv=calcule_h_v_grad_vi(points,connectivite,i)
                A.append(somme(d[0],hvgv[0]))
                A.append(somme(d[1],hvgv[1]))
                A.append(somme(d[2],hvgv[2]))
        for i in range(len(A)):
                M=max([M]+A[i])
        return M

def distrib_reg_moments(points, connectivite):
        x1=[]
        y1=[]
        x2=[]
        y2=[]
        for point in range(len(points)):
                x1.append([])
                y1.append([])
                x2.append([])
                y2.append([])
	#ajout des termes
        C=calc_C(points,connectivite)
        vgv=calcule_h_v_grad_v(points,connectivite)
	for simplex in range(len(connectivite)):
                div_elem=moments(points,connectivite,simplex)
		y1[connectivite[simplex][0]].append(-div_elem[0][0][0]+C[simplex][0][0])
		y1[connectivite[simplex][1]].append(-div_elem[1][0][0]+C[simplex][1][0])
		y1[connectivite[simplex][2]].append(-div_elem[2][0][0]+C[simplex][2][0])
		x1[connectivite[simplex][0]].append(div_elem[0][1][0])
		x1[connectivite[simplex][1]].append(div_elem[1][1][0])
		x1[connectivite[simplex][2]].append(div_elem[2][1][0])

		y2[connectivite[simplex][0]].append(-div_elem[0][0][1]+C[simplex][0][1])
		y2[connectivite[simplex][1]].append(-div_elem[1][0][1]+C[simplex][1][1])
		y2[connectivite[simplex][2]].append(-div_elem[2][0][1]+C[simplex][2][1])
		x2[connectivite[simplex][0]].append(div_elem[0][1][1])
		x2[connectivite[simplex][1]].append(div_elem[1][1][1])
		x2[connectivite[simplex][2]].append(div_elem[2][1][1])
	#moyennage
        for point in range(len(points)):
		x1[point]=moy(x1[point])
		y1[point]=moy(y1[point])
		x2[point]=moy(x2[point])
		y2[point]=moy(y2[point])
        print("x")
        distrib_reg1D(x1,y1)
        print("y")
        distrib_reg1D(x2,y2)
        print("x+y")
        distrib_reg1D(x1+x2,y1+y2)
