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
from operations_vect import *
from donnees import *

g = 9.81
K = 1.

# @brief evalue une fonction forme en un point
# @param gradient le gradient de la fonction
# @param x0 abcisse du point de valeur 1
# @param y0 ordonnée du point de valeur 1
# @param x abscisse du point auquel on évalue
# @param y ordonnée du point auquel on évalue
# @return la valeur de la fonction au point x,y
def fforme_eval (g, x0, y0, x, y) :
	return (1. + g[0] * (x - x0) + g[1] * (y - y0))

##
# @brief Renvoie les gradients des fonctions forme des matrice_points du simplex
# @param matrice_points matrice des matrice_points qui vient de donnees.py
# @param matrice_connectivite matrice de matrice_connectivite qui vient de donnees.py
# @param simplex nous dit la ligne de la matrice de connectivité à lire soit l'indice du simplex
# @return une matrice contenant les gradients
def fforme_gradient (simplex) :
	#les abscisses et coordonnées des matrice_points
	x1 = matrice_points[matrice_connectivite [simplex] [0]] [0]
	x2 = matrice_points[matrice_connectivite [simplex] [1]] [0]
	x3 = matrice_points[matrice_connectivite [simplex] [2]] [0]
	y1 = matrice_points[matrice_connectivite [simplex] [0]] [1]
	y2 = matrice_points[matrice_connectivite [simplex] [1]] [1]
	y3 = matrice_points[matrice_connectivite [simplex] [2]] [1]
	#les normales non normalisées
	u1 = [y3 -y2, x2 - x3]
	u2 = [y3 -y1, x1 - x3]
	u3 = [y1 -y2, x2 - x1]
	#les coefficients
	mu1 =1. / ((x2 - x3) * (y1 - y3) - (x1 - x3) * (y2 - y3))
	mu2 =1. / ((x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3))
	mu3 =1. / ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))
	#les gradients
	g1 = mult (mu1, u1)
	g2 = mult (mu2, u2)
	g3 = mult (mu3, u3)

	return [g1, g2, g3]

##
# @brief vérifie les gradients calculés
# @param matrice_points la matrice de matrice_points définie dans donnees.py
# @param matrice_connectivite la matrice de connectivité définie dans donneer.py
# @return l'erreur du calcul des gradients
# on calcule la somme des gradients sur chaque simplex et on en prend la norme infinie, sachant que cette somme doit être nulle
def verif_gradient () :
	erreur = 0.
	for simplex in range (len (matrice_connectivite)) :
		gradients = fforme_gradient (simplex)
		g1 = gradients[0]
		g2 = gradients[1]
		g3 = gradients[2]
		erreur = max ([erreur, norminf (add (g1, g2, g3))])
	return erreur

##################################################

#Obsolète depuis qu'on a commencé à trier les simplex dès la prise de données, si vous voulez voir ce que faisait
#cette fonction, mettez tous les simplex dans la mémoire
def lieu_stat (matrice_points, matrice_connectivite) :
	m = 1.
	M = 0.
	vals = [20, 40, 60, 80, 100]
	Mm = [[M, m], [M, m], [M, m], [M, m], [M, m]]
	for seuil in range (len (vals)) :
		for simplex in range (len (matrice_connectivite)) :
			if matrice_points [matrice_connectivite[simplex][0]] [0] < vals[seuil] :
				m = min ([m] + div_h_u (simplex))
			if matrice_points [matrice_connectivite[simplex][0]] [0] > vals[len(vals)-seuil-1] :
				M = max ([M] + div_h_u (simplex))
		Mm[j][1] = m
		Mm[len (vals) - seuil - 1][0] = M
	return Mm

################################################### Calculs des termes par simplex

##
# @brief vérifie l'équation de la conservation de la masse sur un simplex
# @param matrice_points la matrice de matrice_points définie dans donnees.py
# @param matrice_connectivite la matrice de connectivité définie dans donneer.py
# @param i l'indice du simplex considéré
# @return les coordonnées de la divergence sur la base des fonctions forme
def div_h_u (simplex) :
	H = (h (simplex, 0) + h (simplex, 1) + h (simplex, 2)) / 3.
	#gradients
	g = fforme_gradient (simplex)
	#les coordonnées de la divergence dans la base des fonctions forme
	div = [0., 0., 0.]
	#la coordonnée selon la iè fonction forme
	for point in range(3) :
		for k in range(3) :
			#dérivée par rapport à x
			div[point] += h (simplex, k) * u (simplex, point) * g[k][0] + h (simplex, point) * u (simplex, k) * g[k][0]
			#dérivée par rapport à y
			div[point] += h (simplex, k) * v (simplex, point) * g[k][1] + h (simplex, point) * v (simplex, k) * g[k][1]
			#on divise par pointa moyenne des h
			div[point] *= 1. / H
        return div;

def moments (simplex) :
	g = fforme_gradient (simplex)
	#les coordonnées de la divergence dans la base des fonctions forme
	div = [[], [], []]
	#la coordonnée selon la iè fonction forme
	for l in range (3) :
                S0 = [0., 0.]
                Sf = [u(simplex, l) * Math.sqrt((u(simplex, l) ** 2) + (v(simplex, l) ** 2)) / (h(simplex, l) ** (4. / 3.)),
				v(simplex, l) * Math.sqrt ((u(simplex, l) ** 2) + (v(simplex, l) ** 2)) / (h(simplex, l) ** (4. / 3.))]
		for k in range (3) :
                        S0[0] += z (simplex, k) * g[k][0]
                        S0[1] += z (simplex, k) * g[k][1]
                        #S0[0] += (z (simplex, k) + h (simplex, k)) * g[k][0]
                        #S0[1] += (z (simplex, k) + h (simplex, k)) * g[k][1]
                div[l] = [S0, Sf]
        return div;

##
# @brief renvoie la matrice gradient de la vitesse du point pt du simplex
# @param matrice_points la matrice de matrice_points définie dans donnees.py
# @param matrice_connectivite la matrice de connectivité définie dans donneer.py
# @param simplex le simplex considéré
# @param pt le point du simplex (entre 0 et 2)
# @return renvoie la matrice gradient de la vitesse
def grad_vitesse (simplex) :
        gradient = [[0., 0.], [0., 0.]]
        grad_fforme = fforme_gradient (simplex)

        for cpt_point in range (3) :
                #du/dx
                gradient[0][0] += u (simplex, cpt_point) * grad_fforme[cpt_point][0]
                #du/dy
                gradient[0][1] += u (simplex, cpt_point) * grad_fforme[cpt_point][1]
                #dv/dx
                gradient[1][0] += v (simplex, cpt_point) * grad_fforme[cpt_point][0]
                #dv/dy
                gradient[1][1] += v (simplex, cpt_point) * grad_fforme[cpt_point][1]
        return csr_matrix (gradient)

def calcule_v_grad_vi (simplex) :
        G = grad_vitesse (simplex)
        V1 = [u (simplex, 0), v (simplex, 0)]
        V2 = [u (simplex, 1), v (simplex, 1)]
        V3 = [u (simplex, 2), v (simplex, 2)]
        return [G.dot (V1), G.dot (V2), G.dot (V3)]

def gradh (simplex) :
        gradient = [0., 0.]
        grad_fforme = fforme_gradient (simplex)
        for cpt_point in range (3) :
                gradient[0] += h (simplex, cpt_point) * grad_fforme[cpt_point][0]
                gradient[1] += h (simplex, cpt_point) * grad_fforme[cpt_point][1]
        return gradient

def calc_g_grad_hi (simplex) :
        G = gradh (simplex)
        return [mult (h (simplex, 0) * g, G), mult (h (simplex, 1) * g, G), mult(h (simplex, 2) * g, G)]

def calc_g_S0Sfi(simplex) :
        moment = moments (simplex)
        S0 = [moment[0][0], moment[1][0], moment[2][0]]
        Sf = [mult (1. / K ** 2, moment[0][1]), mult (1. / K ** 2, moment[1][1]), mult (1. / K ** 2, moment[2][1])]
        S0Sf = [somme (S0[0], Sf[0]),somme (S0[1], Sf[1]), somme (S0[2], Sf[2])]
        for point in range (len (S0Sf)) :
                mult (g, S0Sf[point])
        return S0Sf

def calcule_h_v_grad_vi(simplex) :
        hvgv = calcule_v_grad_vi (simplex)
        hvgv[0] *= h (simplex, 0)
        hvgv[1] *= h (simplex, 1)
        hvgv[2] *= h (simplex, 2)
        return hvgv

def calc_g_grad_h() :
        g_grad_h = []
        for simplex in range (len (matrice_connectivite)) :
                g_grad_h.append (calc_g_grad_hi (simplex))
        returng_grad_h 

def calcule_v_grad_v () :
        v_grad_v = []
        for simplex in range (len (matrice_connectivite)) :
                v_grad_v.append (calcule_v_grad_vi (simplex))
        return v_grad_v

def calcule_h_v_grad_v () :
        hvgv = calcule_v_grad_v ()
        for simplex in range (len (matrice_connectivite)) :
                hvgv[simplex][0] *= matrice_points [matrice_connectivite[simplex][0]] [4]
                hvgv[simplex][1] *= matrice_points [matrice_connectivite[simplex][1]] [4]
                hvgv[simplex][2] *= matrice_points [matrice_connectivite[simplex][2]] [4]
        return hvgv

################################ Traitement des calculs et affichage des résultats

##
# @brief trace un graphe montrant si la solution est stationnaire
# @param matrice_points la matrice de matrice_points définie dans donnees.py
# @param matrice_connectivite la matrice de connectivité définie dans donneer.py
def trace_stationnaire () :
        div = []
        for point in range (len (matrice_points)) :
                div.append ([])
        ax=a3.Axes3D (pl.figure ())
        ax.set_xlim3d (0,120)
        ax.set_ylim3d (-5,5)
        ax.set_zlim3d (0.3,10000)
	#ajout des termes
	for simplex in range (len (matrice_connectivite)) :
                div_elem = div_h_u (simplex)
		div[matrice_connectivite[simplex][0]].append (div_elem[0])
		div[matrice_connectivite[simplex][1]].append (div_elem[1])
		div[matrice_connectivite[simplex][2]].append (div_elem[2])
	#moyennage
        for point in range(len (matrice_points)) :
		div[point] = moy (div[point])
	#affichage
	for simplex in range (len (matrice_connectivite)) :
                sommets = []
                sommets.append ([matrice_points[matrice_connectivite[simplex][0]][0], matrice_points[matrice_connectivite[simplex][0]][1], div[matrice_connectivite[simplex][0]]])
                sommets.append ([matrice_points[matrice_connectivite[simplex][1]][0], matrice_points[matrice_connectivite[simplex][1]][1], div[matrice_connectivite[simplex][1]]])
                sommets.append ([matrice_points[matrice_connectivite[simplex][2]][0], matrice_points[matrice_connectivite[simplex][2]][1], div[matrice_connectivite[simplex][2]]])
                tri = a3.art3d.Poly3DCollection ([sommets])
                tri.set_color (matplotlib.colors.rgb2hex ([1., 0., 0.]))
                tri.set_edgecolor ('k')
                ax.add_collection3d (tri)
        pl.show ()

##
# @brief évalue la véracité de la conservation de la masse
# @param matrice_points la matrice de matrice_points définie dans donnees.py
# @param matrice_connectivite la matrice de connectivité définie dans donneer.py
# @return renvoie l'erreur
def equa_masse () :
	erreur = 0.
	for i in range (len (matrice_connectivite)) :
		erreur = max ([erreur, norminf (div_h_u (i))])
	return erreur

def exe () :
        A = []
        M = 0.
        for i in range (len (matrice_connectivite)) :
                c = calc_g_grad_hi (i)
                d = calc_g_S0Sfi (i)
                hvgv = calcule_h_v_grad_vi (i)
                A.append (somme (d[0], hvgv[0]))
                A.append (somme (d[1], hvgv[1]))
                A.append (somme (d[2], hvgv[2]))
        for i in range (len (A)) :
                M = max ([M] + A[i])
        return M

def reg_moments () :
        x1 = []
        y1 = []
        x2 = []
        y2 = []
        for point in range (len (matrice_points)) :
                x1.append ([])
                y1.append ([])
                x2.append ([])
                y2.append ([])
	#ajout des termes
	for simplex in range (len (matrice_connectivite)) :
                div_elem = moments (simplex)
                C_elem = calc_g_grad_hi (simplex)
                vgv_elem = calcule_v_grad_vi (simplex)
		y1[matrice_connectivite[simplex][0]].append(-div_elem[0][0][0]+C_elem[0][0]+vgv_elem[0][0])
		y1[matrice_connectivite[simplex][1]].append(-div_elem[1][0][0]+C_elem[1][0]+vgv_elem[1][0])
		y1[matrice_connectivite[simplex][2]].append(-div_elem[2][0][0]+C_elem[2][0]+vgv_elem[2][0])
		x1[matrice_connectivite[simplex][0]].append(div_elem[0][1][0])
		x1[matrice_connectivite[simplex][1]].append(div_elem[1][1][0])
		x1[matrice_connectivite[simplex][2]].append(div_elem[2][1][0])

		y2[matrice_connectivite[simplex][0]].append(-div_elem[0][0][1]+C_elem[0][1]+vgv_elem[0][1])
		y2[matrice_connectivite[simplex][1]].append(-div_elem[1][0][1]+C_elem[1][1]+vgv_elem[1][1])
		y2[matrice_connectivite[simplex][2]].append(-div_elem[2][0][1]+C_elem[2][1]+vgv_elem[2][1])
		x2[matrice_connectivite[simplex][0]].append(div_elem[0][1][1])
		x2[matrice_connectivite[simplex][1]].append(div_elem[1][1][1])
		x2[matrice_connectivite[simplex][2]].append(div_elem[2][1][1])
	#moyennage
        for point in range(len(matrice_points)) :
		x1[point]=moy(x1[point])
		y1[point]=moy(y1[point])
		x2[point]=moy(x2[point])
		y2[point]=moy(y2[point])
        print("selon x")
        regx=reg1D(x1,y1)[1]
        if regx < 0  :
                print("Erreur de signe")
        else  :
                print("K="+str(Math.sqrt(1./regx)))
        print("selon y")
        regy=reg1D(x2,y2)[1]
        if regy < 0  :
                print("Erreur de signe")
        else  :
                print("K="+str(Math.sqrt(1./regy)))
        print("selon x et y")
        regxpy=reg1D(x1+x2,y1+y2)[1]
        if regxpy<0 :
                print("Erreur de signe")
        else :
                print("K="+str(Math.sqrt(1./regxpy)))

def distrib_reg_moments () :
        x1 = []
        y1 = []
        x2 = []
        y2 = []
        for point in range (len (matrice_points)) :
                x1.append ([])
                y1.append ([])
                x2.append ([])
                y2.append ([])
	#ajout des termes
	for simplex in range (len (matrice_connectivite)) :
                C = calc_g_grad_hi (simplex)
                vgv = calcule_h_v_grad_vi (simplex)
                div_elem = moments (simplex)
                ligne_point = [         matrice_connectivite[simplex][0],
                                        matrice_connectivite[simplex][1],
                                        matrice_connectivite[simplex][2]]
		y1[ligne_point[0]].append(-div_elem[0][0][0]+C[0][0])
		y1[ligne_point[1]].append(-div_elem[1][0][0]+C[1][0])
		y1[ligne_point[2]].append(-div_elem[2][0][0]+C[2][0])
		x1[ligne_point[0]].append(div_elem[0][1][0])
		x1[ligne_point[1]].append(div_elem[1][1][0])
		x1[ligne_point[2]].append(div_elem[2][1][0])

		y2[ligne_point[0]].append(-div_elem[0][0][1]+C[0][1])
		y2[ligne_point[1]].append(-div_elem[1][0][1]+C[1][1])
		y2[ligne_point[2]].append(-div_elem[2][0][1]+C[2][1])
		x2[ligne_point[0]].append(div_elem[0][1][1])
		x2[ligne_point[1]].append(div_elem[1][1][1])
		x2[ligne_point[2]].append(div_elem[2][1][1])
	#moyennage
        for point in range(len(matrice_points)) :
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
