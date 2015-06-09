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
K = 30.


################################ Traitement des calculs et affichage des résultats

##
# @brief trace un graphe montrant si la solution est stationnaire
# @param matrice_points la matrice de matrice_points définie dans donnees.py
# @param matrice_connectivite la matrice de connectivité définie dans donneer.py
def trace_stationnaire () :
        div = []
        for point in range (len (matrice_points)) :
                div.append ([])
	#ajout des termes
	for simplex in range (len (matrice_connectivite)) :
                div_elem = div_h_u (simplex)
		div[matrice_connectivite[simplex][0]].append (div_elem[0])
		div[matrice_connectivite[simplex][1]].append (div_elem[1])
		div[matrice_connectivite[simplex][2]].append (div_elem[2])
	#moyennage
        moy_par_point(div)
	#affichage
        trace_surface(div)

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
        moments = equa_moments()
        x1 = moments[0]
        y1 = moments[1]
        x2 = moments[2]
        y2 = moments[3]
        Mx = 0.
        My = 0.
        for point in range(len(x1)) :
                Mx = max ([Mx] + [abs(x1[point] - (y1[point] / (K ** 2)))])
                My = max ([My] + [abs(x2[point] - (y2[point] / (K ** 2)))])
        return [Mx, My]

def reg_moments () :
        moments = equa_moments()
        x1 = moments[0]
        y1 = moments[1]
        x2 = moments[2]
        y2 = moments[3]
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
        moments = equa_moments()
        x1 = moments[0]
        y1 = moments[1]
        x2 = moments[2]
        y2 = moments[3]
        print("x")
        distrib_reg1D(x1,y1)
        print("y")
        distrib_reg1D(x2,y2)
        print("x+y")
        distrib_reg1D(x1+x2,y1+y2)

#####################partie graphique

def trace_surface(profondeurPoints) :
        ax=a3.Axes3D (pl.figure ())
        ax.set_xlim3d (0,120)
        ax.set_ylim3d (-5,5)
        ax.set_zlim3d (0.3,10000)
	for simplex in range(len(matrice_connectivite)) :
                sommets = []
                sommets.append( [x(simplex, 0), y(simplex, 0), profondeurPoints[ligne_point(simplex, 0)]])
                sommets.append( [x(simplex, 1), y(simplex, 1), profondeurPoints[ligne_point(simplex, 1)]])
                sommets.append( [x(simplex, 2), y(simplex, 2), profondeurPoints[ligne_point(simplex, 2)]])
                tri = a3.art3d.Poly3DCollection ([sommets])
                tri.set_color (matplotlib.colors.rgb2hex ([1., 0., 0.]))
                tri.set_edgecolor ('k')
                ax.add_collection3d (tri)
        pl.show ()

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
			div[point] += h(simplex, k) * u(simplex, point) * g[k][0] + h(simplex, point) * u(simplex, k) * g[k][0]
			#dérivée par rapport à y
			div[point] += h(simplex, k) * v(simplex, point) * g[k][1] + h(simplex, point) * v(simplex, k) * g[k][1]
			#on divise par pointa moyenne des h
			div[point] *= 1. / H
        return div;

def norm_vit(simplex, point) :
        return Math.sqrt((u(simplex, point) ** 2) + (v(simplex, point) ** 2))

def calcule_S0(simplex) :
	g = fforme_gradient (simplex)
        S0 = [[0., 0.], [0., 0.], [0., 0.]]
	#la coordonnée selon la iè fonction forme
	for point in range (3) :
		for k in range (3) :
                        S0[point][0] += z (simplex, k) * g[k][0]
                        S0[point][1] += z (simplex, k) * g[k][1]
        return S0;

def calcule_Sf(simplex) :
	g = fforme_gradient (simplex)
        Sf=[[], [], []]
	#la coordonnée selon la iè fonction forme
	for point in range (3) :
                Sf[point] = [u(simplex, point) * norm_vit(simplex, point) / (h(simplex, point) ** (4. / 3.)),
			v(simplex, point) * norm_vit(simplex, point) / (h(simplex, point) ** (4. / 3.))]
        return Sf;

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

def div_u(simplex) :
        grad_u = grad_vitesse(simplex)
        return grad_u[0][0] + grad_u[1][1]

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
        return [mult(g, G),
                mult(g, G),
                mult(g, G)]

def calc_g_S0Sfi(simplex) :
        S0 = calcule_S0(simplex)
        Sf = calcule_Sf(simplex)
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

##################### Opérations générales sur tout le maillage

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

def equa_moments() :
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
                S0 = calcule_S0(simplex)
                Sf = calcule_Sf(simplex)
                g_grad_h_elem = calc_g_grad_hi (simplex)
                vgv_elem = calcule_v_grad_vi (simplex)
                ligne_point = [matrice_connectivite[simplex][0],
                                matrice_connectivite[simplex][1],
                                matrice_connectivite[simplex][2]]
                for point in range(3) :
                        y1[ligne_point[point]].append(g * S0[point][0] + g_grad_h_elem[point][0] + vgv_elem[point][0])
                        x1[ligne_point[point]].append(-g * Sf[point][0])

                        y2[ligne_point[point]].append(g * S0[point][1] + g_grad_h_elem[point][1] + vgv_elem[point][1])
                        x2[ligne_point[point]].append(-g * Sf[point][1])
	#moyennage
        moy_par_point(x1)
        moy_par_point(y1)
        moy_par_point(x2)
        moy_par_point(y2)

        return [x1, y1, x2, y2]

########################## Opérations pour les fonctions forme

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
