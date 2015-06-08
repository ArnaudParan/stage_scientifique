#!/usr/bin/python3.4
#-*-coding:utf-8-*

##
# @brief crée la moyenne d'un tableau
# @param tab le tableau qu'on moyenne
# @return la moyenne
def moy (tab) :
	taille_tab = len (tab)
	moyenne = 0.
	for elem in tab :
		moyenne += elem / taille_tab
	return moyenne

##	
# @brief multiplie un vecteur par un scalaire
# @param vect le vecteur
# @param scal le scalaire
# @return le produit
def mult (scal, vect) :
	produit = []
	for i in range (len (vect)) :
		produit.append (vect[i] * scal)
	return produit

def somme (vect1, vect2) :
        somme = vect1
        for cpt in range (len(somme)) :
                somme[cpt] += vect2[cpt]
        return somme

##	
# @brief additionne trois vecteurs
# @param x1 le premier vecteur
# @param x2 le deuxième vecteur
# @param x3 le troisième vecteur
# @return la somme
def add (vect1 ,vect2 ,vect3) :
	somme = []
	#envoie une evectception si les vecteurs n'ont pas la meme taille
	if len (vect1) != len (vect2) or len (vect2) != len (vect3) or len (vect3) != len (vect1) :
		raise "trois vecteurs de taille différente additionnés"
	for i in range (len (vect1)) :
		somme.append (vect1[i] + vect2[i] + vect3[i])
	return somme

##	
# @brief la norme infinie d'un vecteur
# @param x le vecteur
# @return la norme
def norminf (x) :
	sup = max (x)
	inf = min (x)
	return max ([abs (sup), abs (inf)])

