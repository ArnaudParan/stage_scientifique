#!/usr/bin/python3.4
#-*-coding:utf-8-*

import numpy as np
#Algorithme pris sur http://georgioudakis.com/blog/categories/python/cg.html
#Solve [A]{x} = {b} linear equation system with the conjugate gradient method
# http://en.wikipedia.org/wiki/Conjugate_gradient_method

def conjugate_gradient(A, b, x0, MAX_ITERATIONS = 100, TOLERANCE = 1.0e-10):
#   Initializations
	x = x0
	r0 = b - np.dot(A, x)
	p = r0

	#   Start iterations    
	for i in xrange(MAX_ITERATIONS):
		a = float(np.dot(r0.T, r0)/np.dot(np.dot(p.T, A), p))
		x = x + p*a
		ri = r0 - np.dot(A*a, p)
		#print i, np.linalg.norm(ri)

		if np.linalg.norm(ri) < TOLERANCE:
			return x

		b = float(np.dot(ri.T, ri)/np.dot(r0.T, r0))
		p = ri + b * p
		r0 = ri
	print("attention, matrice non inversée avec le gradient conjugué")
	return x


def GC(A,b,MaxIterations=100,erreur=1e-10):
	#préconditionnement
	x=b

	#gradient conjugué
	#gradient
	g=np.dot(A,x)-b
	gAncien=g
	#direction de descente
	h=-g
	for i in xrange(MaxIterations):
		#pas de descente
		rho=-float((np.dot(g,h)/np.dot(np.dot(A,h),h)))
		#calcul de xi+1
		x=x+rho*h
		g=np.dot(A,x)-b
		#gestion de l'erreur
		m_erreur=np.linalg.norm(g)
		if m_erreur < erreur:
			return x
		#calcul du futur pas de descente
		gamma=np.dot(g,g)/np.dot(gAncien,gAncien)
		d=-g+gamma*d
		gAncien=g
	print('Erreur, matrice non inversée')
	return x

