#!/usr/bin/python3.4
#-*-coding:utf-8-*

##
# @file
# @brief la régression linéaire
# @author Arnaud Paran
# @version 0.1
# @date 18 mai 2015

import numpy as np
import pylab as pl
import scipy as scipy
from scipy import special
import math as Math
from pylab import *
from scipy.sparse import csr_matrix
from scipy.sparse.linalg.isolve.iterative import cg
import matplotlib.pyplot as plt
import numpy as np
from operations_vect import*

def reg2D(x,y,z):
	X = 0.
	Y = 0.
	XX = 0.
	YY = 0.
	XY = 0.
	Z = 0.
	ZX = 0.
	ZY = 0.
	for i in range(len(x)):
		X += x[i]
		Y += y[i]
		XX += x[i] * x[i]
		YY += y[i] * y[i]
		XY += x[i] * y[i]
		Z += z[i]
		ZX += z[i] * x[i]
		ZY += z[i] * y[i]
	#La matrice associée au problème
	M=csr_matrix([[len(x), X, Y],
			[X, XX, XY],
			[Y, XY, YY]])
	b=[Z, ZX, ZY]
	res = scipy.sparse.linalg.isolve.iterative.cg(M, b, tol = 1e-100)
	coeffs = res[0]
	print('conv achieved : ' + str(res[1]))
	err_quad=0.
	for i in range(len(x)):
		err_quad+=(z[i]-coeffs[0]-coeffs[1]*x[i]-coeffs[2]*y[i])**2
	err_quad=Math.sqrt(err_quad/len(x))
	print('Erreur de la régression :'+str(err_quad))
	return coeffs

def reg1D(x,z):
	X = 0.
	XX = 0.
	Z = 0.
	ZX = 0.
	for i in range(len(x)):
		X += x[i]
		XX += x[i] * x[i]
		Z += z[i]
		ZX += z[i] * x[i]
	#La matrice associée au problème
	M=csr_matrix([[len(x), X],
			[X, XX]])
	b=[Z, ZX]
	res = scipy.sparse.linalg.isolve.iterative.cg(M, b, tol = 1e-100)
	coeffs = res[0]
	print('conv achieved : ' + str(res[1]))
	err_quad=0.
	for i in range(len(x)):
		err_quad+=(z[i]-coeffs[0]-coeffs[1]*x[i])**2
	err_quad=Math.sqrt(err_quad/len(x))
	print('Erreur de la régression :'+str(err_quad))
	return coeffs

def distrib_reg1D(x,z):
	#les coefficients qui rentrent dans la matrice
	l=[0.,0.,0.,0.,0.,0.,0.,0.,0.]
	for i in range(len(x)):
		l[0]+=1
		l[1]+=x[i]
		l[3]+=x[i]*x[i]
		l[6]+=z[i]
		l[7]+=z[i]*x[i]
	#La matrice associée au problème
	M=csr_matrix([[l[0],l[1],l[2]],[l[1],l[3],l[5]],[l[2],l[5],l[4]]])
	b=[l[6],l[7],l[8]]
	coeffs=scipy.sparse.linalg.isolve.iterative.bicg(M,b)[0]
	err=[]
	for i in range(len(x)):
		err.append((z[i]-coeffs[0]-coeffs[1]*x[i]))
	err_quad=0.
	for i in range(len(x)):
		err_quad+=(z[i]-coeffs[0]-coeffs[1]*x[i])**2
	err_quad=err_quad/len(x)
	err_quad=Math.sqrt(err_quad)
	print(str(err_quad))
	plt.plot(sorted(err),mult(1./len(err),range(len(err))),marker='o')
	plt.plot(sorted(err),sommeScalVect(0.5,mult(0.5,scipy.special.erf(mult(1./(Math.sqrt(2.)*err_quad),sorted(err))))),marker='v')
	plt.show()
