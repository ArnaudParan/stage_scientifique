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

def mult(l,x):
	produit=[]
	for i in range(len(x)):
		produit.append(x[i]*l)
	return produit

def somme(l,x):
	somme=[]
	for i in range(len(x)):
		somme.append(x[i]+l)
	return somme

def reg2D(x,y,z):
	#les coefficients qui rentrent dans la matrice
	l=[0.,0.,0.,0.,0.,0.,0.,0.,0.]
	for i in range(len(x)):
		l[0]+=1
		l[1]+=x[i]
		l[2]+=y[i]
		l[3]+=x[i]*x[i]
		l[4]+=y[i]*y[i]
		l[5]+=x[i]*y[i]
		l[6]+=z[i]
		l[7]+=z[i]*x[i]
		l[8]+=z[i]*y[i]
	#La matrice associée au problème
	M=csr_matrix([[l[0],l[1],l[2]],[l[1],l[3],l[5]],[l[2],l[5],l[4]]])
	b=[l[6],l[7],l[8]]
	coeffs=scipy.sparse.linalg.isolve.iterative.bicg(M,b)[0]
	err_quad=0.
	for i in range(len(x)):
		err_quad+=(z[i]-coeffs[0]-coeffs[1]*x[i]-coeffs[2]*y[i])**2
	err_quad=Math.sqrt(err_quad/len(x))
	print('Erreur de la régression :'+str(err_quad))
	return coeffs

def reg1D(x,z):
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
	err_quad=0.
	for i in range(len(x)):
		err_quad+=(z[i]-coeffs[0]-coeffs[1]*x[i])**2
	err_quad=Math.sqrt(err_quad/len(x))
	print('Erreur de la régression :'+str(err_quad))
	print(coeffs)
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
	plt.plot(sorted(err),somme(0.5,mult(0.5,scipy.special.erf(mult(1./(Math.sqrt(2.)*err_quad),sorted(err))))),marker='v')
	plt.show()
