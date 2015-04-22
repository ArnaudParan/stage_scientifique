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
		return x
