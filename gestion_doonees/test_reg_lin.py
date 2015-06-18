#!/usr/bin/python3.4
#-*-coding:utf-8-*

import unittest
from reg_lin import*

class regLin_testCase(unittest.TestCase) :
	def test_reg2D(self) :
                a0 = 10.
                a1 = 30.
                a2 = 40.
		x = [1.26, 984., 2145., 0.26, 42.02]
		y = [546., 2., 115., 2.1013, 2035.]
		a1x = mult(a1, x)
		a2y = mult(a2, y)
		z = sommeScalVect(a0, somme(a1x, a2y))
		actualMexpected = somme([-a0, -a1, -a2], reg2D(x,y,z))
                error = norminf(actualMexpected)
		self.assertTrue(error <= 1., str(error))

	def test_reg1D(self) :
                a0 = 0.
                a1 = 30.
		x = [1.26, 984., 2145., 0.26, 42.02]
		z = sommeScalVect(a0, mult(a1, x))
		actualMexpected = somme([-a0, -a1], reg1D(x,z))
                error = norminf(actualMexpected)
		self.assertTrue(error <= 1., str(error))

if __name__ == "__main__" :
	unittest.main()
