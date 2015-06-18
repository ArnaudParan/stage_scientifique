#!/usr/bin/python3.4
#-*-coding:utf-8-*

from operations_vect import*
import unittest

class opVect_testCase(unittest.TestCase) :
	def test_moy(self) :
		testTab = [1,2,3,4,5,6,7,8,9,5]
		self.assertEqual(moy(testTab), 5, 'incorrect tab moy' + str(moy(testTab)))

	def test_moy_par_point(self) :
		testTab = [[1,2,3,4,6,7,8,9],
				[2,2,4,4],
				[0]]
		moy_par_point(testTab)
		self.assertEqual(testTab, [5, 3, 0], 'incorrect moy_par_point')

	def test_mult(self) :
		testVect = [1.25, 8.25, 100.25]
		testScal = 4
		self.assertEqual(mult(testScal, testVect), [5, 33, 401])

	def test_somme(self) :
		testVect1 = [0, 1, 2, 3, 4, 6, 7, 8, 9, 10]
		testVect2 = [10, 9, 8, 7, 6, 4, 3, 2, 1, 0]
		self.assertEqual(somme(testVect1, testVect2), [10]*10)

	def test_add(self) :
		testVect1 = [3, 3, 3, 0, 0, 0, 9, 9, 9]
		testVect2 = [3, 3, 3, 6, 6, 6, 0, 0, 0]
		testVect3 = [3, 3, 3, 3, 3, 3, 0, 0, 0]
		self.assertEqual(add(testVect1, testVect2, testVect3), [9]*9)

	def test_norminf(self) :
		testTab = [-10, 5, 2, 100, -50, 30, -500]
		self.assertEqual(norminf(testTab), 500)

	def test_sommeScalVect(self) :
		testVect = [0, 1, 2, 3, 4, 6, 7, 8, 9, 10]
		testScal = 1
		self.assertEqual(sommeScalVect(testScal, testVect), [1, 2, 3, 4, 5, 7, 8, 9, 10, 11])

	def test_variance(self) :
		testVect = [5, 5, 5, 5, 5, 5]
		self.assertEqual(variance(testVect), 0)

if __name__ == "__main__" :
	unittest.main()
