#!/usr/bin/python3.4
#-*-coding:utf-8-*

import numpy as np
import matplotlib.pyplot as plt


class Tens4 :
	"""classe codant un tenseur d'ordre 4 afin de réaliser des
	opérations linéaires sur le champ eulérien à deux dimensions
	"""
	 def _init_(self, n):
	 	"""Tenseur pour des matrices d'ordre ixj
		"""
	 	self.n=n
		self.tab=np.array([[0]*n*n]*n*n)
	
	@staticmethod
	def ident(n):
		"""méthode statique qui renvoie Id(n,n)"""
		retour=Tens4(n,n)
		retour.tab=np.identity(n*n)
		return retour
	
	def _getitem(self,i,j,k,l):
		"""renvoie l'élément i,j,k,l du tenseur"""
		return self.tab[i+self.n*j,k+self.n*l]

	def _setitem(self,nbre,i,j,k,l):
		"""change l'élément i,j,k,l"""
		self.tab[i+self.n*j,k+self.n*l]=nbre
	
	def _add(self,autre)
		if autre.n==self.n
			somme=Tens4(self.n)
			for i in range (0,self.n)
				
