#!/usr/bin/python3.4
#-*-coding:utf-8-*

import numpy as np
import matplotlib.pyplot as plt

def B(n):
	B=np.diag(np.ones(n-1),1)
	B+=B.transpose()
	B+=(-2)*np.identity(n)
	return (-1)*B

def I(n):
	return np.identity(n)
