# Implementation of RANSAC 3D algorithm for plane detection
# Written By Mohammad Haghir Ebrahimabadi


import numpy as np
from functools import reduce
from math import gcd
# from copy import deepcopy

from genPlane import genPlane

class RANSAC3d:

	"""
	3D RANSAC algorithm from the paper:
	https://www.isprs.org/proceedings/xxxvi/3-w52/final_papers/Tarsha-Kurdi_2007.pdf

	Input:
		- matrix of points in x-y-z space
		- threshold of distance t between the chosen plane and the other points
		- forseeble_support is the maximum probable number of points belonging
		  to the same plane
		- minimum probability of finding at least one good set of observations
		  in N trials --  0.9 < alpha < 0.99

	Output:
		- parameters of a plane?
	"""
	
	def __init__(self, t, fs, alpha):
		# self.points = points
		self.t = t
		self.alpha = alpha
		self.fs = fs
		# self.eps = 1 - fs/np.shape(points)[0]
		# self.N = int(np.log10(1 - alpha)/np.log10(1 - np.power(1 - self.eps, 3)))
		self.bestSupport = 0
		self.bestPlane = [0, 0, 0, 0]
		self.bestStd = float('inf')

	def pick3pts(self, points):
		rand_pts_ids = np.random.choice(np.shape(points)[0], 3, replace = False)
		return points[rand_pts_ids, :]

	def pts2plane(self, pts3):
		# pts3 = self.points[pts_row_id, :]
		A, B, C = pts3[0, :], pts3[1, :], pts3[2, :]
		v1, v2 = B - A, C - A 
		n = np.cross(v1, v2)
		# gd = reduce(gcd, n)
		# n = n / gd
		k = - n[0] * pts3[0, 0] - n[1] * pts3[0, 1] - n[2] * pts3[0, 2]
		return np.concatenate((n, np.asarray([k])), axis = 0)

	def dist2plane(self, pl, points):
		d1, d2, d3 = points[:, 0] * pl[0], points[:, 1] * pl[1],  points[:, 2] * pl[2]
		d = d1 + d2 + d3 + pl[3]
		return d

	def find_plane(self, points):
		eps = 1 - self.fs/np.shape(points)[0]
		N = int(np.log10(1 - self.alpha)/np.log10(1 - np.power(1 - eps, 3)))
		for i in range(N):
			rand_pts = self.pick3pts(points)
			pl = self.pts2plane(rand_pts)
			d = self.dist2plane(pl, points)
			s = d[np.where(d <= self.t)[0]]
			st = np.std(s)
			if (np.shape(s)[0] > self.bestSupport) or (np.shape(s)[0] == self.bestSupport and st < self.bestStd):
				self.bestSupport, self.bestPlane, self.bestStd = np.shape(s)[0], pl, st
		return self.bestPlane, self.bestSupport, self.bestStd, N



# a, b, c, d = 1, -2, 1, 5
# pl1 = genPlane(a, b, c, d, [-1, 2], [-1, 2], [-1, 2], 2000)
# points = pl1.plane()

# ransac3d = RANSAC3d(points, 2, 200, 0.95)

# print(ransac3d.bestStd > 5)
# print(ransac3d.eps)
# print(ransac3d.N)
# pts_row_id = ransac3d.pick3pts()
# pl = ransac3d.pts2plane(pts_row_id)
# d = ransac3d.dist2plane(pl)
# print(pl[3])
# print(np.shape(d))
# # print(reduce(gcd, [9, 3, -3]))
# pl, bestS, bestSt = ransac3d.find_plane()
# print(pl)
# print(bestSt)
# print(bestS)
# print(np.dot(np.asarray([1, -2, 1]), pl[0:3]))








