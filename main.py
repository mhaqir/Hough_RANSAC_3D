# Evaluating the implementation of Hough and RANSAC 3D algorithms 
# for plane detection
# Written by Mohammad Haghir Ebrahimabadi


from PIL import Image
from glob import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cv2 as cv


from util import project_image_to_rect, project_depth_to_points
# from objLoader import objloader
from hough3d import Hough3d
from ransac3d import RANSAC3d
from genPlane import genPlane
from util import compare_planes


params = {'plane_params': [2, 1, -3, 3], 'x_lim': [-1, 2], 'y_lim': [-1, 2], 'z_lim': [-1, 2], 'n': 6000, \
			'eps': 0.05, 'theta_step': 1, 'phi_step': 1, 'rho_step': 0.1, \
			't': 1, 'fs': 1000, 'alpha': 0.99}

class evaluation(genPlane, Hough3d, RANSAC3d):
	"""

	Input
		- a list of plane parameters for generating a data set
		- limits on the coordinate values
		- number of 3d points
		- limit for adding noise to 3d points


	Output
		- detected plane parameters
		- angle between the normal of original and detected planes 

	"""

	def __init__(self, **kwargs):
		for k, v in kwargs['kwargs'].items():
			setattr(self, k, v)
		genPlane.__init__(self, self.x_lim, self.y_lim, self.z_lim, self.n, self.eps)
		self.points = self.plane(self.plane_params)
		self.noisy_points = self.noisy_plane(self.points)
		Hough3d.__init__(self, self.points, self.theta_step, self.phi_step, self.rho_step)
		RANSAC3d.__init__(self, self.t, self.fs, self.alpha)

	def evaluate_Hough(self):
		params, conf = self.a_b_c_d()
		a, b, c, d = self.plane_params[0], self.plane_params[1], self.plane_params[2], self.plane_params[3]
		with open('Hough_a_{a}_b_{b}_c_{c}_d_{d}.txt'.format(a = a,b = b, c = c, d = d), 'w') as f:
			f.write('Output planes from 3D Hough for the original plane with parameters:\n')
			f.write('a = {a}, b = {b}, c = {c}, d = {d}\n'.format(a = a, b = b, c = c, d = d))
			f.write('*******************************************************************\n')
			f.write('a    b    c    d    parallel    angle/distance    conf\n')
			for r, abcd in params.items():
				for item in abcd:
					p, a_d = compare_planes(np.asarray(item), np.asarray([a, b, c, d]))
					if p == True:
						p_o = 'True'
						a_d_o = '0/{}'.format(round(a_d, 2))
					else:
						p_o = 'False'
						a_d_o = '{}/-'.format(round(a_d, 2))
					f.write('{a}	{b}    {c}    {d}    {p}    {a_d}    {conf}\n'.format(a = round(item[0], 2), b = round(item[1], 2),\
											 c = round(item[2], 2), d = round(item[3], 2), p = p_o, a_d = a_d_o, conf = conf[r]))


	def evaluate_RANSAC(self):
		bestPl, bestS, bestSt, N = self.find_plane(self.points)
		bestPl_n, bestS_n, bestSt_n, N = self.find_plane(self.noisy_points)
		p_a_d_o_nonn = compare_planes(self.plane_params, bestPl)
		p_a_d_o_n = compare_planes(self.plane_params, bestPl_n)


		a, b, c, d = self.plane_params[0], self.plane_params[1], self.plane_params[2], self.plane_params[3]
		with open('ransac_a_{a}_b_{b}_c_{c}_d_{d}_eps_{eps}.txt'.format(a = a, b = b, c = c, d = d, eps = self.eps), 'w') as f:
			f.write('Original plane: a = {a}, b = {b}, c = {c}, d = {d}\n'.format(a = a, b = b, c = c, d = d))
			f.write('\n')
			f.write('RANSAC parameters:\n')
			f.write('\n')
			f.write('t: {}\n'.format(self.t))
			f.write('fs: {}\n'.format(self.fs))
			f.write('alpha: {}\n'.format(self.alpha))
			f.write('*******************************************************************\n')
			f.write('Detected planes:\n')
			f.write('\n')
			# f.write('N for data without noise: {}\n'.format(ransac3d.N))
			f.write('Best plane: a = {a}, b = {b}, c = {c}, d = {d}\n'.format(a = round(bestPl[0], 4), b = round(bestPl[1], 4), c = round(bestPl[2], 4), d = round(bestPl[3], 4)))
			f.write('Best std: {}\n'.format(round(bestSt, 4)))
			f.write('Best support: {}\n'.format(round(bestS, 4)))
			f.write('\n')
			f.write('eps: {}\n'.format(self.eps))	
			# f.write('N for noisy data: {}\n'.format(ransac3d_noisy.N))
			f.write('Best plane: a = {a}, b = {b}, c = {c}, d = {d}\n'.format(a = round(bestPl_n[0], 4), b = round(bestPl_n[1], 4), c = round(bestPl_n[2], 4), d = round(bestPl_n[3], 4)))
			f.write('Best std: {}\n'.format(round(bestSt_n, 4)))
			f.write('Best support: {}\n'.format(round(bestS_n, 4)))
			f.write('*******************************************************************\n')
			f.write('Comparison:\n')
			f.write('			parallel	angle/distance\n')
			if p_a_d_o_nonn is not None:
				if p_a_d_o_nonn[0] == 'True':
					f.write('original vs. non-noisy    True        -/{}\n'.format(round(p_a_d_o_nonn[1], 4)))
				else:
					f.write('original vs. non-noisy    False        {}/-\n'.format(round(p_a_d_o_nonn[1], 4)))
			if p_a_d_o_n is not None:
				if p_a_d_o_n[0] == 'True':
					f.write('original vs. noisy    True        -/{}\n'.format(round(p_a_d_o_n[1], 4)))
				else:
					f.write('original vs. noisy        False        {}/-\n'.format(round(p_a_d_o_n[1], 4)))



e = evaluation(kwargs = params)
e.evaluate_Hough()
e.evaluate_RANSAC()

