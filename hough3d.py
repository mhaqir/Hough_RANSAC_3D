# Implementation of Hough 3D algorithm for 3D plane detection
# Written By Mohammad Haghir Ebrahimabadi


import numpy as np
import cv2 as cv

class Hough3d:

	"""
	3D Hough transform from the paper:
	https://www.isprs.org/proceedings/xxxvi/3-w52/final_papers/Tarsha-Kurdi_2007.pdf

	Input:
		- points in x-y-z space
		- theta_step
		- phi_step
		- rho_step

	Output:
		- matrix H
		- paramaters of the output planes based on the peaks in the matrix H
	"""
	
	def __init__(self, points, theta_step, phi_step, rho_step):
		self.points = points
		self.X_min, self.X_max = np.min(points[:, 0]), np.max(points[:, 0])
		self.Y_min, self.Y_max = np.min(points[:, 1]), np.max(points[:, 1])
		self.Z_min, self.Z_max = np.min(points[:, 2]), np.max(points[:, 2])
		self.n_theta = int(360/theta_step)
		self.n_phi = int(180/phi_step)
		self.Dis_min = np.sqrt(self.X_min**2 + self.Y_min**2 + self.Z_min**2)
		self.Dis_max = np.sqrt(self.X_max**2 + self.Y_max**2 + self.Z_max**2)
		print('self.Dis_min: ', self.Dis_min)
		print('self.Dis_max: ', self.Dis_max)
		self.n_rho = 2*int(np.abs(self.Dis_max - self.Dis_min)/rho_step)
		print('n_rho:', self.n_rho)
		self.theta = np.expand_dims(np.linspace(0, 360, self.n_theta) * np.pi/180, 1)
		self.phi = np.expand_dims(np.linspace(-90, 90, self.n_phi) * np.pi/180, 1)
		self.rho = np.expand_dims(np.linspace(self.Dis_min, self.Dis_max, self.n_rho), 1)
		print('rho:', self.rho)
		self.H = np.zeros((3*self.n_rho, self.n_theta, self.n_phi), dtype = np.uint16)
		print('H', self.H.shape)
		self.ratio = (self.n_rho - 1)/(self.rho[self.n_rho - 1] - self.rho[0])
		print('ratio: ', self.ratio)

	# def rho_index(self, q): # q is a point (x, y, z)
	# 	rho_mat = np.matmul(np.cos(self.theta), np.cos(self.phi.T)) * q[0] + \
	# 			np.matmul(np.cos(self.theta), np.sin(self.phi.T)) * q[1] + \
	# 			np.sin(self.phi.T) * q[2]
	# 	return np.rint(self.ratio*(rho_mat - self.rho[0] + 1)) #  + int(self.n_rho/2)

	def rho_index(self, q): # q is a point (x, y, z)
		rho_mat = np.matmul(np.cos(self.theta), np.sin(self.phi.T)) * q[0] + \
				np.matmul(np.sin(self.theta), np.sin(self.phi.T)) * q[1] + \
				np.cos(self.phi.T) * q[2]
		return np.rint(self.ratio*(rho_mat - self.rho[0] + 1)) #  + int(self.n_rho/2)

	# Given a data point, updates H for that
	def update_H(self, q):  
		temp_H = np.zeros_like(self.H)
		rho_indices = self.rho_index(q).astype(int)
		unique_rho_indices = np.unique(rho_indices)
		# print('# unique:', len(unique_rho_indices))
		# print(unique_rho_indices)
		for i in unique_rho_indices:
			w = np.where(rho_indices == i)
			temp_H[i, w[0], w[1]] = temp_H[i, w[0], w[1]] + 1
		return temp_H

	# Goes through all of the data points and builds H based on that
	def build_H(self):
		for i in range(np.shape(self.points)[0]):
			if i % 1000 == 0:
				print('point: ', i)
			temp_H = self.update_H(self.points[i,:])
			self.H = self.H + temp_H
		return self.H

	# Finds the peaks of the matrix H
	def H_max(self):
		H = self.build_H()
		H_m = np.zeros_like(H, dtype = np.uint8)
		H_m_values = []
		for i in range(np.shape(H)[0]):
			h = np.zeros_like(H[i, :, :], dtype = np.uint8)
			m = np.max(H[i, :, :])
			w = np.where(H[i, :, :] == m)
			h[w[0], w[1]] = 1
			H_m[i, :, :] = h
			H_m_values.append(m)
		return H_m, H_m_values

	# Finds rho, theta, and phi for the peaks
	def rho_theta_phi(self):
		H_m, H_m_values = self.H_max()
		r_t_p = {}
		conf = {}
		for i in range(np.shape(H_m)[0]):
			contours, hierarchy = cv.findContours(H_m[i, :, :], cv.RETR_TREE, cv.CHAIN_APPROX_SIMPLE)
			lst_coordinates = []
			if len(contours) > 0:
				for j in range(len(contours)):
					cimg = np.zeros_like(H_m[i, :, :], dtype = np.uint8)
					cv.drawContours(cimg, contours, j, color = 255, thickness = -1)
					pts = np.where(cimg == 255)
					lst_coordinates.append((np.mean(pts[0]), np.mean(pts[1])))
				r_t_p[i] = lst_coordinates
				conf[i] = H_m_values[i]
		return r_t_p, conf

	# Converts the peaks from phi-theta-rho space to their 
	# corresponding a, b, c, d parameters in x-y-z space 
	def a_b_c_d(self):
		r_t_p, conf = self.rho_theta_phi()
		a_b_c_d_values = {}
		for r, lst_t_p in r_t_p.items():
			lst_a_b_c_d = []
			for t_p in lst_t_p:
				a = np.cos(t_p[0]*np.pi/180)*np.sin(t_p[1]*np.pi/180)# -np.cos(t_p[0]*np.pi/180)*np.cos(t_p[1]*np.pi/180)/np.sin(t_p[1]*np.pi/180) 
				b = np.sin(t_p[0]*np.pi/180)*np.sin(t_p[1]*np.pi/180)# -np.sin(t_p[0]*np.pi/180)*np.cos(t_p[1]*np.pi/180)/np.sin(t_p[1]*np.pi/180)
				c = np.cos(t_p[1]*np.pi/180)#(r/self.ratio + self.rho[0] - 1)/np.sin(t_p[1]*np.pi/180)
				d = - r/self.ratio - self.rho[0] + 1
				lst_a_b_c_d.append([a, b, c, d[0]])
			a_b_c_d_values[r] = lst_a_b_c_d 
		return a_b_c_d_values, conf

	# Same implementation ** without confidence score **
	# # Finds the peaks of the matrix H
	# def H_max(self):
	# 	H = self.build_H()
	# 	H_m = np.zeros_like(H, dtype = np.uint8)
	# 	for i in range(np.shape(H)[0]):
	# 		h = np.zeros_like(H[i, :, :], dtype = np.uint8)
	# 		w = np.where(H[i, :, :] == np.max(H[i, :, :]))
	# 		h[w[0], w[1]] = 1
	# 		H_m[i, :, :] = h
	# 	return H_m

	# # Finds rho, theta, and phi for the peaks
	# def rho_theta_phi(self):
	# 	H_m = self.H_max()
	# 	r_t_p = {}
	# 	for i in range(np.shape(H_m)[0]):
	# 		contours, hierarchy = cv.findContours(H_m[i, :, :], cv.RETR_TREE, cv.CHAIN_APPROX_SIMPLE)
	# 		lst_intensities = []
	# 		if len(contours) > 0:
	# 			for j in range(len(contours)):
	# 				cimg = np.zeros_like(H_m[i, :, :], dtype = np.uint8)
	# 				cv.drawContours(cimg, contours, j, color = 255, thickness = -1)
	# 				pts = np.where(cimg == 255)
	# 				lst_intensities.append((np.mean(pts[0]), np.mean(pts[1])))
	# 			r_t_p[i] = lst_intensities
	# 	return r_t_p, H_m

	# Converts the peaks from phi-theta-rho space to their 
	# corresponding a, b, c, d parameters in x-y-z space 
	# def a_b_c_d(self):
	# 	r_t_p, H_m = self.rho_theta_phi()
	# 	a_b_c_d_values = {}
	# 	for r, lst_t_p in r_t_p.items():
	# 		lst_a_b_c_d = []
	# 		for t_p in lst_t_p:
	# 			a = np.cos(t_p[0]*np.pi/180)*np.cos(t_p[1]*np.pi/180)# -np.cos(t_p[0]*np.pi/180)*np.cos(t_p[1]*np.pi/180)/np.sin(t_p[1]*np.pi/180) 
	# 			b = np.sin(t_p[0]*np.pi/180)*np.cos(t_p[1]*np.pi/180)# -np.sin(t_p[0]*np.pi/180)*np.cos(t_p[1]*np.pi/180)/np.sin(t_p[1]*np.pi/180)
	# 			c = np.sin(t_p[1]*np.pi/180)#(r/self.ratio + self.rho[0] - 1)/np.sin(t_p[1]*np.pi/180)
	# 			d = - r/self.ratio - self.rho[0] + 1
	# 			lst_a_b_c_d.append([a, b, c, d[0]])
	# 		a_b_c_d_values[r] = lst_a_b_c_d
	# 	return a_b_c_d_values












