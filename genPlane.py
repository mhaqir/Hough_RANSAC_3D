# Implementation of a class for generating 3D planes
# Written By Mohammad Haghir Ebrahimabadi


import numpy as np
from itertools import product, count, islice


class genPlane:

	"""
	This class is for generating 3d points on a plane.

	Input:
		- plane parameters: a, b, c, d
		- limits on the coordinate values
		- number of 3d points
		- limit for adding noise to 3d points

	Output:
		- 3d points
		- noisy 3d points
	"""

	def __init__(self, x_lim, y_lim, z_lim, n, eps):
		self.x_lim = x_lim
		self.y_lim = y_lim
		self.z_lim = z_lim
		self.n = n
		self.eps = eps

	# from https://stackoverflow.com/questions/28057307/factoring-a-number-into-roughly-equal-factors
	def factors(self):
		i = int(self.n**0.5 + 0.5)
		while self.n % i != 0:
			i -= 1
		return i, self.n/i

	# from https://stackoverflow.com/questions/4114167/checking-if-a-number-is-a-prime-number-in-python
	def is_prime(self):
		return self.n > 1 and all(self.n % i for i in islice(count(2), int(np.sqrt(self.n) - 1)))

	# def plane(self):
	# 	if self.is_prime():
	# 		print('n is a prime number!')
	# 		# X_Y = []
	# 	else:
	# 		Y_n, X_n = self.factors()
	# 		X = np.linspace(self.x_lim[0], self.x_lim[1], int(X_n))
	# 		Y = np.linspace(self.y_lim[0], self.y_lim[1], int(Y_n))
	# 		X_Y = np.round(np.asarray(list(product(X, Y))), 4)
	# 		Z = np.expand_dims(-(self.a * X_Y[:, 0] + self.b * X_Y[:, 1] + self.d)/self.c, 1)
	# 		points = np.concatenate((X_Y, Z), axis = 1)
	# 		return points

	def draw_2d_points(self, n1, n2, lim1, lim2):
		return np.linspace(lim1[0], lim1[1], n1), np.linspace(lim2[0], lim2[1], n2)

	# def third_dimension(self, a1, a2, a3, lim1, lim2):
	# 	P, Q = self.draw_points(int(n1), int(n2), self.x_lim, self.z_lim)
	# 	X_Z = np.round(np.asarray(list(product(X, Z))), 4)
	# 	Y = np.expand_dims(-(self.a * X_Z[:, 0] + self.c * X_Z[:, 1] + self.d)/self.b, 1)

	def plane(self, plane_params):
		points = None
		a, b, c, d = plane_params[0], plane_params[1], plane_params[2], plane_params[3] 
		if self.is_prime():
			print('n is a prime number!')
			# X_Y = []
		else:
			n1, n2 = self.factors()
			if c != 0:
				X, Y = self.draw_2d_points(int(n1), int(n2), self.x_lim, self.y_lim)
				X_Y = np.round(np.asarray(list(product(X, Y))), 4)
				Z = np.expand_dims(-(a * X_Y[:, 0] + b * X_Y[:, 1] + d)/c, 1)
				points = np.concatenate((X_Y, Z), axis = 1)
			elif b != 0:
				X, Z = self.draw_2d_points(int(n1), int(n2), self.x_lim, self.z_lim)
				X_Z = np.round(np.asarray(list(product(X, Z))), 4)
				Y = np.expand_dims(-(a * X_Z[:, 0] + c * X_Z[:, 1] + d)/b, 1)
				points = np.concatenate((np.expand_dims(X_Z[:, 0], axis = 1), Y, \
					np.expand_dims(X_Z[:, 1], axis = 1)), axis = 1)
			elif a != 0:
				Y, Z = self.draw_2d_points(int(n1), int(n2), self.y_lim, self.z_lim)
				Y_Z = np.round(np.asarray(list(product(Y, Z))), 4)
				X = np.expand_dims(-(b * Y_Z[:, 0] + c * Y_Z[:, 1] + d)/a, 1)
				points = np.concatenate((X, np.expand_dims(Y_Z[:, 0], axis = 1), \
					             np.expand_dims(Y_Z[:, 1], axis = 1)), axis = 1)
			# else:
			# 	print('normal is a zero-vector!')
		if points is not None:
			return points
		else:
			print('normal is a zero-vector!')

	def noisy_plane(self, points):
		# points = self.plane()
		noise_values = np.random.uniform(-self.eps, self.eps, (self.n, 3))
		# print(noise_values)
		# print('min: ', np.min(noise_values))
		# print('max: ', np.max(noise_values))
		if points is not None:
			return points + noise_values
		else:
			print('normal is a zero-vector!')





