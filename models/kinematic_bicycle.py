"""	
Kinematic single track model.
x is a 5x1 vector: [x, y, psi, v, 𝛿]^T
u is a 2x1 vector: [acc (or ΔT),Δ𝛿]^T
""" 

import numpy as np
import math
from params import vehicle
import casadi as cs
import os
import sys
import time

class Kinematic:
	def __init__(self, params = None):
		if params == None:
			vehicle_params = vehicle.AV21()
		else:
			vehicle_params = params
		self.g      = vehicle_params['g']
		self.lf     = vehicle_params['lf']
		self.lr     = vehicle_params['lr']
		self.mass   = vehicle_params['mass']
		self.Iz     = vehicle_params['Iz']
		self.hcog   = vehicle_params['hcog']
		self.mu     = vehicle_params['mu']
		self.min_v  = vehicle_params['min_v']
		self.max_v  = vehicle_params['max_v']
		self.h_aero = vehicle_params['haero']
		self.tire_p = vehicle_params['l_pressure']
		# Physical constraints
		self.max_acc     = vehicle_params['max_acc']
		self.min_acc     = vehicle_params['min_acc']
		self.max_steer   = vehicle_params['max_steer']
		self.min_steer   = vehicle_params['min_steer']
		self.max_steer_v = vehicle_params['max_steer_vel']
		self.min_steer_v = vehicle_params['min_steer_vel']
		self.width       = vehicle_params['width']
		self.length      = vehicle_params['length']
		'''
		self.max_inputs  = vehicle_params['max_inputs']
		self.min_inputs  = vehicle_params['min_inputs']
		self.max_rates   = vehicle_params['max_rates']
		self.min_rates   = vehicle_params['min_rates']
		'''
		self.n_states    = 5
		self.n_inputs    = 2

	def sim_continuous(self, x0, u, t, wheel_v=None, roll=None):
		"""	simulates the nonlinear continuous model with given input vector
			by numerical integration using 6th order Runge Kutta method
			x0 is the initial state of size 5x1 [x, y, psi, v, 𝛿]^T
			u is the input vector of size 2x1 [acc (or ΔT),Δ𝛿]^T
			t is the time vector of size 2x1 [0, Ts]^T
		"""
		n_steps    = len(t)-1 # 1
		x          = np.zeros([n_steps+1, self.n_states]) # size: 2x5
		dxdt       = np.zeros([n_steps+1, self.n_states]) # size: 2x5
		dxdt[0, :] = self.derivative_eqs(None, x0, u[:,0])
		x[0, :]    = x0
		for ids in range(1, n_steps+1):
			x[ids, :]  = self.odeintRK4(x[ids-1, :], [t[ids-1],t[ids]], u[:,0])
			# if using odeintRK6 method: x[ids, :]  = self.odeintRK6(x[ids-1, :], t[ids], u)
			dxdt[ids, :] = self.derivative_eqs(None, x[ids, :], u[:,0])
		
		return x, dxdt

	def sim_continuous_multistep(self, x0, u, t, step, step_ls, wheel_v = None, roll=None):
		"""	simulates the nonlinear continuous model with given input vector
			by numerical integration using Runge Kutta method
			propagating sim for n steps
		"""
		x          = np.zeros([step+1, self.n_states]) # size: (step+1)x6
		dxdt       = np.zeros([step+1, self.n_states]) # size: (step+1)x6
		x[0, :]    = x0
		for ids in range(1, step+1):
			# x[ids, :]  = self.odeintScipy(x[ids-1, :], [t[0],t[1]], u[ids-1, :])
			x[ids, :]  = self.odeintRK4(x[ids-1, :], [t[0],t[1]], u[ids-1, :])
			dxdt[ids, :] = self.derivative_eqs(None, x[ids-1, :], u[ids-1, :])
			step_ls.append(step)
		return x, dxdt, step_ls

	def derivative_eqs(self, t, x, u):
		"""	
		write dynamics as first order ODE: dxdt = f(x(t))
		x is a 7x1 vector: [x, y, psi, v, 𝛿]^T
		u is a 2x1 vector: [acc(or T), Δ𝛿]^T
		""" 
		steer_v = u[1]
		acc     = u[0]
		psi     = x[2]
		v       = x[3] 
		𝛿       = x[4]
		β 		= np.arctan(self.lr*np.tan(𝛿)/(self.lf + self.lr))
		
		size    = len(x)
		dxdt    = np.zeros(size)
		dxdt[0] = v * np.cos(psi + β)
		dxdt[1] = v * np.sin(psi + β)
		dxdt[2] = np.sin(β) * v / self.lr
		dxdt[3] = acc
		dxdt[4] = steer_v 
		
		return dxdt
		
	# Numerical integration methods 
	# 1. Runge-Kutta 6th Order Method
	def odeintRK6(self, y0, h, u):
		A = np.asarray([[1/3],
						[0, 2/3],
						[1/12, 1/3, -1/12],
						[25/48, -55/24, 35/48, 15/8],
						[3/20, -11/20, -1/8, 1/2, 1/10],
						[-261/260, 33/13, 43/156, -118/39, 32/195, 80/39]], dtype=object) # dependence of the stages on derivatives
		B      = np.asarray([13/200, 0, 11/40, 11/40, 4/25, 4/25, 13/200])  # quadrature  weights
		C      = np.asarray([1/3, 2/3, 1/3, 5/6, 1/6, 1]) 					# nodes weights within the step
		fun    = self.derivative_eqs
		y_next = np.zeros([1, len(y0)])
		K      = np.zeros((len(B), len(B)))
		K[0]   = h * fun(0, y0, u)
		K[1]   = h * fun(C[0]*h, y0+A[0][0]*K[0], u)
		K[2]   = h * fun(C[1]*h, y0+A[1][0]*K[0]+A[1][1]*K[1], u)
		K[3]   = h * fun(C[2]*h, y0+A[2][0]*K[0]+A[2][1]*K[1]+A[2][2]*K[2], u)
		K[4]   = h * fun(C[3]*h, y0+A[3][0]*K[0]+A[3][1]*K[1]+A[3][2]*K[2]+A[3][3]*K[3], u)
		K[5]   = h * fun(C[4]*h, y0+A[4][0]*K[0]+A[4][1]*K[1]+A[4][2]*K[2]+A[4][3]*K[3]+A[4][4]*K[4], u)
		K[6]   = h * fun(C[5]*h, y0+A[5][0]*K[0]+A[5][1]*K[1]+A[5][2]*K[2]+A[5][3]*K[3]+A[5][4]*K[4]+A[5][5]*K[5], u)
		y_next = y0 + B@K
		
		return y_next

	# 2. 3/8-rule, variation of Runge-Kutta 4th Order Method
	def odeintRK4(self, y0, t, u):
		A   = np.asarray([ [1/3],
						   [-1/3, 1],
						   [1, -1, 1]], dtype=object)
		B   = np.asanyarray([1/8, 3/8, 3/8, 1/8])
		C   = np.asarray([1/3, 2/3, 1])
		fun = self.derivative_eqs
		y_next = np.zeros([len(t)-1, len(y0)])
		K      = np.zeros((len(B), len(y0)))
		for i in range(len(t)-1):
			h      = t[i+1] - t[i]
			K[0]   = h * fun(t[i], y0, u)
			K[1]   = h * fun(C[0] * h, y0+A[0][0]*K[0], u)
			K[2]   = h * fun(C[1] * h, y0+A[1][0]*K[0]+A[1][1]*K[1], u)
			K[3]   = h * fun(C[2] * h, y0+A[2][0]*K[0]+A[2][1]*K[1]+A[2][2]*K[2], u)
			y_next[i, :] = y0 + B@K
			y0 = y_next[i, :]
		return y_next

	# 3. Heun's method
	def odeintHeun2(self, y0, t, u):
		'''
		x: [x, y, vx, vy, phi, 𝛿, ω]^T
		u: [acc (or ΔT), Δ𝛿]^T
		'''
		A      = 2/3 							# dependence of the stages on derivatives
		B      = np.asarray([1/4, 3/4]) 		# quadrature  weights
		C      = 2/3 							# nodes weights within the step
		fun    = self.derivative_eqs
		y_next = np.zeros([len(t)-1, len(y0)])
		K      = np.zeros((len(B), len(B)))
		for i in range(len(t)-1):
			h = t[i+1] - t[i]
			K[0]   = h * fun(0, y0, u)
			K[1]   = h * fun(C*h, y0*A*K[0], u)
			y_next[i, :] = y0 + B@K
			y0 = y_next[i, :]

		return y_next

	def casadi(self, x, u, dxdt):
		"""	write dynamics as first order ODE: dxdt = f(x(t))
			x is a 5x1 vector: [x, y, psi, v, 𝛿]^T
			u is a 1x1 vector: Δ𝛿 (was [acc, Δ𝛿])
			dxdt is a casadi.SX variable
		"""
		steer_v = u   # was u[1]
		acc = 0       # was u[0]
		𝛿   = x[4] 
		v   = x[3]		
		psi = x[2]
		β   = np.arctan(self.lr*np.tan(𝛿)/(self.lf + self.lr))

		v = cs.if_else(v < self.min_v, self.min_v, v)
		𝛿 = cs.if_else(v < self.min_v, 0, 𝛿)
		
		dxdt[0] = v * cs.cos(psi + β)
		dxdt[1] = v * cs.sin(psi + β)
		dxdt[2] = np.sin(β) * v / self.lr
		dxdt[3] = acc
		dxdt[4] = steer_v
		return dxdt
