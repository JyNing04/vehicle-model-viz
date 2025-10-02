"""	
Dynamic single track model.
x is a 6x1 vector: [x, y, vx, vy, psi, ω, delta]^T
u is a 2x1 vector: [delta_v, ax]^T
""" 

import numpy as np
import math
from params import vehicle
import os
import sys
import time


class Dynamic:
	def __init__(self, params = None, tire_params = None):
		if params == None:
			vehicle_params = vehicle.AV24()
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
		self.switch_v = 4.0
		self.camber = vehicle_params['camber']
		self.roll_scale = vehicle_params['roll_scale']
		self.steer_ratio = vehicle_params['steer_ratio']
		self.steer_offset = vehicle_params['steer_offset'] * (math.pi / 180) * (1 / self.steer_ratio)
		# Physical constraints
		self.max_acc     = vehicle_params['max_acc']
		self.min_acc     = vehicle_params['min_acc']
		self.max_steer   = vehicle_params['max_steer']
		self.min_steer   = vehicle_params['min_steer']
		self.max_steer_v = vehicle_params['max_steer_vel']
		self.min_steer_v = vehicle_params['min_steer_vel']
		self.width       = vehicle_params['width']
		self.length      = vehicle_params['length']
		# self.MF_long     = vehicle_params['MF_long']
		self.MF_lat	 	 = vehicle_params['MF_lat']
		self.tire_params  = False
		if tire_params:
			self.tire_params   = True
			self.MF_lat['Bf'] = tire_params['Bf']
			self.MF_lat['Cf'] = tire_params['Cf']
			self.MF_lat['Df'] = tire_params['Df']
			self.MF_lat['Ef'] = tire_params['Ef']
			self.MF_lat['Br'] = tire_params['Br']
			self.MF_lat['Cr'] = tire_params['Cr']
			self.MF_lat['Dr'] = tire_params['Dr']
			self.MF_lat['Er'] = tire_params['Er']
			self.MF_lat['Shfy']= tire_params['Shfy']
			self.MF_lat['Shry']= tire_params['Shry']
			self.MF_lat['Svfy']= tire_params['Svfy']
			self.MF_lat['Svry']= tire_params['Svry']
		'''
		self.max_inputs  = vehicle_params['max_inputs']
		self.min_inputs  = vehicle_params['min_inputs']
		self.max_rates   = vehicle_params['max_rates']
		self.min_rates   = vehicle_params['min_rates']
		'''
		self.n_states    = 7
		self.n_inputs    = 2
		# Force coefficients
		self.rho   = 1.225							     # Mass density of air [kg/m^3]
		self.Af    = 1.6 + 0.00056 * (self.mass - 765)   # Frontal area  (ref: Book VDC Eq 4.3)
		self.Cd    = 0.725 								 # Aerodynamic drag coefficient: 0.7 to 1.1 typical values for Formula Once car
		self.theta = 0.21								 # Banked curve	
		# self.usemd = 111
		# Tire model

	def sim_continuous(self, x0, u, wheel_v, t, roll=None):
		"""	simulates the nonlinear continuous model with given input vector
			by numerical integration using 4th order Runge Kutta method
		"""
		n_steps    = len(t)-1 # 1
		x          = np.zeros([n_steps+1, self.n_states]) # size: 2x7
		dxdt       = np.zeros([n_steps+1, self.n_states]) # size: 2x7
		dxdt[0, :] = self.derivative_eqs(None, x0, wheel_v, u[:,0], roll)
		x[0, :]    = x0
		for ids in range(1, n_steps+1):
			x[ids, :]  = self.odeintRK4(x[ids-1, :], [t[ids-1],t[ids]], u[:,ids-1], wheel_v, roll)
			# if using odeintRK6 method: x[ids, :]  = self.odeintRK6(x[ids-1, :], t[ids], u)
			dxdt[ids, :] = self.derivative_eqs(None, x[ids, :], u[:,ids-1], wheel_v, roll)
		
		return x, dxdt
	
	def sim_continuous_multistep(self, x0, u, t, step, step_ls, wheel_v = None, roll=None):
		"""	simulates the nonlinear continuous model with given input vector
			by numerical integration using Runge Kutta method
			propagating sim for n steps
		"""
		x       = np.zeros([step+1, self.n_states]) # size: (step+1)x6
		dxdt    = np.zeros([step+1, self.n_states]) # size: (step+1)x6
		x[0, :] = x0
		for ids in range(1, step+1):
			x[ids, :]  = self.odeintRK4(x[ids-1, :], [t[0],t[1]], u[ids-1, :], wheel_v[ids-1], roll[ids-1])
			dxdt[ids, :] = self.derivative_eqs(None, x[ids-1, :], u[ids-1, :], wheel_v[ids-1], roll[ids-1])
			step_ls.append(step)
		return x, dxdt, step_ls

	def calc_forces(self, x, u, wheel_v, roll):
		acc     = u[-1] 
		vx      = x[2] 
		vy      = x[3]
		psi     = x[4]
		ω       = x[5] 
		delta   = x[6] - self.steer_offset
		# Normal forces on tires
		froll  = 0.029					# rolling resistance coefficient f varies in the range 0.01 to 0.04. (Wong, 2001)
		'''
		Rolling resistance force = Cr * Normal tire force
		Cr = coefficient of rolling resistance,
		Normal tire force = Weight of the vehicle * (g + a * sin(θ)) / number of tires
		g = gravitational acceleration (9.8 m/s^2 on the surface of the Earth)
		a = lateral acceleration due to the centripetal force
		θ = banking angle
		'''
		Fznorm = self.mass * (self.g + acc * math.cos(-abs(roll))) #/ 4 # This should be for one axle 0.42 front, 0.58 rear
		Rx     = froll * Fznorm
		#######################################
		# Calculate forces on tires
		#######################################
		MF_params_lat = self.MF_lat
		Bf = MF_params_lat['Bf']
		Cf = MF_params_lat['Cf']
		Df = MF_params_lat['Df']
		Ef = MF_params_lat['Ef']
		Br = MF_params_lat['Br']
		Cr = MF_params_lat['Cr']
		Dr = MF_params_lat['Dr']
		Er = MF_params_lat['Er']
		Shfy = MF_params_lat['Shfy']
		Shry = MF_params_lat['Shry']
		Svfy = MF_params_lat['Svfy']
		Svry = MF_params_lat['Svry']
		# MF_params_long = self.MF_long
		# long_B = MF_params_long['B']
		# long_C = MF_params_long['C']
		# long_D = MF_params_long['D']
		# long_E = MF_params_long['E']
		# Slip angles & slip ratio
		kappa, alpha_f, alpha_r = self.calc_slips(x, u, wheel_v)
		alpha_f += Shfy
		alpha_r += Shry
		# Lateral tire forces
		Fyf = Svfy + Df * math.sin(Cf * math.atan(Bf * alpha_f - Ef * (Bf * alpha_f - math.atan(Bf * alpha_f))))
		Fyr = Svry + Dr * math.sin(Cr * math.atan(Br * alpha_r - Er * (Br * alpha_r - math.atan(Br * alpha_r))))
		if self.tire_params:
			Fyf = Fyf
			Fyr = Fyr
		# Fyf = Svfy + Df * math.sin(Cf * math.atan(Bf * alpha_f))
		# Fyr = Svry + Dr * math.sin(Cr * math.atan(Br * alpha_r))
		# Longitudinal tire forces
		# Fxr = long_D * np.sin(long_C * np.arctan(long_B * kappa - long_E * (long_B * kappa - np.arctan(long_B * kappa))))
		Fcx = self.mass * vy * ω 						 		 # Centripetal force
		Fbx = self.mass * self.g * math.sin(math.radians(roll)) * math.sin(psi)
		# Aerodynamic drag force
		Faero = 1/2 * self.rho * self.Cd * self.Af * vx ** 2
		# Longitudinal total forces
		Fxr = self.mass * acc
		# Lateral force total
		Fby = self.mass * self.g * math.sin(-abs(roll+self.camber))       # Lateral forces due to banked curve
		Fcy = self.mass * vx * ω							  # Centripetal force
		Fyl = Fyr + Fyf * math.cos(delta) - Fby - Fcy
		return Fyr, Fyf, Fby, Fcy, Fxr

	def calc_slips(self, x, u, wheel_v):
		acc   = u[-1] 
		delta = x[-1] - self.steer_offset
		vx    = x[2] 
		vy    = x[3]
		ω     = x[5] 
		# Adding a small constant (e.g., 1e-6) to avoid division by zero
		epsilon = 1e-6
		wheel_v += epsilon
		ACCELERATION = True if acc > 0 else False
		# Longitudinal slip ratio is defined differently during acceleration and braking
		if ACCELERATION:
			kappa = (wheel_v - vx) / wheel_v    			# Eq. 4.11
		else:
			kappa = (wheel_v - vx) / vx						# Eq. 4.12
		# Slip angles of front wheel and rear wheel
		β = math.atan(self.lr*np.tan(delta)/(self.lf + self.lr))
		theta_vf = math.atan((vy + self.lf * ω) / vx)          # Eq. 2.29  ; β + self.lf*ω/vx	math.atan((vy + self.lf * ω) / vx) 	
		theta_vr = math.atan((vy - self.lr * ω) /  vx)			# Eq. 2.30  ; β - self.lr*ω/vx   math.atan((vy - self.lr * ω) / vx)
		alpha_f  = delta - theta_vf					 		         # Eq. 2.23
		alpha_r  = -1 * theta_vr								    # Eq. 2.24
		if vx <= self.switch_v:
			alpha_f = 0.0
			alpha_r = 0.0
			kappa   = 0.0
		# print("front slip: {:.4f}, rear slip: {:.4f}".format(alpha_f, alpha_r))
		return kappa, alpha_f, alpha_r

	def derivative_eqs(self, t, x, u, wheel_v, roll):
		"""	
		write dynamics as first order ODE: dxdt = f(x(t))
		x is a 6x1 vector: [x, y, vx, vy, psi,ω]^T
		u is a 4x1 vector: [T, B, delta, ax]^T
		""" 
		delta_v = u[0]
		ax 		= u[1]
		vx      = x[2] 
		vy      = x[3]
		psi     = x[4]
		ω       = x[5] 
		delta   = x[6] - self.steer_offset

		Fyr, Fyf, Fby, Fcy, Fxr = self.calc_forces(x, u, wheel_v, roll)
		# Fby = self.mass*(vx*ω)*math.tan(roll)
		size    = len(x)
		dxdt    = np.zeros(size)
		dxdt[0] = vx * math.cos(psi) - vy * math.sin(psi)
		dxdt[1] = vx * math.sin(psi) + vy * math.cos(psi)
		dxdt[2] = 1./self.mass * (Fxr - Fyf * math.sin(delta)) + vy * ω
		dxdt[3] = 1./self.mass * (Fyr + Fyf * math.cos(delta) - Fby*vx/self.roll_scale - Fcy)
		dxdt[4] = ω 
		dxdt[5] = 1/self.Iz * (self.lf * Fyf * math.cos(delta) - self.lr * Fyr ) 
		if vx <= self.switch_v: # switch to ekin model 
			if vx <= -0.01:
				ax = 0
			dxdt[0] = vx * math.cos(psi) - vy * math.sin(psi)
			dxdt[1] = vx * math.sin(psi) + vy * math.cos(psi)
			dxdt[2] = ax
			dxdt[3] = ax * np.tan(delta)
			dxdt[4] = ω 
			dxdt[5] = (1 / (self.lf + self.lr)) * ax * ω 
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
	def odeintRK4(self, y0, t, u, wheel_v, roll):
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
			K[0]   = h * fun(t[i], y0, u, wheel_v, roll)
			K[1]   = h * fun(C[0] * h, y0+A[0][0]*K[0], u, wheel_v, roll)
			K[2]   = h * fun(C[1] * h, y0+A[1][0]*K[0]+A[1][1]*K[1], u, wheel_v, roll)
			K[3]   = h * fun(C[2] * h, y0+A[2][0]*K[0]+A[2][1]*K[1]+A[2][2]*K[2], u, wheel_v, roll)
			y_next[i, :] = y0 + B@K
			y0 = y_next[i, :]
		return y_next

	# 3. Heun's method
	def odeintHeun2(self, y0, t, wheel_v, u, roll):
		'''
		x: [x, y, vx, vy, psi, delta, ω]^T
		u: [acc (or ΔT), Δdelta]^T
		'''
		A      = 2/3 							# dependence of the stages on derivatives
		B      = np.asarray([1/4, 3/4]) 		# quadrature  weights
		C      = 2/3 							# nodes weights within the step
		fun    = self.derivative_eqs
		y_next = np.zeros([len(t)-1, len(y0)])
		K      = np.zeros((len(B), len(B)))
		for i in range(len(t)-1):
			h = t[i+1] - t[i]
			K[0]   = h * fun(0, y0, wheel_v, u, roll)
			K[1]   = h * fun(C*h, y0*A*K[0], wheel_v, u, roll)
			y_next[i, :] = y0 + B@K
			y0 = y_next[i, :]

		return y_next
