
def AV24():

	mu = 1.0489 		    # Surface friction coefficient
	lf = 1.7238 			# Distance from center of gravity to front axle [m]
	lr = 1.248 			    # Distance from center of gravity to rear axle [m]
	tw_f = 1.638			# Track width of front [m]
	tw_r = 1.523			# Track width of rear [m]
	hcog = 0.275			# Height of center of gravity [m]
	haero = 0.15			# height of the location at which the equivalent aerodynamic force acts 
	mass = 787.2858		   	# Total mass of the vehicle [kg] 787.2858
	distribut_weight_kg = [154.22, 164.65, 208.65, 230.42] # FL, FR, RL, RR 
	Iz = 1000	   			# Moment of inertial of the entire vehicle about the z axis [kgm^2]
	width = 1.5815			# width of the vehicle [m]
	length = 4.921          # length of the vehicle [m]
	g = 9.81  				# Gravitational acceleration [m/s^2]
	p = 151685              # Tire pressure of left front and rear tires (Pa)

	min_steer = -0.279 	    # Minimum steering angle constraint [rad]
	max_steer = 0.279		# Maximum steering angle constraint [rad]
	max_steer_vel = 3.2 	# Maximum steering velocity constraint [rad/s]
	min_steer_vel = -3.2 	# Minimum steering velocity constraint [rad/s]
	min_v = -5.0			# Minimum longitudinal velocity [m/s]
	max_v = 89.408			# Maximum longitudinal velocity [m/s]
	# switch_v =7.319		# Switching velocity (velocity at which the acceleration is no longer able to create wheel spin) [m/s]
	max_acc = 9.51 			# max acceleration [m/s^2]
	min_acc = -13.26 		# max deceleration [m/s^2]
	roll_scale = 10.0
	camber = 0.0
	steer_offset = -2.5     # Steering offset [deg]
	steer_ratio  = 15
	#######################
	# Tire model parameters
	#######################
	tire_FL = {
		'PCx1' : 2.0,									
		'PDx1' : 1.7168, 
		'PDx2' : -0.289, 								 
		'PEx1' : -5.04761e-06, 
		'PKx1' : 63.75,
		'PKx3' : 0.2891,
		'LMUX' : 0.93,
		'PCy1' : 1.5496,
		'PDy1' : 1.70839,
		'PDy2' : -0.354062,
		'PEy1' : -1.98642,
		'PKy1' : -53.0946 ,
		'PKy2' : 2.23333,
		'LMUY' : 1,
		'N0': 3112
	}
	tire_FR = {
		'PCx1' : 1.53051,									
		'PDx1' : 1.30026 , 
		'PDx2' : -0.13349, 								 
		'PEx1' :-0.0104179, 
		'PKx1' : 45.8603,
		'PKx3' : -0.0713791,
		'LMUX' : 0.93,
		'PCy1' : 1.37872,
		'PDy1' : 1.1501,
		'PDy2' : -0.267779,
		'PEy1' : -1.39187,
		'PKy1' : -31.2787,
		'PKy2' : 2.17339,
		'LMUY' : 1,
		'N0': 6224
	}
	tire_RL = {
		'PCx1' : 1.60634,									
		'PDx1' : 1.57278, 
		'PDx2' : -0.307988, 								 
		'PEx1' : 0., 
		'PKx1' : 61.035,
		'PKx3' : 0.137978,
		'LMUX' : 0.93,
		'PCy1' : 1.53471,
		'PDy1' : 1.55128,
		'PDy2' : -0.388936,
		'PEy1' : -2.46857,
		'PKy1' : -44.4747,
		'PKy2' : 2.10251,
		'LMUY' : 1,
		'N0': 5347
	}
	tire_RR = {
		'PCx1' : 1.49561,									
		'PDx1' : 1.36527, 
		'PDx2' : -0.24548, 								 
		'PEx1' : 0., 
		'PKx1' : 46.5687,
		'PKx3' : 0.0272879,
		'LMUX' : 0.93,
		'PCy1' : 1.35624,
		'PDy1' : 1.14004,
		'PDy2' : -0.316978,
		'PEy1' : -1.20379,
		'PKy1' : -45.8825,
		'PKy2' : 3.20514,
		'LMUY' : 1,
		'N0': 7109
	}
	
	Long_params = {
		'mu': 0.4,        # friction coefficient of brake pad (typical range= 0.3 to 0.45)
		'Rcf': 0.163,     # frontbrake pad radius [m]
		'Rcr': 0.150,     # rear brake pad radius [m]
		'Rw': 0.3118,     # wheel radius [m]
		'thr_x': [0.0, 20.0, 40.0, 60.0, 100.0],
		'thr_y': [0.0, 12.0, 25.0, 70.0, 100.0],
		'trans_efficiency': 0.9938593859385938,
		'engine_poly': [-12.62710387327076, 0.1684030849164447, 1.0310252238667894e-05, -5.808995267727239e-09],
		'gear_ratios': [0, 2.9170, 1.875, 1.3810, 1.1150, 0.96, 0.8889],
		'final_drive': 3.0,
		'fr': 0.060796723472803696,
		'Af': 1.0,
		'Ap': 1.0,
		'Cd': 0.3,
		'Cl': 1.
	}

	MF_Lat = {
		'Bf': 17.197514461629723,
		'Cf': 1.6717081677023067,
		'Df': 2651.105946201609,
		'Ef': 0.1826, 
		'Br': 17.472250073164812,
		'Cr': 1.7576987654036749,
		'Dr': 3451.8729290576975,
		'Er': -0.0826,
		'Shfy': 0.007627728874956934,
		'Shry': 0.010287248963999114,
		'Svfy': 134.73955296770842,
		'Svry': 194.0794300197951,
	}
	MF_Lat_ksc = {
		'Bf': 17.197514461629723,
		'Cf': 1.6717081677023067,
		'Df': 2651.105946201609,
		'Ef': 0.1826, 
		'Br': 17.472250073164812,
		'Cr': 1.7576987654036749,
		'Dr': 3451.8729290576975,
		'Er': -0.0826,
		'Shfy': 0.0,
		'Shry': 0.0,
		'Svfy': 0.0,
		'Svry': 0.0,
	}

	max_inputs = [max_acc, max_steer]
	min_inputs = [min_acc, min_steer]

	max_rates = [None, max_steer_vel]
	min_rates = [None, min_steer_vel]	

	params_dict = {
		'lf': lf,
		'lr': lr,
		'mass': mass,
		'distribut_weight_kg':distribut_weight_kg,
		'Iz': Iz,
		'hcog': hcog,
		'haero':haero,
		'mu': mu,
		'min_v': min_v,
		'max_v': max_v,
		# 'switch_v': switch_v,
		'max_acc': max_acc,
		'min_acc': min_acc,
		'max_steer': max_steer,
		'min_steer': min_steer,
		'max_steer_vel': max_steer_vel,
		'min_steer_vel': min_steer_vel,
		'width': width,
		'length': length,
		'max_inputs': max_inputs,
		'min_inputs': min_inputs,
		'max_rates': max_rates,
		'min_rates': min_rates,
		'g': g,
		'l_pressure': p, 
		'long_params': Long_params,
		'MF_lat': MF_Lat,
		'MF_lat_ksc': MF_Lat_ksc,
		'tire_FL':tire_FL,
		'tire_FR':tire_FR,
		'tire_RL':tire_RL,
		'tire_RR':tire_RR,
		'roll_scale': roll_scale,
		'camber': camber,
		'steer_offset': steer_offset,
		'steer_ratio': steer_ratio,
		'tw_f': tw_f,
		'tw_r': tw_r
		}
	return params_dict
