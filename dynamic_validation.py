"""	Generate data by simulating (kinematic/dynamic) bicycle model.
"""
import math
import os
import sys
path_ = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../'))
sys.path.append(path_)
import numpy as np
from params import vehicle
from models.dynamic_bicycle import Dynamic
from models.kinematic_bicycle import Kinematic
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, Circle
from matplotlib.lines import Line2D
from matplotlib.animation import FuncAnimation
from collections import deque
import time

class Vehicle:
    def __init__(self, params):
        self.params = params
        self.true_states = None
        self.pred_states = None
        self.inputs      = None
        self.roll_data   = None
        self.wheel_data  = None
        self.states      = None

    def load_data(self, file_path, idx):
        data = np.loadtxt(file_path, delimiter=',', comments='#', max_rows=None if idx == -1 else idx)
        wheel_v = data[:,10:14] #* 0.277778
        states  = data[:,1:8]
        n_inputs = 2
        inputs  = np.zeros((data.shape[0], n_inputs))
        steer = data[:,7]
        inputs[1:,0] = np.diff(steer)
        inputs[:,1] = data[:,8]
        roll   = data[:,-1] #* -1	
        times  = data[:,0]
        print(states.shape)
        # x(m),y(m),vx(m/s),vy(m/s),phi(rad),delta(rad),omega(rad/s),ax(m/s^2),deltadelta(rad/s),wheel_fl/fr/ll/lr
        states = np.stack((states[:,0],states[:,1],states[:,2],states[:,3],states[:,4],states[:,5],states[:,6],wheel_v[:,0],wheel_v[:,1],wheel_v[:,2],wheel_v[:,3],roll),axis=1)
        return times, states, inputs 
    
    def gen_states(self, file_path, models):
        data_limit = -1
        times, states, self.inputs = self.load_data(file_path, idx = data_limit)
        # ground truth vehicle states
        self.true_states = states[:, :7]   # omega(rad/s)
        # model prediction vehicle states
        self.pred_states = np.zeros(self.true_states.shape)
        self.pred_states[0, :] = self.true_states[0, :]
        self.wheel_data = states[:, -5:-1]
        self.roll_data  = states[:, -1]
        self.states = states

class Simulation:
    def __init__(self, vehicle, model, states=None, inputs=None):
        self.vehicle = vehicle
        self.SAMPLING_TIME = 0.04
        self.states = vehicle.states
        self.inputs = inputs if inputs is not None else vehicle.inputs
        self.true_states = states if states is not None else vehicle.true_states
        self.N_SAMPLES = len(self.inputs)
        self.states_pred = np.zeros(self.true_states.shape)
        self.wheel_data = vehicle.wheel_data
        # if model == 'single_track':
        #     self.wheel_data = np.mean(vehicle.wheel_data, axis=1, keepdims=True)
        self.roll_data  = vehicle.roll_data
        self.dxdt_ek = np.zeros(self.true_states.shape)
        self.base_horizon = 10
        self.steer_inputs = self.true_states[:, 6]
        self.veh_model = model

    def calc_inrange(self, arr, a, b):
        lower_diff = arr - a
        lower_diff[lower_diff < 0] = 0

        upper_diff = b - arr
        upper_diff[upper_diff<0] = 0
        in_range_diff = np.minimum(lower_diff, upper_diff)
        return in_range_diff
    
    # Function to calculate out-range difference
    def calc_steer_diff(self, steering, std_dev):
        # Check if the current steering angle is out of range based on standard deviation
        if abs(steering) > std_dev:
            # Calculate the out range difference between the current steering angle and standard deviation
            out_range_difference = abs(steering) - std_dev
            in_range_difference  = 0
        else:
            # Calculate the in range difference between the current steering angle and standard deviation
            in_range_difference  = std_dev - abs(steering) 
            out_range_difference = 0
        return in_range_difference, out_range_difference
    
    def open_sim(self, multi_step, step_base, c1, c2, horizon_fixed=True):
        # Simulation logic here
        self.base_horizon = step_base
        step = step_base
        if self.veh_model == 'dynamic':
            model = Dynamic()
        elif self.veh_model == 'kinematic':
            model = Kinematic()
            # x: [x, y, psi, v, 𝛿]^T
			# u: [acc, Δ𝛿]^T
            kin_states = np.zeros((len(self.true_states), model.n_states))
            kin_inputs = np.zeros((len(self.steer_inputs), model.n_inputs))
            kin_states[:, :2] = self.states[:,:2]
            kin_states[:, 2] = self.states[:,4]
            kin_states[:, 3] = self.states[:,2]
            kin_states[:, 4] = self.states[:,5]
            kin_inputs[:, 0] = self.states[:,7]
            kin_inputs[1:, 1] = np.diff(self.steer_inputs)
            self.true_states = kin_states
            self.inputs = kin_inputs
            self.states_pred = np.zeros(kin_states.shape)
            self.states_pred[0,:] = kin_states[0, :]
            self.dxdt_ek = np.zeros(kin_states.shape)
        std_steer = np.std(self.steer_inputs)
        c1_scaler = c1
        c2_scaler = c2
        self.step_ls = [step]
        idn = 0
        if horizon_fixed:
            while idn < (self.N_SAMPLES-1):
                if idn > self.N_SAMPLES - 1 - step:
                    step = self.N_SAMPLES - 1 - idn
                inputs = np.vstack((self.inputs[idn:idn+step, :]))
                x_next, dxdt_next, self.step_ls = model.sim_continuous_multistep(self.true_states[idn, :], inputs, [0, self.SAMPLING_TIME], step, self.step_ls, wheel_v=self.wheel_data[idn:idn+step, :] , roll=self.roll_data[idn:idn+step])
                self.states_pred[idn+1:idn+step+1, :] = x_next[1:, :]
                self.dxdt_ek[idn+1:idn+step+1, :] = dxdt_next[1:, :]
                idn += step # step idx update
        else:
            while idn < (self.N_SAMPLES-1):
                step  = step_base
                curr_steer = self.inputs[idn, 2]
                in_range_diff, out_range_diff = self.calc_steer_diff(curr_steer, std_steer)
                step += math.floor(in_range_diff*c1_scaler*step) - math.ceil(out_range_diff*c2_scaler*step)
                if idn > self.N_SAMPLES - 1 - step:
                    step = self.N_SAMPLES - 1 - idn
                step = 5 if step < 5 else step # Lower bound of step
                # step_ls.append(step)
                inputs = np.vstack((self.inputs[idn:idn+step, :]))
                x_next, dxdt_next, self.step_ls = model.sim_continuous_multistep(self.true_states[idn, :], inputs, [0, self.SAMPLING_TIME], step, self.step_ls, wheel_v=self.wheel_data[idn:idn+step,:] , roll=self.roll_data[idn:idn+step])
                self.states_pred[idn+1:idn+step+1, :] = x_next[1:, :]
                self.dxdt_ek[idn+1:idn+step+1, :] = dxdt_next[1:, :]
                idn += step
    
class Plotting:
    def __init__(self, simulation, vx_start, left_np, right_np, pit_np, Whole_track = False):
        self.simulation = simulation
        model = simulation.veh_model
        horizon = simulation.base_horizon + 2
        self.vx = simulation.true_states[:, 2]
        self.steering_angle = simulation.steer_inputs
        self.step_ls = simulation.step_ls
        self.id = np.where(self.vx>vx_start)[0][0]

        self.x = simulation.true_states[self.id:, 0]
        self.y = simulation.true_states[self.id:, 1]
        self.vx = self.vx[self.id:]
        self.steering_angle = self.steering_angle[self.id:]
        self.x_ekin = simulation.states_pred[self.id:, 0]
        self.y_ekin = simulation.states_pred[self.id:, 1]
        self.angles = simulation.true_states[self.id:, 4]
        if model == 'kinematic':
            self.angles_ekin = simulation.states_pred[self.id:, 2]
        else:
            self.angles_ekin = simulation.states_pred[self.id:, 4]

        # Initialize the figure and axis
        self.fig, self.ax = plt.subplots(figsize=(12,8))
        self.ax.plot(left_np[:,0], left_np[:,1], label='left bound')
        self.ax.plot(right_np[:, 0], right_np[:,1], label='right bound')
        if pit_np.any():
            self.ax.plot(pit_np[:, 0], pit_np[:,1], label='pit')
        self.ax.set_xlabel("$X [m]$", fontsize=14)
        self.ax.set_ylabel("$Y [m]$", fontsize=14)
        self.ax.set_title("Animated Car Movement: Ground Truth vs. Model [{}]".format(model), fontsize=18)
        # Initialize a text label
        # text_label = ax.text(0.1, 0.1, 'Car Speed: ${:.2f} [m/s]$'.format(vx[0]), transform=ax.transAxes)

        self.focus_on = False if Whole_track else True
        if Whole_track:
            self.ax.set_xlim(min(self.x)-10, max(self.x)+10)
            self.ax.set_ylim(min(self.y)-10, max(self.y)+10)
            self.ax.set_aspect('equal', 'box')
        # Initialize arrow storage and maximum number of arrows to keep
        self.arrow1_storage = deque(maxlen=horizon*2)  # keep only arrows during horizons * 2
        self.arrow2_storage = deque(maxlen=horizon*2)


        # Initialize arrow (vehicle)
        self.length = 0.5
        arrow = plt.Arrow(self.x[0], self.y[0], self.length * np.cos(self.angles[0]), self.length * np.sin(self.angles[0]), width=0.1, color='green', label='True position')
        arrow_ekin = plt.Arrow(self.x_ekin[0], self.y_ekin[0], self.length * np.cos(self.angles_ekin[0]), self.length * np.sin(self.angles_ekin[0]), width=0.1, color='red', label='Pred position')
        self.ax.add_patch(arrow)
        self.ax.add_patch(arrow_ekin)
        self.ax.legend()

        # Initialize car patches
        self.true_patches = self.draw_car(self.x[0], self.y[0], self.angles[0], carID='true')
        self.ekin_patches = self.draw_car(self.x[0], self.y[0], self.angles[0], carID='pred')
        for true_patch, ekin_patch in zip(self.true_patches,self.ekin_patches):
            self.ax.add_patch(ekin_patch)
            self.ax.add_patch(true_patch)
            
        self.last_time = None
        # Initialize steering wheel
        self.steering_wheel_radius = 15.0  # Adjust as needed
        self.steering_wheel_offset_x = 0  # Horizontal offset
        self.steering_wheel_offset_y = 35  # 10 units above the car
        # Define the offset from the steering wheel's position
        self.text_offset_x = 0  # No horizontal offset
        self.text_offset_y = -self.steering_wheel_radius - 10  # 5 units below the steering wheel
        self.draw_steering_wheel(self.ax, self.x[0] + self.steering_wheel_offset_x, self.y[0] + self.steering_wheel_offset_y, self.steering_wheel_radius, self.angles[0])
        self.steering_angle_text = self.ax.text(self.x[0] + self.steering_wheel_offset_x, self.y[0] + self.steering_wheel_offset_y - self.steering_wheel_radius - 10, '', ha='center')


    def draw_car(self, x, y, angle, carID):
        patches = []
    
        # Define the basic body of the car as a triangle
        body_x = [-2, -2, 2.5]
        body_y = [-1., 1., 0]
        
        # Define the front wing
        front_wing_x = [1.5, 2, 2, 1.5]
        front_wing_y = [-1.2, -1.2, 1.2, 1.2]
        
        # Define the rear wing
        rear_wing_x = [-2.5, -2, -2, -2.5]
        rear_wing_y = [-0.7, -0.7, 0.7, 0.7]

        # Tires (drawn as circles; you may replace with rectangles or other shapes)
        tire_positions = [(-1.5, 0.9), (-1.5, -0.9), (1.2, 0.4), (1.2, -0.4)]

        # Rotate and translate
        angle_rad = angle
        rot_matrix = np.array([
            [np.cos(angle_rad), -np.sin(angle_rad)],
            [np.sin(angle_rad), np.cos(angle_rad)]
        ])
        
        body = np.dot(rot_matrix, np.array([body_x, body_y]))
        front_wing = np.dot(rot_matrix, np.array([front_wing_x, front_wing_y]))
        rear_wing = np.dot(rot_matrix, np.array([rear_wing_x, rear_wing_y]))
        
        # Translate
        body[0, :] += x
        body[1, :] += y
        front_wing[0, :] += x
        front_wing[1, :] += y
        rear_wing[0, :] += x
        rear_wing[1, :] += y

        tire_positions = np.dot(rot_matrix, np.array(tire_positions).T)
        tire_positions[0, :] += x
        tire_positions[1, :] += y
        
        if carID == 'true':
            patches.append(Polygon(body.T, closed=True, facecolor='green', edgecolor='black'))
            patches.append(Polygon(front_wing.T, closed=True, facecolor='orange', edgecolor='black'))
            patches.append(Polygon(rear_wing.T, closed=True, facecolor='orange', edgecolor='black'))
        elif carID == 'pred':
            patches.append(Polygon(body.T, closed=True, facecolor='red', edgecolor='black'))
            patches.append(Polygon(front_wing.T, closed=True, facecolor='blue', edgecolor='black'))
            patches.append(Polygon(rear_wing.T, closed=True, facecolor='blue', edgecolor='black'))
        
        for pos in tire_positions.T:
            patches.append(Circle((pos[0], pos[1]), radius=0.2, facecolor='black'))

        return patches

    def draw_steering_wheel(self, ax, center_x, center_y, radius, angle):
        # Draw outer circle
        outer_circle = Circle((center_x, center_y), radius, fill=False, color='blue')
        outer_circle.is_steering_wheel_component = True  # Custom attribute
        ax.add_patch(outer_circle)
        
        # Draw inner circle
        inner_circle = Circle((center_x, center_y), radius * 0.2, fill=True, color='blue')
        inner_circle.is_steering_wheel_component = True  # Custom attribute
        ax.add_patch(inner_circle)
        
        # Draw horizontal line (spoke), rotated by the steering angle
        angle_rad = np.radians(angle)
        dx = radius * np.cos(angle_rad)
        dy = radius * np.sin(angle_rad)
        spoke = Line2D([center_x - dx, center_x + dx], [center_y - dy, center_y + dy], color='blue')
        spoke.is_steering_wheel_component = True  # Custom attribute
        ax.add_line(spoke)

    def update(self, frame):
        # Remove the oldest path for ground truth
        if len(self.arrow1_storage) == self.arrow1_storage.maxlen:
            oldest_arrow1 = self.arrow1_storage.popleft()  
            oldest_arrow1.remove()
        # Remove the oldest path for ekin
        if len(self.arrow2_storage) == self.arrow2_storage.maxlen:
            oldest_arrow2 = self.arrow2_storage.popleft()  
            oldest_arrow2.remove()
        # Path of ground truth
        dx = self.length * np.cos(self.angles[frame])
        dy = self.length * np.sin(self.angles[frame])
        arrow = plt.Arrow(self.x[frame], self.y[frame], dx, dy, width=0.1, color='green')
        self.ax.add_patch(arrow)
        # Path of ekin prediction
        dx_ekin = self.length * np.cos(self.angles_ekin[frame])
        dy_ekin = self.length * np.sin(self.angles_ekin[frame])
        arrow_ekin = plt.Arrow(self.x_ekin[frame], self.y_ekin[frame], dx_ekin, dy_ekin, width=0.1, color='red')
        self.ax.add_patch(arrow_ekin)
        self.arrow1_storage.append(arrow)
        self.arrow2_storage.append(arrow_ekin)
        # Create and add new car position
        new_patches_ekin = self.draw_car(self.x_ekin[frame], self.y_ekin[frame], self.angles_ekin[frame], carID='pred')
        new_patches_true = self.draw_car(self.x[frame], self.y[frame], self.angles[frame], carID='true')

        for old_patch, new_patch in zip(self.ekin_patches, new_patches_ekin):
            if isinstance(old_patch, Circle) and isinstance(new_patch, Circle):
                old_patch.set_center(new_patch.center)
            elif isinstance(old_patch, Polygon) and isinstance(new_patch, Polygon):
                old_patch.set_xy(new_patch.get_xy())

        for old_patch, new_patch in zip(self.true_patches, new_patches_true):
            if isinstance(old_patch, Circle) and isinstance(new_patch, Circle):
                old_patch.set_center(new_patch.center)
            elif isinstance(old_patch, Polygon) and isinstance(new_patch, Polygon):
                old_patch.set_xy(new_patch.get_xy())
        speed_text = f'Car Speed: {self.vx[frame]:.2f} [m/s]'
        horizon_text = f'Horizon: {self.step_ls[frame+self.id]} [step]'

        # Remove existing steering wheel components
        patches_to_remove = [patch for patch in self.ax.patches if hasattr(patch, 'is_steering_wheel_component') and patch.is_steering_wheel_component]
        lines_to_remove = [line for line in self.ax.lines if hasattr(line, 'is_steering_wheel_component') and line.is_steering_wheel_component]

        for patch in patches_to_remove:
            patch.remove()
        for line in lines_to_remove:
            line.remove()

        # Zoomed in mode
        if self.focus_on:
            # Update the plot limits to focus on the arrow
            window_size = 50
            x_center = self.x[frame]
            y_center = self.y[frame]
            x_min = x_center - window_size * 3 / 4
            x_max = x_center + window_size * 3 / 4
            y_min = y_center - window_size * 1 / 3 
            y_max = y_center + window_size * 1 / 3
            self.ax.set_xlim(x_min, x_max)
            self.ax.set_ylim(y_min, y_max)
            self.ax.set_aspect('equal', 'box')
            # Set steering wheel icon
            self.steering_wheel_radius = window_size * 0.05
            steering_wheel_center_x = self.x[frame] + window_size/10
            steering_wheel_center_y = self.y[frame] + window_size/5
            text_x = steering_wheel_center_x 
            text_y = steering_wheel_center_y - window_size/10
            self.steering_angle_text.set_position((text_x, text_y))
        else:
            # Draw new steering wheel based on the current steering angle
            steering_wheel_center_x = self.x[frame] + self.steering_wheel_offset_x
            steering_wheel_center_y = self.y[frame] + self.steering_wheel_offset_y
            text_x = steering_wheel_center_x + self.text_offset_x
            text_y = steering_wheel_center_y + self.text_offset_y
            self.steering_angle_text.set_position((text_x, text_y))
        if pit_np.any():
            self.ax.legend(['left bound', 'right bound', 'pitlane', 'True position', 'Pred position', speed_text, horizon_text])
        else:
            self.ax.legend(['left bound', 'right bound', 'True position', 'Pred position', speed_text, horizon_text])
        current_steering_angle = self.steering_angle[frame] / 0.0666 / (math.pi/180.0) # Steering  wheel angle data
        self.draw_steering_wheel(self.ax, steering_wheel_center_x, steering_wheel_center_y, self.steering_wheel_radius, current_steering_angle)
        self.steering_angle_text.set_text(f'Steering wheel Angle: {current_steering_angle:.2f}°')
        # Calculate and display FPS
        current_time = time.time()
        if self.last_time is not None:
            fps = 1.0 / (current_time - self.last_time)
            # print(f"FPS: {fps}")
        self.last_time = current_time

    def animate(self):
        # self.fig, self.ax = plt.subplots(figsize=(12, 8))
        ani = FuncAnimation(self.fig, self.update, frames=range(len(self.x)), interval=50)
        plt.show()

if __name__ == "__main__":
    # Load racetrack boundaries
    desk = True
    usr = 'ning' if desk else 'jyning'

    trackID = 0
    racetracks = ['lvms', 'laguna']

    # Centralized config for each track
    track_config = {
        'lvms': {
            'left': 'inner_bound.csv',
            'right': 'outer_bound.csv',
            'pit': 'pitlane.csv',
            'dataset': 'LVMS_2024_01_11_B',
            'load_pit': True
        },
        'laguna': {
            'left': 'inner_bound.csv',
            'right': 'outer_bound.csv',
            'pit': 'pitlane_center.csv',
            'dataset': 'Laguna_no_lidar',
            'load_pit': True
        },
    }

    # Pick active track
    track_name = racetracks[trackID]
    cfg = track_config[track_name]
    bounds_path = f'data/maps/{track_name}/'

    # Load map files
    left_np  = np.loadtxt(bounds_path + cfg['left'],  delimiter=',', dtype=np.float64)
    right_np = np.loadtxt(bounds_path + cfg['right'], delimiter=',', dtype=np.float64)
    pit_np   = (np.loadtxt(bounds_path + cfg['pit'], delimiter=',', dtype=np.float64)
                if cfg['load_pit'] else np.zeros(1))

    # Vehicle model & params
    params = vehicle.AV24()
    vehicle = Vehicle(params)
    models = ['dynamic', 'kinematic']
    modelID = 0
    dataset = cfg['dataset']
    file_path = f'data/veh_data/{dataset}.csv'
    vehicle.gen_states(file_path, models)
    simulation = Simulation(vehicle, models[modelID])
    multi_step = True
    base_step  = 40 
    c1_scaler  = 9
    c2_scaler  = 2
    simulation.open_sim(multi_step, base_step, c1_scaler, c2_scaler, horizon_fixed=True)
    vx_start = 20 # start from the point vx >= [] m/s
    # Dynamic plot
    plotting = Plotting(simulation, vx_start, left_np, right_np, pit_np, Whole_track = False)
    plotting.animate()

