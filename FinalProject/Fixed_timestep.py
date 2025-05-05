import numpy as np
import time 
from skyfield.api import EarthSatellite, load, utc
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import csv


# For interpreting two-line element (TLE) data
from sgp4.api import Satrec
from sgp4.api import jday # Is this needed?

from PIL import Image




"""Global Constants"""
G = 6.674e-20 # Gravity Constant
DAYS = 1
TIME_PERIOD = DAYS * 24 * 3600 

NUM_STEPS = int(1e4)

Δt = TIME_PERIOD / NUM_STEPS

# Defines spherical body that satelites orbit
# In the use cases here, only will be the Earth, the Moon, and the Sun
class CelestialBody:
    
	def __init__(self,name,radius, mass=None, position=None,velocity=None):
		self.name = name
		self.radius = radius
		self.mass = mass
		if position is not None:
			self.d = np.array(position)
		else:
			self.d = np.array([0.0,0.0,0.0]) # Unless otherwise specificed, centered at origin

		if velocity is not None:
			self.v = np.array(velocity)
		else:
			self.v = np.array([0.0,0.0,0.0]) # Unless otherwise specified, 0 velocity
		self.a = np.zeros(3)
		self.positionArray = [self.d.copy()]


# Defines satelite that is orbiting the larger celestial body
class Satellite(CelestialBody):	
	def __init__(self, data):
	
		self.name = data[0]
		self.radius = float(data[1])
		mass = float(data[2])
		name = data[0]
		radius = float(data[1])
		self.mass = float(data[2])
		super().__init__(name,radius,mass)
		TLE1 = data[3]
		TLE2 = data[4]
		self.convertTwoLine(TLE1, TLE2)
		self.errors = []
    

	def convertTwoLine(self,TLE1,TLE2):
		
		satellite = EarthSatellite(TLE1, TLE2, self.name)
		t = satellite.epoch  # Use epoch embedded in TLE
		geocentric = satellite.at(t)
		position_km = geocentric.position.km
		velocity_kms = geocentric.velocity.km_per_s

		self.d = np.array(position_km, dtype=float)
		self.v = np.array(velocity_kms, dtype=float)

		self.positionArray.clear()
		self.positionArray.append(self.d.copy())
		# Add more or automate this process later
		
class OrbitalSimulation:
    # Specifying a setuptype will override all successive variables in the simulation.
	def __init__(self,bodies,method):

		self.bodies = bodies
		self.import_satelite_csv() # Add the satelites to the simulation
		self.current_time = 0
		self.times = [self.current_time]
		self.method = method

	def import_satelite_csv(self):
		with open('FinalProject/satelites.csv', newline='') as csvfile:
			reader = csv.reader(csvfile)
			for row in reader:
				satellite = Satellite(row)
				self.bodies.append(satellite)
		

	"""Run Simulation Function"""
	def run_orbital_simulation(self):
		# For knowing how long the simulation is taking
		start = time.perf_counter()
		for step in range(NUM_STEPS):
			if step % int(.01 * NUM_STEPS) == 0:
				print(f"{int(step / NUM_STEPS * 100)} Percent Completed!")
			self.advanceTimestep()
			# self.sanityCheck()
		end = time.perf_counter()
		print("Completed!")
		print(f"Took {end-start} seconds.")


	# Kept from N-Body Simulation
	def compute_acceleration(self, body):
		acceleration = np.zeros(3)
		for other_body in self.bodies:
			if other_body is body:
				continue
			diff = other_body.d - body.d  # Vector pointing from body to other_body
			distance_squared = np.dot(diff, diff)  # r^2
			distance = np.sqrt(distance_squared)  
			distance_cubed = distance_squared * distance  

			# Gravitational acceleration formula: G * m * xj-xi / |xj-xi|^3 
			acceleration += G * other_body.mass * diff / distance_cubed  
		return acceleration
	

	# Checks that no objects are colliding with each other, or flying off into space.
	def sanityCheck(self):
		# Write method here
		return 

	"""Numerical Approximation Methods"""

	def explicit_euler(self):
		new_accelerations = [self.compute_acceleration(body) for body in self.bodies]
		for body, a_new in zip(self.bodies, new_accelerations):
			body.d += body.v * Δt
			body.v += a_new * Δt
			body.positionArray.append(body.d.copy())


	# Uses the newton method for solving the implicit equation
	def implicit_euler(self):
		# Extrapolates arrays of all positions, velocities, and masses
		positions = np.array([body.d for body in self.bodies])
		velocities = np.array([body.v for body in self.bodies])
		masses = np.array([body.mass for body in self.bodies])

		r_next = positions.copy()

		for iteration in range(len(self.bodies)):
			return

	def verlet_integration(self):
		for body in self.bodies:
			body.v += 0.5 * body.a * Δt

		# Full-step for position
		for body in self.bodies:
			body.d += body.v * Δt
			body.positionArray.append(body.d.copy())

		# Compute all accelerations using their updated position
		new_accelerations = [self.compute_acceleration(body) for body in self.bodies]

		# Second half-step for velocity using the new accelerations
		for body, a_new in zip(self.bodies, new_accelerations):
			body.v += 0.5 * a_new * Δt
			body.a = a_new 

	# Other approximation methods here


	def computePositionError(self,b5_position, b4_position,):
		return abs(np.linalg.norm(b5_position - b4_position))

	def computeVelocityError(self,b5_velocity, b4_velocity):
		return np.linalg.norm(b5_velocity - b4_velocity)

	def dopri45(self):
		"""Butcher Tableau for dopri45"""
		A = np.array([
		[0, 0, 0, 0, 0, 0, 0],
		[1/5, 0, 0, 0, 0, 0, 0],
		[3/40, 9/40, 0, 0, 0, 0, 0],
		[44/45, -56/15, 32/9, 0, 0, 0, 0],
		[19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0, 0],
		[9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0, 0],
		[35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0]
	])

		C = np.array([0, 1/5, 3/10, 4/5, 8/9, 1, 1])
  
		"""4th and 5th Order Coefficients"""
		b5 = np.array([35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0])
		b4 = np.array([5179/57600,0,7571/16695,393/640,-92097/339200,187/2100,1/40])
	
		return A,b4,b5,C,7 # Number of Steps

	def bogacki_shampine(self):
		"""Butcher Tableau for Bogacki–Shampine 3(2) method"""
		A = np.array([
			[0,    0,   0],
			[1/2,  0,   0],
			[0,   3/4,  0],
			[2/9, 1/3, 4/9]
		])

		C = np.array([0, 1/2, 3/4, 1])  # Nodes (c)

		"""2nd and 3rd Order Coefficients"""
		b3 = np.array([2/9, 1/3, 4/9, 0])  # 3rd order (main solution)
		b2 = np.array([7/24, 1/4, 1/3, 1/8])  # 2nd order (embedded)

		return A, b2, b3, C, 4  # Number of stages

	def heun_euler(self):
		"""Butcher Tableau for Heun-Euler 2(1) method"""
		A = np.array([
			[0,   0],
			[1,   0],
			[1/2, 1/2]  # This row is just for Heun's final combination
		])

		C = np.array([0, 1, 1])  # Nodes (c)

		"""1st and 2nd Order Coefficients"""
		b2 = np.array([1/2, 1/2, 0])  # 2nd order (Heun)
		b1 = np.array([1,   0,   0])  # 1st order (Euler)

		return A, b1, b2, C, 3  # Number of stages
			
	def runge_kutta_stepper(self, method):
		A, b4, b5, C, steps = method()
		initial_positions = [body.d.copy() for body in self.bodies]
		initial_velocities = [body.v.copy() for body in self.bodies]
		k = [[] for _ in range(steps)]

		for i in range(steps):
			for b_idx, body in enumerate(self.bodies):
				body.d = initial_positions[b_idx].copy()
				body.v = initial_velocities[b_idx].copy()
				for j in range(i):
					body.d += Δt * A[i, j] * k[j][b_idx][1]
					body.v += Δt * A[i, j] * k[j][b_idx][0]

			accelerations = [self.compute_acceleration(body) for body in self.bodies]
			k[i] = [(a.copy(), body.v.copy()) for a, body in zip(accelerations, self.bodies)]

		
		for b_idx, body in enumerate(self.bodies):
			acc_sum = sum(b5[i] * k[i][b_idx][0] for i in range(steps))
			vel_sum = sum(b5[i] * k[i][b_idx][1] for i in range(steps))
			b5_velocity = initial_velocities[b_idx] + Δt * acc_sum
			body.v = b5_velocity
			b5_position = initial_positions[b_idx] + Δt * vel_sum
			body.d = b5_position
			body.positionArray.append(body.d.copy())

			if isinstance(body, Satellite):
				# For error purposes, compute b4 positions and velocities
				acc_sum = sum(b4[i] * k[i][b_idx][0] for i in range(steps))
				vel_sum = sum(b4[i] * k[i][b_idx][1] for i in range(steps))
				b4_velocity = initial_velocities[b_idx] + Δt * acc_sum
				b4_position = initial_positions[b_idx] + Δt * vel_sum
				

				body.errors.append(self.computePositionError(b5_position,b4_position))

	
	
	# advances timestep for all bodies in the system
	def advanceTimestep(self):
		if self.method== 'Explicit':
			self.explicit_euler()
		elif self.method == 'Implicit':
			self.implicit_euler()
		elif self.method == 'Verlet':
			self.verlet_integration()
		elif self.method == "Dopri":
			self.runge_kutta_stepper(self.dopri45)
		else:
			raise ValueError("Unsupported method")
		# Advance timestep
		self.current_time += Δt
		self.times.append(self.current_time)


	"""Graphing Methods"""

	def graph_errors(self):
		# Graph individual errors for each satellite
		for body in self.bodies:
			if isinstance(body, Satellite):
				plt.plot(self.times[1:],body.errors,label="Times Vs Errors")
				plt.xlabel("Time (s)")
				plt.ylabel("Error (Embedded)")
				plt.title(f"Embedded Error vs Time for {body.name}")
				plt.grid(True)
				plt.show()


	def graph_positions(self):
		fig = go.Figure()
			
		total_points = 1000


		for idx, body in enumerate(self.bodies):
			if body.name == "Earth":
				num_positions = len(body.positionArray)
				x_vals = []
				y_vals = []
				z_vals = []
				for t in range(0,num_positions):
					x_vals.append(body.positionArray[t][0])
					y_vals.append(body.positionArray[t][1])
					z_vals.append(body.positionArray[t][2])
				earth_img = Image.open("FinalProject/textures/Earth.jpg")

				# Resize the image if necessary, to reduce performance overhead
				earth_img = earth_img.resize((360, 180))  # Resize based on performance needs

				# Convert image to a numpy array (this will hold RGB values)
				earth_array = np.asarray(earth_img)

				# Generate spherical coordinates (longitude and latitude)
				R = 6371  # Earth radius in km
				lon = np.linspace(-np.pi, np.pi, earth_array.shape[1])  # Longitude range: -pi to pi
				lat = np.linspace(-np.pi/2, np.pi/2, earth_array.shape[0])  # Latitude range: -pi/2 to pi/2

				# Create meshgrid for latitude and longitude
				lon, lat = np.meshgrid(lon, lat)

				# Convert spherical coordinates to Cartesian coordinates (for 3D surface)
				x = R * np.cos(lat) * np.cos(lon)
				y = R * np.cos(lat) * np.sin(lon)
				z = R * np.sin(lat)

				x_shifted = x + x_vals[-1]
				y_shifted = y + y_vals[-1]
				z_shifted = z + z_vals[-1]

				# Map the RGB texture onto the sphere
				# We need to flatten the image and map it onto the sphere

				# Normalize the image texture
				r = earth_array[:, :, 0] / 255.0  # Red channel, normalized to [0, 1]
				g = earth_array[:, :, 1] / 255.0  # Green channel, normalized to [0, 1]
				b = earth_array[:, :, 2] / 255.0  # Blue channel, normalized to [0, 1]

				# Stack the RGB channels for color mapping (this will be the surfacecolor)
				rgb_texture = np.stack([r, g, b], axis=-1)

				# Add 3d line for body trajectory
				fig.add_trace(go.Scatter3d(
					x=x_vals,
					y=y_vals,
					z=z_vals,
					mode='lines',
					name=f'{body.name}',
					line=dict(width=2)
				))

				# Plot the Earth with the image as a mapped texture
				fig = go.Figure(data=[
					go.Surface(
						x=x_shifted, y=y_shifted, z=z_shifted,
						surfacecolor=rgb_texture,  # Use full RGB color mapping
						showscale=False,  # Disable the color scale (it's not needed for textures)
						lighting=dict(ambient=0.6),  # Lighting properties for better appearance
						opacity=1  # Full opacity for the Earth surface
					)
				])
				

				
			elif isinstance(body,Satellite):
    
				num_positions = len(body.positionArray)
				x_vals = []
				y_vals = []
				z_vals = []
				for t in range(0,num_positions):
					x_vals.append(body.positionArray[t][0])
					y_vals.append(body.positionArray[t][1])
					z_vals.append(body.positionArray[t][2])

				# Add 3d line for body trajectory
				fig.add_trace(go.Scatter3d(
					x=x_vals,
					y=y_vals,
					z=z_vals,
					mode='lines',
					name=f'{body.name}',
					line=dict(width=2)
				))
			
				# Add a grey sphere to represent the final location of the satellite
				fig.add_trace(go.Scatter3d(
					x=[x_vals[-1]],  # Final x position
					y=[y_vals[-1]],  # Final y position
					z=[z_vals[-1]],  # Final z position
					mode='markers',
					marker=dict(
						size=5,  # Size of the sphere (adjust to fit your scale)
						color='grey',  # Grey color for the sphere
						opacity=1  # Full opacity for clear visibility
					),
					name=f'{body.name} (Final Location)'  # Name of the body for trace identification
				))
					
				
				# Set up axes and layout 
				fig.update_layout(
					title=f'Satelite Trajectories for {int(DAYS * 24)} Hours - Timestep: Approx {Δt} Seconds',
					scene=dict(
						xaxis_title='X',
						yaxis_title='Y',
						zaxis_title='Z'
					),
					width=800,
					height=800,
					showlegend=True
				)
		fig.show()




def main():
	# CHANGE METHODS HERE FOR DIFFERENT THINGS
    method = "Dopri"
    graphErrors = True

    earth = CelestialBody("Earth",6371,5.972e24)
    sun = CelestialBody("Sun", 696340,1.989e30,[150790000.0,0.0,0.0])
    moon = CelestialBody("Moon", 1737.4,7.348e22,[0.0,384400.0,0.0])

    bodies = [earth, sun, moon] 


    system = OrbitalSimulation(bodies,method)
    system.run_orbital_simulation()
    if method == "Dopri" and graphErrors:
        system.graph_errors()
    system.graph_positions()


if __name__ == "__main__":
	main()





			