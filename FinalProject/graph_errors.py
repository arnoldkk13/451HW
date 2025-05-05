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
	def run_orbital_simulation(self, NUM_STEPS,Δt):
		# For knowing how long the simulation is taking
		start = time.perf_counter()
		for step in range(NUM_STEPS):
			if NUM_STEPS > 100:
				if step % int(.01 * NUM_STEPS) == 0:
					print(f"{int(step / NUM_STEPS * 100)} Percent Completed!")
			self.advanceTimestep(Δt)
			# self.sanityCheck()
		end = time.perf_counter()
		print(f"{Δt} Timstep Completed!")
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


	def computePositionError(self,b5_position, b4_position,):
		return abs(np.linalg.norm(b5_position - b4_position))

	def computeVelocityError(self,b5_velocity, b4_velocity):
		return abs(np.linalg.norm(b5_velocity - b4_velocity))

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
			
	def compute_acceleration_from_position(self, position, mass, epsilon=1e-10):
		acceleration = np.zeros(3)
		for other_body in self.bodies:
			if np.array_equal(other_body.d, position):
				continue  # Skip the self-interaction (body can't interact with itself)

			diff = other_body.d - position
			distance_squared = np.dot(diff, diff)
			distance = np.sqrt(distance_squared)
			
			# Prevent division by zero or very small numbers
			if distance < epsilon:
				continue  # Skip this pair or handle separately if needed

			distance_cubed = distance_squared * distance
			acceleration += G * other_body.mass * diff / distance_cubed
		return acceleration

	def runge_kutta_stepper(self, method, Δt):
		# Get the Butcher tableau parameters (A, b4, b5, C, steps)
		A, b4, b5, C, steps = method()

		# Store initial positions and velocities
		initial_positions = [body.d.copy() for body in self.bodies]
		initial_velocities = [body.v.copy() for body in self.bodies]

		# Initialize k array to store intermediate results
		k = [[None] * len(self.bodies) for _ in range(steps)]

		# Loop over the stages of the Runge-Kutta method
		for i in range(steps):
			# Temporary storage for stage positions and velocities
			stage_positions = []
			stage_velocities = []

			# Prepare the stage positions and velocities for each body
			for b_idx in range(len(self.bodies)):
				d_temp = initial_positions[b_idx].copy()
				v_temp = initial_velocities[b_idx].copy()

				# Accumulate the results from previous stages
				for j in range(i):
					d_temp += Δt * A[i, j] * k[j][b_idx][1]  # Velocity part
					v_temp += Δt * A[i, j] * k[j][b_idx][0]  # Acceleration part

				stage_positions.append(d_temp)
				stage_velocities.append(v_temp)

			# Compute accelerations for each body at the current stage positions
			accelerations = [
				self.compute_acceleration_from_position(stage_positions[b_idx], self.bodies[b_idx].mass)
				for b_idx in range(len(self.bodies))
			]

			# Store the accelerations and velocities at each stage
			for b_idx, acc in enumerate(accelerations):
				k[i][b_idx] = (acc.copy(), stage_velocities[b_idx].copy())

		# Update positions and velocities using b4 and b5 coefficients
		for b_idx, body in enumerate(self.bodies):
			# Compute the total accelerations and velocities based on b5 coefficients
			acc_sum = sum(b5[i] * k[i][b_idx][0] for i in range(steps))
			vel_sum = sum(b5[i] * k[i][b_idx][1] for i in range(steps))

			# Final update for position and velocity
			body.v = initial_velocities[b_idx] + Δt * acc_sum
			body.d = initial_positions[b_idx] + Δt * vel_sum
			body.positionArray.append(body.d.copy())

			# For error checking: compute b4 positions and velocities (error estimation)
			acc_sum_b4 = sum(b4[i] * k[i][b_idx][0] for i in range(steps))
			vel_sum_b4 = sum(b4[i] * k[i][b_idx][1] for i in range(steps))

			# Compute b4 velocities and positions for error calculation
			b4_velocity = initial_velocities[b_idx] + Δt * acc_sum_b4
			b4_position = initial_positions[b_idx] + Δt * vel_sum_b4

			# Error computation if the body is a satellite
			if isinstance(body, Satellite):
				# Compute the error between b5 and b4
				body.errors.append(self.computePositionError(body.d, b4_position))
			
	
	# advances timestep for all bodies in the system
	def advanceTimestep(self,Δt):
		if self.method== 'Bogacki':
			self.runge_kutta_stepper(self.bogacki_shampine,Δt)
		elif self.method == 'Heun':
			self.runge_kutta_stepper(self.heun_euler,Δt)
		elif self.method == "Dopri":
			self.runge_kutta_stepper(self.dopri45,Δt)
		else:
			raise ValueError("Unsupported method")
		# Advance timestep
		self.current_time += Δt
		self.times.append(self.current_time)


	def getMaxSatelliteError(self):
		max_error = -1 * float('inf')
		for body in self.bodies:
			if isinstance(body,Satellite):
				error = max(body.errors)
				if error > max_error:
					max_error = error
		if max_error == -1 * float('inf'):
			print("Error not updating")
			exit(1)
		return max_error,self.method
			



"""Graphing Methods"""

def graph_errors():
	methods = ['Heun', 'Bogacki', "Dopri"]

	earth = CelestialBody("Earth",6371,5.972e24)
	sun = CelestialBody("Sun", 696340,1.989e30,[150790000.0,0.0,0.0])
	moon = CelestialBody("Moon", 1737.4,7.348e22,[0.0,384400.0,0.0])

	



	# Graph the maximum error accumulated across each method
	# First, must run each method for various Δt values and store their maximum error. This will take a while for sure
	Δt = TIME_PERIOD / 100 # Start with this
	NUM_STEPS = 100
	heun_error = []
	bogacki_error = []
	dopri_error = []
	Δts = []
	while Δt > 1:
		bodies = [earth] 
		for method in methods:
			system = OrbitalSimulation(bodies,method)
			system.run_orbital_simulation(NUM_STEPS, Δt)
			error,method = system.getMaxSatelliteError()
			if method == 'Heun':
				heun_error.append(error)
			elif method == 'Bogacki':
				bogacki_error.append(error)
			elif method == 'Dopri':
				dopri_error.append(error)
			else:
				print("Unsupported method. Check spelling and Capitalization")
			bodies = [earth] # Reset here

		Δts.append(Δt)
		Δt /= 2
		NUM_STEPS *= 2
		
		
	
	plt.plot(Δts,heun_error, label="First Order (Euler Method)")
	plt.plot(Δts,bogacki_error, label="Second Order (Verlet Method)")
	plt.plot(Δts,dopri_error, label="Fourth Order (RK-4 Method)")
	plt.xlabel("Time Step Size Δt (s)")
	plt.semilogx()
	plt.semilogy()
	plt.ylabel("Error (Embedded) km}")
	plt.legend()
	plt.title(f"Embedded Errors vs Δt (Logarithmic)")
	plt.grid(True)
	plt.show()






def main():
	# CHANGE METHODS HERE FOR DIFFERENT THINGS
    
    graphErrors = True


    if graphErrors:
        graph_errors()


if __name__ == "__main__":
	main()





			