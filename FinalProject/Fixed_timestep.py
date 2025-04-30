import numpy as np

# For interpreting two-line element data
from sgp4.api import Satrec
from sgp4.api import jday


import plotly.graph_objects as go



"""Global Constants"""
G = 6.674e-11 # Gravity Constant
DAYS = 1
TIME_PERIOD = DAYS * 24 * 3600 

NUM_STEPS = 1e5

Δt = TIME_PERIOD / NUM_STEPS

# TLE of the ISS Zarya, in NEO
# taken at 3:00PM 4/30/2025
ISS_ZARYA_TLE1 = "1 25544U 98067A   25119.19035294  .00013779  00000+0  25440-3 0  9995"
ISS_ZARYA_TLE2 = "2 25544  51.6352 189.7367 0002491  81.0639 279.0631 15.49383308507563"

# TLE of ARKTIKA-M 1, an arctic weather monitoring satellite in HEO
# Taken at 3:05PM 4/30/2025
ARKTIKA_TLE1 = "1 47719U 21016A   25118.78057943  .00000078  00000-0  00000-0 0  9993"
ARKTIKA_TLE2 = "2 47719  63.1448 113.9682 7084235 268.5740  16.7891  2.00604244 30475"



# Defines spherical body that satelites orbit
# In the use cases here, only will be the Earth, the Moon, and the Sun
class CelestialBody:
	def __init__(self,name,radius, mass, position=None,velocity=None):
		self.name = name
		self.radius = radius
		self.mass = mass
		if position is not None:
			self.d = np.array(position)
		else:
			self.d = np.array([0,0,0]) # Unless otherwise specificed, centered at origin

		if velocity is not None:
			self.v = np.array(velocity)
		else:
			self.v = np.array([0,0,0]) # Unless otherwise specified, 0 velocity

		self.positionArray = [self.d.copy()]


		
    
# Defines satelite that is orbiting the larger celestial body
class Satellite(CelestialBody):	
	def __init__(self, name):
		super().__init__(self,name)
		self.convertTwoLine()
    

	def convertTwoLine(self):
		if self.name.lower() == "zarya":
			self.mass = 19323 # mass in kg
			self.radius = 12.56 # radius in kg

			satellite = Satrec.twoline2rv(ISS_ZARYA_TLE1, ISS_ZARYA_TLE2)
			epoch_jd = satellite.jdsatepoch  # Julian date (days)
			epoch_fr = satellite.jdsatepochF # Fractional part of day
			error_code, r, v = satellite.sgp4(epoch_jd, epoch_fr)
			if error_code != 0:
				print("Error:", error_code)
			else:
				self.d = r
				self.v = v

		elif self.name.lower() == "arktika":
			self.mass = 2100 # mass in kg
			self.radius = 5 # estiamted radius in kg
			satellite = Satrec.twoline2rv(ARKTIKA_TLE1, ARKTIKA_TLE2)
			epoch_jd = satellite.jdsatepoch  # Julian date (days)
			epoch_fr = satellite.jdsatepochF # Fractional part of day
			error_code, r, v = satellite.sgp4(epoch_jd, epoch_fr)
			if error_code != 0:
				print("Error:", error_code)
			else:
				self.d = r
				self.v = v

		# Add more or automate this process later
		
class OrbitalSimulation:
    # Specifying a setuptype will override all successive variables in the simulation.
	def __init__(self,bodies):
		self.bodies = bodies
		self.current_time = 0
		self.times = [self.current_time]

	"""Run Simulation Function"""
	def orbital_simulation(self):
		# For knowing how long the simulation is taking
		for step in range(NUM_STEPS):
			if step % int(.01 * NUM_STEPS) == 0:
				print(f"{int(step / NUM_STEPS * 100)} Percent Completed!")
			self.advanceTimestep()
			self.sanityCheck()

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
			acceleration += G * other_body.m * diff / distance_cubed  
		return acceleration
	
	
	# advances timestep for all bodies in the system
	def advanceTimestep(self):
		if self.method== 'explicit':
			self.explicit_euler()
		elif self.method == 'implicit':
			self.implicit_euler()
		elif self.method == 'verlet':
			self.verlet_integration()
		else:
			raise ValueError("Unsupported method")
		# Advance timestep
		self.current_time += Δt
		self.times.append(self.current_time)

	# Checks that no objects are colliding with each other, or flying off into space.
	def sanityCheck(self):
		# Write method here
		return 


	"""Numerical Approximation Methods"""

	def explicit_euler(self):

		for body in self.bodies:
			new_a = body.compute_acceleration(body)
			body.v += new_a * Δt
			body.d += body.v * Δt
			body.positionArray.append(body.d.copy())

	def implicit_euler():
		return 0
	
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

		
