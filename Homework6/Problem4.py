import numpy as np
import random 
import time
import plotly.graph_objects as go

  
"""Simulation Variables"""
N = 5 # Number of stars in the cluster. Change to simulate a different number of stars


SimulationYears = 1e10 # Number of total years

"""Different Timesteps depending on number of stars"""
# Δtyears = SimulationYears / 1e3 # Number of years for the timestep. 1000 Iterations guarenteed (~ 9 seconds for 10 stars)
Δtyears = SimulationYears / 1e5 # Number of years for the timestep. 100000 Iterations guarenteed (~ 1 minute for 10 stars)
# Δtyears = SimulationYears / 1e6 # Number of years for the timestep. 1,000,000 Iterations guarenteed (~15 minutes for 10 stars)
Δt = Δtyears * 3.154e+7
SimulationSeconds = SimulationYears * 3.154e+7

NumSteps = int(np.ceil(SimulationSeconds / Δt))
BINARY_DISTANCE = 1e18

"""Global Constants"""
m = 1.989e30 # mass of star. Assumed to all be one solar mass
G = 6.674e-11 # Gravity Constant
σx = 50 * 9.4607e15 # Position upper bound for gaussian distribution
σv = 1.2 * np.sqrt(N)  # Velocities given as assumption from the number of stars considered

"""Classes"""
class Star:
    # Each star object in the cluster
    # position in the three dimensional vector, d, that dictates the x, y, and z position of the star. 
    # Same with velocity, v
	def __init__(self,position=None, velocity=None):
		self.m = m
		# If position and velocity are provided, use them
		if position is not None:
			self.d = np.array(position)
		else:
			self.d =  np.array([random.gauss(0,σx), random.gauss(0,σx), random.gauss(0,σx)])
		if velocity is not None:
			self.v = np.array(velocity)
		else:
			self.v = np.array([random.gauss(0,σv), random.gauss(0,σv), random.gauss(0,σv)])
		self.a = np.zeros(3) # acceleration is computed at each step
		self.positionArray = [self.d.copy()]

class StarSystem:
	def __init__(self,setup_type):
		if setup_type == "binary":
			print("Running binary N Body Simulation...")
			self.stars = self.setup_binary_system()

		else:
			print("Running N Body Simulation...")
			self.stars = self.createRandomStars()
			
		self.current_time = 0
		self.times = [0]

	def setup_binary_system(self):
		stars = []
		distance =  BINARY_DISTANCE
		
		r1 = np.array([-0.5 * distance, 0.0, 0.0])
		r2 = np.array([+0.5 * distance, 0.0, 0.0])
		VELOCITY = np.sqrt(G * m / (distance*.5)) /2

		v1 = np.array([0.0, +VELOCITY, 0.0])
		v2 = np.array([0.0, -VELOCITY, 0.0])

		star1 = Star(position=r1, velocity=v1)
		star2 = Star(position=r2, velocity=v2)
  
		stars.append(star1)
		stars.append(star2)

		return stars

	def createRandomStars(self):
		stars = []
		for _ in range(N):
			star = Star()
			stars.append(star)
		return stars
	
	def compute_acceleration(self, star):
		acceleration = np.zeros(3)
		for other_star in self.stars:
			if other_star is star:
				continue
			diff = other_star.d - star.d  # Vector pointing from star to other_star
			distance_squared = np.dot(diff, diff)  # r^2
			distance = np.sqrt(distance_squared)  
			distance_cubed = distance_squared * distance  

			# Gravitational acceleration formula: G * m * xj-xi / |xj-xi|^3 
			acceleration += G * other_star.m * diff / distance_cubed  

		return acceleration
		
	# advances timestep for all stars in the system
	def advanceTimestep(self):
		# First half-step for velocity
		for star in self.stars:
			star.v += 0.5 * star.a * Δt

		# Full-step for position
		for star in self.stars:
			star.d += star.v * Δt
			star.positionArray.append(star.d.copy())

		# Compute all accelerations using their updated position
		new_accelerations = [self.compute_acceleration(star) for star in self.stars]

		# Second half-step for velocity using the new accelerations
		for star, a_new in zip(self.stars, new_accelerations):
			star.v += 0.5 * a_new * Δt
			star.a = a_new 

		# Advance time once
		self.current_time += Δt
		self.times.append(self.current_time)
				
	"""Run Simulation Function"""
	def N_Body_Simulation(self):
		for step in range(NumSteps):
			if step % int(.01 * NumSteps) == 0 and step != 0:
				print(f"{int(step / NumSteps * 100)} Percent Completed!")
			self.advanceTimestep()
      
	"""Graphing Functions"""
	def graph_positions(self):
		fig = go.Figure()

		total_points = 1000
		timestep_interval = NumSteps // total_points

		for idx, star in enumerate(self.stars):
			num_positions = len(star.positionArray)
			x_vals = []
			y_vals = []
			z_vals = []
			for t in range(0,num_positions,timestep_interval):
				x_vals.append(star.positionArray[t][0])
				y_vals.append(star.positionArray[t][1])
				z_vals.append(star.positionArray[t][2])

			# Add 3d line for star trajectory
			fig.add_trace(go.Scatter3d(
				x=x_vals,
				y=y_vals,
				z=z_vals,
				mode='lines',
				name=f'Star {idx + 1}',
				line=dict(width=2)
			))
			
			# Add dot that represents the final location of the star
			fig.add_trace(go.Scatter3d(
				x=[x_vals[-1]],
				y=[y_vals[-1]],
				z=[z_vals[-1]],
				mode='markers',
				name=f'Star {idx + 1} Current Location',
				marker=dict(size=6, symbol='circle')
			))

			# Set up axes and layout 
			fig.update_layout(
				title=f' {N} Star System Trajectories for {int(SimulationYears):.0e} Years - Timestep: {int(Δtyears)} Years',
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
	"""Uncomment to run random stars or binary perfect star system"""
	system = StarSystem(setup_type=None)
	# system = StarSystem(setup_type="binary")
 
	start = time.perf_counter()
	system.N_Body_Simulation()
	end = time.perf_counter()
	print("Completed!")
	print(f"Took {end-start} seconds.")
 
	system.graph_positions()


if __name__ == "__main__":
	main()
