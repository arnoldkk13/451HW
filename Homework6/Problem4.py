import numpy as np
import random 
  
  
"""Simulation Variables"""
N = 2 # Number of stars in the cluster. Change to simulate a different number of stars
timestep = 100 # Number of years for the timestep

simulationYears = 10e6 # Number of total years
simulationSeconds = 3.15 * 10e13


"""Global Constants"""
mi = 1.989 * 10e30 # mass of star. Assumed to all be one solar mass
G = 6.674 * 10e-11 # Gravity Constant
σx = 50 * 9.4607 * 10e15 # Position upper bound for gaussian distribution
σv = 5 * np.sqrt(N) # Velocities given as assumption from the number of stars considered

"""Classes"""
class Star:
    # Each star object in the cluster
    # position in the three dimensional vector that dictates the x,y, and z position of the star
    # Same with velocity
	def __init__(self):
		self.mass = mi
		self.position = np.array([random.gauss(0,σx), random.gauss(0,σx), random.gauss(0,σx)])
		self.velocity = np.array([random.gauss(0,σv), random.gauss(0,σv), random.gauss(0,σv)])
		self.positionArray = [self.position.copy()]

	# Updates the euler position and writes to position array
	# Satisfies the euler step xi+1 = xi + f(ti,xi)Δt
	def updateEuler(self,Δt,acceleration):
		self.position += self.velocity * Δt
		self.velocity += acceleration * Δt
		self.positionArray.append(self.position.copy())

class StarSystem:
	def __init__(self, stars):
		self.stars = stars
		self.G = G # Gravity constant

	# Compute the accelerations of each star based on the position of other stars
	# Acceleration is a 3d vector that has an x,y,z component
	def compute_accelerations(self):
		accelerations = [np.zeros(3) for _ in self.stars]
		for i, star_i in enumerate(self.stars):
			for j, star_j in enumerate(self.stars):
				if i == j:
					continue 
				accelerations[i] += (G * star_i.mass * star_j.mass) / abs(star_j.position - star_i.position)**2 * (star_j.position - star_i.position) / abs(star_j.position - star_i.position) 
				
	

"""Approximation Functions"""
def f(t,x):
    return

# The explicit euler method is defined as xi+1 = xi + f(ti,xi)Δt
def explicitEuler(stars,years,Δt):
    ts = [0]
    for step in range(years):
        for star in stars:
            acceleration = 0 # Define acceleration on paper here
            star.updateEuler(Δt, acceleration)
            
    
    for i in range(N):
        f_vector = f(t0,x0)
        x1 = x0 + f_vector * Δt
        x0 = x1
        t0 = t0 + Δt
        ts.append(t0)
        xs.append(x1)
    return xs,ts

"""Graphing Functions"""


"""Simulation Function"""
def N_Body_Simulation(stars,years,Δt):
    explicitEuler(stars,years,Δt)


def main():
    stars = [Star(n) for n in N]
    N_Body_Simulation(stars,simulationYears,timestep)
    
    
# The explicit euler method is defined as xi+1 = xi + f(ti,xi)Δt
def explicitEuler(stars,years,Δt):
    N = int(36000/Δt)
    xs = [x0]
    ts = [t0]
    for i in range(N):
        f_vector = f(t0,x0)
        x1 = x0 + f_vector * Δt
        x0 = x1
        t0 = t0 + Δt
        ts.append(t0)
        xs.append(x1)
    return xs,ts