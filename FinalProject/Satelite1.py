import numpy as np
import plotly.graph_objects as go


# Defines spherical body that satelites orbit
# In the use cases here, only will be the Earth and the Moon
class CelestialBody:
	def __init__(self,name,radius,position=None,velocity=None):
		self.name = name
		self.radius = radius
		if position is not None:
			self.d = np.array(position)
		else:
			self.d = np.array([0,0,0]) # Unless otherwise specificed, centered at origin

		if velocity is not None:
			self.v = np.array(velocity)
		else:
			self.v = np.array([0,0,0]) # Unless otherwise specified, 0 velocity
    
# Defines satelite that is orbiting the larger celestial body
# Here, theta defines the polar angle [0,pi] and phi defines the azumuthal angle [0,2pi]
class Satellite:	
	def __init__(self,name,altitude,theta,phi,CelestialBody,position=None,velocity=None):
		self.name = name
		self.altitude = altitude
		self.theta = theta
		self.phi = phi
		self.CelestialBody = CelestialBody
		if position is not None:
			self.d = np.array(position)
		else:
			self.d = np.array([0,0,0]) # Unless otherwise specificed, centered at origin

		if velocity is not None:
			self.v = np.array(velocity)
		else:
			self.v = np.array([0,0,0]) # Unless otherwise specified, 0 velocity
    
    
class OrbitalSimulation:
    # Specifying a setuptype will override all successive variables in the simulation.
	def __init__(self,setup_type,):
		if setup_type.lower() == "placeholder...":
			print(f"Running {setup_type.lower()} simulation...")
			
		else:
			print("Running unique simulation...")
			self.stars = self.createRandomStars()