import numpy as np
import time
import matplotlib.pyplot as plt


def f(t,x_vect): 
	d,v,h,m = x_vect
	
	"""Internal Function Definitions"""
	DV = 1.5 if v < 8 and h <=0 else -1/10*v**2
	DM = 0 if v < 8 and h <= 0 else -4 + v

	f_d = v
	f_v = DV
	f_h = m
	f_m = DM
	return np.array([f_d,f_v,f_h,f_m])
	
# The explicit euler method is defined as xi+1 = xi + f(ti,xi)Δt
def explicitEuler(Δt,x0,t0):
    N = int(25/Δt)
    xs = [x0]
    ts = [t0]
	# Change dynamically to stop once the turkey reaches the ground again
    while not (x0[2] <= 0 and x0[3] < 0): # This is the case where the turkey's vertical velocity is down but the position hits 0
        f_vector = f(t0,x0)
        x1 = x0 + f_vector * Δt
        x0 = x1
        t0 = t0 + Δt
        ts.append(t0)
        xs.append(x0)
    print("Completed!")
    return xs,ts

def displayResults(Δt,x0,t0):
	start = time.perf_counter()
	xs,ts = explicitEuler(Δt,x0,t0)

	ds = [(position[0]) for position in xs]
	vs = [(velocity[1]) for velocity in xs]
	hs = [(vert_position[2]) for vert_position in xs]
	ms = [(vert_velocity[3]) for vert_velocity in xs]

	end = time.perf_counter()
	print(f"Finished calculating. Took {end-start} seconds.")


	"""Horizontal Distance vs Time"""
	plt.plot(ts,ds,color="blue")
	plt.xlabel("Time (s)")
	plt.ylabel("Horizontal Position (m)")
	plt.title(f"Horizontal Position vs Time (Time Step: {Δt})")
	plt.xlim(0, 12)
	plt.ylim(0, 50)
	plt.grid(True)
	plt.show()

	"""Vertical Distance vs Time"""
	plt.plot(ts, hs, color="green")
	plt.xlabel("Time (s)")
	plt.ylabel("Vertical Position (m)")
	plt.title(f"Vertical Position vs Time (Time Step: {Δt})")
	plt.xlim(0, 12)
	plt.ylim(-1, 6)
	plt.grid(True)
	plt.show()

	"""Trajectory Plot"""
	plt.plot(ds,hs,color="red")
	plt.xlabel("Position (m)")
	plt.ylabel("Height (m)")
	plt.title(f"Turkey Flight Trajectory (Time Step: {Δt})")
	plt.xlim(0, 40)
	plt.ylim(-5, 35)
	plt.grid(True)
	plt.show()
	



def main():
	d0 = 0
	v0 = 0
	h0 = 0
	m0 = 0
	x0 = np.array([d0,v0,h0,m0])
	displayResults(.00001,x0,t0=0)



if __name__ == "__main__":
	main()