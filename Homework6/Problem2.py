import numpy as np
import time
import matplotlib.pyplot as plt
# In this case, x is a vector with components(v(t), (T(t)+G(t))/m(t))
# Returns f_d,f_v as a vector
def f(t,x): 
	d,v = x
	f_d = v

	"""Internal Function Definitions"""
	T = 12 * 10e6 if t < 600 else 0
	M = (1 - ((.9*t)/600)) * 10e6 if t < 600 else .1 * 10e6
	G = -(6371000)**2 * ((10 * M) /d**2)

	f_v = ((T + G) / M)
	return np.array([f_d,f_v])
    

# The explicit euler method is defined as xi+1 = xi + f(ti,xi)Δt
def explicitEuler(Δt,x0,t0):
    N = int(36000/Δt)
    xs = [x0]
    ts = [t0]
    for i in range(N):
        if i % int(.25 * N) == 0 and i != 0:
            print(f"{int(i / N * 100)} Percent Completed!")
        f_vector = f(t0,x0)
        x1 = x0 + f_vector * Δt
        x0 = x1
        t0 = t0 + Δt
        ts.append(t0)
        xs.append(x1)
    print("Completed!")
    return xs,ts



def displayResults(Δt,x0,t0):
	start = time.perf_counter()
	xs,ts = explicitEuler(Δt,x0,t0)
	ds = [(position[0]/1000) for position in xs]
	vs = [(velocity[1]/1000) for velocity in xs]
	end = time.perf_counter()
	print(f"Height for steps size {Δt}s: {ds[-1]} km. Took {end-start} seconds.")

	plt.plot(ts,ds,label="Rocket Position")
	plt.xlabel("Time (seconds)")
	plt.ylabel("Position (km)")
	plt.title(f"Rocket Trajectory: Position vs Time (time step {Δt})")
	plt.legend()
	plt.grid(True)
	plt.show()
 
	plt.plot(ts,vs,color="red",label="Rocket Velocity")
	plt.xlabel("Time (seconds)")
	plt.ylabel("Velocity (km)")
	plt.title(f"Rocket Trajectory: Velocity vs Time (time step {Δt})")
	plt.legend()
	plt.grid(True)
	plt.show()

	

def main():
	d0 = 6371000
	v0 = 0
	x0 = np.array([d0,v0])
	displayResults(.001,x0,t0=0)


    
if __name__ == "__main__":
    main()