import numpy as np 


"""Compute approximations to x(4) = e4 for the scalar 
ordinary differential equation x'(t) = x(t), x(0) = 1 using
the explicit Euler method, implicit Euler method, and trapezoidal
(Crank-Nicolson) method with the steps sizes Δt = 2,1,1/2,...1/32"""

"""Function Definitions"""
def f(t,x):
    return x

def f_prime(t,x):
    return f(t,x)

def error(xn):
    return abs(np.e ** 4 - xn)

"""Helper Methods for generation"""
def generateStepSizes():
    i = 2
    steps = []
    while i >=1/32:
        steps.append(i)
        i = i / 2
    return steps

"""Methods"""
# The explicit euler method is defined as xi+1 = xi + f(ti,xi)Δt
def explicitEuler(f,Δt,t0,x0):
    N = int(4/Δt)
    for _ in range(N):
        x1 = x0 + f(t0,x0) * Δt
        x0 = x1
        t0 = t0 + Δt
    return x1

# The implicit euler method is defined as xi+1 = xi + f(ti+1,xi+1)Δt
def implicitEuler(Δt,x0,t0):
    # Calculated closed form expression on paper
	def xi_1(xi):
		return xi / (1 - Δt)
	N = int(4/Δt)
	t = t0
	x = x0
	for _ in range(N):
		x_next = xi_1(x)
		t += Δt
		x = x_next
	return x_next

def trapezoidal(Δt,x0,t0):
    # Calculated closed form expression on paper	
	def xi_1(xi):
		return xi * ((1 + Δt/2) / (1 - Δt/2))
	N = int(4/Δt)
	t = t0
	x = x0
	for _ in range(N):
		x_next = xi_1(x)
		t += Δt
		x = x_next
	return x_next


"""Results generation and visualization"""

def printResults(steps):
    for Δt in steps:
        print(f"Explicit Euler for steps size {Δt}:", explicitEuler(f,Δt,t0=0,x0=1))
        if (Δt != 1):
             print(f"Implicit Euler for steps size {Δt}:", implicitEuler(Δt,t0=0,x0=1))
        if (Δt != 2):
             print(f"Trapezoidal method for steps size {Δt}:", trapezoidal(Δt,t0=0,x0=1))
       
# def graphResults(steps,t0,x0):
#     ts = [t0]
#     xs_explicit = [x0]
#     xs_implicit = [x0]
#     xs_trapezoidal = [x0]
#     for Δt in steps:
        
        
    
def main():
    steps = generateStepSizes()
    printResults(steps)

    
if __name__ == "__main__":
    main()




