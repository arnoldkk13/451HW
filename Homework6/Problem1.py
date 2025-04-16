import numpy as np 
import matplotlib.pyplot as plt

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
def generateStepSizes(stepSize):
    i = 2
    steps = []
    while i >= stepSize:
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

# Trendlines are found here based off of visualization 
def trendLineTrapezoidal(x):
    C = 270
    N = 2
    return C * (x ** (-1 * N))

def trendLineEuler(x):
    C = 430
    N = 1
    return C * (x ** (-1 * N))
     

def graphResults(t_exp,t_imp,t_trap, x_exp,x_imp,x_trap, stepSize):
    plt.plot(t_exp,x_exp, label="Explicit Euler")
    plt.plot(t_imp,x_imp, label="Implicit Euler")
    plt.plot(t_trap,x_trap,label="Trapezoidal Rule")
    
    # Emplace trendlines
    trend_x = np.linspace(2,t_exp[-1],300)
    euler_trend = [trendLineEuler(x) for x in trend_x]
    trap_trend = [trendLineTrapezoidal(x) for x in trend_x]

    plt.plot(trend_x,euler_trend, '--',label="TrendLine = 430x^-1")
    plt.plot(trend_x,trap_trend,'--',label="TrendLine = 270x^-2")

    plt.title(f"ODE Solvers' Function Evaluations vs Error (Step Size {stepSize:2f})")
    plt.xlabel("Function Evaluations")
    plt.ylabel("Error")
    plt.loglog()
    plt.legend()
    plt.grid(True)
    plt.show()

def displayResults(steps):
    # Different time steps are required to be stored as there are values where implicit and trapezoidal fails
    t_explicit = []
    t_implicit = []
    t_trapezoid = []

    x_explicit = []
    x_implicit = []
    x_trapezoid = []

    for Δt in steps:
        N = int(4/Δt)

        # Print to terminal all the errors for evaluation
        print(f"Explicit Euler for steps size {Δt}:", error(explicitEuler(f,Δt,t0=0,x0=1)))
        if (Δt != 1):
             print(f"Implicit Euler for steps size {Δt}:", error(implicitEuler(Δt,t0=0,x0=1)))
        if (Δt != 2):
             print(f"Trapezoidal method for steps size {Δt}:", error(trapezoidal(Δt,t0=0,x0=1)))
        

        t_explicit.append(N)
        x_explicit.append(error(explicitEuler(f,Δt,t0=0,x0=1)))
        if (Δt != 1):
            t_implicit.append(N)
            x_implicit.append(error(implicitEuler(Δt,t0=0,x0=1)))
        if (Δt != 2):
            t_trapezoid.append(N)
            x_trapezoid.append(error(trapezoidal(Δt,t0=0,x0=1)))

        # Add tabulation if time
        # tabulateResults()

        # Only graph when there are sufficient results
        if Δt < 1/16:
            graphResults(t_explicit, t_implicit, t_trapezoid, x_explicit, x_implicit, x_trapezoid,Δt)
       
        
        
    
def main():
    steps = generateStepSizes(10**-6)
    displayResults(steps)

    
if __name__ == "__main__":
    main()




