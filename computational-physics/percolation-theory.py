import numpy as np
import matplotlib.pyplot as plt
import math

"""
Determining the coefficients of the Hermite polynomials by recursion relations
"""



#highest order of the Hermite polynomial as a two-dimensional array of coefficients 
NDMN = 11
NARY = np.zeros((NDMN+1,NDMN+1))

# initializing the coefficients of the Hermite polynomial
#define coefficient a_0,0
NARY[0][0] = 1
#define coefficient a_1,0
NARY[1][0] = 0
#define coefficient a_1,1
NARY[1][1] = 2

#define the recursion relationship for determining the coefficients based on the initial
def coeff_NARY():
    coeff = 2*NARY[n][k-1] - 2*n*NARY[n-1][k]
    return coeff

#insert the value of the coefficient into the initalized array using the recursion function coeff_NARY  
for k in np.arange(0, NDMN+1):
    for n in np.arange(0, NDMN):
        NARY[n+1][k] = coeff_NARY()
        
print(NARY)

#define a function for hermite polynomials
def hermite_POLY(n,rho):
        H_n = 0 #initialize the ground state energy level n = 0
        for k in np.arange(0, n+1):
            H_n = H_n + NARY[n][k]*rho**k
        return H_n

#create values of rho, a dimensionless quantity that replaces constants multiplied by a distance varaible x
rho = np.linspace(-5,5,100)

#determine the output of the hermite polynomials based on rho
#to visually examine the relationship between them, divide by n^2, the energy state
H_0, H_1, H_2, H_3, H_4 = ([] for i in np.arange(0,5))
for i in rho:
    H_0.append(hermite_POLY(0,i)/(1**3))
    H_1.append(hermite_POLY(1,i)/(1**3))
    H_2.append(hermite_POLY(2,i)/(2**3))
    H_3.append(hermite_POLY(3,i)/(3**3))
    H_4.append(hermite_POLY(4,i)/(4**3))
    
#print(H_0)
#print(H_1)
#print(H_2)
#print(H_3)
#print(H_4)

plt.figure()    
#plt.plot(rho, H_0, label="$H_0/{⍴}^3$")
plt.plot(rho, H_1, label="$H_1/{1}^3$")
plt.plot(rho, H_2, label="$H_2/{2}^3$")
plt.plot(rho, H_3, label="$H_3/{3}^3$")
plt.plot(rho, H_4, label="$H_4/{4}^3$")
plt.xlim(0,3)
plt.ylim(-1,6)
plt.xlabel("$⍴$ [Unitless]")
plt.ylabel("$H_n/⍴ [Unitless]$")
plt.title("Hermite Polynomials $H_n(⍴)/n^3$ for n = 1, 2, 3, 4")
plt.grid()
plt.legend()
plt.show()



"""
Assessing theoretical limitations of the Hermite polynomial in relation to the harmonic oscillator
Compare and assess the behavior of a particle confined to a well depending on its energy state n
We see that a particle with greater values of n become highly excited, thus are more likely to escape the
potential well confinement. This concept is called the Correspondence Principle. This is in agreement with
the theoretical model. Based on the Uncertainty Principle, we know that a particle is more likely to be found
outside of the potential well as its energy increases. 
"""



#define the wave function for a simple harmonic oscillator, where n represents the energy state
def psi(n,rho):
    return (1 / (np.sqrt(2**n*math.factorial(n)*np.sqrt(np.pi)))) * np.exp(-(rho**2) / 2)

psi_0, psi_1, psi_2, psi_3, psi_4 = ([] for i in np.arange(0,5))
for x in rho:
    psi_0.append(psi(0,x)*(hermite_POLY(0,x)/(1**3)))
    psi_1.append(psi(1,x)*(hermite_POLY(1,x)/(1**3)))
    psi_2.append(psi(2,x)*(hermite_POLY(2,x)/(2**3)))
    psi_3.append(psi(3,x)*(hermite_POLY(3,x)/(3**3)))
    psi_4.append(psi(4,x)*(hermite_POLY(4,x)/(4**3)))
    
#print(psi_0)
#print(psi_1)
#print(psi_2)
#print(psi_3)
#print(psi_4)

plt.figure()
plt.plot(rho, psi_0, label="n=0")
plt.plot(rho, psi_1, label="n=1")
plt.plot(rho, psi_2, label="n=2")
plt.plot(rho, psi_3, label="n=3")
plt.plot(rho, psi_4, label="n=4")
plt.xlabel("$⍴$ [Unitless]")
plt.ylabel("$Ψ_n(⍴)$ [Energy]")
plt.title("1-D Harmonic Oscillator Wave Function $Ψ_n(⍴)$ for n = 0, 1, 2, 3, 4")
plt.grid()
plt.legend()
plt.show()

H_4_compare = []
H_4_compare.append(hermite_POLY(4,x))
H_10 = []
H_10.append(hermite_POLY(10,x)/(10**3))
psi_4_compare = []
psi_10 = []
for x in rho:
    psi_4_compare.append(psi(4,x)*hermite_POLY(4,x))
    psi_10.append(psi(10,x)*(hermite_POLY(10,x)))
    
plt.figure()
plt.plot(rho, psi_4_compare, label="n=4")
plt.plot(rho, psi_10, label="n=10")
plt.xlabel("$⍴$ [Unitless]")
plt.ylabel("$Ψ_n(⍴)$ [Energy]")
plt.title("1-D Harmonic Oscillator Wave Function $Ψ_n(⍴)$ for n = 4 compared to n = 10")
plt.grid()
plt.legend()
plt.show()



