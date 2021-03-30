import math
import matplotlib.pyplot as plt
import numpy as np
from scipy import constants as cnst
from scipy import optimize
from scipy import integrate
from decimal import Decimal, getcontext
plt.style.use('Solarize_Light2') 

""" Solving the Shroedinger's Equation for a particle in a finite 1-dimensional 
rectangular potential well where L = 7 fm and u_0 = 10 MeV, find all possible
energy levels and their corresponding wave functions, then compare to the 
infinite potential well.
"""



#DETERMINING ENERGY LEVELS

#define constants
u_0 = (10*10**6) / (6.242*10**18)  #eV to J
L = 7 * 10**-15 # m
STEP = 10**4 #steps between an interval

#define k, a, cotangent
def k(E):
    return ((2*cnst.proton_mass*E)/(cnst.hbar**2))**(1/2) 
def a(E):
    return ((2*cnst.proton_mass*(u_0-E))/(cnst.hbar**2))**(1/2)
def cot_func(E):
    return 2/(math.tan(k(E)*L))

#find energy E approx 0 = 2cot(k*L) - k/a + a/k, for values of 0 > E > u_0
def transc(E):
    return(cot_func(E) - k(E)/a(E) + a(E)/k(E))
energy = np.linspace(u_0/10**3,u_0-(u_0/10**3),100)
transc_out = []
for x in energy:
    transc_out.append(transc(x))

#visually determine where the expected roots should occur to determine ranges for optimize.root_scalar
"""transc_plot = plt.figure()
plt.suptitle("Energy vs Transcendental Output")
plt.plot(energy,  transc_out, label = "Transcendental Output")
#plt.vlines(r1.root, -150, 100, label = "Energy = $3.288 * 10^-13$ J", color = "r")
#plt.vlines(r2.root, -150, 100, label = "Energy = $1.195 * 10^-12$ J", color = "r")
plt.xlabel("Energy E [J]")
plt.ylabel("Transcendental Output [Unitless]")
plt.grid()
plt.legend(loc="best")
plt.show(transc_plot)"""

#determine roots (where the transcendental output is = 0) numerically
r1 = optimize.root_scalar(transc, bracket=[.2*10**(-12),.6*10**(-12)],method='bisect', xtol=10**(-16))
print("First Energy [J] = " , r1)
r2 = optimize.root_scalar(transc, bracket=[1.0*10**(-12),1.4*10**(-12)],method='bisect', xtol=10**(-16))
print("Second Energy [J] = " , r2)

#assess to what extent root values satisfy the transcendental equation
print("First Energy Root satisfies the transcendental equation up to ", transc(r1.root), " from 0")
print("Second Energy Root satisfies the transcendental equation up to ", transc(r2.root), "from 0")

#confirm and visualize energy roots in relation to transcendental output
transc_plot = plt.figure()
plt.suptitle("Energy vs Transcendental Output")
plt.plot(energy,  transc_out, label = "Transcendental Output")
plt.vlines(r1.root, -150, 100, label = "Energy = $3.288 * 10^-13$ J", color = "r")
plt.vlines(r2.root, -150, 100, label = "Energy = $1.195 * 10^-12$ J", color = "g")
plt.xlabel("Energy E [J]")
plt.ylabel("Transcendental Output [Unitless]")
plt.legend(loc="best")
plt.show(transc_plot)



#DETERMINING CORRESPONDING WAVE FUNCTIONS FOR ENERGY LEVELS

#defining precision
getcontext().prec = 30
getcontext().Emin = -999999999999999999
getcontext().Emax = 999999999999999999

#determining constants of finite potential well wave function 
#constant C

#use previous energies to establish corresponding values of alpha and kappa from functions a and k
energies = [r1.root, r2.root]
alpha = [a(r1.root), a(r2.root)]
kappa = [k(r1.root), k(r2.root)]

#redefine all wave function equations in terms of C, eliminating other constants A, B, and G
#utilize the normalization condition to determine constant C
def constant_C(E,alpha,kappa):
        psi1_C = lambda x: (math.e**(alpha*x)) ** 2
        psi1_C_integrate, error1 = integrate.quad(psi1_C, -np.inf, 0)
        psi2_C = lambda x: ((alpha / kappa) * math.sin(kappa * x) + math.cos(kappa * x)) ** 2
        psi2_C_integrate, error2 = integrate.quad(psi2_C, 0, L)
        psi3_C = lambda x: (((alpha / kappa) * math.sin(kappa * L) + math.cos(kappa * L)) * math.e**(alpha * L) * math.e**(-alpha * x)) ** 2
        psi3_C_integrate, error3 = integrate.quad(psi3_C, L, np.inf)
        normie_C = (1 / (psi1_C_integrate + psi2_C_integrate + psi3_C_integrate)) ** (1 / 2)
        return normie_C
C = []
C1 = Decimal(constant_C(energies[0],alpha[0],kappa[0]))
C.append(C1)
C2 = Decimal(constant_C(energies[1],alpha[1],kappa[1]))
C.append(C2)

#determine constants A, B, and G based on constant C

#constant A
def constant_A(C,alpha,kappa):
    normie_A = Decimal(alpha)/Decimal(kappa)*C
    return normie_A
A = []
A1 = Decimal(constant_A(C[0],alpha[0],kappa[0]))
A.append(A1)
A2 = Decimal(constant_A(C[1],alpha[1],kappa[1]))
A.append(A2)

#constant B
def constant_B(C):
    normie_B = C
    return normie_B
B = []
B1 = Decimal(constant_B(C[0]))
B.append(B1)
B2 = Decimal(constant_B(C[1]))
B.append(B2)

#constant G
def constant_G(A,B,alpha,kappa):
    normie_G = (A*Decimal(math.sin(Decimal(kappa)*Decimal(L))) + B*Decimal(math.cos(Decimal(kappa)*Decimal(L)))) * Decimal(math.e)**(Decimal(alpha)*Decimal(L))
    return normie_G
G = []
G1 = Decimal(constant_G(A[0],B[0],alpha[0],kappa[0]))
G.append(G1)
G2 = Decimal(constant_G(A[1],B[1],alpha[1],kappa[1]))
G.append(G2)

print("A1,A2 = :", A)    
print("B1,B2 = :", B)
print("C1,C2 = :", C)
print("G1,G2 = :", G)

#plug in A, B, C, and G to wave function equation and plot for range of lengths
# spanning slightly beyond the potential well
#region1
psi1_x = np.linspace(-L,0,STEP)
psi1_e1 = []
for x in psi1_x:
    x = Decimal(x)
    psi1_e1.append(C1*Decimal(math.e)**(Decimal(alpha[0])*x))
psi1_e2 = []
for x in psi1_x:
    x = Decimal(x)
    psi1_e2.append(C2*Decimal(math.e)**(Decimal(alpha[1])*x))

#region 2
psi2_x = np.linspace(0,L,STEP)
psi2_e1 = []
for x in psi2_x:
    x = Decimal(x)
    psi2_e1.append(A1*Decimal(math.sin(Decimal(kappa[0])*x)) + B1*Decimal(math.cos(Decimal(kappa[0])*x)))
psi2_e2 = []
for x in psi2_x:
    x = Decimal(x)
    psi2_e2.append(A2*Decimal(math.sin(Decimal(kappa[1])*x)) + B2*Decimal(math.cos(Decimal(kappa[1])*x)))
        
#region 3
psi3_x = np.linspace(L,2*L,STEP)
psi3_e1 = []
for x in psi3_x:
    x = Decimal(x)
    psi3_e1.append(G1*Decimal(math.e)**(Decimal(-alpha[0])*x))
psi3_e2 = []
for x in psi3_x:
    x = Decimal(x)
    psi3_e2.append(G2*Decimal(math.e)**(Decimal(-alpha[1])*x))

fin_well = plt.figure()
plt.suptitle("Finite Well Wave Function")
plt.hlines(u_0,-L,0)
plt.hlines(u_0,L,2*L)
plt.vlines(0,0,u_0)
plt.vlines(L,0,u_0)
plt.plot(psi1_x, psi1_e1,'tab:red', label = "Wave Corresponding to Energy = $3.288 * 10 ** -13$ J")
plt.plot(psi1_x, psi1_e2,'tab:blue', label = "Wave Corresponding to Energy = $1.195 * 10 ** -12$ J")
plt.plot(psi2_x, psi2_e1,'tab:red')
plt.plot(psi2_x, psi2_e2,'tab:blue')
plt.plot(psi3_x, psi3_e1,'tab:red')
plt.plot(psi3_x, psi3_e2,'tab:blue')
plt.xlabel("Length [km]")
plt.ylabel("Energy [J]")
plt.legend()
plt.show(fin_well)

#infinite well comparison
psi_length = np.linspace(0,L,STEP)

def E_inf(n):
    return (n*(np.pi)*cnst.hbar/L)**2*(1/(2*cnst.proton_mass))

A_inf = math.sqrt(2/L)
def psi_inf_eq(n,l):
    return A_inf*math.sin(n*l*np.pi/L)

psi_inf_1 = []
psi_inf_2 = []
psi_inf_3 = []
psi_inf_4 = []

for l in psi_length:
   psi_inf_1.append(psi_inf_eq(1,l))
   psi_inf_2.append(psi_inf_eq(2,l))
   psi_inf_3.append(psi_inf_eq(3,l))
   psi_inf_4.append(psi_inf_eq(4,l))

"""inf_well = plt.figure()
plt.suptitle("Infinite Well Wave Function")
fig, subpl = plt.subplots(5,1)
subpl[3].plot(psi_length, psi_inf_1,'tab:orange', label = "n = 1")
subpl[2].plot(psi_length, psi_inf_2,'tab:green', label = "n = 2")
subpl[1].plot(psi_length, psi_inf_3,'tab:blue', label = "n = 3")
subpl[0].plot(psi_length, psi_inf_4,'tab:purple', label = "n = 4")
plt.xlabel("Length [km]")
plt.ylabel("Energy, for n = [1,2,3,4] [J]")
subpl[3].legend(loc = "upper left")
subpl[2].legend(loc = "upper left")
subpl[1].legend(loc = "upper left")
subpl[0].legend(loc = "upper left")
plt.show(inf_well)"""

fin_and_inf_well = plt.figure()
plt.suptitle("Finite and Infinite Well Wave Function")
plt.plot(psi1_x, psi1_e1, 'tab:red', label = "Finite Wave Function for Energy = $3.288 * 10 ^ {-13}$ J")
plt.plot(psi1_x, psi1_e2, 'tab:pink', label = "Finite Wave Function for Energy = $1.195 * 10 ^ {-12}$ J")
plt.plot(psi2_x, psi2_e1, 'tab:red')
plt.plot(psi2_x, psi2_e2, 'tab:pink')
plt.plot(psi3_x, psi3_e1, 'tab:red')
plt.plot(psi3_x, psi3_e2, 'tab:pink')
plt.plot(psi_length, psi_inf_1,'tab:orange', label = "Infinite Wave Function for Energy n = 1")
plt.plot(psi_length, psi_inf_2,'tab:green', label = "Infinite Wave Function for Energy n = 2")
plt.plot(psi_length, psi_inf_3,'tab:blue', label = "Infinite Wave Function for Energy n = 3")
plt.plot(psi_length, psi_inf_4,'tab:purple', label = "Infinite Wave Function for Energy n = 4")
plt.ylim(-7.0*10**7,2*10**7)
plt.xlabel("Length [km]")
plt.ylabel("Energy [J]")
plt.legend(loc = "best")
plt.show(fin_and_inf_well)


    
