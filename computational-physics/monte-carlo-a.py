"""
PHY 104B Project 3 Part A
Bottreypich Cassey Chap
March 5, 2021
Winter Quarter 2021
"""

import random
import matplotlib.pyplot as plt
import numpy as np
plt.style.use("bmh")

TOLERANCE = 10**-5
INITIAL_N = int(40)
N = int(40)
SEED = int(1)
random.seed(SEED)
N_list = []

A_1 = 0.6826895

MAX_IT = 15
ITERATIONS = np.arange(1, MAX_IT+1)

#defining the function we want to integrate, which is known as the normal probability function
def NPF(x):
    return np.sqrt(2/np.pi) * np.e**(-1/2*x**2)

#Want to put the values from NPF(x_rand) into an array to be summed eventually
def NPF_rand(rand, NPF_x):
    for i in rand:
        NPF_x.append(NPF(i))
        
NPF_MEANS = []        
for i in ITERATIONS:
    rand_x_0 = []
    rand_x_1 = []
    
    for j in range(N):
        rand_x_0.append(random.random())
        rand_x_1.append(random.random())
        
    NPF_0 = []
    NPF_1 = []    
    NPF_rand(rand_x_0, NPF_0)
    NPF_rand(rand_x_1, NPF_1)
    INTEGRANDS = [np.mean(NPF_0), np.mean(NPF_1)]
    
    print("Iteration: {}".format(str(i)))
    print("N: {}".format(str(N)))
    print("NPF 1: {}".format(str(INTEGRANDS[0])))
    print("NPF 2: {}".format(str(INTEGRANDS[1])))
     

    NPF_MEANS.append(np.mean(INTEGRANDS))
    
    if np.abs(INTEGRANDS[1]-INTEGRANDS[0]) > TOLERANCE:
       print("Mean of NPF Integrals: {}".format(NPF_MEANS[-1]))
       N_list.append(N)
       N *= 2
    else:
        print("A(x = 1): {}".format(NPF_MEANS[-1]))
        break
        
PCT_ERR = abs((NPF_MEANS[-1]-A_1)/A_1)
print("Percent Error: {}%".format(PCT_ERR*100))


monte_carlo_a = plt.figure(figsize = (10,6))
plt.plot(ITERATIONS[0:len(NPF_MEANS)], NPF_MEANS, label='A(x = 1) vs Iterations')
plt.axhline(A_1,color='k',label='Exact A(x = 1): {}'.format(A_1))
plt.suptitle('A(x = 1) vs Iterations\nSeed: {}, Initial N: {}'.format(str(SEED), str(INITIAL_N)))
plt.legend()
plt.xlabel('Iteration')
plt.ylabel('A(x = 1)')
plt.show()

"""
plt.figure()
plt.plot(np.arange(655360),rand_x_0)
plt.show
"""