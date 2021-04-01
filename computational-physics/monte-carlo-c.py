"""
PHY 104B Project 3 Part C
Bottreypich Cassey Chap
March 5, 2021
Winter Quarter 2021
"""
import random
import matplotlib.pyplot as plt
import numpy as np
plt.style.use("bmh")

SEED = int(1)
random.seed(SEED)
INITIAL_N = int(50)
N = int(50)
N_list = []

B = float(1.0)

TOLERANCE = 10**-5  
DIFFERENCE = 10**-8

MAX_IT = 15
ITERATIONS = np.arange(1, MAX_IT+1)
THROWAWAY = 0

TRU_VAL = np.arcsin(B)

def asin(x):    
    return B*1/np.sqrt(1-x**2)

def asin_rand(rand, asin_x):
    for i in rand:
        asin_x.append(asin(i))

ASIN_MEANS = []

for i in ITERATIONS:
    rand_x_0 = []
    rand_x_1 = []
    
    for j in range(N):
        rand_x_0.append(random.uniform(0,B))
        
        if np.abs(rand_x_0[j]-1) <= DIFFERENCE: 
            rand_x_0.append(random.uniform(0,B))
            THROWAWAY += 1
            if np.abs(rand_x_0[j]-1) > DIFFERENCE:
                break
            
        rand_x_1.append(random.uniform(0,B))   
        
        if np.abs(rand_x_1[j]-1) <= DIFFERENCE: 
            rand_x_1.append(random.uniform(0,B))
            THROWAWAY += 1
            if np.abs(rand_x_1[j]-1) > DIFFERENCE:
                break  
            
    ASIN_0 = []
    ASIN_1 = []    
    asin_rand(rand_x_0, ASIN_0)
    asin_rand(rand_x_1, ASIN_1)
    INTEGRANDS = [np.mean(ASIN_0), np.mean(ASIN_1)]
    
    print("Iteration: {}".format(str(i)))
    print("N: {}".format(str(N)))
    print("I 1: {}".format(str(INTEGRANDS[0])))
    print("I 2: {}".format(str(INTEGRANDS[1])))
     

    ASIN_MEANS.append(np.mean(INTEGRANDS))
    
    if np.abs(INTEGRANDS[1]-INTEGRANDS[0]) > TOLERANCE:
       print("Mean of I(B): {}".format(ASIN_MEANS[-1]))
       N_list.append(N)
       N *= 3
    else:
        print("I(B): {}".format(ASIN_MEANS[-1]))
        break

print("# of Values Thrown Away: {}".format(THROWAWAY))
PCT_ERR = abs((ASIN_MEANS[-1]-TRU_VAL)/TRU_VAL)
print("Percent Error: {}%".format(PCT_ERR*100))

monte_carlo_c = plt.figure()
plt.plot(ITERATIONS[0:len(ASIN_MEANS)], ASIN_MEANS, label='I vs Iterations')
plt.axhline(np.arcsin(B),color='k', label='I({}) = arcsin({})'.format(B,B))
plt.suptitle('I vs Iterations\nSeed: {} & initial N: {}'.format(str(SEED), str(INITIAL_N)))
plt.legend()
plt.xlabel('Iteration')
plt.ylabel('I')
plt.show()

