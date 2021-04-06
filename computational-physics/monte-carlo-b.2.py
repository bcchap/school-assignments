"""
PHY 104B Project 3 Part B Limit
Bottreypich Cassey Chap
March 5, 2021
Winter Quarter 2021
"""
import random
import matplotlib.pyplot as plt
import numpy as np
plt.style.use("bmh")

#Initialize the random number generator
TOLERANCE = 10**-3

#Let N be starting number of poits for integration
INITIAL_N = int(100)
N = int(100)
N_list = []
#Input seed for random number generator
SEED = int(1)
random.seed(SEED)

Q_ACCEPTED = 1

#Zero the iteration counter
MAX_IT = 15
ITERATIONS = np.arange(1, MAX_IT+1)

#defining the function we want to integrate, which is known as the Yukawa charge distribution
#We are using symmetry to just do 1/8 of the sphere
def Q(a, b, c):
    r = np.sqrt(a**2+b**2+c**2)
    return 1/(np.pi*8) * np.e**(-r) * 8000

#Want to put the values from NPF(x_rand) into an array to be summed eventually
def Q_rand(rand_a, rand_b, rand_c, Q_r):
    for i in range(N):
        Q_r.append(Q(rand_a[i], rand_b[i], rand_c[i]))
        
Q_MEANS = []   
     
for i in ITERATIONS:
    rand_x_0 = []
    rand_y_0 = []
    rand_z_0 = []
    rand_x_1 = []
    rand_y_1 = []
    rand_z_1 = []
    
    for j in range(N):
        rand_x_0.append(random.uniform(-10,10))
        rand_y_0.append(random.uniform(-10,10))
        rand_z_0.append(random.uniform(-10,10))
        rand_x_1.append(random.uniform(-10,10))
        rand_y_1.append(random.uniform(-10,10))
        rand_z_1.append(random.uniform(-10,10))
              
    YUKAWA_0 = []
    YUKAWA_1 = []
    Q_rand(rand_x_0, rand_y_0, rand_z_0, YUKAWA_0)
    Q_rand(rand_x_1, rand_y_1, rand_z_1, YUKAWA_1)
    INTEGRANDS = [np.mean(YUKAWA_0), np.mean(YUKAWA_1)]
    

    print("Iteration: {}".format(str(i)))
    print("N: {}".format(str(N)))
    print("Q 1: {}".format(str(INTEGRANDS[0])))
    print("Q 2: {}".format(str(INTEGRANDS[1])))
    Q_MEANS.append(np.mean(INTEGRANDS))
    
#    print(YUKAWA_0)
#    print(YUKAWA_1)
     
    if np.abs(INTEGRANDS[1]-INTEGRANDS[0]) > TOLERANCE:
       print("Mean of Q: {}".format(Q_MEANS[-1]))
       N_list.append(N)
       N *= 2
    else:
        print("Q: {}".format(Q_MEANS[-1]))
        break

PCT_ERR = abs((Q_MEANS[-1]-Q_ACCEPTED)/Q_ACCEPTED)
print("Percent Error: {}%".format(PCT_ERR*100))

monte_carlo_b = plt.figure(figsize = (10,6))
plt.plot(ITERATIONS[0:len(Q_MEANS)], Q_MEANS, label='Q vs Iterations')
plt.axhline((Q_ACCEPTED),color='k',label='Exact Q: {}'.format(Q_ACCEPTED))
plt.suptitle('Q vs Iterations\nSeed: {}, Initial N: {}'.format(str(SEED), str(N_list[0])))
plt.legend()
plt.xlabel('Iteration')
plt.ylabel('Q')
plt.show()
