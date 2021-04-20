#%%
""" 
Code Author: B. C. Chap
Course: UCD PHY 104B Computational Methods of Mathematical Physics
Instructor: D. Ferenc
Date: Winter Quarter 2021
Textbook: Computational Methods in Physics and Engineering by Samuel S. M. Wong

Topic:
    Methods of Least Squares - The Central Limit Theorem
References:
    Chapter 6 Methods of Least Squares
    Problem 6-3:
        Use a random number generator with an even distribution in the range [—1, +1] 
        to produce n = 6 values and store the sum as x. Collect 1000 such sums and 
        plot their distribution. Compare the results with a normal distribution of the 
        same mean and variance as the x collected. Calculate the chi^2 value. Repeat the 
        calculations with n = 50. Compare the two chi^2 obtained.
"""

#%% IMPORTING MODULES:
    
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
plt.style.use('Solarize_Light2')

#%% INITIALIZING VALUES:
    
n = 1 #number of values to sum, changing this value determines the distribution
lower_bound = -1
upper_bound = 1
observations = 1000 #number of times the value is measured/observed 
BINS = 25
mu = 0

#%% DEFINING SAMPLE DISTRIBUTIONS:
    
def even(events): #even is taken to mean uniform    
    ray = []
    for i in range(events):
        ray.append(np.random.uniform(lower_bound, upper_bound))
    sumray = np.sum(ray)
    return sumray

evenray = []
for i in range(observations):
    evenray.append(even(n))    
 
def gauss(events):
    ray = []
    for i in range(events):
        ray.append(np.random.normal(loc = mu, scale = .5))
    sumray = np.sum(ray)
    return sumray

normray = []
for i in range(observations):
    normray.append(gauss(n))

#%% DETERMINING VARIANCE, STD DEVIATION, AND CHI^2 / REDUCED CHI^2:

heven, binseven = np.histogram(evenray, bins = BINS, range = (lower_bound, upper_bound))
hnorm, binsnorm = np.histogram(normray, bins = BINS, range = (lower_bound, upper_bound))
binsceven = (binseven[1:] + binseven[:-1])/2 #centering bins

def variance(data):
    mean = np.mean(data)
    deviations = [(x - mean) ** 2 for x in data]
    variance = np.sum(deviations) / len(data)
    return variance

evenvar = variance(evenray)
sigma = np.sqrt(evenvar)

def chisq(nobs, nexp):
    return sum((nobs-nexp)**2/nexp)
chisqnormeve = chisq(heven, hnorm)

def redchisq(chi_sq, n, DOF):
    return chi_sq/(len(n)- DOF)
redchisqnormeve = redchisq(chisqnormeve, heven, 0)

print("Chi Squared for n = {}: {}".format(n, chisqnormeve))
print("Reduced Chi Squared for n = {}: {}".format(n, redchisqnormeve))

#%% VISUALIZING SAMPLE DISTRIBUTIONS: 
    #IDEAL VS RANDOMLY GENERATED
    
fx = np.linspace(np.floor(np.min(evenray)), np.ceil(np.max(evenray)), 10**3)
fy = norm.pdf(fx, mu, sigma)

plt.figure(figsize=(10,8))
plt.hist(evenray, BINS, label = "Uniform Distribution \n $\chi^2$ = {}".format(chisqnormeve), edgecolor='black', color = "skyblue", density = True)
#plt.hist(normray, BINS, label = "Normal Distribution", edgecolor='black', color = "pink", density = True)
plt.plot(fx,fy, linewidth = 2, label = "Ideal Normal Distribution with $\sigma$ = {:.2f}".format(sigma))
plt.xlabel('Sum of {} Random Generated Number(s)'.format(n), fontsize=20)
plt.ylabel('Number of Occurences',fontsize=20)
plt.suptitle('Comparison of Normalized Uniform and Normal Distribution'.format(),fontsize=22)
plt.title('for n = {} value(s) generated between [{},{}], \n and nsum = {} observations'.format(n, lower_bound, upper_bound, observations), fontsize = 16)
plt.legend(loc='upper right')
plt.show()
