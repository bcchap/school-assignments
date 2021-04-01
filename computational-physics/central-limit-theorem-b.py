import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import scipy.stats as stats

plt.style.use('Solarize_Light2')

""" Use a random number generator with an even distribution in the range [—1, +1] 
to produce n = 6 values and store the sum as x. Collect 1000 such sums and 
plot their distribution. Compare the results with a normal distribution of the 
same mean and variance as the x collected. Calculate the chi^2 value. Repeat the 
calculations with n = 50. Compare the two chi^2 obtained.""" 

n = 1
lower_bound = -1
upper_bound = 1
observations = 1000
BINS = 25
mu = 0

def even(events):
    ray = []
    for i in range(events):
        ray.append(np.random.uniform(lower_bound, upper_bound))
    sumray = np.sum(ray)
    return sumray

evenray = []
for i in range(observations):
    evenray.append(even(n))
    
def variance(data):
    mean = np.mean(data)
    deviations = [(x - mean) ** 2 for x in data]
    variance = np.sum(deviations) / len(data)
    return variance

evenvar = variance(evenray)
sigma = np.sqrt(evenvar)
 
def gauss(events):
    ray = []
    for i in range(events):
        ray.append(np.random.normal(loc = mu, scale = .5))
    sumray = np.sum(ray)
    return sumray

normray = []
for i in range(observations):
    normray.append(gauss(n))

   
heven, binseven = np.histogram(evenray, bins = BINS, range = (lower_bound, upper_bound))
hnorm, binsnorm = np.histogram(normray, bins = BINS, range = (lower_bound, upper_bound))
binsceven = (binseven[1:] + binseven[:-1])/2

def chisq(nobs, nexp):
    return sum((nobs-nexp)**2/nexp)
chisqnormeve = chisq(heven, hnorm)

def redchisq(chi_sq, n, DOF):
    return chi_sq/(len(n)- DOF)
redchisqnormeve = redchisq(chisqnormeve, heven, 0)

print("Chi Squared for n = {}: {}".format(n, chisqnormeve))
print("Reduced Chi Squared for n = {}: {}".format(n, redchisqnormeve))


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
