import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import scipy.stats as stats
import math

plt.style.use('Solarize_Light2') 

""" Use a random number generator with an even distribution in the range [—1, +1] 
to produce n = 6 values and store the sum as x. Collect 1000 such sums and 
plot their distribution. Compare the results with a normal distribution of the 
same mean and variance as the x collected. Calculate the chi^2 value. Repeat the 
calculations with n = 50. Compare the two chi^2 obtained.

A straight sloped line of your choice.""" 

n = 1
lower_bound = -1
upper_bound = 1
observations = 1000
BINS = 25
mu = 0
slope, yint = -10, 6


def lin_func(m, x, b):
    return np.random.normal(slope, 1/(upper_bound-lower_bound)) * x + np.random.normal(b, 1/(upper_bound-lower_bound))

xuniray = []
xnormray = []
linuniray = []
linnormray = []

for i in range(observations):
    xuniray.append(np.sum(np.random.uniform(lower_bound, upper_bound, n)))
    linuniray.append(lin_func(slope, xuniray[i], yint))
    xnormray.append(np.sum(np.random.normal(mu, 1/(upper_bound-lower_bound), n)))
    linnormray.append(lin_func(slope, xnormray[i], yint))
    
linvar = np.var(linuniray)
linsig = np.sqrt(linvar)

huni, binuni = np.histogram(linuniray, bins = BINS, range = (np.min(linuniray), np.max(linuniray)))
hnorm, binnorm = np.histogram(linnormray, bins = BINS, range = (np.min(linnormray), np.max(linnormray)))
binsc = (binuni[1:] + binuni[:-1])/2

sigmai = np.full(len(huni), linsig)
    
def chisq(nobs, nexp,sigsq):
    return sum((nobs-nexp)**2/sigsq)
chisqlin = chisq(huni, hnorm, sigmai)
print("Chi Squared for n = {}: {}".format(n, chisqlin))

def redchisq(chi_sq, n, DOF):
    return chi_sq/(len(n)- DOF)
redchilin = redchisq(chisqlin, huni, 0)
print("Reduced Chi Squared for n = {}: {}".format(n,redchilin))

fx = np.linspace(np.floor(np.min(linuniray)), np.ceil(np.max(linuniray)), 10**3)
fy = norm.pdf(fx, yint, linsig)

xline = np.linspace(np.min(xuniray), np.max(xuniray), observations)
yline = slope * xline + yint
plt.figure(figsize=(10,8))
plt.plot(xuniray, linuniray, ".", label = "Uniform Distribution Along  y = ({}$\pm .5$)*x + ({}$\pm .5$)".format(slope, yint))
plt.plot(xline, yline, linewidth = 2, label = "y = {}*x + {}".format(slope,yint))
plt.xlabel('Random Generated Numbers between [{},{}]'.format(lower_bound, upper_bound), fontsize=20)
plt.ylabel('Corresponding Random Linear Data',fontsize=20)
plt.suptitle('Randomly Generated x and y values for y = {}*x + {}'.format(slope, yint), fontsize=22)
plt.legend(loc='upper right')
plt.show()

plt.figure(figsize=(10,8))
plt.hist(linuniray, BINS, label = "Uniform Distribution \n $\chi^2$ = {:.2f}".format(chisqlin), edgecolor='black', color = "skyblue", density = True)
#plt.hist(linnormray, BINS, label = "Normal Distribution", edgecolor='black', color = "pink", density = True)
plt.plot(fx,fy, linewidth = 2, label = "Normal Distribution with \n $\mu = {:.2f}$ $\sigma$ = {:.2f}".format(yint, linsig))
plt.xlabel('Sum of {} Random Generated Number(s)'.format(n), fontsize=20)
plt.ylabel('Number of Occurences',fontsize=20)
plt.suptitle('Normalized Uniform vs Normal Distributions Along y = ({}$\pm .5$)*x + ({}$\pm .5$)'.format(slope, yint),fontsize=16)
plt.title('for n = {} value(s) generated between [{},{}], \n and nsum = {} observations'.format(n, lower_bound, upper_bound, observations), fontsize = 12)
plt.legend(loc='upper right')
plt.show()