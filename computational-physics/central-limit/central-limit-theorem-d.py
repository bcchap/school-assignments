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
Goal:
    Repeat the same process for Problem 6-3 with a Mickey-mouse shaped distribution instead of even
"""

#%% IMPORTING MODULES:

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
plt.style.use('Solarize_Light2') 

#%% INITIALIZING VALUES:
    
lower_bound = -2
upper_bound = 2
n = 1 #number of values to sum, changing this value determines the distribution
nrand = 10**4
nsum = 10**4
BINS = 20

#%% DEFINING THE SAMPLE DISTRIBUTION:
    
x = np.random.uniform(lower_bound, upper_bound, nrand)
y = np.random.uniform(lower_bound, upper_bound, nrand)

head = x**2 + y**2
inhead = (head <= 1)

lear = ((x-1)**2 + (y-1)**2)
inlear = (lear <=.25)
rear = ((x+1)**2 + (y-1)**2)
inrear = (rear <=.25)
out = np.logical_not(inhead) & np.logical_not(inlear) & np.logical_not(inrear)

#%% VISUALIZING THE SAMPLE DISTRIBUTION:
    #FROM UNIFORM AND MICKEY CONSTRAINED SAMPLES
plt.figure(figsize =(10,7))
plt.plot(x[out], y[out], ".", label = "outside")
plt.plot(x[inhead], y[inhead], "c.", label = "inside")
plt.plot(x[inlear], y[inlear], "c.")
plt.plot(x[inrear], y[inrear], "c.")
plt.legend(loc="best")
plt.xlabel('Random Generated Numbers between [{},{}]'.format(lower_bound, upper_bound), fontsize=20)
plt.ylabel('Random Generated Numbers between [{},{}]'.format(lower_bound, upper_bound), fontsize=20)
plt.suptitle('Randomly Generated x and y Values\n within a Mickey Mouse Distribution'.format(), fontsize=22)
plt.show()

#%% REPEATEDLY RANDOMLY GENERATING THE SAMPLE DISTRIBUTION:
    
xray, yray, heads, inheads, lears, inlears, rears, inrears, out, inside = ([] for i in range(10))
for i in range(nsum):
    xray.append(np.random.uniform(lower_bound, upper_bound, n))
    yray.append(np.random.uniform(lower_bound, upper_bound, n))
    heads.append(xray[i]**2 + yray[i]**2)
    inheads.append(heads[i]<=1)
    lears.append((xray[i]-1)**2 + (yray[i]-1)**2)
    inlears.append(lears[i] <=.25)
    rears.append((xray[i]+1)**2 + (yray[i]-1)**2)
    inrears.append(rears[i] <=.25)
    out.append(np.mean(np.logical_not(inheads[i]) & np.logical_not(inlears[i]) & np.logical_not(inrears[i])))
    inside.append(nsum - out[i])

mu = np.mean(inside)
sigma = np.sqrt(np.var(inside))
fx = np.linspace((np.min(inside)), (np.max(inside)), 10**3)
fy = norm.pdf(fx, mu, sigma)

#%% VISUALIZING SUMMED SAMPLES: 
    #FROM UNIFORM AND MICKEY CONSTRAINED SAMPLES
    
plt.figure(figsize=(10,7))
plt.hist(inside, BINS, label = "Sum of n = {} values within\n Mickey Mouse Distribution".format(n), edgecolor='black', color = "skyblue", density = True)
plt.plot(fx,fy, linewidth = 2, label = "Ideal Normal Distribution\n with $\mu$ {:.2f} and $\sigma$ = {:.2f}".format(mu, sigma))
plt.xlabel('Sum of {} Random Generated Number(s)'.format(n), fontsize=20)
plt.ylabel('Number of Occurences',fontsize=20)
plt.suptitle('Comparison of Normalized Uniform and Normal Distribution '.format(),fontsize=22)
plt.title('for value(s) generated within a Mickey Mouse Distribution'.format(), fontsize = 16)
plt.legend(loc='upper right')
plt.show()

