import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage
import pylab

plt.style.use('Solarize_Light2') 

L = 100
probray =  [0.5,0.525,0.550,0.575,0.6,0.625]
lattice = np.random.rand(L,L)
for i in range(len(probray)):
    prob = probray[i]
    matrix_i = prob > lattice
    clust_i, nconnect = scipy.ndimage.measurements.label(matrix_i)
    area = scipy.ndimage.measurements.sum(matrix_i, clust_i, index=np.arange(clust_i.max() + 1))
    cluster = area[clust_i]
    plt.subplot(2, 3, i+1).set_title('p = {:.3f}'.format(prob), fontsize = 12)
    plt.imshow(cluster, cmap = 'summer', origin='lower')
    plt.grid()
    plt.suptitle("Clusters for Lattice Size $100x100$ with Various Occupation Probabilities")

startprob = .5
stopprob = .7
probsize = 500
sampsize = 1000

probs = np.linspace(startprob, stopprob, probsize)

def Pi(L):
    perctimes = np.zeros(probsize)
    for i in range(sampsize):
        lattice = np.random.rand(L, L)
        for p_i in range(probsize):
            matrix_i = probs[p_i] > lattice
            clust_i, nconnect = scipy.ndimage.measurements.label(matrix_i)
            perc_i = pylab.intersect1d(clust_i[0,:], clust_i[-1,:])
            perc = perc_i[np.where(perc_i>0)]
            if (len(perc)>0):
                perctimes[p_i] += 1
    percprobpi = perctimes/sampsize
    return percprobpi


plt.figure(figsize=(10,7))
#plt.plot(probs, Pi(10), "-", label = "Lattice Size = $10x10$")
plt.plot(probs, Pi(25), "-", label = "Lattice Size = $25x25$")
plt.plot(probs, Pi(50), "-", label = "Lattice Size = $50x50$")
plt.plot(probs, Pi(75), "-", label = "Lattice Size = $75x75$")
plt.plot(probs, Pi(100), "-", label = "Lattice Size = $100x100$")
plt.xlim(startprob,stopprob)
plt.ylim(0,1)
plt.vlines(0.5927, 0 , 1)
plt.xlabel('Occupation Probability $p$'.format(), fontsize=20)
plt.ylabel('Percolation Probability $\Pi$',fontsize=20)
plt.suptitle('Percolation Probability as a Function of Occupation Probability'.format(),fontsize=18)
plt.title('for various lattice sizes'.format(), fontsize = 16)
plt.legend(loc = "best")
plt.show()

