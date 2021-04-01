import numpy as np
import matplotlib.pyplot as plt

plt.style.use('Solarize_Light2')

""" Box 6.1 """

def lin_func(m, x, b):
    return m*x+b

x_i = np.array([.25, 1.05, 2.25, 2.88, 2.97, 3.64, 3.92, 4.94, 5.92])
y_i = np.array([.86, 2.18, 4.84, 5.80, 6.99, 8.84, 8.71, 11.98, 12.40])
sigma_i = np.array([.27, 1.16, 1.14, .93, .31, .66, .98, .93, .60])

alpha = np.sum(1/sigma_i**2)
beta = np.sum(x_i/sigma_i**2)
gamma = beta
delta = np.sum(x_i**2/sigma_i**2)
theta = np.sum(y_i/sigma_i**2)
phi = np.sum(x_i*y_i/sigma_i**2) 


matrix = np.matrix([[alpha, beta], [gamma, delta]])
matrixinvflip = np.flip(np.linalg.inv(matrix))
print(matrix)
print(matrixinvflip)

print("alpha = {:.4f}".format(alpha))
print("beta = {:.4f}".format(beta))
print("gamma = {:.4f}".format(gamma))
print("delta = {:.4f}".format(delta))
print("theta = {:.4f}".format(theta))
print("phi = {:.4f}".format(phi))

det = (alpha*delta - beta*gamma)
print("determinant = {:.4f}".format(det))

fit_m= 1/det * (alpha*phi - theta*gamma)
fit_b = 1/det * (theta*delta - beta*phi)

print("best fit value of m (theta): {:.4f}".format(fit_m))
print("best fit value of b (phi): {:.4f}".format(fit_b))

sigma_m = np.sqrt(np.diag(matrixinvflip))[0]
sigma_b = np.sqrt(np.diag(matrixinvflip))[1]

print("sigma_m = {:.4f}".format(sigma_m))
print("sigma_b = {:.4f}".format(sigma_b))

sigma_mbsq = matrixinvflip[1,0]

print("sigma_m,b ** 2 = {:.4f}".format(sigma_mbsq))

fx = np.linspace(np.floor(x_i[0]), np.ceil(x_i[-1]), 10**2)
fy = lin_func(fit_m, fx, fit_b)
plt.figure()
plt.errorbar(x_i, y_i, yerr = sigma_i, label = "Measured x and y values \n with Associated Uncertainty", fmt="o")
plt.plot(fx, fy, "--", label = "Linear Fit y = mx + b \n m = {:.2f} b= {:.2f}".format(fit_m, fit_b))
plt.xlabel("x")
plt.ylabel("y")
plt.title("Linear Least Squares Fit for Measured x and y values \n from Table 6.2 in Wong")
plt.legend()