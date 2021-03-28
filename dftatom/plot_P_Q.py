import numpy as np
import matplotlib.pyplot as plt
import matplotlib

rmesh = np.loadtxt("fort.1000")
Nr = rmesh.shape[0]

P = np.loadtxt("fort.1003")
Q = np.loadtxt("fort.1004")

plt.clf()
plt.plot(rmesh, P, label="P")
plt.plot(rmesh, Q, label="Q")
plt.xlim(0.0, 10.0)
plt.grid(True)
plt.legend()
plt.savefig("IMG_P_Q.pdf")
