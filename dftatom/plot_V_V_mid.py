import numpy as np
import matplotlib.pyplot as plt

rmesh = np.loadtxt("fort.1000")
Nr = rmesh.shape[0]
V = np.loadtxt("fort.1001")

rmid = 0.5*(rmesh[1:Nr] + rmesh[0:Nr-1])
Vmid = np.loadtxt("fort.1002")

idx_plot = range(500,510)
plt.clf()
plt.plot(rmesh[idx_plot], V[idx_plot], label="V", marker="o", alpha=0.5)
plt.plot(rmid[idx_plot], Vmid[idx_plot], label="Vmid", marker="o", linewidth=0)
plt.legend()
plt.savefig("IMG_V_Vmid.png", dpi=150)
