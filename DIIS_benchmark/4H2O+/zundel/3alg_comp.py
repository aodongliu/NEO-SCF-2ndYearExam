import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from matplotlib import ticker
import h5py

# Basis: hatz/pb4d
# Load Data
zundel_stw_e = np.loadtxt("/Users/aodongliu/LiGroup/gaussian/DIIS_benchmark/4H2O+/zundel/zundel_stw_e.dat", dtype=None)
zundel_stw_pe = np.loadtxt("/Users/aodongliu/LiGroup/gaussian/DIIS_benchmark/4H2O+/zundel/zundel_stw_pe.dat", dtype=None)
zundel_stw_pp = np.loadtxt("/Users/aodongliu/LiGroup/gaussian/DIIS_benchmark/4H2O+/zundel/zundel_stw_pp.dat", dtype=None)

zundel_smtnd_e = np.loadtxt("/Users/aodongliu/LiGroup/gaussian/DIIS_benchmark/4H2O+/zundel/zundel_smtnd_e.dat", dtype=None)
zundel_smtnd_pe = np.loadtxt("/Users/aodongliu/LiGroup/gaussian/DIIS_benchmark/4H2O+/zundel/zundel_smtnd_pe.dat", dtype=None)
zundel_smtnd_pp = np.loadtxt("/Users/aodongliu/LiGroup/gaussian/DIIS_benchmark/4H2O+/zundel/zundel_smtnd_pp.dat", dtype=None)

zundel_smtd_e = np.loadtxt("/Users/aodongliu/LiGroup/gaussian/DIIS_benchmark/4H2O+/zundel/zundel_smtd_e.dat", dtype=None)
zundel_smtd_pe = np.loadtxt("/Users/aodongliu/LiGroup/gaussian/DIIS_benchmark/4H2O+/zundel/zundel_smtd_pe.dat", dtype=None)
zundel_smtd_pp = np.loadtxt("/Users/aodongliu/LiGroup/gaussian/DIIS_benchmark/4H2O+/zundel/zundel_smtd_pp.dat", dtype=None)

# Create x axis for plotting
zundel_stw_iter = np.arange(1,zundel_stw_e.size,1)
zundel_stw_iter2 = np.arange(1,zundel_stw_pe.size,1)
zundel_smtnd_iter = np.arange(1,zundel_smtnd_e.size,1)
zundel_smtd_iter = np.arange(1,zundel_smtd_e.size,1)

# Subtract the final energy for energy at each iteration
zundel_stw_e -= zundel_stw_e[-1]
zundel_smtnd_e -= zundel_smtnd_e[-1]
zundel_smtd_e -= zundel_smtd_e[-1]

plt.style.use('default')
fig = plt.figure(figsize=(16,9))
gs = fig.add_gridspec(3, hspace=0)
ax = gs.subplots(sharex=True)
fig.suptitle('NEO-SCF with Different Algorithms, hatz/pb4d, $4H_2O+$ zundel Isomer')
ax[0].plot(zundel_smtnd_iter,zundel_smtnd_e[1:],'ob-',ms=2,label="Simul NODIIS")
ax[0].plot(zundel_smtd_iter,zundel_smtd_e[1:],'or-',ms=2,label="Simul DIIS")
ax[0].plot(zundel_stw_iter,zundel_stw_e[1:],'og-',ms=2,label="STW DIIS")
ax[0].axhline(y=10**(-8),xmin=0,xmax=zundel_stw_iter.size+1,c="black",linewidth=2,zorder=0,label='Threshold')
ax[0].set_yscale("log")
ax[0].set_ylabel("$E$ Diff against $E_{final}$")
ax[0].legend(loc='upper right')


ax[1].plot(zundel_smtnd_iter,zundel_smtnd_pe[1:],'ob',ms=2,label="Simul NODIIS")
ax[1].plot(zundel_smtd_iter,zundel_smtd_pe[1:],'or',ms=2,label="Simul DIIS")
ax[1].plot(zundel_stw_iter2,zundel_stw_pe[1:],'og',ms=2,label="STW DIIS")
ax[1].axhline(y=10**(-8),xmin=0,xmax=zundel_stw_iter.size+1,c="black",linewidth=2,zorder=0,label='Threshold')
ax[1].set_yscale("log")
ax[1].set_ylabel("$RMSD P_E$")
ax[1].legend(loc='upper right')


ax[2].plot(zundel_smtnd_iter,zundel_smtnd_pp[1:],'ob',ms=2,label="Simul NODIIS")
ax[2].plot(zundel_smtd_iter,zundel_smtd_pp[1:],'or',ms=2,label="Simul DIIS")
ax[2].plot(zundel_stw_iter2,zundel_stw_pp[1:],'og',ms=2,label="STW DIIS")
ax[2].axhline(y=10**(-8),xmin=0,xmax=zundel_stw_iter.size+1,c="black",linewidth=2,zorder=0,label='Threshold')
ax[2].set_yscale("log")
ax[2].set_ylabel("$RMSD P_P$")
ax[2].legend(loc='upper right')
ax[2].set_xlabel("# of Iterations")
ax[2].set_xticks(np.arange(0, zundel_stw_iter.size+1, step=200))


# Hide x labels and tick labels for all but bottom plot.
for axes in ax:
    axes.label_outer()


fig.set_facecolor('w')

plt.tight_layout()
plt.savefig('/Users/aodongliu/LiGroup/gaussian/DIIS_benchmark/4H2O+/plots/zundel_hatz_pb4d_3AlgComp.png', dpi=300)
plt.show()

