import matplotlib.pyplot as plt
import numpy as np

data = np.genfromtxt('/home/fmq/MESA/work/my/parsec/grid/Z0.001Y0.25/Z0.001Y0.25OUTA1.74_F7_M0.497.HB.DAT', names=True)
data['L_GRAV'] = data['LOG_L'] * data['L_GRAV']
# print(data['L_GRAV'])
def plot_HR():
    plt.figure(dpi=100, figsize=(4, 3))
    ax = plt.gca()
    # ax.plot(data['LOG_TE'], data['LOG_L'], color='dodgerblue')
    ax.plot(data['MODELL'], data['LOG_L'], color='r')
    ax.plot(data['MODELL'], data['LY'], color='b')
    data['L_GRAV'] = data['LOG_L'] * data['LY']
    ax.plot(data['MODELL'], data['L_GRAV'], color='y')
    ax.set(xlabel=r'$\log Teff_{\rm eff}~{[\rm K]}$',
           ylabel=r'$\log L~{[\rm L_\odot}]$')
    ax.invert_xaxis()
    plt.tight_layout()
    plt.savefig('HR.png',dpi=200)
    plt.show()

plot_HR()
# def plot_HR():
#     plt.figure(dpi=100, figsize=(4, 3))
#     ax = plt.gca()
#     ax.plot(data['LOG_TE'], data['L_GRAV'], color='dodgerblue')
#     ax.set(xlabel=r'$\log Teff_{\rm eff}~{[\rm K]}$',
#            ylabel=r'$\log L~{[\rm L_\odot}]$')
#     # ax.invert_xaxis()
#     plt.tight_layout()
#     plt.savefig('HR.png',dpi=200)
#     plt.show()

# plot_HR()