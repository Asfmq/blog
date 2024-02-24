import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits

data_culpan = pd.read_csv('/home/fmq/MESA/work/my/observe_data/knownhsd.dat',header=None, sep='\t')

Teff = data_culpan[0].apply(lambda x: x[305:311])
logg = data_culpan[0].apply(lambda x: x[318:322])


def convert_Teff_to_float(x):
    try:
        return np.log10(float(x))
    except ValueError:
        return -10**3
def convert_logg_to_float(x):
    try:
        return float(x)
    except ValueError:
        return -10**3
log_Teff = Teff.apply(convert_Teff_to_float)
log_g = logg.apply(convert_logg_to_float)

# print(log_g)
# print(log_Teff)

num_rows = data_culpan.shape[0]
print(num_rows)
n=0
plt.figure(dpi=100, figsize=(6, 4))
ax = plt.gca()
for i in range(num_rows):
    if log_Teff[i] != -10**3 and log_g[i] != -10**3:
        ax.scatter(log_Teff[i], log_g[i],marker="x",color='gray')
    else:
        n=n+1
    
# print(n)



data_lei = fits.open('/home/fmq/MESA/work/my/observe_data/sd_mass/sd_mass_0p20.fits')
num_rows = data_lei[1].data.shape[0]
print(f"The number of rows is {num_rows}")
Teff=data_lei[1].data['teff']
logg=data_lei[1].data['logg']
loghe=data_lei[1].data['loghe']
he_class = data_lei[1].data['sp_class']
mass = data_lei[1].data['mass_median']
L=data_lei[1].data['l_div_lsun_median']

log_Teff = np.log10(Teff)
log_L = np.log10(L)
ax = plt.gca()
# 创建一个字典来映射 he_class 的值到颜色和标记
class_to_marker = {
    'sdB': ('k', 'o'),  # 黑色圆圈
    'sdO': ('g', 's'),  # 绿色正方形
    'sdOB': ('m', '^'),  # 紫色正三角形
    'He-sdB': ('pink', '*'),  # 粉红色五角星
    'He-sdO': ('b', '<'),  # 蓝色左三角
    'He-sdOB': ('r', 'D')  # 红色菱形
}

# 绘制散点图
for cls, (color, marker) in class_to_marker.items():
    mask = he_class == cls
    ax.scatter(log_Teff[mask], logg[mask], color=color, marker=marker, facecolors='none', label=cls)

# 添加图例
ax.legend()
ax.invert_xaxis()
ax.set(xlabel=r'$\log T_{\rm eff}/K$', ylabel=r'$\log g$ /cm s$^{-2}$')
plt.savefig('culpan&lei.png', dpi=200)
plt.show()