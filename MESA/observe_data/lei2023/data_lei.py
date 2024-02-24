import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits


data_lei = fits.open('/home/fmq/MESA/work/my/observe_data/sd_mass/sd_mass_0p20.fits')
num_rows = data_lei[1].data.shape[0]
print(f"The number of rows is {num_rows}")
Teff=data_lei[1].data['teff']
logg=data_lei[1].data['logg']
loghe=data_lei[1].data['loghe']
he_class = data_lei[1].data['sp_class']
mass = data_lei[1].data['mass_median']
L=data_lei[1].data['l_div_lsun_median']
# 将 loghe 数组中的每一个元素转换为字符串，然后删除 '>'
loghe = np.array([str(x).replace('>', '') for x in loghe])

# 将 loghe 数组中的每一个元素转换为浮点数
loghe = np.array([float(x) for x in loghe])
Y_surf = 0.98*10**loghe/(1/4+10**loghe)


log_Teff = np.log10(Teff)
log_L = np.log10(L)

plt.figure(dpi=100, figsize=(6, 4))
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
# for i in range(num_rows):
#     if he_class == 'He-sdB' or he_class == 'He-sdO' or he_class == 'He-sdOB':
#         ax.scatter(loghe[i])
# 绘制散点图
for cls, (color, marker) in class_to_marker.items():
    mask = he_class == cls
    ax.scatter(log_Teff[mask], Y_surf[mask], color=color, marker=marker, facecolors='none', label=cls)
# for cls, (color, marker) in class_to_marker.items():
#     if cls in ['He-sdB', 'He-sdO', 'He-sdOB']:
#         mask = he_class == cls
#         ax.scatter(log_Teff[mask], loghe[mask], color=color, marker=marker, facecolors='none', label=cls)
# 添加图例
plt.axhline(y=0.28, color='black', linestyle='dashed')
ax.legend()
ax.invert_xaxis()
# ax.set(xlabel=r'$\log T_{\rm eff}/K$', ylabel=r'$\log L$ /$L_\odot$')
ax.set(xlabel=r'$\log T_{\rm eff}/K$', ylabel=r'$Y_{surf}$')
# start, end = ax.get_ylim()
# ax.yaxis.set_ticks(np.arange(start, end, step=50))
plt.savefig('lei2.png', dpi=200)
plt.show()

