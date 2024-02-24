import os
import numpy as np
import emcee
import getdist
from getdist import plots, MCSamples
import matplotlib.pyplot as plt
import corner
from multiprocessing import Pool
import time
from scipy.stats import norm
import pandas as pd
from astropy.io import fits
from astropy.table import Table
import re

# 模型导入
def load_method():
    path_method = "/home/fmq/MESA/work/my/parsec/grid/Z0.001Y0.25"
    # path_method = "/home/fmq/MESA/work/my/MCMC/hot_sd/temp"
    global files
    global data
    global all_data
    files = sorted(os.listdir(path_method))
    data = {}
    all_data = []
    for i in range(len(files)):
        # if i > 1:continue
        data[i] = np.genfromtxt(path_method + '/' + files[i], names=True)
        if_se = 0
        if if_se != 0:
            try:
                    start_index = np.where(data[i]['center_he4'] < 0.93)[0][0]
            except IndexError:
                print(f"Error ZAHB: {files[i]}")
                continue
            try:
                    end_index=np.where(data[i]['center_he4']<0.001)[0][0]
            except IndexError:
                print(f"Error TAHB: {files[i]}")
                end_index=-1
            index_range=slice(start_index,end_index)
            all_data.append(data[i][index_range])
            data[i] = data[i][index_range]
    # return files, data, all_data

def plot():
    load_method()
    plt.figure(dpi=100, figsize=(4, 3))
    ax = plt.gca()
    for i in range(len(files)):
        ax.plot(data[i]['LOG_TE'], data[i]['LOG_L'], color='dodgerblue')
    ax.set(xlabel=r'$\log Teff_{\rm eff}~{[\rm K]}$',
           ylabel=r'$\log L~{[\rm L_\odot}]$')
    ax.invert_xaxis()
    plt.tight_layout()
    plt.savefig('HR_all.png',dpi=200)
    plt.show()

# 观测数据
def convert_to_logfloat(x):
    try:
        return np.log10(float(x))
    except ValueError:
        return -10**3
def convert_to_float(x):
    try:
        return float(x)
    except ValueError:
        return -10**3
def load_lei():

    data_lei = fits.open('/home/fmq/MESA/work/my/observe_data/lei2023/sd_mass_0p20.fits')
    num_rows = data_lei[1].data.shape[0]

    star_name=data_lei[1].data['spec_name']
    
    he_class = data_lei[1].data['sp_class']
    loghe=data_lei[1].data['loghe']
    # 将 loghe 数组中的每一个元素转换为字符串，然后删除 '>'
    loghe = np.array([str(x).replace('>', '') for x in loghe])
    # 将 loghe 数组中的每一个元素转换为浮点数
    log_he = np.array([float(x) for x in loghe])
    log_he_err = data_lei[1].data['loghe_err']
    # Y_surf = 0.98*10**log_he/(1/4+10**loghe)
    log_Teff = np.log10(data_lei[1].data['teff'])
    log_Teff_err = np.log10(1+data_lei[1].data['teff_err']/log_Teff)
    log_g=data_lei[1].data['logg']
    log_g_err=data_lei[1].data['logg']

    mass = data_lei[1].data['mass_median']
    log_L = np.log10(data_lei[1].data['l_div_lsun_median'])
    radius = data_lei[1].data['radius_median']

    data_lei_list = []
    for i in range(num_rows):
        data_lei_list.append({
            'star_name': star_name[i],
            'he_class': he_class[i],
            'log_Teff': log_Teff[i],
            'log_Teff_err': log_Teff_err[i],
            'log_g': log_g[i],
            'log_g_err': log_g_err[i],
            'log_he': log_he[i],
            'log_he_err': log_he_err[i],
            'mass_lei': mass[i],
            'log_L_lei':log_L[i],
            'radius_lei':radius[i]
        })
    return data_lei_list

def loda_huo():
    data_huo = pd.read_csv('/home/fmq/MESA/work/my/observe_data/huo2023/huo.csv')
    print(data_huo)
    star_name = data_huo['Destigation']
    log_Teff = np.log10(data_huo['Teff']*1000)
    log_g = data_huo['log g']
    log_he = data_huo['[M/H]']
    num_rows = len(data_huo)
    data_huo_list = []
    for i in range(num_rows):
        data_huo_list.append({
            'star_name': star_name[i],
            'log_Teff': log_Teff[i],
            'log_Teff_err': 0.26,
            'log_g': log_g[i],
            'log_g_err': 0.25,
            'log_he': log_he[i],
            'log_he_err': log_he_err[i],
            'mass_lei': mass[i],
            'log_L_lei':log_L[i],
            'radius_lei':radius[i]
        })
    return data_lei_list

loda_huo()
#先验函数
def log_prior(p0):
    star_mass, log_L, radius, age = p0
    if 0.5 < log_L < 3.5 and 0.3 < star_mass < 0.9 and 0<radius<0.5 and 1e8<age<3e8:
        return 0.0
    return -np.inf

# 似然函数
def chi2_component(observed, predicted, observed_err):
    return ((observed - predicted) / observed_err) ** 2

def log_likelihood(p0, observed_log_Teff, observed_log_Teff_err, observed_log_g, observed_log_g_err, observed_log_he, observed_log_he_err):
    star_mass, log_L, radius, age = p0
    chi2 = []
    step = [1e-2,1e-1,1e-2,2e7]
    components = ['star_mass','log_L',  'radius', 'star_age']
    for i in range(len(files)):
        match = re.match(r'(\d+\.\d+)-(\d+\.\d+)-(\d+\.\d+).data', files[i])
        if match:
            Y_surf = float(match.group(3))
        chi2_0 = 0
        for component in components:
            # if p0[components.index(component)]>max(data[i][component]) or p0[components.index(component)]<min(data[i][component]):
            #     continue
            nearest_index = np.abs(data[i][component] - p0[components.index(component)]).argmin()
            predicted_log_teff = data[i]['log_Teff'][nearest_index]
            predicted_log_g = data[i]['log_g'][nearest_index]
            predicted_component = data[i][component][nearest_index]
            predicted_log_he = np.log10(Y_surf/4/(1-Y_surf))
            
            chi2_0 += chi2_component(observed_log_Teff, predicted_log_teff, observed_log_Teff_err) + chi2_component(observed_log_g, predicted_log_g, observed_log_g_err) + chi2_component(observed_log_he, predicted_log_he, observed_log_he_err) + chi2_component(p0[components.index(component)],predicted_component,step[components.index(component)])
        chi2.append(chi2_0)
        # print(chi2)
    c0 = min(chi2)
    # print(c0)
    return -0.5 * c0

#概率函数
def log_probability(p0, observed_log_Teff, observed_log_Teff_err, observed_log_g, observed_log_g_err, observed_log_he, observed_log_he_err):
    lp = log_prior(p0)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(p0, observed_log_Teff, observed_log_Teff_err, observed_log_g, observed_log_g_err, observed_log_he, observed_log_he_err)

#高斯拟合
def gaussian_fit(samples, ax, label, color, truth):
    mu, std = norm.fit(samples)
    x = np.linspace(np.min(samples), np.max(samples), 100)
    p = norm.pdf(x, mu, std)
    
    print(samples)# # 计算直方图的值
    hist_values, bin_edges = np.histogram(samples, bins=30, density=True)
    max_hist_value = max(hist_values)
    print(max_hist_value)
    
    # 缩放 PDF 的高度
    p_scaled = p * max_hist_value
    # ax.plot(x, p_scaled, color=color, linewidth=2)
    ax.hist(samples, bins=30, density=True, alpha=0.5, color=color)
    ax.text(0.95, 0.95, f'$True={truth:.3f}$\n$\mu={mu:.3f}$\n$\sigma={std:.3f}$', transform=ax.transAxes, ha='right', va='top')
    

# MCMC参数
nwalkers = 32
ndim = 4
nsteps = 1000

def MCMC(observed_data):
    # 初始化步行者的起始位置
    initial_guess = [0.45, 1.4, 0.2, 2e8]
    step = [1e-2,1e-1,1e-2,2e7]
    p0 = initial_guess + step * np.random.randn(nwalkers, ndim)
    # 创建MCMC采样器并运行
    pool_if = 0
    if pool_if == 0:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(observed_data['log_Teff'], observed_data['log_Teff_err'], observed_data['log_g'], observed_data['log_g_err'], observed_data['log_he'], observed_data['log_he_err']))
        sampler.run_mcmc(p0, nsteps, progress=True)
    else:
        with Pool() as pool:
            sampler = emcee.EnsemR. AndraebleSampler(nwalkers, ndim, log_probability, args=(all_data[0]['log_Teff'], all_data[0]['log_L'], observed_data['log_Teff'], observed_data['log_L']), pool=pool)
            start = time.time()
            sampler.run_mcmc(p0, nsteps, progress=True)
            end = time.time()
            multi_time = end - start
            print("Multiprocessing took {0:.1f} seconds".format(multi_time))
    # 获取MCMC采样结果
    samples = sampler.chain[:, :, :].reshape((-1, ndim))
    samples_chain = sampler.get_chain()
    # np.save('samples2.npy', samples)
    # np.save('samples2_chain.npy',samples_chain)
    star_name = observed_data['star_name']
    np.save(f'./samples/samples_{star_name}.npy', samples)
    np.save(f'./samples_chain/samples_chain_{star_name}.npy', samples_chain)



# 绘制角图
def plot_corner(observed_data):
    star_name = observed_data['star_name']
    samples = np.load(f'./samples/samples_{star_name}.npy')
    # 使用getdist创建MCSamples对象
    names = ['Star Mass', 'log_L', 'radius', 'age']
    labels = ['M', 'log_L', 'radius', 'age']
    # samples_gd = MCSamples(samples=samples, names=names, labels=labels)
    samples_gd = getdist.MCSamples(samples=samples, names=names, labels=labels)
    fig = corner.corner(samples, labels=labels,  sharex=True, sharey=True,show_titles=True)
    axes = np.array(fig.axes).reshape((ndim, ndim))
    # 获取 'Star Mass' 参数的最可能值、下限和上限
   
    # param_info = samples_gd.getParams()
    # print('info:')
    # print(param_info)
    # 对样本进行高斯拟合并绘制在角图上
    # gaussian_fit(samples[:, 0].flatten(), axes[0][0], '星体质量', 'blue', truth=star_mass_true)
    # gaussian_fit(samples[:, 1].flatten(), axes[1][1], 'log_L', 'green', truth=log_L_true)

    # 添加图例
    fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.9))
    # 添加图例，设置文字大小为 10
    plt.legend(fontsize=1)
    # 保存图像
    plt.savefig(f'./corner/corner_{star_name}.png', dpi=200)
    plt.show()

#绘制MC链
def plot_chain(observed_data):
    star_name = observed_data['star_name']
    samples_chain = np.load(f'./samples_chain/samples_chain_{star_name}.npy')
    ig, axes = plt.subplots(4, figsize=(10, 7), sharex=True)
    labels = ['M', 'log_L', 'radius', 'age']
    for i in range(ndim):
        ax = axes[i]
        ax.plot(samples_chain[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples_chain))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)

    axes[-1].set_xlabel("step number")
    plt.savefig(f'./chain/chain_{star_name}.png',dpi=200)

    plt.show()

def save_fits(observed_data):
    star_name = observed_data['star_name']
    samples = np.load(f'./samples/samples_{star_name}.npy')
    # 获取恒星质量和光度的估计值及误差
    estimated_star_mass = np.median(samples[:, 0])
    estimated_log_L = np.median(samples[:, 1])

    star_mass_error = np.std(samples[:, 0])
    log_L_error = np.std(samples[:, 1])

    # 将标量值和误差转换为一维数组
    estimated_star_mass = np.array([estimated_star_mass])
    estimated_log_L = np.array([estimated_log_L])
    star_mass_error = np.array([star_mass_error])
    log_L_error = np.array([log_L_error])

    # 创建Astropy表格
    data = Table({
        'NAME': [star_name],
        'STAR_MASS': estimated_star_mass,
        'STAR_MASS_ERROR': star_mass_error,
        'LOG_L': estimated_log_L,
        'LOG_L_ERROR': log_L_error
    })

    # 创建FITS文件
    hdu_list = fits.HDUList([
        fits.PrimaryHDU(),
        fits.BinTableHDU(data, name='DATA')
    ])

    # 保存到文件
    hdu_list.writeto('feng.fits', overwrite=True)


def calculate_parameter_statistics(samples, param_index):
    best_fit = np.median(samples[:, param_index])
    lower_limit = np.percentile(samples[:, param_index], 16)
    upper_limit = np.percentile(samples[:, param_index], 84)
    return best_fit, lower_limit, upper_limit

def save_dat(observed_data):
    # 获取恒星的相关数据
    star_name = observed_data['star_name']
    samples = np.load(f'./samples/samples_{star_name}.npy')
    # 获取恒星质量和光度的估计值及误差
    
    names = ['Star Mass', 'log_L', 'radius', 'age']
    star_mass_index = names.index('Star Mass')
    mass_median, mass_low, mass_up = calculate_parameter_statistics(samples, 0)
    log_L_median, log_L_low, log_L_up = calculate_parameter_statistics(samples, 1)
    radius_median, radius_low, radius_up = calculate_parameter_statistics(samples, 2)
    age_median, age_low, age_up = calculate_parameter_statistics(samples, 3)

    mass_lei = observed_data['mass_lei']
    log_L_lei = observed_data['log_L_lei']
    radius_lei = observed_data['radius_lei']
    he_class = observed_data['he_class']
    log_Teff = observed_data['log_Teff']
    log_Teff_err = observed_data['log_Teff_err']
    log_g = observed_data['log_g']
    log_g_err = observed_data['log_g_err']
    log_he = observed_data['log_he']
    log_he_err = observed_data['log_he_err']
    # 将结果以追加模式写入 .dat 文件
    with open('feng.dat', 'a', newline='') as f:
        # 如果文件为空，写入表头
        if os.path.getsize('feng.dat') == 0:
            f.write("NAME MASS_lei log_L_lei radius_lei he_class log_he log_he_err log_Teff log_Teff_err log_g log_g_err MASS_LOW MASS_UP LOG_L LOG_L_LOW LOG_L_UP RADIUS RADIUS_LOW RADIUS_UP AGE AGE_LOW AGE_UP\n")
        # 写入每颗恒星的结果到文件
        f.write(f"{star_name} {mass_lei} {log_L_lei} {radius_lei} {he_class} {log_he} {log_he_err} {log_Teff} {log_Teff_err} {log_g} {log_g_err} {mass_median} {mass_low} {mass_up} {log_L_median} {log_L_low} {log_L_up} {radius_median} {radius_low} {radius_up} {age_median} {age_low} {age_up}\n")


observed_data = {
    'log_Teff': 4.512,
    'log_Teff_err': 0.013,
    'log_g': 5.758,
    'log_g_err': 0.1,
    'log_he': 0,
    'log_he_err': 0.1
}


# 导入模型
load_method()
observed_data =load_lei()

# 准确值
star_mass_true = 0.4982
log_L_true = 1.381
radius_true = 0.154
age_true = 2.55e8

for star_index, observed_data_star in enumerate(observed_data):
    # MCMC(observed_data_star)
    # print(star_index)
    # print(observed_data_star)
    # plot_corner(observed_data_star)
    # plot_chain(observed_data_star)
    # save_fits(observed_data_star)
    # save_dat(observed_data_star)


# MCMC(observed_data)
# # save_fits()
# plot_corner(observed_data)
# plot_chain()


# fig = corner.corner(samples, labels=labels, show_titles=True, quantiles=[0.16, 0.5, 0.84],truths=[star_mass_true, log_L_true])