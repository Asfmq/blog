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
    # path_method = "/home/fmq/MESA/work/my/parsec/grid/Z0.001Y0.25"
    path_method = "/home/fmq/MESA/work/my/MCMC/hot_sd/temp"
    global files
    global data
    global all_data
    files = sorted(os.listdir(path_method))
    data = {}
    for i in range(len(files)):
        # if i > 1:continue
        data[i] = np.genfromtxt(path_method + '/' + files[i], names=True, skip_header=5)
        if_se = 1
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
            data[i] = data[i][index_range]
    # return files, data, all_data



#先验函数
def log_prior(p0):
    star_mass, log_L, radius, age = p0
    if 0.5 < log_L < 3.5 and 0.3 < star_mass < 0.9 and 0<radius<0.5 and 5e7<age<15e8:
        return 0.0
    return -np.inf

def chi2_component(observed, predicted, observed_err):
    return ((observed - predicted) / observed_err) ** 2

def log_likelihood(p0, observed_log_Teff, observed_log_Teff_err, observed_log_g, observed_log_g_err, observed_log_he, observed_log_he_err):
    star_mass, log_L, radius, age = p0
    chi2 = []
    step = [1e-2,1e-1,1e-2,5e6]
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



# MCMC参数
nwalkers = 32
ndim = 4 
nsteps = 1000

def MCMC(observed_data):
    # 初始化步行者的起始位置
    initial_guess = [0.45, 1.4, 0.2, 2e8]
    step = [1e-2,1e-1,1e-2,1e7]
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
    np.save('samples2.npy', samples)
    np.save('samples2_chain.npy',samples_chain)
    # star_name = observed_data['star_name']
    # np.save(f'./samples/samples_{star_name}.npy', samples)
    # np.save(f'./samples_chain/samples_chain_{star_name}.npy', samples_chain)


# 绘制角图
def plot_corner(observed_data):
    # star_name = observed_data['star_name']
    samples = np.load('samples2.npy')
    samples = samples[int(50*nwalkers):]
    # 使用getdist创建MCSamples对象
    names = ['Star Mass', 'log_L', 'radius', 'age']
    labels = ['M', 'log_L', 'radius', 'age']
    samples_gd = MCSamples(samples=samples, names=names, labels=labels)
    # fig = corner.corner(samples, labels=labels, truths=[star_mass_true, log_L_true])
    fig = corner.corner(samples, labels=labels,  sharex=True, sharey=True,show_titles=True,truths=[star_mass_true, log_L_true,radius_true,age_true])
    axes = np.array(fig.axes).reshape((ndim, ndim))

    # 添加图例
    fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.9))
    # 添加图例，设置文字大小为 10
    plt.legend(fontsize=1)
    # 保存图像
    # plt.savefig(f'./corner/corner_{star_name}.png', dpi=200)
    plt.savefig('t.png', dpi=200)
    plt.show()

#绘制MC链
def plot_chain(observed_data):
    # star_name = observed_data['star_name']
    samples_chain = np.load(f'samples2_chain.npy')
    ig, axes = plt.subplots(4, figsize=(10, 7), sharex=True)
    labels = ['M', 'log_L', 'radius', 'age']
    for i in range(ndim):
        ax = axes[i]
        ax.plot(samples_chain[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples_chain))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)


    axes[-1].set_xlabel("step number")
    # plt.savefig(f'./chain/chain_{star_name}.png',dpi=200)

    plt.show()

observed_data = {
    'log_Teff': 4.512,
    'log_Teff_err': 0.013,
    'log_g': 5.758,
    'log_g_err': 0.1,
    'log_he': 0,
    'log_he_err': 0.1
}


# 准确值
star_mass_true = 0.4982
log_L_true = 1.381
radius_true = 0.154
age_true = 2.55e8

# load_method()
# MCMC(observed_data)
# # save_fits()
plot_corner(observed_data)
plot_chain(observed_data)


# fig = corner.corner(samples, labels=labels, show_titles=True, quantiles=[0.16, 0.5, 0.84],truths=[star_mass_true, log_L_true])