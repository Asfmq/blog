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
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
import random
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
    if 0.5 < log_L < 3.5 and 0.3 < star_mass < 0.9 and 0<radius<0.5 and 1e8<age<3e8:
        return 0.0
    return -np.inf

def interp_data(x, y, y0):
    intersection_points = []
    dx = np.diff(x)
    dy = np.diff(y)
    
    change_indices = np.where((dx[:-1] * dx[1:] < 0) | (dy[:-1] * dy[1:] < 0))[0] + 1
    x_subs = np.split(x, change_indices)
    y_subs = np.split(y, change_indices)
    for x_sub, y_sub in zip(x_subs, y_subs):
        if y0 >= np.min(y_sub) and y0 <= np.max(y_sub):
            intersection_point = np.interp(y0, y_sub, x_sub)
            intersection_points.append(intersection_point)
    return intersection_points
def chi2_sol(observed, predicted, observed_err):
    return ((observed - predicted) / observed_err) ** 2

def log_likelihood(p0, observed_log_Teff, observed_log_Teff_err, observed_log_g, observed_log_g_err, observed_log_he, observed_log_he_err):
    star_mass, log_L, radius, age = p0
    chi2_component = []
    components = ['star_mass','log_L', 'radius', 'star_age']
    step = [1e-2,1e-1,1e-2,5e6]
    for i in range(len(files)):
        match = re.match(r'(\d+\.\d+)-(\d+\.\d+)-(\d+\.\d+).data', files[i])
        # Y_surf = None
        chi2_component0 = 0
        if match:
            Y_surf = float(match.group(3))
            # print(Y_surf)
        for component in components:
            x0=p0[components.index(component)]
            predicted_log_teff = interp_data(data[i]['log_Teff'],data[i][component],x0)
            predicted_log_g = interp_data(data[i]['log_g'],data[i][component],x0)
            predicted_log_he = np.log10(Y_surf/4/(1-Y_surf))
            # print('n1=')
            # print(predicted_log_teff)
            # print('n2=')
            # print(predicted_log_g)
            if len(predicted_log_teff) == 0:
                predicted_component = data[i][component][0]
                chi2 = chi2_sol(p0[components.index(component)],predicted_component,step[components.index(component)])
            else:
                chi2_0 = [0]
                for n in range(min(len(predicted_log_teff),len(predicted_log_g))):
                    chi2_1=chi2_sol(observed_log_Teff, predicted_log_teff[n], observed_log_Teff_err) + chi2_sol(observed_log_g, predicted_log_g[n], observed_log_g_err) + chi2_sol(observed_log_he, predicted_log_he, observed_log_he_err)
                    chi2_0.append(chi2_1)
                chi2=min(chi2_0)
            chi2_component0 += chi2
        chi2_component.append(chi2_component0)
    c0 = min(chi2_component)
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
nsteps = 100

def MCMC(observed_data):
    # 初始化步行者的起始位置
    initial_guess = [0.45, 1.4, 0.2, 2e8]  
    p0 = initial_guess + [1e-2,1e-1,1e-2,1e7] * np.random.randn(nwalkers, ndim)
    # 创建MCMC采样器并运行
    pool_if = 1
    if pool_if == 0:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(observed_data['log_Teff'], observed_data['log_Teff_err'], observed_data['log_g'], observed_data['log_g_err'], observed_data['log_he'], observed_data['log_he_err']))
        sampler.run_mcmc(p0, nsteps, progress=True)
    else:
        with Pool() as pool:
            sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(observed_data['log_Teff'], observed_data['log_Teff_err'], observed_data['log_g'], observed_data['log_g_err'], observed_data['log_he'], observed_data['log_he_err']), pool=pool)
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
def load_test():
    path_test_star = "/home/fmq/MESA/work/my/MCMC/code/lei/test_star"
    files = sorted(os.listdir(path_test_star))
    # 指定要使用的列
    usecols = ["star_mass", "surface_he4", "star_age", "log_Teff", "log_g", "log_L", "radius", "center_he4"]
    # 用于存储所有数据的列表
    data_list = []
    # 循环处理每个文件
    for i, file in enumerate(files):
        # 构建文件路径
        file_path = os.path.join(path_test_star, file)
        df = pd.read_csv(file_path, skiprows=5, delim_whitespace=True, usecols=usecols)
        df = df.apply(pd.to_numeric, errors='coerce')
        try:
            start_index = np.where(df['center_he4'] < 0.93)[0][0]
        except IndexError:
            # print(f"Error ZAHB: {files[i]}")
            continue
        df = df[start_index:-1]
        df.drop('center_he4', axis=1, inplace=True)
        df['star_age'] = np.log10(df['star_age'])
        df['surface_he4'] = np.log10(df['surface_he4']/4/(1-df['surface_he4']))
        # 在 start_index 和 df 的长度减一之间生成一个随机整数
        random_index = random.randint(0, len(df) - 1)
        data_list.append({
            'star_name': f'star{i}',
            'log_Teff': df["log_Teff"].iloc[random_index],
            'log_Teff_err': 0.013,
            'log_g': df["log_g"].iloc[random_index],
            'log_g_err': 0.1,
            'log_he': df["surface_he4"].iloc[random_index],
            'log_he_err': 0.1,
            'mass_true':df["star_mass"].iloc[random_index],
            'log_L_true':df["log_L"].iloc[random_index],
            'radius_true':df["radius"].iloc[random_index],
            'age_true':df["star_age"].iloc[random_index]
        })
    # 将所有的 DataFrame 连接成一个单一的 DataFrame
    return data_list
def main():
    # index = next(i for i, data in enumerate(observed_data) if data['star_name'] == 'SDSSJ011506.17+140513.5')
    for star_index, observed_data_star in enumerate(observed_data):
        # if star_index <=index:continue
        print(star_index)
        print(observed_data_star)
        # MCMC(observed_data_star)
        plot_corner(observed_data_star)
        plot_chain(observed_data_star)
        # save_dat(observed_data_star)
        break

# load_method()
# MCMC(observed_data)
# save_fits()
plot_corner(observed_data)
plot_chain(observed_data)


# fig = corner.corner(samples, labels=labels, show_titles=True, quantiles=[0.16, 0.5, 0.84],truths=[star_mass_true, log_L_true])