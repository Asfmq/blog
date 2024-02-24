import os
import numpy as np
import pandas as pd
# from isochrones.interp import DFInterpolator
import emcee
import matplotlib.pyplot as plt
from multiprocessing import Pool
import time
from scipy.stats import norm
from astropy.io import fits
from astropy.table import Table
import re
from scipy.interpolate import interp1d
from scipy.stats import gaussian_kde
import random

from src.load_method import load_method
from src.load_observed import load_lei, load_culpan, load_fontaine, load_test
from src.save_result import plot_corner, plot_chain, save_dat

#先验函数
def log_prior(p0):
    star_mass, log_L, radius, age = p0
    if 0.30 < star_mass < 0.85 and 0 < log_L < 3.5 and 0<radius<1.0 and 0<age<3.714*10**10:
        return 0.0
    return -np.inf

# 似然函数

def log_likelihood(p0, method_data, observed_log_Teff, observed_log_Teff_err, observed_log_g, observed_log_g_err, observed_log_he, observed_log_he_err):
    
    def interp1_data(x, y, x0, zero_indices):
        intersection_points = []
        for i in range(len(zero_indices)):
            x_sub = x.iloc[zero_indices[i]:zero_indices[i]+2]
            y_sub = y.iloc[zero_indices[i]:zero_indices[i]+2]
            interp_func = interp1d(y_sub, x_sub)
            intersection_point = interp_func(x0)
            intersection_points.append(intersection_point)
        return np.array(intersection_points)
    def chi2_sol(observed, predicted, observed_err):
        return ((observed - predicted) / observed_err) ** 2
    
    components = ["star_mass", "log_L", "radius", "star_age"]
    component = components[0]
    x0 =p0[components.index(component)]
    interp_indices = np.where(abs(method_data[component]-x0)<0.01)[0]
    method_data0 = method_data.iloc[interp_indices]
    lens = len(interp_indices)
    if lens == 0 :
        chi2_interp = np.inf
    else:
        chi2_log_Teff = chi2_sol(observed_log_Teff, method_data0['log_Teff'].values, observed_log_Teff_err)
        chi2_log_g = chi2_sol(observed_log_g, method_data0['log_g'].values, observed_log_g_err)
        chi2_log_he = chi2_sol(observed_log_he, method_data0['log_he'].values, observed_log_he_err)

        chi2_mass = chi2_sol(p0[0], method_data0['star_mass'].values, 0.01)
        chi2_interp = np.amin(chi2_log_Teff + chi2_log_g + chi2_log_he + chi2_mass)
    chi2_component = chi2_interp
    for component in ["log_L", "radius", "star_age"]:
        x0=p0[components.index(component)]
        zero_indices = np.where((method_data[component].iloc[:-1].values-x0) * (method_data[component].iloc[1:].values-x0) <0)[0]
        change_indices = np.where((method_data["mass"].iloc[:-1].values != method_data["mass"].iloc[1:].values) | (method_data["he"].iloc[:-1].values != method_data["he"].iloc[1:].values))[0]
        interp_indices = np.setdiff1d(zero_indices, change_indices)
        lens = len(interp_indices)
        if lens == 0 :
            chi2_interp = np.inf
        else:
            predicted_log_teff = interp1_data(method_data['log_Teff'],method_data[component], x0, interp_indices)
            predicted_log_g = interp1_data(method_data['log_g'], method_data[component], x0, interp_indices)
            predicted_log_he = method_data['log_he'].iloc[interp_indices].values
            chi2_log_Teff = chi2_sol(observed_log_Teff, predicted_log_teff, observed_log_Teff_err)
            chi2_log_g = chi2_sol(observed_log_g, predicted_log_g, observed_log_g_err)
            chi2_log_he = chi2_sol(observed_log_he, predicted_log_he, observed_log_he_err)
            chi2_interp = np.amin(chi2_log_Teff + chi2_log_g + chi2_log_he)  # 所有插值点的最小卡方值
        chi2_component += chi2_interp
    chi2 = chi2_component
    return -0.5 * chi2

#概率函数
def log_probability(p0, method_data, observed_log_Teff, observed_log_Teff_err, observed_log_g, observed_log_g_err, observed_log_he, observed_log_he_err):
    lp = log_prior(p0)
    if not np.isfinite(lp):
        return -np.inf
    return log_likelihood(p0, method_data, observed_log_Teff, observed_log_Teff_err, observed_log_g, observed_log_g_err, observed_log_he, observed_log_he_err)

def MCMC(observed_data, method_data):
    # 初始化步行者的起始位置
    initial_guess = [0.45, 1.1, 0.10, 1e7]
    step = [0.1, 0.5, 0.20, 1e7]
    # step = [1e-2,1e-1,1e-2,0.1]
    p0 = initial_guess + step * np.random.rand(nwalkers, ndim)
    # p0 =  np.random.randn(nwalkers, ndim) * [np.std(method_data['star_mass']), np.std(method_data['log_L']), np.std(method_data['radius']), np.std(method_data['star_age'])] + [np.mean(method_data['star_mass']), np.mean(method_data['log_L']), np.mean(method_data['radius']), np.mean(method_data['star_age'])]

    # 创建MCMC采样器并运行
    pool_if = 1
    if pool_if == 0:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(method_data, observed_data['log_Teff'], observed_data['log_Teff_err'], observed_data['log_g'], observed_data['log_g_err'], observed_data['log_he'], observed_data['log_he_err']))
        sampler.run_mcmc(p0, nsteps, progress=True)
    else:
        with Pool() as pool:
            sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(method_data, observed_data['log_Teff'], observed_data['log_Teff_err'], observed_data['log_g'], observed_data['log_g_err'], observed_data['log_he'], observed_data['log_he_err']), pool=pool)
            start = time.time()
            sampler.run_mcmc(p0, nsteps, progress=True)
            end = time.time()
            multi_time = end - start
            print("Multiprocessing took {0:.1f} seconds".format(multi_time))
    # 获取MCMC采样结果
    samples_chain = sampler.get_chain()
    star_name = observed_data['star_name']
    np.save(f'{sample_path}/sample_{star_name}', samples_chain)


def main():
    # index = next(i for i, data in enumerate(observed_data) if data['star_name'] == 'star8')
    for star_index, observed_data_star in enumerate(observed_data):
        # if star_index <=index:continue
        print(star_index)
        print(observed_data_star)
        # nerr =10
        # method_data = all_method_data[
        #     (observed_data_star['log_Teff'] - nerr*observed_data_star['log_Teff_err']  < all_method_data['log_Teff']) & 
        #     (all_method_data['log_Teff'] < observed_data_star['log_Teff'] + nerr*observed_data_star['log_Teff_err']) &
        #     (observed_data_star['log_g']-nerr*observed_data_star['log_g_err'] < all_method_data['log_g']) & 
        #     (all_method_data['log_g'] < observed_data_star['log_g'] + nerr*observed_data_star['log_g_err']) &
        #     (observed_data_star['log_he'] - nerr*observed_data_star['log_he_err'] < all_method_data['log_he']) & 
        #     (all_method_data['log_he'] < observed_data_star['log_he'] + nerr*observed_data_star['log_he_err'])
        # ]
        MCMC(observed_data_star, all_method_data)
        plot_corner(observed_data_star, sample_path, corner_path, ndim, nsteps, nwalkers)
        plot_chain(observed_data_star, sample_path, chian_path, ndim, nsteps, nwalkers)
        save_dat(observed_data_star, sample_path, data_files, ndim, nsteps, nwalkers)
        # break


sample_path = './sample/sample_smooth_lei'
corner_path = './corner/corner_smooth_lei'
chian_path = './chain/chain_smooth_lei'
os.makedirs(sample_path, exist_ok=True)
os.makedirs(corner_path, exist_ok=True)
os.makedirs(chian_path, exist_ok=True)
data_files = 'lei_smooth.dat'

# 导入模型
# path_methods = ["/home/zxlei/pfiles/fmq/sdb/data_hb", "/home/zxlei/pfiles/fmq/sdb/data_wd"]
# method_data = load_method(path_methods, 'all_data.csv')
all_method_data = pd.read_csv('/home/zxlei/pfiles/fmq/mcmc/all_data.csv')

# observed_data = load_test()
# observed_data = pd.read_csv('/home/zxlei/pfiles/fmq/mcmc/test_star.csv').to_dict('records')
# observed_data = load_test(method_data)
observed_data = load_lei()

nwalkers = 128
ndim = 4
nsteps = 100


main()
