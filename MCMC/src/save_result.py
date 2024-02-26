import os
import numpy as np
import corner
import getdist
import matplotlib.pyplot as plt

# 绘制角图
def plot_corner(observed_data, sample_path, corner_path, ndim, nsteps, nwalkers):
    star_name = observed_data['star_name']
    samples_chain = np.load(f'{sample_path}/sample_{star_name}.npy')
    samples = samples_chain[:, :, :].reshape((-1, ndim))[int(nsteps/2*nwalkers):]
    # samples = np.load(f'./samples/samples_{nwalkers}_{nsteps}.npy')
    # 使用getdist创建MCSamples对象
    names = ['core Mass', 'env mass' 'log_L', 'radius', 'age']
    labels = ['core_mass', 'env_mass', 'log_L', 'radius', 'age']
    # samples_gd = MCSamples(samples=samples, names=names, labels=labels)
    samples_gd = getdist.MCSamples(samples=samples, names=names, labels=labels)
    fig = corner.corner(samples, labels=labels,  sharex=True, sharey=True,show_titles=True) 
    #truths=[observed_data['mass_true'], observed_data['log_L_true'], observed_data['radius_true'],observed_data['age_true']])
    axes = np.array(fig.axes).reshape((ndim, ndim))
    plt.savefig(f'{corner_path}/corner_{star_name}.png', dpi=200)
    # plt.show()

#绘制MC链
def plot_chain(observed_data, sample_path, chian_path, ndim, nsteps, nwalkers):
    star_name = observed_data['star_name']
    samples_chain = np.load(f'{sample_path}/sample_{star_name}.npy')
    # samples_chain = np.load(f'./samples_chain/samples_chain_{nwalkers}_{nsteps}.npy')
    ig, axes = plt.subplots(4, figsize=(10, 7), sharex=True)
    labels = ['core_mass', 'env_mass', 'log_L', 'radius', 'age']
    for i in range(ndim):
        ax = axes[i]
        ax.plot(samples_chain[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples_chain))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)

    axes[-1].set_xlabel("step number")
    os.makedirs(chian_path, exist_ok=True)
    plt.savefig(f'{chian_path}/chain_{star_name}.png', dpi=200)


def save_dat(observed_data, sample_path, data_files, ndim, nsteps, nwalkers):
    # 获取恒星的相关数据
    star_name = observed_data.get('star_name', np.nan)
    samples_chain = np.load(f'{sample_path}/sample_{star_name}.npy')
    samples = samples_chain[:, :, :].reshape((-1, ndim))[int(nsteps/2*nwalkers):]
    # samples = np.load(f'./samples/samples_{nwalkers}_{nsteps}.npy')


    def calculate_parameter_statistics(samples, param_index):
        percentiles = np.percentile(samples[:, param_index], [25, 50, 75])
        return tuple(percentiles)
    
    # 使用字典来存储参数及其索引
    parameters = {'core_mass': 0, 'env_mass':1, 'log_L': 2, 'radius': 3, 'age': 4}
    
    # 计算参数的统计信息
    parameter_stats = {name: calculate_parameter_statistics(samples, index) for name, index in parameters.items()}
    
    # 解包字典，以便更容易地构建写入的字符串
    core_mass_stats, env_mass_stats, log_L_stats, radius_stats, age_stats = [parameter_stats[name] for name in parameters.keys()]
    mass_stats = core_mass_stats + env_mass_stats
    
    mass_lei = observed_data.get('mass_lei', np.nan)
    log_L_lei = observed_data.get('log_L_lei', np.nan)
    radius_lei = observed_data.get('radius_lei', np.nan)
    mass_true = observed_data.get('mass_true', np.nan)
    log_L_true = observed_data.get('log_L_true', np.nan)
    radius_true = observed_data.get('radius_true', np.nan)
    age_true = observed_data.get('age_true', np.nan)
    he_class = observed_data.get('he_class', np.nan)
    log_Teff = observed_data.get('log_Teff', np.nan)
    log_Teff_err = observed_data.get('log_Teff_err', np.nan)
    log_g = observed_data.get('log_g', np.nan)
    log_g_err = observed_data.get('log_g_err', np.nan)
    log_he = observed_data.get('log_he', np.nan)
    log_he_err = observed_data.get('log_he_err', np.nan)

    
    # 将结果以追加模式写入 .dat 文件
    with open(data_files, 'a', newline='') as f:
        # 如果文件为空，写入表头
        if os.path.getsize(data_files) == 0:
            header = "name mass_lei log_L_lei radius_lei mass_true log_L_true radius_true age_true he_class log_he log_he_err log_Teff log_Teff_err log_g log_g_err "
            header += "mass_low mass mass_up log_L_low log_L log_L_up radius_low radius radius_up age_low age age_up "
            header += "mass_max_density log_L_max_density radius_max_density age_max_density\n"  # 添加新的表头
            f.write(header)

        # 构建数据字符串
        data_str = f"{star_name} {mass_lei} {log_L_lei} {radius_lei} {mass_true} {log_L_true} {radius_true} {age_true} {he_class} {log_he} {log_he_err} {log_Teff} {log_Teff_err} {log_g} {log_g_err} "
        data_str += f"{mass_stats[0]} {mass_stats[1]} {mass_stats[2]} {log_L_stats[0]} {log_L_stats[1]} {log_L_stats[2]} {radius_stats[0]} {radius_stats[1]} {radius_stats[2]} {age_stats[0]} {age_stats[1]} {age_stats[2]}\n"  # 添加新的数据
        f.write(data_str)
