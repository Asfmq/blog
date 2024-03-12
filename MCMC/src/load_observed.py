
import pandas as pd
import numpy as np
from astropy.io import fits

def load_lei():
    data_lei = fits.open('/home/fmq/MESA/work/my/observe_data/lei2023/sd_mass_0p20.fits')
    df = pd.DataFrame(data_lei[1].data)

    df['log_he'] = df['loghe'].replace('>', '', regex=True).astype(float)
    df['log_Teff'] = np.log10(df['teff'])
    df['log_Teff_err'] = np.log10(1 + df['teff_err'] / df['teff'])
    df['log_L'] = np.log10(df['l_div_lsun_median'])

    selected_columns = ['spec_name', 'sp_class', 'log_Teff', 'log_Teff_err', 'logg', 'logg_err', 'log_he', 'loghe_err', 'mass_median', 'log_L', 'radius_median']
    df = df[selected_columns]

    df.columns = ['star_name', 'he_class', 'log_Teff', 'log_Teff_err', 'log_g', 'log_g_err', 'log_he', 'log_he_err', 'mass_lei', 'log_L_lei', 'radius_lei']

    # 过滤数据
    df = df[(df['he_class'] != 'sdO') & (df['he_class'] != 'He-sdO')]

    return df.to_dict(orient='records')

def load_culpan():
    # 读取数据
    data_culpan = pd.read_csv('~/pfiles/fmq/observe_data/culpan2022/knownhsd.dat', header=None, sep='\t')

    # 将字符串转换为浮点数的辅助函数
    def convert_to_float(x):
        try:
            return float(x)
        except ValueError:
            return -10**3

    # 处理每一行的辅助函数
    def process_row(row):
        star_name = row[0][0:27].replace(" ", "").replace('[?/]', '')
        Teff = convert_to_float(row[0][305:311])
        Teff_err = convert_to_float(row[0][312:317])
        logg = convert_to_float(row[0][318:322])
        logg_err = convert_to_float(row[0][323:327])
        loghe = convert_to_float(row[0][328:333])
        loghe_err = convert_to_float(row[0][334:339])

        # 如果有无效数据，则返回 None
        if Teff == -10**3 or logg == -10**3 or loghe == -10**3 or loghe_err == -10**3 or logg_err == -10**3:
            return None

        # 处理缺失的错误值
        if Teff_err == -10**3 or Teff_err == 0:
            Teff_err = 0.02
        if logg_err == -10**3 or logg_err == 0:
            logg_err = 0.1
        if loghe_err == -10**3 or loghe_err == 0:
            loghe_err = 0.1

        # 返回处理后的数据
        return {
            'star_name': star_name,
            'log_Teff': np.log10(Teff),
            'log_Teff_err': np.log10(1 + Teff_err / Teff),
            'log_g': logg,
            'log_g_err': logg_err,
            'log_he': loghe,
            'log_he_err': loghe_err
        }

    # 使用 apply 函数处理每一行，并丢弃包含无效数据的行，最后转换为列表
    data_culpan_list = data_culpan.apply(process_row, axis=1).dropna().tolist()

    return data_culpan_list

def load_fontaine():
    data_fontaine = pd.read_csv('/home/zxlei/pfiles/fmq/observe_data/fontaine2012/asteroseismology_26.txt', header=0, sep=',')
    star_name = data_fontaine['rawid']
    log_Teff = np.log10(data_fontaine['Teff'])
    log_Teff_err = np.log10(1+data_fontaine['Teff-err']/data_fontaine['Teff'])
    log_g = data_fontaine['logg']
    log_g_err = data_fontaine['logg-err']
    log_he = data_fontaine['loghe']
    log_he_err = data_fontaine['loghe-err']
    mass_true = data_fontaine['mass']

    new_data = pd.DataFrame({
        'star_name': star_name,
        'log_Teff': log_Teff,
        'log_Teff_err': log_Teff_err,
        'log_g': log_g,
        'log_g_err': log_g_err,
        'log_he': log_he,
        'log_he_err': log_he_err,
        'mass_true': mass_true
    })

    data_dict = new_data.to_dict(orient='records')
    return data_dict

def load_test(all_data):
    bins = np.arange(0.35, 0.95, 0.05)  # 定义区间
    all_data['bins'] = pd.cut(all_data['star_mass'], bins)  # 将star_mass列分为指定的区间

    data_list = []
    for bin in bins[:-1]:  # 遍历每个区间
        df_bin = all_data[(all_data['mass'] > bin) & (all_data['mass'] < bin + 0.05)]

        random_index = np.random.choice(df_bin.index, 10, replace=False)  # 随机选择10个数据
        for index in random_index:
            data_list.append({
                'star_name': f'star{len(data_list)}',
                'log_Teff': df_bin["log_Teff"].loc[index],
                'log_Teff_err': 0.013,
                'log_g': df_bin["log_g"].loc[index],
                'log_g_err': 0.1,
                'log_he': df_bin["log_he"].loc[index],
                'log_he_err': 0.1,
                'mass_true': df_bin["star_mass"].loc[index],
                'log_L_true': df_bin["log_L"].loc[index],
                'radius_true': df_bin["radius"].loc[index],
                'age_true': df_bin["star_age"].loc[index]
            })
    # data_list = pd.DataFrame(data_list)
    # data_list.to_csv('test_star.csv', index=False)
    return data_list