import os
import re
import pandas as pd
import numpy as np
from tqdm import tqdm
def load_method(paths, file_name):
    # 指定要使用的列
    usecols = ["star_mass", "surface_he4", "surface_h1", "star_age", "log_Teff", "log_g", "log_L", "radius", "center_he4", "he_core_mass"]
    # 用于存储所有数据的列表
    data_list = []

    # 循环处理每个路径
    for path_method in paths:
        # 获取文件列表
        files = sorted(os.listdir(path_method))

        # 循环处理每个文件
        for file in tqdm(files):
            match = re.match(r'(\d+\.\d+)-(\d+\.\d+)-(\d+\.\d+).data', file)
            if match:
                mass = match.group(1)
                env = match.group(2)
                Y_surf = match.group(3)
            if float(mass) < 2.4:continue
            # 使用 Pandas 读取数据文件
            file_path = os.path.join(path_method, file)
            df = pd.read_csv(file_path, skiprows=5, delim_whitespace=True, usecols=usecols)
            df = df.apply(pd.to_numeric, errors='coerce')
            try:
                start_index = np.where(df['center_he4'] < 0.93)[0][0]
            except IndexError:
                # print(f"Error ZAHB: {files[i]}")
                continue
            df = df[start_index:-1]
            df.drop('center_he4', axis=1, inplace=True)
            df['he'] = Y_surf
            df['log_he'] = df['log_he'] = np.log10(df['surface_he4'] / 4 / df['surface_h1'])
            if not df['star_mass'].empty:
                df['mass'] = df['star_mass'].iloc[0]
                df['star_age'] = df['star_age']-df['star_age'].iloc[0]
                if df['star_age'].iloc[0] < 0:
                    print('age', file)
            else:
                print(file)
                continue
            # df['file'] = file
            df.drop('surface_he4', axis=1, inplace=True)
            df.drop('surface_h1', axis=1, inplace=True)
            # print(df)
            if df.empty:
                continue
            data_list.append(df)
    # 将所有的 DataFrame 连接成一个单一的 DataFrame
    data = pd.concat(data_list, ignore_index=True)
    
    data.to_csv(file_name, index=False)
    return data