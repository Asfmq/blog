import matplotlib.pyplot as plt
import numpy as np
import os
import re
import sys

# 获取命令行参数
new_mass = sys.argv[1]
new_env_mass = sys.argv[2]
Zbase = sys.argv[3]
files = "{}-{}-{}.data".format(new_mass, new_env_mass, Zbase)
# path1+'/'+new_mass+'-'+new_env_mass+'-'+Zbase+'.data'
path1 = "./data_remove_env"
data = []
data = np.genfromtxt(path1+'/'+files, names=True, skip_header=5)


key = ['star_mass', 'he_core_mass']

tip_index=-1
star_mass = data[[key[0]][-1]]
he_core_mass = data[[key[1]][-1]]
q = he_core_mass / star_mass
print(q)




