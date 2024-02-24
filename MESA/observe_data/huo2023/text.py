import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d



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

def interp1_data(x, y, y0):
    intersection_points = []
    y=y-y0
    zero_indices = np.where(y[:-1] * y[1:] < 0)[0]
    print(zero_indices)
    for i in range(len(zero_indices)):
        x_sub = x[zero_indices[i]:zero_indices[i]+2]
        y_sub = y[zero_indices[i]:zero_indices[i]+2]
        print(x_sub)
        print(y_sub)
        interp_func = interp1d(y_sub, x_sub)
        intersection_point = interp_func(0)
        intersection_points.append(intersection_point)
    return intersection_points

# 示例数据
# x_data = np.array([1, 2, 3, 4, 5, 4, 3, 2, 1])
# y_data = np.array([4, 2, 1, 3, 6, 5, 3, 1, 4])
# # 示例数据
# # x_data = np.array([1, 2, 3, 4, 5])
# # y_data = np.array([4, 2, 1, 3, 6])

# # 输入参数
# y0 = 2.5
# # interp_data(x_data, y_data, y0)
# # 调用函数
# intersection_points = interp1_data(x_data, y_data, y0)

# # 打印结果
# print("Intersection points (x):", intersection_points)

# # # 绘制结果
# plt.scatter(x_data, y_data, color='r')
# plt.plot(x_data, y_data, '-', label='Data Points')
# plt.axhline(y=y0, color='r', linestyle='--', label=f'y={y0}')
# plt.scatter(intersection_points, [y0]*len(intersection_points), color='g', label='Intersection Points')
# plt.legend()
# plt.show()
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# 生成随机数
x = np.random.randn(10000)

# 计算直方图
hist, bin_edges = np.histogram(x, bins='auto', density=True)

# 计算密度函数
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
pdf = norm.pdf(bin_centers)

# 绘制直方图和密度函数
plt.hist(x, bins='auto', density=True, alpha=0.6, color='g')
plt.plot(bin_centers, pdf, 'r')
plt.show()