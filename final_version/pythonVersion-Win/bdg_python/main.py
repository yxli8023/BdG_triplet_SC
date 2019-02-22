# 函数库导入
from numpy import *
import numpy as np # 矩阵
import random # 随机数模块
import os # 系统信息获取
from param import *
from matRebulit import *
from check import *
from visualShow import *

# 文件路径
current_path = os.path.abspath(__file__) # 读取当前文件的位置
father_path = os.path.abspath(os.path.dirname(current_path) + os.path.sep + ".") # 获取当前文件的父目录


#数据存储
f = open("param.txt","w")
f.write(" T = ")
f.write(repr(T)) #repr 将数值转换为字符串输出
f.write(" eps = ")
f.write(repr(eps))
f.write(" h = ")
f.write(repr(h))
f.close()

#======================== 自洽过程 ===============
def self_consis():
    global ham,w
    w,ham_diag = np.linalg.eig(ham)
    del_loop = np.zeros((xn,yn))
    del_err = np.zeros((xn,yn))
    del_cal(w,ham_diag)#计算出一组delta值
    num = 0
    ref = 3
    while ref > 1:
        f1 = open(father_path+"\\"+"count.dat","w")
        f2 = open(father_path+"\\"+"check.dat","w")
        ref = 0
        for m in range(xn):
            for l in range(yn):
                del_err[m,l] = delta[m,l] - del_loop[m,l]
                del_loop[m,l] = delta[m,l]
                if abs(real(del_err[m,l])) > eps:
                    ref += 1
                if abs(imag(del_err[m,l])) > eps:
                    ref += 1
        f1.write(str(num))
        f2.write(str(ref))
        num += 1
        matrix_construct(12,w,ham_diag)
        w,ham_diag = np.linalg.eig(ham)
        del_cal(w,ham_diag)
        print(ref)

#============= Start calculation ==============
matrix_construct(0,0,0)
# 矩阵本征值求解
w,ham_diag = np.linalg.eig(ham)
self_consis()
energy_plot(w)
delta_plot(delta)