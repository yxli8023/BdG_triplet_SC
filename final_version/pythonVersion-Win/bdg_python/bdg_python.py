# 函数库导入
from numpy import *
import numpy as np # 矩阵
import random # 随机数模块
import os # 系统信息获取
import seaborn as sns # 热图
import matplotlib.pyplot as plt #作图
from param import *
from diagterm import *
from pariterm import *
from spin_orbital import *
from Hop import *

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

#================== set value complete =========================
def hermitian():
    cc = 0
    f = open(father_path+"\\"+"hermitian.dat","w")
    for m in range(N):
        for l in range(N):
            if ham[m,l] != conj(ham[l,m]):
                string = "("+str(m)+","+str(l)+")"+":"
                f.write(string)
                f.write(str(ham[m,l]))
                f.write("\n")
                cc += 1
    f.close
    print("non-hermitian number:%d"% cc) 
#============================================================
def count():
    cc = 0
    for m in range(N):
        for l in range(N):
            if ham[m,l] != 0:
                cc += 1
    print("non-zero number:%d"% cc) 
#============================================================
def matrix_construct(x):
    diagele()
    pair_Init(x)
    hopping()
    couple()
#==============================================================
def delta_plot(x):
    sns.set()
    ax = sns.heatmap(x,fmt="d",cmap="cool",square=True,xticklabels=2,yticklabels=2)
    ax.invert_yaxis() # y轴反序
    plt.xlabel("x")
    plt.ylabel("y")
    plt.xticks(rotation=0)    # 将字体进行旋转
    #plt.show()
    plt.pause(1.5)
    plt.close()
#==============================================
def energy_plot(w):
    # 本征值作图(python中需要先对本征值排序然后作图)
    f = open(father_path+"\\"+"eigvalue.dat","w")
    for m in range(len(w)):
        a = sort(w)
        temp = real(a[m])
        f.write(str(temp))
        f.write("\n")
    nMin = round(len(a)/2-50)
    nMax = round(len(a)/2+50)
    plt.plot(real(a[nMin:nMax]),"r.")
    plt.pause(1.5)
    plt.close()
#============= Start calculation ==============
matrix_construct(0)

#print(ham)

# 矩阵本征值求解
w,ham_diag = np.linalg.eig(ham)

# 格点位置上序参量的计算
for i in range(xn):
    for j in range(yn):
        m = i*j
        delta[i,j] = abs(pair_energy(m,m,w,ham_diag))
energy_plot(w)
delta_plot(delta)