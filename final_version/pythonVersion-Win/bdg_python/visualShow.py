import os
from numpy import *
import seaborn as sns # 热图
import matplotlib.pyplot as plt #作图
# 文件路径
current_path = os.path.abspath(__file__) # 读取当前文件的位置
father_path = os.path.abspath(os.path.dirname(current_path) + os.path.sep + ".") # 获取当前文件的父目录

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