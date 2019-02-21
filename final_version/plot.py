# 在运行该文件作图是，请保持该文件与数据文件在相同的文件夹下
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
#预定义变量
xn = 48
yn = 48
# 获取文件相关的路径信息
current_path = os.path.abspath(__file__) # 读取当前文件的位置
father_path = os.path.abspath(os.path.dirname(current_path) + os.path.sep + ".") # 获取当前文件的父目录
def order_plot(file_name):
    data_path = father_path + "\\" + file_name #数据路径
    if os.path.exists(data_path):
        f = open(data_path,"r") # 数据读入
        order = f.read()
        order1 = order.split('\n')
        order1.pop() # 删除文件中最后的空白行
        cc = [] # 收集数值型的order
        for m in range(len(order1)):
            cc.append(eval(order1[m])) #将字符串数据转换为数值
        aa = np.array(cc)
        dd = np.reshape(aa,(xn,yn))
        sns.set()
        ax = sns.heatmap(dd,fmt="d",cmap="cool",square=True,xticklabels=10,yticklabels=10)
        ax.invert_yaxis() # y轴反序
        plt.xlabel("x")
        plt.ylabel("y")
        plt.xticks(rotation=0)    # 将字体进行旋转
        plt.savefig(father_path+"\\"+"order.png",dpi = 800)
        plt.pause(1.5)
        plt.close()
    else:
        print("Your data file isn't exit")

def energy_plot(file_name):
    data_path = father_path + "\\" + file_name #数据路径
    if os.path.exists(data_path):
        f = open(data_path,"r") # 数据读入
        order = f.read()
        order1 = order.split('\n')
        order1.pop() # 删除文件中最后的空白行
        cc = [] # 收集数值型的order
        for m in range(len(order1)):
            cc.append(eval(order1[m])) #将字符串数据转换为数值
        aa = np.array(cc)
        plt.plot(aa,'r.')
        plt.xlim(xmax=4660,xmin=4550)
        plt.ylim(ymax=0.4,ymin=-0.4)
        plt.xlabel("n")
        plt.ylabel("En")
        plt.savefig(father_path+"\\"+"energy.png",dpi = 800)
        plt.pause(1.5)
        plt.close()
    else:
        print("Your data file isn't exit")

def Majorana_plot():
    for i in range(1,3):
        data_path = father_path + "\\" + "gamma" + str(i) + ".dat" #数据路径
        if os.path.exists(data_path):
            f = open(data_path,"r") # 数据读入
            order = f.read()
            order1 = order.split('\n')
            order1.pop() # 删除文件中最后的空白行
            cc = [] # 收集数值型的order
            for m in range(len(order1)):
                cc.append(eval(order1[m])) #将字符串数据转换为数值
            aa = np.array(cc)
            dd = np.reshape(aa,(xn,yn))
            sns.set()
            ax = sns.heatmap(dd,fmt="d",cmap="cool",square=True,xticklabels=10,yticklabels=10)
            ax.invert_yaxis() # y轴反序
            plt.xlabel("x")
            plt.ylabel("y")
            plt.xticks(rotation=0)    # 将字体进行旋转
            plt.savefig(father_path+"\\"+"Majorana"+str(i)+".png",dpi = 800)
            plt.pause(1.5)
            plt.close()
        else:
            print("Your data file isn't exit")

order_plot("order_module.dat")
energy_plot("energy.dat")
Majorana_plot()
print("All of output data and pictures are in "+father_path)
