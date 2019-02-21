# 函数库导入
from numpy import *
import numpy as np # 矩阵
import math # 数学函数
import random # 随机数模块
import os # 系统信息获取
import seaborn as sns # 热图
import matplotlib.pyplot as plt #作图
# 变量定义
xn=4
yn=4
N = xn*yn*4
len2 = xn*yn
len1 = xn
xp = 1
yp = 1
ZeroPoint = xn*yn*2
eps = 1e-5 # error
# 参数定义
h = 0.6 # h is external magnetic field
V = 0 # Point potential
t0 = 1 # hoppping strength
lam = 0.5*t0 # spin-couple strength
mu = -4 # chemical potential
kb = 1 ## Boltzmann constant
a = 1 # lattice constant 
T = 1e-4 #Temperature
#pi = 3.14159265359
beta = 1/(kb*T)
phi0 = pi
B = 2*phi0/(xn*yn) #Zeeman Magnetic
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
# 矩阵初始化
ham = np.zeros((N,N),dtype=complex) # Hamilition matrix is complex number
# 在python矩阵的索引从0开始，所以在对ham进行索引取值时，最大不超过N-1
ham_diag = ham
delta = np.zeros((xn,yn)) # order parameter
# hopping过程中位相考虑
def phi(y):
    return exp(-1j*pi/phi0*B*y)
# 超导配对能(自洽计算)
def pair_energy(x,y,eig_val,eig_vec):
    s = 0+0*1j
    for m in range(N):
        s = s + eig_vec[x,m]*conj(eig_vec[y+len2*2,m])*tanh(eig_val[m]/2*beta)
    return 5.0/2.0*s
# 杂质势
def potential(x,y):
    xy = (y-1)*xn + x
#       Positive energy (正能位置的杂质拖动正能对应的态)
    ham[xy,xy] = ham[xy,xy] + V # position at spin-up oriention   (1,1)block
    ham[xy+len2,xy+len2] = ham[xy+len2,xy+len2] + V  # position at spin-down oriention   (2,2) block
#   Negative energy (负能对应的杂质拖动负能对应的态)
    ham[xy+len2*2,xy+len2*2] = ham[xy+len2*2,xy+len2*2] - V     # (3,3)block
    ham[xy+len2*3,xy+len2*3] = ham[xy+len2*3,xy+len2*3] - V       #(4,4)block
    return

# 对角线元素填充(tested)
def diagele():
	for m in range(xn*yn):
            ham[m,m] = h - mu
            ham[len2+m,len2+m] = -h - mu
            ham[len2*2+m,len2*2+m] = -(-h-mu)
            ham[len2*3+m,len2*3+m] = -(h-mu)
#Pariing energy(tested)
def pair_Init(input):
    if input == 0:
        # 最初用随机数对超导配对能进行赋值
        s = random.random()
        #s = 0.5
        for m in range(xn*yn):
            ham[len2*2+m,m] = s # C-up*C-down    (1,3)block
            ham[m,len2*2+m] = s #  C-down*C-up    (3,1)block 
            ham[len2+m,len2*3+m] = -s  # C-up*C-down   (2,4)block
            ham[len2*3+m,len2+m] = -s  #  E-up*E-down    (4,2)block
    elif input == 100:
        pass
    else:
        for m in range(xn*yn):
            ham[len2*2+m,m] = pair_energy(m,m,w,ham_diag) # C-up*C-down    (1,3)block
            ham[m,len2*2+m] = conj(pair_energy(m,m,w,ham_diag)) #  C-down*C-up    (3,1)block 
            ham[len2+m,len2*3+m] = -pair_energy(m,m,w,ham_diag)  # C-up*C-down   (2,4)block
            ham[len2*3+m,len2+m] = -conj(pair_energy(m,m,w,ham_diag))  #  E-up*E-down    (4,2)block

#========================= hopping term (tested) ===================================================
def hopping():
#------------------------  x-direction right hopping -------------------
    for l in range(yn):
        for m in range(1+l*xn,l*xn+xn-1):
            ham[m,m+1] = -t0*phi(l+1)# (1,1) block
            ham[m,m-1] = -t0*conj(phi(l+1))
            #  (2,2)block
            ham[len2+m,len2+m+1] = -t0*phi(l+1)
            ham[len2+m,len2+m-1] = -t0*conj(phi(l+1))  #(19,18)
#     --------------------------negative energy ------------------------
            #  (3,3)  block
            ham[len2*2+m,len2*2+m+1] = t0*conj(phi(l+1))
            ham[len2*2+m,len2*2+m-1] = t0*(phi(l+1))
            #  (4,4)   block
            ham[len2*3+m,len2*3+m+1] = t0*conj(phi(l+1))
            ham[len2*3+m,len2*3+m-1] = t0*(phi(l+1))
#----------------------------------------- y-direction up hooping -------------------------------------
#----------------------------------------- 最后一行向上跳跃并未考虑------------------------------------
# 由朗道规范可以知道，在y方向的跳跃不包含位相信息(除了最后一行)
# 所以在此处没有考虑最后一行在y方向的跳跃情况
    for m in range(xn,xn*(yn-1)):
        ham[m,m+len1] = -t0 # up(positive) y方向上的向上跳跃
        ham[m,m-len1] = -t0 # up(negative) y 方向上的向下跳跃
        # (2,2)  block
        ham[len2+m,len2+m+len1] = -t0
        ham[len2+m,len2+m-len1] = -t0
        # (3,3) block
        ham[len2*2+m,len2*2+m+len1] = t0
        ham[len2*2+m,len2*2+m-len1] = t0
        # (4,4) block
        ham[len2*3+m,len2*3+m+len1] = t0
        ham[len2*3+m,len2*3+m-len1] = t0
# ------------ first row down-direction hopping ----------------------------------------------------
# 最后一行向第一行跳跃的过程中，人为的加上一个位相，保证整个波函数跳跃一圈的位相是2pi
#     up hopping is positive               down hopping is negative
    for m in range(xn):
# (1,1)block
            ham[m,m+len1] = -t0  #  first row hopping to up direction
            ham[m,xn*(yn-1)+m] = -t0*(phi((m+1)*yn))  #
            ham[xn*(yn-1)+m,m] = -t0*conj(phi((m+1)*yn))  # last row hopping t first row
            ham[xn*(yn-1)+m,xn*(yn-1)+m-len1] = -t0 #last row hopping to down direction
# (2,2)block
            ham[len2+m,len2+m+len1] = -t0  #  first row hopping to up direction
            ham[len2+m,len2+xn*(yn-1)+m] = -t0*(phi((m+1)*yn))  # 
            ham[len2+xn*(yn-1)+m,len2+m] = -t0*conj(phi((m+1)*yn))  # last row hopping t first row
            ham[len2+xn*(yn-1)+m,len2+xn*(yn-1)+m-len1] = -t0 #last row hopping to down direction
# (3,3)block
            ham[len2*2+m,len2*2+m+len1] = t0  #  first row hopping to up direction
            ham[len2*2+m,len2*2+xn*(yn-1)+m] = t0*conj(phi((m+1)*yn))  #
            ham[len2*2+xn*(yn-1)+m,len2*2+m] = t0*(phi((m+1)*yn))  # last row hopping t first row
            ham[len2*2+xn*(yn-1)+m,len2*2+xn*(yn-1)+m-len1] = t0 #last row hopping to down direction 
# (4,4)block
            ham[len2*3+m,len2*3+m+len1] = t0  #  first row hopping to up direction
            ham[len2*3+m,len2*3+xn*(yn-1)+m] = t0*conj(phi((m+1)*yn))  # 
            ham[len2*3+xn*(yn-1)+m,len2*3+m] = t0*(phi((m+1)*yn))  # last row hopping t first row
            ham[len2*3+xn*(yn-1)+m,len2*3+xn*(yn-1)+m-len1] = t0 #last row hopping to down direction
# first coloum left-direction hopping (include last coloum hopping towards first cloum cause of matrix symmetry consider)----
    for m in range(yn):
# (1,1) block
        ham[m*xn,m*xn+1] = -t0*phi(m+1)#第一列上的点向右跳跃  
        ham[m*xn,(m+1)*xn-1] = -t0*conj(phi(m+1))#第一列上的点向右跳跃(周期性边界条件决定它跳跃到最后一列)
        ham[(m+1)*xn-1,m*xn] = -t0*phi(m+1)#最后一列的点向右跳跃
        ham[(m+1)*xn-1,(m+1)*xn-2] = -t0*conj(phi(m+1))#最后一列的点向左跳跃
# (2,2) block
        ham[len2+m*xn,len2+m*xn+1] = -t0*phi(m+1)#第一列上的点向右跳跃   
        ham[len2+m*xn,len2+(m+1)*xn-1] = -t0*conj(phi(m+1))#第一列上的点向右跳跃(周期性边界条件决定它跳跃到最后一列)
        ham[len2+(m+1)*xn-1,len2+m*xn] = -t0*phi(m+1)#最后一列的点向右跳跃
        ham[len2+(m+1)*xn-1,len2+(m+1)*xn-2] = -t0*conj(phi(m+1))#最后一列的点向左跳跃
# (3,3) block
        ham[len2*2+m*xn,len2*2+m*xn+1] = t0*conj(phi(m+1))#第一列上的点向右跳跃   
        ham[len2*2+m*xn,len2*2+(m+1)*xn-1] = t0*(phi(m+1))#第一列上的点向右跳跃(周期性边界条件决定它跳跃到最后一列)
        ham[len2*2+(m+1)*xn-1,len2*2+m*xn] = t0*conj(phi(m+1))#最后一列的点向右跳跃
        ham[len2*2+(m+1)*xn-1,len2*2+(m+1)*xn-2] = t0*(phi(m+1))#最后一列的点向左跳跃
# (4,4) block
        ham[len2*3+m*xn,len2*3+m*xn+1] = t0*conj(phi(m+1))#第一列上的点向右跳跃   
        ham[len2*3+m*xn,len2*3+(m+1)*xn-1] = t0*(phi(m+1))#第一列上的点向右跳跃(周期性边界条件决定它跳跃到最后一列)
        ham[len2*3+(m+1)*xn-1,len2*3+m*xn] = t0*conj(phi(m+1))#最后一列的点向右跳跃
        ham[len2*3+(m+1)*xn-1,len2*3+(m+1)*xn-2] = t0*(phi(m+1))#最后一列的点向左跳跃
            #     spin-orbital couple energy  
#XXXX                              - x direction copule (right and left couple) --------
#----------------------------------positive energy couple (1,2)block----------------------
def couple():
    for l in range(yn):
        for m in range(1+l*xn,xn*(l+1)-1):
            ham[m,len2+m+1] = 1j*lam*phi(l+1) # spin_up-spin_down couple right  正向耦合(1--->2)
            ham[m,len2+m-1] = -1j*lam*conj(phi(l+1))  # left反向耦合(2--->1)
# -------------------------------- positive energy copule (2,1)block----------------
      #(2,1)block
            ham[len2+m,m+1] = 1j*lam*phi(l+1) # spin_up-spin_down couple  right
            ham[len2+m,m-1] = -1j*lam*conj(phi(l+1))  # left
# ------------------------------- negative energy couple (3,4) block-----------------
            ham[len2*2+m,len2*3+m+1] = 1j*lam*conj(phi(l+1))
            ham[len2*2+m,len2*3+m-1] = -1j*lam*(phi(l+1))
# ------------------------------- negative energy couple (4,3) block-----------------------
            ham[len2*3+m,len2*2+m+1] = 1j*lam*conj(phi(l+1)) # spin_up-spin_down couple  right
            ham[len2*3+m,len2*2+m-1] = -1j*lam*(phi(l+1))  # left
# XXXX---------------------------------------------------------------------------
# ----------------------  positive energy couple (1,2)block -----------------------------
    for m in range(xn,xn*(yn-1)):
        ham[m,len2+m+len1] = lam
        ham[m,len2+m-len1] = -lam
# ----------------------------------------- (2,1) block ------------------------------------------------
        ham[len2+m,m+len1] = -lam
        ham[len2+m,m-len1] = lam
# ----------------------------------------- (3,4) block -------------------------------------------------
        ham[len2*2+m,len2*3+m+len1] = lam
        ham[len2*2+m,len2*3+m-len1] = -lam
# -----------------------------------------(4,3) block -----------------------------------------------
        ham[len2*3+m,len2*2+m+len1] = -lam
        ham[len2*3+m,len2*2+m-len1] = lam
## ------------------------ right and left boundary---------------------------------
#------------------------ positive energy ---------------------------------------        
    for m in range(yn):
    #  (1,2)block
        ham[m*xn,len2+m*xn+1] = 1j*lam*phi(m+1)
        ham[m*xn,len2+(m+1)*xn-1] = -1j*lam*conj(phi(m+1)) # 增加一层会使自旋取向改变一次
        ham[(m+1)*xn-1,len2+m*xn] = 1j*lam*phi(m+1)
        ham[(m+1)*xn-1,len2+(m+1)*xn-2] = -1j*lam*conj(phi(m+1))
      #  (2,1)block
        ham[len2+m*xn,m*xn+1] = 1j*lam*phi(m+1)
        ham[len2+m*xn,(m+1)*xn-1] = -1j*lam*conj(phi(m+1)) # 增加一层会使自旋取向改变一次
        ham[len2+(m+1)*xn-1,m*xn] = 1j*lam*phi(m+1)
        ham[len2+(m+1)*xn-1,(m+1)*xn-2] = -1j*lam*conj(phi(m+1))
          
#---------------------------------------negative energy---------------------------------
      #  (3,4)block
        ham[len2*2+m*xn,len2*3+m*xn+1] = -conj(1j*lam*phi(m+1))
        ham[len2*2+m*xn,len2*3+(m+1)*xn-1] = -conj(-1j*lam*conj(phi(m+1))) # 增加一层会使自旋取向改变一次
        ham[len2*2+(m+1)*xn-1,len2*3+m*xn] = -conj(1j*lam*phi(m+1))
        ham[len2*2+(m+1)*xn-1,len2*3+(m+1)*xn-2] = -conj(-1j*lam*conj(phi(m+1)))
        #  (4,3)block
        ham[len2*3+m*xn,len2*2+m*xn+1] = -conj(1j*lam*phi(m+1))
        ham[len2*3+m*xn,len2*2+(m+1)*xn-1] = -conj(-1j*lam*conj(phi(m+1))) # 增加一层会使自旋取向改变一次
        ham[len2*3+(m+1)*xn-1,len2*2+m*xn] = -conj(1j*lam*phi(m+1))
        ham[len2*3+(m+1)*xn-1,len2*2+(m+1)*xn-2] = -conj(-1j*lam*conj(phi(m+1)))
#----------------------- up and down boundary -------------------------------
# --------------------- positive energy ------------------------------------
# -------------向下耦合，第一行会和最后一行进行周期性边界条件的转换----------------            
    for m in range(xn):
        # (1,2)block
        ham[m,len2+m+len1] = lam   # up
        ham[m,len2+xn*(yn-1)+m] = -lam*(phi((m+1)*yn)) # down  
        ham[xn*(yn-1)+m,len2+m] = lam*conj(phi((m+1)*yn))#up  
        ham[xn*(yn-1)+m,len2+xn*(yn-1)-len1+m] = -lam  #down
       # (2,1)block
        ham[len2+m,m+len1] = -lam   #up
        ham[len2+m,xn*(yn-1)+m] = lam*(phi((m+1)*yn))  # down
        ham[len2+m+xn*(yn-1),m] = -lam*conj(phi((m+1)*yn))
        ham[len2+m+xn*(yn-1),m+xn*(yn-1)-len1] = lam

# ------------------------------------ negative energy ----------------------------
      # (3,4)block
        ham[len2*2+m,len2*3+m+len1] = lam
        ham[len2*2+m,len2*3+xn*(yn-1)+m] = -lam*conj(phi((m+1)*yn))
        ham[len2*2+m+xn*(yn-1),len2*3+m] = lam*(phi((m+1)*yn))
        ham[len2*2+m+xn*(yn-1),len2*3+m+xn*(yn-1)-len1] = -lam
      # (4,3)block
        ham[len2*3+m,len2*2+m+len1] = -lam
        ham[len2*3+m,len2*2+xn*(yn-1)+m] = lam*conj(phi((m+1)*yn))
        ham[len2*3+m+xn*(yn-1),len2*2+m] = -lam*(phi((m+1)*yn))
        ham[len2*3+m+xn*(yn-1),len2*2+m+xn*(yn-1)-len1] = lam
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
def count():
    cc = 0
    for m in range(N):
        for l in range(N):
            if ham[m,l] != 0:
                cc += 1
    print("non-zero number:%d"% cc) 
#============================================================
def matrix_construct(x):
    if x == 0 :
        pass
    else:
        ham = 0
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
    plt.pause(3)
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
    plt.pause(3)
    plt.close()
#===============================================
def del_cal():
    # 格点位置上序参量的计算
    for i in range(xn):
        for j in range(yn):
            m = i*j
            delta[i,j] = abs(pair_energy(m,m,w,ham_diag))
#======================== 自洽过程 ===============
def self_consis():
    del_loop = np.zeros((xn,yn))
    del_err = np.zeros((xn,yn))
    del_cal()#计算出一组delta值
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
        print(del_err)
        f1.write(str(num))
        f2.write(str(ref))
        num += 1
        matrix_construct(12)
        w,ham_diag = np.linalg.eig(ham)
        del_cal()
        print(ref)
        #energy_plot(w)
#============= Start calculation ==============
matrix_construct(0)
# 矩阵本征值求解
w,ham_diag = np.linalg.eig(ham)
self_consis()
# 序参量计算
#del_cal()
#energy_plot(w)
#delta_plot(delta)
