from numpy import *
import numpy as np # 矩阵
# 变量定义
xn=12
yn=12
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
# 矩阵初始化
ham = np.zeros((N,N),dtype=complex) # Hamilition matrix is complex number
# 在python矩阵的索引从0开始，所以在对ham进行索引取值时，最大不超过N-1
ham_diag = ham
delta = np.zeros((xn,yn)) # order parameter