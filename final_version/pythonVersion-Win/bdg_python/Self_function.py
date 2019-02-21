from param import *
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