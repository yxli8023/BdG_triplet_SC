from param import *
from Self_function import *



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