from param import *
from subfunction import *
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


# 对角线元素填充(tested)
def diagele():
	for m in range(xn*yn):
            ham[m,m] = h - mu
            ham[len2+m,len2+m] = -h - mu
            ham[len2*2+m,len2*2+m] = -(-h-mu)
            ham[len2*3+m,len2*3+m] = -(h-mu)


#Pariing energy(tested)
def pair_Init(input,eig_val,eig_vec):
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
            ham[len2*2+m,m] = pair_energy(m,m,eig_val,eig_vec) # C-up*C-down    (1,3)block
            ham[m,len2*2+m] = conj(pair_energy(m,m,eig_val,eig_vec)) #  C-down*C-up    (3,1)block 
            ham[len2+m,len2*3+m] = -pair_energy(m,m,eig_val,eig_vec)  # C-up*C-down   (2,4)block
            ham[len2*3+m,len2+m] = -conj(pair_energy(m,m,eig_val,eig_vec))  #  E-up*E-down    (4,2)block
#     spin-orbital couple energy  
def couple():
#XXXX   x direction copule (right and left couple) --------
#----------------------------------positive energy couple (1,2)block----------------------
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

#========================= 矩阵重构 =======================
def matrix_construct(x,eig_val,eig_vec):
    if x == 0 :
        pass
    else:
        ham = np.zeros((N,N),dtype=complex)
    diagele()
    pair_Init(x,eig_val,eig_vec)
    hopping()
    couple()