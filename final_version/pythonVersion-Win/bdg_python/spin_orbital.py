from param import *
from Self_function import *
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