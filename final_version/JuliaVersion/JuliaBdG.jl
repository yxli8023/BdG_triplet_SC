
using LinearAlgebra
# 变量定义
xn=24
yn=24
N = xn*yn*4
len2 = xn*yn
len1 = xn
xp = 1
yp = 1
ZeroPoint = xn*yn*2
err = 1e-5 # error
# 参数定义
h = 0.6 # h is external magnetic field
V = 0 # Point potential
t0 = 1 # hoppping strength
lambda = 0.5*t0 # spin-couple strength
mu = -4 # chemical potential
kb = 1 ## Boltzmann constant
a = 1 # lattice constant 
T = 1e-4 #Temperature
#pi = 3.14159265359
beta = 1/(kb*T)
phi0 = pi
B = 2*phi0/(xn*yn); #Zeeman Magnetic
image = 1im

ham = zeros(ComplexF32,N,N);
ham_diag = ham;
delta = zeros(Float32,N,N);

father_path = pwd()

#数据存储
f = open("param.dat","w")
write(f," T = ")
write(f,T) #repr 将数值转换为字符串输出
write(f," eps = ")
write(f,err)
write(f," h = ")
write(f,h)
close(f)

function phi(y)
    exp(-1im*pi/phi0*B*y)
end

function pair_energy(x,y,eig_val,eig_vec)
    s = 0+0*1im
    for m=1:N
        s = s + eig_vec[x,m]*conj(eig_vec[y+len2*2,m])*tanh(eig_val[m]/2*beta)
    end
    return 5.0/2.0*s
end

# 杂质势
function potential(x,y)
    xy = (y-1)*xn + x
#       Positive energy (正能位置的杂质拖动正能对应的态)
    ham[xy,xy] = ham[xy,xy] + V # position at spin-up oriention   (1,1)block
    ham[xy+len2,xy+len2] = ham[xy+len2,xy+len2] + V  # position at spin-down oriention   (2,2) block
#   Negative energy (负能对应的杂质拖动负能对应的态)
    ham[xy+len2*2,xy+len2*2] = ham[xy+len2*2,xy+len2*2] - V     # (3,3)block
    ham[xy+len2*3,xy+len2*3] = ham[xy+len2*3,xy+len2*3] - V       #(4,4)block
end

# 对角线元素填充(tested)
function diagele()
	for m=1:xn*yn
            ham[m,m] = h - mu
            ham[len2+m,len2+m] = -h - mu
            ham[len2*2+m,len2*2+m] = -(-h-mu)
            ham[len2*3+m,len2*3+m] = -(h-mu)
    end
end

#Pariing energy(tested)
function pair_Init(input,eig_val,eig_vec)
    if input == 0
        # 最初用随机数对超导配对能进行赋值
        s = rand()
        #s = 0.5
        for m=1:xn*yn
            ham[len2*2+m,m] = s # C-up*C-down    (1,3)block
            ham[m,len2*2+m] = s #  C-down*C-up    (3,1)block 
            ham[len2+m,len2*3+m] = -s  # C-up*C-down   (2,4)block
            ham[len2*3+m,len2+m] = -s  #  E-up*E-down    (4,2)block
        end
    elseif input == 100
        
    else
        for m =1:xn*yn
            ham[len2*2+m,m] = pair_energy(m,m,eig_val,eig_vec) # C-up*C-down    (1,3)block
            ham[m,len2*2+m] = conj(pair_energy(m,m,eig_val,eig_vec)) #  C-down*C-up    (3,1)block 
            ham[len2+m,len2*3+m] = -pair_energy(m,m,eig_val,eig_vec)  # C-up*C-down   (2,4)block
            ham[len2*3+m,len2+m] = -conj(pair_energy(m,m,eig_val,eig_vec))  #  E-up*E-down    (4,2)block
        end
    end
end


#                   hopping  term  
function hopping()
#     ------------------------  x-direction right hopping -------------------
#     ---------------------- positive energy  -----------------------------------------
#     right hopping is positive       left hopping is negatie 
    for l = 1:yn                 # 最后一列向右跳跃没有考虑到，但是向左跳跃考虑率到了
        for m = 2+(l-1)*xn:xn*l-1    # (1,1)block   spin-up
            ham[m,m+1] = -t0*phi(l)  # (i--->i+1)  正向跳跃
            ham[m,m-1] = -t0*conj(phi(l))  # (i+1--->i) 相当于考虑了反向跳跃
    # (2,2)block   spin-forwn
            ham[len2+m,len2+m+1] = -t0*phi(l)
            ham[len2+m,len2+m-1] = -t0*conj(phi(l))
 
#     --------------------------negative energy ------------------------
# (3,3) block
            ham[len2*2+m,len2*2+m+1] = t0*conj(phi(l))
            ham[len2*2+m,len2*2+m-1] = t0*(phi(l))
# (4,4) block
            ham[len2*3+m,len2*3+m+1] = t0*conj(phi(l))
            ham[len2*3+m,len2*3+m-1] = t0*(phi(l))
        end 
    end 
#----------------------------------------- y-direction up hooping -------------------------------------
#----------------------------------------- 最后一行向上跳跃并未考虑------------------------------------
# 由朗道规范可以知道，在y方向的跳跃不包含位相信息(除了最后一行)
# 所以在此处没有考虑最后一行在y方向的跳跃情况

    for m = (xn+1):xn*(yn-1)   # (1,1) block
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
    end
# ------------ first row forwn-direction hopping ----------------------------------------------------
# 最后一行向第一行跳跃的过程中，人为的加上一个位相，保证整个波函数跳跃一圈的位相是2pi
#     up hopping is positive               forwn hopping is negative
    for m = 1:xn   # (1,1)block
        ham[m,m+len1] = -t0  #  first row hopping to up direction
        ham[m,xn*(yn-1)+m] = -t0*(phi(m*yn))  #
        ham[xn*(yn-1)+m,m] = -t0*conj(phi(m*yn))  # last row hopping t first row
        ham[xn*(yn-1)+m,xn*(yn-1)+m-len1] = -t0 #last row hopping to forwn direction
# (2,2)block
        ham[len2+m,len2+m+len1] = -t0  #  first row hopping to up direction
        ham[len2+m,len2+xn*(yn-1)+m] = -t0*(phi(m*yn))  # 
        ham[len2+xn*(yn-1)+m,len2+m] = -t0*conj(phi(m*yn))  # last row hopping t first row
        ham[len2+xn*(yn-1)+m,len2+xn*(yn-1)+m-len1] = -t0 #last row hopping to forwn direction
# (3,3)block
        ham[len2*2+m,len2*2+m+len1] = t0  #  first row hopping to up direction
        ham[len2*2+m,len2*2+xn*(yn-1)+m] = t0*conj(phi(m*yn))  #
        ham[len2*2+xn*(yn-1)+m,len2*2+m] = t0*(phi(m*yn))  # last row hopping t first row
        ham[len2*2+xn*(yn-1)+m,len2*2+xn*(yn-1)+m-len1] = t0 #last row hopping to forwn direction 
# (4,4)block
        ham[len2*3+m,len2*3+m+len1] = t0  #  first row hopping to up direction
        ham[len2*3+m,len2*3+xn*(yn-1)+m] = t0*conj(phi(m*yn))  # 
        ham[len2*3+xn*(yn-1)+m,len2*3+m] = t0*(phi(m*yn))  # last row hopping t first row
        ham[len2*3+xn*(yn-1)+m,len2*3+xn*(yn-1)+m-len1] = t0 #last row hopping to forwn direction
    end      
#------------------------ first coloum left-direction hopping (include last coloum hopping towards first cloum cause of matrix symmetry consider)----
    for m = 1:yn   # (1,1) block
        ham[m*xn,m*xn-(xn-1)] = -t0*phi(m)   # last point hopping to right
        ham[m*xn,m*xn-1] = -t0*conj(phi(m))# last point hopping to left
        ham[1+(m-1)*xn,1+(m-1)*xn+1] = -t0*phi(m)# first point hopping to right
        ham[1+(m-1)*xn,m*xn] = -t0*conj(phi(m))# first point hopping to left
# (2,2) block
        ham[len2+m*xn,len2+m*xn-(xn-1)] = -t0*phi(m)
        ham[len2+m*xn,len2+m*xn-1] = -t0*conj(phi(m))
        ham[len2+1+(m-1)*xn,len2+1+(m-1)*xn+1] = -t0*phi(m)  # first point hopping to right
        ham[len2+1+(m-1)*xn,len2+m*xn] = -t0*conj(phi(m))  # first point hopping to left
# (3,3) block
        ham[len2*2+m*xn,len2*2+m*xn-(xn-1)] = t0*conj(phi(m))
        ham[len2*2+m*xn,len2*2+m*xn-1] = t0*(phi(m))
        ham[len2*2+1+(m-1)*xn,len2*2+1+(m-1)*xn+1] = t0*conj(phi(m))# first point hopping to right
        ham[len2*2+1+(m-1)*xn,len2*2+m*xn] = t0*(phi(m))# first point hopping to left
# (4,4) block
        ham[len2*3+m*xn,len2*3+m*xn-(xn-1)] = t0*conj(phi(m))
        ham[len2*3+m*xn,len2*3+m*xn-1] = t0*(phi(m))
        ham[len2*3+1+(m-1)*xn,len2*3+1+(m-1)*xn+1] =t0*conj(phi(m))# first point hopping to right
        ham[len2*3+1+(m-1)*xn,len2*3+m*xn] = t0*(phi(m))# first point hopping to left
    end
end
                                                                            

#           spin-orbital couple 
function couple()
 
#XXXX                              - x direction copule (right and left couple) --------
#----------------------------------positive energy couple (1,2)block----------------------
    for l = 1:yn                  # 用l代表层数，每一层的最右边的位置上向右边的耦合没有考虑
        for m = 2+(l-1)*xn:xn*l-1     #(1,2)block    
            ham[m,len2+m+1] = image*lambda*phi(l) # spin_up-spin_forwn couple right  正向耦合(1--->2)
            ham[m,len2+m-1] = -image*lambda*conj(phi(l))  # left反向耦合(2--->1)
# -------------------------------- positive energy copule (2,1)block----------------
#(2,1)block
            ham[len2+m,m+1] = image*lambda*phi(l) # spin_up-spin_forwn couple  right
            ham[len2+m,m-1] = -image*lambda*conj(phi(l))  # left
# ------------------------------- negative energy couple (3,4) block-----------------
            ham[len2*2+m,len2*3+m+1] = image*lambda*conj(phi(l))
            ham[len2*2+m,len2*3+m-1] =  -image*lambda*(phi(l))
# ------------------------------- negative energy couple (4,3) block-----------------------
            ham[len2*3+m,len2*2+m+1] = image*lambda*conj(phi(l)) # spin_up-spin_forwn couple  right
            ham[len2*3+m,len2*2+m-1] = -image*lambda*(phi(l))  # left
        end
    end
# XXXX--------------------------------------------------------------------------------------------------------------


#YYYY ------------------------------------ y direction copule (up and forwn couple)加一层---------------------------
# -------------- 不包括最后一行向上耦合的行为 ---------------------------------
# ----------------------  positive energy couple (1,2)block ---------------------------------------------
    for m = xn+1:xn*(yn-1)
        ham[m,len2+m+len1] = lambda
        ham[m,len2+m-len1] = -lambda
# ----------------------------------------- (2,1) block ------------------------------------------------
        ham[len2+m,m+len1] = -lambda
        ham[len2+m,m-len1] = lambda
# ----------------------------------------- (3,4) block -------------------------------------------------
        ham[len2*2+m,len2*3+m+len1] = lambda
        ham[len2*2+m,len2*3+m-len1] = -lambda
# -----------------------------------------(4,3) block -----------------------------------------------
        ham[len2*3+m,len2*2+m+len1] = -lambda
        ham[len2*3+m,len2*2+m-len1] = lambda
    end

# ------------------------ right and left boundary---------------------------------
#------------------------ positive energy ---------------------------------------
    for m = 1:yn   #  (1,2)block
        ham[m*xn,len2+m*xn-(xn-1)] = image*lambda*phi(m)
        ham[m*xn,len2+m*xn-1] = -image*lambda*conj(phi(m)) # 增加一层会使自旋取向改变一
        ham[m*xn-(xn-1),len2+m*xn-(xn-1)+1] = image*lambda*phi(m)
        ham[m*xn-(xn-1),len2+m*xn] = -image*lambda*conj(phi(m))
#  (2,1)block
        ham[len2+m*xn,m*xn-(xn-1)] = image*lambda*phi(m)
        ham[len2+m*xn,m*xn-1] = -image*lambda*conj(phi(m)) # 增加一层会使自旋取向改变一次(1------>48)
        ham[len2+m*xn-(xn-1),m*xn-(xn-1)+1] = image*lambda*phi(m)
        ham[len2+m*xn-(xn-1),m*xn] = -image*lambda*conj(phi(m))  
#---------------------------------------negative energy---------------------------------
#  (3,4)block
        ham[len2*2+m*xn,len2*3+m*xn-(xn-1)] = image*lambda*conj(phi(m))
        ham[len2*2+m*xn,len2*3+m*xn-1] = -image*lambda*(phi(m))# 增加一层会使自旋取向改变一次
        ham[len2*2+m*xn-(xn-1),len2*3+m*xn-(xn-1)+1] = image*lambda*conj(phi(m))
        ham[len2*2+m*xn-(xn-1),len2*3+m*xn] = -image*lambda*(phi(m))
#  (4,3)block
        ham[len2*3+m*xn,len2*2+m*xn-(xn-1)] = image*lambda*conj(phi(m))
        ham[len2*3+m*xn,len2*2+m*xn-1] = -image*lambda*(phi(m)) # 增加一层会使自旋取向改变一次
        ham[len2*3+m*xn-(xn-1),len2*2+m*xn-(xn-1)+1] = image*lambda*conj(phi(m))
        ham[len2*3+m*xn-(xn-1),len2*2+m*xn] = -image*lambda*(phi(m))
    end           

#----------------------- up and forwn boundary -------------------------------
# --------------------- positive energy ------------------------------------
# -------------向下耦合，第一行会和最后一行进行周期性边界条件的转换----------------
# 
    for m = 1:xn       # (1,2)block
        ham[m,len2+m+len1] = lambda   # up
        ham[m,len2+xn*(yn-1)+m] = -lambda*(phi(m*yn)) # forwn  
        ham[xn*(yn-1)+m,len2+m] = lambda*conj(phi(m*yn))#up  
        ham[xn*(yn-1)+m,len2+xn*(yn-1)-len1+m] = -lambda  #forwn
# (2,1)block
        ham[len2+m,m+len1] = -lambda   #up
        ham[len2+m,xn*(yn-1)+m] = lambda*(phi(m*yn))  # forwn
        ham[len2+m+xn*(yn-1),m] = -lambda*conj(phi(m*yn))
        ham[len2+m+xn*(yn-1),m+xn*(yn-1)-len1] = lambda

# ------------------------------------ negative energy ----------------------------
# (3,4)block
        ham[len2*2+m,len2*3+m+len1] = lambda
        ham[len2*2+m,len2*3+xn*(yn-1)+m] = -lambda*conj(phi(m*yn))
        ham[len2*2+m+xn*(yn-1),len2*3+m] =  lambda*(phi(m*yn))
        ham[len2*2+m+xn*(yn-1),len2*3+m+xn*(yn-1)-len1] = -lambda
# (4,3)block
        ham[len2*3+m,len2*2+m+len1] = -lambda
        ham[len2*3+m,len2*2+xn*(yn-1)+m] = lambda*conj(phi(m*yn))
        ham[len2*3+m+xn*(yn-1),len2*2+m] = -lambda*(phi(m*yn))
        ham[len2*3+m+xn*(yn-1),len2*2+m+xn*(yn-1)-len1] = lambda
    end
#YYYY--------------------------------------------------------------------
end

function hermitian()
    cc = 0
    father_path = pwd()
    f = open(father_path+"\\"+"hermitian.dat","w")
    for m = 1:N
        for l = 1:N
            if ham[m,l] != conj(ham[l,m])
                s = "("+string(m)+","+string(l)+")"+":"
                write(f,string)
                write(f,ham[m,l])
                write(f,"\n")
                cc += 1
            end
        end
    end
    close(f)
    print("non-hermitian number:%d"% cc)
end

function countZeros()
    cc = 0
    for m = 1:N
        for l = 1:N
            if ham[m,l] != 0
                cc += 1
            end
        end
    end
    print("non-zero number:%d"% cc) 
end

function matrix_construct(x,eig_val,eig_vec)
    if x == 0 
        
    else
        ham = zeros(ComplexF32,N,N)
    end
    diagele()
    pair_Init(x,eig_val,eig_vec)
    hopping()
    couple()
end

function del_cal(eig_val,eig_vec)
#    global delta
    # 格点位置上序参量的计算
    for i = 1:xn
        for j = 1:yn
            m = i*j
            delta[i,j] = abs(pair_energy(m,m,eig_val,eig_vec))
        end
    end
end

# ======================== 自洽过程 ===============
function self_consis()
#    global ham,w
    w,ham_diag = eigen(ham)
    del_loop =zeros(ComplexF32,xn,yn)
    del_err = zeros(ComplexF32,xn,yn)
    del_cal(w,ham_diag)#计算出一组delta值
    num = 0
    ref = 3
    while ref > 1
        f1 = open("count.dat","w")
        f2 = open("check.dat","w")
        ref = 0
        for m = 1:xn
            for l = 1:yn
                del_err[m,l] = delta[m,l] - del_loop[m,l]
                del_loop[m,l] = delta[m,l]
                if abs(real(del_err[m,l])) > err
                    ref += 1
                end
                if abs(imag(del_err[m,l])) > err
                    ref += 1
                end
            end
        end
        write(f1,num)
        write(f2,ref)
        num += 1
        matrix_construct(12,w,ham_diag)
        w,ham_diag = eigen(ham)
        del_cal(w,ham_diag)
        println(ref)
        print("\t")
    end
end

# ========= 开始计算  ========================
matrix_construct(0,0,0)

#Pkg.add("CPUTime") # 只需要在第一次使用时执行该命令，保证你能有这个拓展包
using CPUTime # 提供程序计算的CPU用时
# u,v=eigen(ham);
@time @CPUtime u,v=eigen(ham);

#import Pkg;
#Pkg.add("GR");# 只需要在第一次使用时执行该命令，保证你能有这个拓展包
using GR;
xMin=xn*yn*2-50;
xMax=xn*yn*2+50;
plot(sort(real(u[xMin:xMax])))

#aa=1:length(real(u[xMin:xMax]));
#scatter(aa,real(u[xMin:xMax]))

self_consis()
