from param import *
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
        #delread()
        pass
    else:
        for m in range(xn*yn):
            ham[len2*2+m,m] = pair_energy(m,m,w,ham_diag) # C-up*C-down    (1,3)block
            ham[m,len2*2+m] = conj(pair_energy(m,m,w,ham_diag)) #  C-down*C-up    (3,1)block 
            ham[len2+m,len2*3+m] = -pair_energy(m,m,w,ham_diag)  # C-up*C-down   (2,4)block
            ham[len2*3+m,len2+m] = -conj(pair_energy(m,m,w,ham_diag))  #  E-up*E-down    (4,2)block