from param import *
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
#============================================================
def count():
    cc = 0
    for m in range(N):
        for l in range(N):
            if ham[m,l] != 0:
                cc += 1
    print("non-zero number:%d"% cc) 