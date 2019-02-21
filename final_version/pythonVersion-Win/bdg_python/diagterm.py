from param import *
# 对角线元素填充(tested)
def diagele():
	for m in range(xn*yn):
            ham[m,m] = h - mu
            ham[len2+m,len2+m] = -h - mu
            ham[len2*2+m,len2*2+m] = -(-h-mu)
            ham[len2*3+m,len2*3+m] = -(h-mu)
