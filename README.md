# BdG equation solve

在超导的问题中，通常会涉及到BdG方程的求解。本代码是用来重复[Vortex pinning by the point potential in topological superconductors:A scheme for braiding Majorana bound states](https://github.com/yxli8023/BdG/blob/master/article/PhysRevB.96.184508.pdf)论文中关于Majorana费米子在拓扑超导中存在的一些问题。

# [BdG equation](https://www.springer.com/cda/content/document/cda_downloaddocument/9783319313122-c2.pdf?SGWID=0-0-45-1576973-p179895452)

哈密顿量本征值问题求解的关键问题在于哈密顿量的对角化，而BdG方程的关键点在于引入[Bogoliubov变换](https://www.wikiwand.com/en/Bogoliubov_transformation)，引 进了准粒子算符来对非对角化的哈密顿量进行对角化处理，该准粒子算符具有粒子空穴对称性。

# Important

关于如何将论文中的模型一步步的转换成Fortran程序，在这里提供一份尚未完善的帮助手册，读者可以借助这份手册同时对照提供的程序，可以帮助读者更好的明白程序中每一步的含义，以及如何去重复该论文的结果，同时也可以在该程序的基础上进行进一步的数值计算。

文件地址为[https://github.com/yxli8023/BdG/blob/master/article/Matrix%20construct.pdf](https://github.com/yxli8023/BdG/blob/master/article/Matrix%20construct.pdf)

# 程序介绍

> 整个项目中有多个文件夹，按英文序号进行命名，记录了作者一步一步对代码的编写，同时里面也包括了对论文中公式到矩阵形式进行转换的想法。
>
> final_version中有两个文件夹Linux 和Windows，分别是程序的Linux版本和windows版本
>
> 每个文件夹中都有程序源代码，已经通过了检验，同时也附上了计算得到的最终数据。
>
> 同时作者还提供了一个由Mathematic编写的简单的绘图程序，可以用来及时的展示计算得到的结果（visual.nb）
>
> 在运行Visual.nb文件时，只需将该文件和数据结果放在同一个文件夹之中即可
>
> final_version中同时有origin的绘图文件pic，用来比较清晰的展示计算结果

![Image text](https://github.com/yxli8023/BdG/blob/master/Pictures/Mathematica_6N0HcvNprH.jpg)

![](https://github.com/yxli8023/BdG/blob/master/Pictures/Origin93_64_hWVZ1U68Dp.png)

> 原始文章在[artile](https://github.com/yxli8023/BdG/blob/master/article/PhysRevB.96.184508.pdf)中

## 蹭热度的我

最近python比较火，也比较好用，开源工具包也比较多，功能齐全，所以在final_version文件夹中，提供了一个利用python的代码来展示计算结果的，不过要注意的是，请将此代码和你的数据文件放在同一个人文件夹下面，否则运行结果可能不正确。

## 参数设置

* xn: x方向格点数
* yn: y方向格点数
* t0: hopping energy
* lambda: couple energy
* T: 系统温度
* eps: 自恰精度

![Image text](https://github.com/yxli8023/BdG/blob/master/Pictures/sublime_text_M0uKPPBlOy.jpg)

# 函数介绍

### diag

>diag 子过程用来对矩阵中的对角线上的元素进行赋值

### hopping 

>hopping 子过程是对哈密顿量中的hopping项进行赋值

### couple 

>copule 子过程代表自旋轨道耦合项

### pair_Init 

>pari_Init 是对超导的配对能进行赋值，在最开始，首先通过随机数对每个格点位置的配对能进行赋值，然后通过自洽求解逐渐求解出每个位置的配对能，程序中的eps就是用来衡量每一次自洽过程之后，每个格点位置的配对能$\Delta$与上一次之间的差，当每个位置上的差都小于eps时，程序终止。

### potential(xp,yp)

> 该子过程用来在晶格点阵中增加杂质势，(xp,yp)则代表了杂质势在二维晶格中的位置，同时参数V则代表了杂质势的大小

### eig

> 哈密顿矩阵本征值与本征矢量的求解

### delRead

> 该子过程用来读取一个数据文件，将读取出来的值作为每个格点位置上的超导配对能

### delta

> 该函数是自洽过程中很重要的一步，通过求解得到的本征值和本征矢量来计算出每个格点位置的超导配对能

### Majorana

> 分离展示每个Majorana零能模在晶格点阵中的位置

### matrix_verify

> 验证矩阵是否为厄密矩阵，用作程序调试检验

### count_number

> 统计矩阵中非零元素的个数，用作程序调试检验

# 输出文件

### param.dat

> 用来保存程序运行时一些物理参数的值

### eigenvalue.dat

> 程序结束之前会输出的最终哈密顿量的本征值

### order.dat

> 计算得到的每个格点位置的$\Delta$，即就是超导序参量

### order_module.dat

> 计算得到的每个格点位置的$\Delta$的绝对值，即就是超导序参量的绝对值

### energy.dat

> 用来记录每一次自洽过程中求解得到的本征值，每一次自洽都会将上一次的数据重写

### gamma1.dat

> Majorana零能模之一

### gamma2.dat

> Majorana零能模之一

# 结果展示

![Image text](https://github.com/yxli8023/BdG/blob/master/Pictures/del.jpg)

![Image text](https://github.com/yxli8023/BdG/blob/master/Pictures/energy.jpg)



# 运行平台

### Linux Platform

编译命令

> ifort -mkl code.f -o a.x

> -mkl 用来链接对角化的函数库

> -mcmodel 可以作为备用选项，函数编写过程中，如果不采用动态数组声明的方式，可以直接通过使用

> ifort -mkl -mcmodel=large code.f -o a.x 直接编译，可以避免由于数组过大引起的编译器报错

### Windows Platform

> windows平台下利用Visual Studio 2017 + intel visual fortran 2019 进行编译运行

> 由于计算过程中会涉及到较大矩阵的计算，所以需要对编译器进行一定的设置

1.先在VS中打开项目，并右键单击如图所示位置

![](https://github.com/yxli8023/BdG/blob/master/Pictures/devenv_Afw3Vxj1VI.jpg) 

点击属性选项，进入下面的界面，相应的设置如图所示

![](https://github.com/yxli8023/BdG/blob/master/Pictures/devenv_72vHIPRTBw.jpg)

![devenv_tDihBL2BXM](https://github.com/yxli8023/BdG/blob/master/Pictures/devenv_tDihBL2BXM.jpg)

### [补充内容(VS+IVF安装)]()



# 物理内容

> 在本篇论文中，哈密顿量是厄密，对应着[厄密矩阵](https://www.wikiwand.com/en/Hermitian_matrix)，其本征值是实数。

> 正则变换后，存在Partical-Hole symmetric，本征值应该是正负对称

# Tips

> 在原始论文中，晶格大小为48*48，转换成矩阵形式后为48 * 48 * 4 的方阵，矩阵大小偏大，计算过程会非常缓慢，可以通过调节晶格的大小，利用小晶格先验证程序的正确与否，比如通过本征值去检验粒子空穴对称性，当对称性满足后可以通过序参量的绝对值，利用绘图的方式去观察零能点出现的位置。

> 在程序使用过程中如果有任何的问题，或者您对程序的某些不足的地方进行了修改，都请联系我
>
> email:  yxli@m.scnu.edu.cn