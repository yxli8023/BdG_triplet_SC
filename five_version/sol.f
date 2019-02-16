      module param
      implicit none
      ! two direction x and y,the index are i and j
      ! i=48  j=48  
      ! include two spin and create and elimination operator so have 4
      ! types
      integer(kind=8),parameter::N=48*48*4
      real*8,parameter::pi=3.1415926535
      real h ! h is external magnetic field
      real V ! V is point potential
      real a ! Lettice constant
      real kb ! Boltzmann constant
      real T ! Temperature
      real t0 ! t0 is the couple strength of Hubbrd model
      real lambda ! lambda is the spin-orbital coupling strength
      complex(kind=8)::image=(0,1.0) ! Image number
      real mu ! mu is the chemcial potential
      real(kind=8) Hamiltonian(N,N)
      integer l,m ! Both of them are loop variation
      integer i,j ! i and j is index of x and y
      complex::phi_ij=(2.0,1.0)
      integer,parameter::x_len = 48
      integer,parameter::y_len = 48
      integer,parameter::len2 = 48*48
      !complex::Ht = (0.0,0.0)  ! Simplify parameter
      !complex::Hso = (0.0,0.0) ! Simplify parameter
      !complex::del = (0.3,0.4)
      real::Ht=12.1  ! Simplify parameter
      real::Hso=12.8 ! Simplify parameter
      real::del=4.9
      integer counter ! loop number counter
      integer num ! test number
      ! lapack parameter 
      integer::lda=N ! row of matrix
      real::w(N) ! eigen-value
      integer,parameter::lwork=10*N+2*N**2
      integer iwork(N)
      integer::liwork=6*N
      real(kind=8)work(lwork)
      integer info
      end module param
      !call dsyevd('V','U',N,Hamiltonian,lda,w,work,lwork,iwork,liwork,info)      
      !call ssyev('V','U',N,a,lda,w,work,lwork,info)
      program bdg
      ! Author :YuXuanLi
      ! E-mail:yxli406@gmail.com
      ! Article:Phys.Rev.B 96,184508
      use param
      open(15,file='row.dat')! 测试矩阵的赋值情况
c===============================Parameter Setting=======================
      h = 0.6
      mu = -4
      V = 0
      t0 = 1.0
      lambda = 1.0
c==============================Construct Hamiltonian Matrix================
      Hamiltonian = 0   ! Initinal Matrix
c----------------------------- Setting the diagonal element---------------      
      ! spin-up create operator
      do l=1,int(N/4)
            Hamiltonian(l,l)=h-mu
      end do
      ! spin-down create operator
      do l=N/4+1,int(N/4)*2
            Hamiltonian(l,l)=-h-mu
      end do      
      ! spin-down elimination operator
      do l=int(N/4)*2+1,int(N/4)*3
            Hamiltonian(l,l)=-h-mu
      end do
      ! spin-up elimination operator
      do l=int(N/4)*3+1,N
            Hamiltonian(l,l)=h-mu
      end do
c----------------------------- Diagonal element setting complete-----------
      ! add a point potential at (24,24)
      Hamiltonian(24,24)=V ! position at spin-up oriention
      Hamiltonian(24+48,24+48)=V ! position at spin-down oriention



c----------------------------- periodic boundary condition----------------
! Cause of periodic boundary condition was used that point which live in the boundary 
! are special can be condiserd first
! Boundary point include (1,:) (:,1) (48,:) (:,48) when we construct matrix encounter these
! points we should think in the different direction(x or y) (1,:) and (48,:) or (:,1) and (:,48)
! are nearest-neighbor point each other.
c============================  Hamiltonian=(48*48*4)^2=N^2  ====================
c  网格点数为  48*48 每个格点处包括产生算符与消灭算符并区分自旋取向，即每个点都需要四个算符来描述
c     将算符形式转换为矩阵形式的时候，也就是矩阵的大小为(48*48*4)^2
c  在计算过程中选取周期性边界条件：(下面所说的行列是格点的行列，并不是矩阵的行列)
c     1. 第1行的下面为第48行，第48行的上面为第1行
c     2. 第1列的左边为第48列，第48列的右端为第1列
c  将算符形式哈密顿量转换为矩阵形式的时候有如下注意点：
c     在格点的表示下，左右相邻的两个格点，在矩阵中位置上也属于相邻(+1or-1)，但是在格点中上下相邻的点
c     在矩阵中表示的时候，相差为48，这是由于算符形式在向矩阵形式转变时，构造上产生的一种形式变化，且注意
c     这种相邻的行为描述不包括格点边界上的点，因为周期性边界条件的原因，导致边界上点在构造矩阵时有一些不同
c  算符向矩阵构造的过程中，可以形象的将矩阵分成四大块(在行方向)：
c     O={C-up,C-down,E-down,E-up}H{E-up,E-down,C-down,C-up}'      H 则是需要构造的矩阵
c     C-up:自旋向上产生算符 C-down:自旋向下产生算符  E-up:自旋向上消灭算符 E-down:自旋向下消灭算符
c     每一个算符都是一个48*48的列向量(要保证哈密顿量的构造中包含每一个格点)===>  H 是一个 (48*48*4)^2 的矩阵
c     当从一个确定格点的C-up到C-down的过程中，在矩阵的表示上为(i,i)->(i,i+48*48)  i 共有48*48*4个
c     所以矩阵在行方向上可以分成4块，每一块分别和上面算符表达式中的含义对应
c     对于构造完成的矩阵，通过矩阵形式和算符形式的相互转换可以得出，矩阵第1列中的数都对应着格点中(1,1)这个
c     位置上自旋向上的消灭算符前的系数，剩下的对应项依次类推(包括自旋向下消灭，自旋向下产生，自旋向上产生)
c=============================== start down boundary ============================================================                        
      do i = 0,1
c                 hopping energy    跳跃是同一块中的操作
            ! left-right hopping
            Hamiltonian(1+len2*i,1+2) = Ht  ! right hopping  
            Hamiltonian(1+len2*i,48) = Ht  ! left hopping 向右跳跃时利用周期性边界条件跳跃到最右端
            ! up-down hopping
            Hamiltonian(1+len2*i,48*1+1) = Ht ! 向上跳跃时加一层
            Hamiltonian(1+len2*i,48*47+1) = Ht! 向下跳跃时通过周期性边界条件到最上端的一层
c                 spin-orbital couple    耦合作用涉及到不同块之间的操作，一般是跳跃到下一个块去耦合作用  每个块之间间隔大小为48*48
            ! left-right couple
            Hamiltonian(1+len2*i,len2+2) = Hso ! right  代表第2个位置上的自旋向下产生算符  48*48=len2 增加了一个块
            Hamiltonian(1+len2*i,len2+48) = Hso ! left   (1,1)格点位置向左方耦合为最右端的格点  一个len2增加一个块
            ! up-down couple
            Hamiltonian(1+len2*i,len2+48*1+1) = Hso !up  上下跳跃的时候x方向的坐标是相同的，故为+1   而+48*1代表跳跃了一层
            Hamiltonian(1+len2*i,len2+48*47+1) = Hso ! down !cause of periodic first row and last row are nearest-neighbor
c                superconductor pairing energy
            Hamiltonian(1+len2*i,len2*2+1) = del !
            !write(15,*)1+len2*i
c=================以上是计算了(1,1)这个格点位置上E-down前的所有系数====================================================
c=================下面是计算从(2,1)-->(47,1)格点上E-down算符前所有系数==================================================
            do m=2,47
c                 hopping energy
                  ! left-right hopping
                  Hamiltonian(m+len2*i,m+1) = Ht  ! right hopping
                  Hamiltonian(m+len2*i,m-1) = Ht  ! left hopping
                  ! up-down hopping
                  Hamiltonian(m+len2*i,48*1+m) = Ht
                  Hamiltonian(m+len2*i,48*47+m) = Ht
c                 spin-orbital couple
                  ! left-right couple
                  Hamiltonian(m+len2*i,len2+m+1) = Hso ! right
                  Hamiltonian(m+len2*i,len2+m-1) = Hso ! left  对于(2,1)-->(47,1)他们属于中间位置，左右方向都是有确定的格点，不会涉及到周期性边界的问题
                  ! up-down couple
                  Hamiltonian(m+len2*i,len2+48*1+m) = Hso !up
                  Hamiltonian(m+len2*i,len2+48*47+m) = Hso ! down (2,1)-->(47,1)在向下耦合时，遇到的边界问题相同，都是周期性边界条纹
c                superconductor pairing energy
                  Hamiltonian(m+len2*i,len2*2+m) = del ! conjg is complex conjugate function
                  !write(15,*)m+len2*i
            end do
c==============下面来计算(48,1)这个格点，该格点的情况与(1,1)情况几乎相同=============================================
c                 hopping energy
                  ! left-right hopping
                  Hamiltonian(48+len2*i,1) = Ht  ! right hopping  周期性边界说明48右端是(1,1)
                  Hamiltonian(48+len2*i,47) = Ht  ! left hopping  左端是(47,1)格点位
                  ! up-down hopping
                  Hamiltonian(48+len2*i,48*1+48) = Ht ! up   
                  Hamiltonian(48+len2*i,48*47+48) = Ht ! down   向下跳跃时继续遇到周期性边界条件，连接到最上面的一行
c                 spin-orbital couple
                  ! left-right couple
                  Hamiltonian(48+len2*i,len2+1) = Hso ! right 向右跳跃会遇到周期性边界条件，连接到(1,1)
                  Hamiltonian(48+len2*i,len2+48-1) = Hso ! left  一个len2跳跃一块  +48 到跳跃块的对应位置上
                  ! up-down couple
                  Hamiltonian(48+len2*i,len2+48*1+48) = Hso !up
                  Hamiltonian(48+len2*i,len2+48+48*47) = Hso ! down  向下进行耦合时满足周期性边界条件---->(48,48)
c                superconductor pairing energy
                  Hamiltonian(48+len2*i,len2*2+1) = del ! conjg is complex conjugate function
                  !write(15,*)48+len2*i
c=============================  down boundary complete=============================================================
                  


c============================  start left boundary=================================================================
c   进行左边界计算时每增加一次位置，矩阵中的位置变化为48
            do m = 2,47
c                 hopping energy
                  ! left-right hopping
                  Hamiltonian(1+(m-1)*48+len2*i,(m-1)*48+1+1) = Ht  ! right hopping
                  Hamiltonian(1+(m-1)*48+len2*i,48+(m-1)*48) = Ht  ! left hopping 向左方向跳跃，直接跳跃到右边界
                  ! up-down hopping
                  Hamiltonian(1+(m-1)*48+len2*i,1+(m-1)*48+48) = Ht ! up 向上跳跃增加一层   层的增加为48  块的增加为48*48
                  Hamiltonian(1+(m-1)*48+len2*i,1+(m-1)*48-48) = Ht ! down 向下跳跃减少一层
c                 spin-orbital couple
                  ! left-right couple
                  Hamiltonian(1+(m-1)*48+len2*i,len2+1+(m-1)*48+1) = Hso! right  对于耦合 首先增加一个块，然后再从对角位置向左或者向右跳跃
                  Hamiltonian(1+(m-1)*48+len2*i,len2+(m-1)*48+48) = Hso! left 对于左边界上的点，向左跳跃碰到周期性边界条件，直接增加到左右边就可以
                  ! up-down couple
                  Hamiltonian(1+(m-1)*48+len2*i,len2+1+(m-1)*48+48*1)=Hso !up 在原来的位置上直接增加一层即可
                  Hamiltonian(1+(m-1)*48+len2*i,len2+1+(m-1)*48-48*1)=Hso !down (2,1)->(47,1)属于中间位置，上下跳跃或者耦合不会出现周期性边界条件的问题 
c                superconductor pairing energy
                  Hamiltonian(1+(m-1)*48+len2*i,len2*2+1+(m-1)*48) = del
                  !write(15,*)1+(m-1)*48+len2*i
            end do
c===========================  (48,1) ============================================================================================
c                 hopping energy
                  ! left-right hopping
                  Hamiltonian(48*47+1+len2*i,48*47+1+1) = Ht  ! right hopping
                  Hamiltonian(48*47+1+len2*i,48*47+48) = Ht  ! left hopping
                  ! up-down hopping
                  Hamiltonian(48*47+1+len2*i,1) = Ht ! up hopping
                  Hamiltonian(48*47+1+len2*i,48*47+1-48) = Ht ! down hopping 
c                 spin-orbital couple
                  ! left-right couple
                  Hamiltonian(48*47+1+len2*i,len2+48*47+1+1) = Hso ! right
                  Hamiltonian(48*47+1+len2*i,len2+48*48) = Hso ! left
                  ! up-down couple
                  Hamiltonian(48*47+1+len2*i,len2+1) = Hso !up
                  Hamiltonian(48*47+1+len2*i,len2+48*47+1-48) = Hso ! down !cause of periodic first row and last row are nearest-neighbor
c                superconductor pairing energy
                  Hamiltonian(48*47+1+len2*i,len2*2+48*47+1) = del ! conjg is complex conjugate function
                  !write(15,*)48*47+1+len2*i
c==========================  left boundary complete ================================================================
                              



c=========================   up boundary start ==================================================================
            do m = 48*47+2,48*48-1
c                 hopping energy
                  ! left-right hopping
                  Hamiltonian(m+len2*i,m+1) = Ht  ! right hopping
                  Hamiltonian(m+len2*i,m-1) = Ht  ! left hopping
                  ! up-down hopping
                  Hamiltonian(m+len2*i,m-48*47) = Ht ! up hopping  上边界向上跳跃是跳跃到下边界
                  Hamiltonian(m+len2*i,m-48*1) = Ht ! down hopping  上边界向下跳跃是减去一层
c                 spin-orbital couple
                  ! left-right couple
                  Hamiltonian(m+len2*i,len2+m+1) = Hso ! right
                  Hamiltonian(m+len2*i,len2+m-1) = Hso ! left
                  ! up-down couple
                  Hamiltonian(m+len2*i,len2-48*47+m) = Hso !up
                  Hamiltonian(m+len2*i,len2-48*1+m) = Hso ! down !cause of periodic first row and last row are nearest-neighbor
c                superconductor pairing energy
                  Hamiltonian(m+len2*i,len2*2+m) = del 
                  !write(15,*)m+len2*i
            end do
c====================================================(48,48)============================================
c                 hopping energy
                  ! left-right hopping
                  Hamiltonian(48*48+len2*i,48*47+1) = Ht  ! right hopping 最右端的点向右跳跃满足周期性边界条件跳跃到最左端
                  Hamiltonian(48*48+len2*i,48*48-1) = Ht  ! left hopping
                  ! up-down hopping
                  Hamiltonian(48*48+len2*i,48) = Ht ! up hopping 向上跳跃到(1,48)
                  Hamiltonian(48*48+len2*i,48*47) = Ht ! down hopping 向下跳跃到(47,48)
c                 spin-orbital couple
                  ! left-right couple
                  Hamiltonian(48*48+len2*i,len2+48*47+1) = Hso ! right
                  Hamiltonian(48*48+len2*i,len2+48*48-1) = Hso ! left
                  ! up-down couple
                  Hamiltonian(48*48+len2*i,len2+48) = Hso !up
                  Hamiltonian(48*48+len2*i,len2+48*47) = Hso ! down !cause of periodic first row and last row are nearest-neighbor
c                superconductor pairing energy
                  Hamiltonian(48*48+len2*i,len2*2+48*48) = del ! conjg is complex conjugate function
                  !write(15,*)48*48+len2*i
c=========================== up boundary complete =========================================================


c=========================== right boundry start ===========================================================
            do m = 2,47
c           hopping energy
                  ! left-right hopping
                  Hamiltonian(48*m+len2*i,(m-1)*48+1) = Ht  ! right hopping  最右端的点通过周期性边界条件跳跃到最左端
                  Hamiltonian(48*m+len2*i,m*48-1) = Ht  ! left hopping
                  ! up-down hopping
                  Hamiltonian(48*m+len2*i,48*(m+1)) = Ht
                  Hamiltonian(48*m+len2*i,48*(m-1)) = Ht
c                 spin-orbital couple    //  While hopping to the next block
                  ! left-right couple
                  Hamiltonian(48*m+len2*i,len2+(m-1)*48+1) = Hso ! right
                  Hamiltonian(48*m+len2*i,len2+48*m-1) = Hso ! left
                  ! up-down couple
                  Hamiltonian(48*m+len2*i,len2+48*(m+1)) = Hso !up
                  Hamiltonian(48*m+len2*i,len2+48*(m-1)) = Hso ! down !cause of periodic first row and last row are nearest-neighbor
c                superconductor pairing energy
                  Hamiltonian(48*m+len2*i,len2*2+48*m) = del ! conjg is complex conjugate function
                  !write(15,*)48*m+len2*i
            end do
      
      end do
c========================== right boundary complete  =======================================================
c========================== Above all have benn check correct ==============================================
! According periodic condiction only boundary point have something trouble,we can handle point without boundary
! In the eare without boundary hopping and couple condition is identical







c====================================== spin-up eliminate operator===========================================
c       Hopping energy
      ! left-right hop energy (x-direction)
      do l = 2,47
      do m=48+2,48*l-1      !
           Hamiltonian(m,m-1) = Ht ! left
           Hamiltonian(m,m+1) = Ht ! right
c    ! up-down hop energy(y-driection)
           Hamiltonian(m,48*1+m) = Ht ! up layer have same x coordinate
           Hamiltonian(m,m-48*1) = Ht ! down layer  last: 48+48*47-1 
c      Spin-orbital couple energy
c     ! left-right couple energy
           Hamiltonian(m,len2+m-1) = Hso ! next block  增加一个block 即增加len2
           Hamiltonian(m,len2+m+1) = Hso ! next block
c     ! up-down couple energy
           Hamiltonian(m,len2-48*1+m) = Hso ! down 增加一层为48
           Hamiltonian(m,len2+48*1+m) = Hso ! up 
c        superconductor piring energy
           Hamiltonian(m,len2*2+m) = del ! last:48*48*2+48*47-1
            !write(15,*)m
      end do
      end do
c========================================  spin-down elinate operator ================================================
c       Hopping energy
      ! left-right hop energy (x-direction)
      do l = 2,47
      do m=len2+48+2,len2+48*l-1 !second block // decrease 1 is ensure avoid boundary point first
            Hamiltonian(m,m-1) = Ht
            Hamiltonian(m,m+1) = Ht 
c     ! up-down hop energy(y-driection)
            Hamiltonian(m,48*1+m) = Ht ! up layer   // up-down hopping points have same x-coordinate
            Hamiltonian(m,m-48*1) = Ht ! down layer
c       Spin-orbital couple energy
c      ! left-right couple energy
            Hamiltonian(m,m+1-len2) = Hso
            Hamiltonian(m,m-1-len2) = Hso
c      ! up-down couple energy
            Hamiltonian(m,m-len2-48*1) = Hso
            Hamiltonian(m,m-len2+48*1) = Hso
c         superconductor piring energy
            Hamiltonian(m,len2*2+m) = del ! In do loop has increase 1 layer
            !write(15,*)m
      end do
      end do
c========================================== spin-down create operator ======================================
c       Hopping energy
      ! left-right hop energy (x-direction)
c      do l=2,47
c      do m=2*len2+48+2,2*len2+48*l-1  ! (4658-->6863)
c            Hamiltonian(m,m-1) = Ht
c            Hamiltonian(m,m+1) = Ht 
c     ! up-down hop energy(y-driection)
c            Hamiltonian(m,48*1+m) = Ht ! up layer
c            Hamiltonian(m,m-48*1) = Ht ! down layer
c       Spin-orbital couple energy
c      ! left-right couple energy
c            Hamiltonian(m,m+len2+1) = Hso
c            Hamiltonian(m,m+len2-1) = Hso
c      ! up-down couple energy
c            Hamiltonian(m,m+len2-48*1) = Hso
c            Hamiltonian(m,m+len2+48*1) = Hso
c         superconductor piring energy
c            Hamiltonian(m,m-2*len2) = del ! In do loop has increase 1 layer
            !write(15,*)m
c      end do
c      end do
c================================================ spin-up create operator =====================================
c       Hopping energy
c      ! left-right hop energy (x-direction)
c      do l = 2,47
c      do m=3*len2+48+2,3*len2+48*l-1 
c            Hamiltonian(m,m-1) = Ht
c            Hamiltonian(m,m+1) = Ht 
c     ! up-down hop energy(y-driection)
c            Hamiltonian(m,48*1+m) = Ht ! up layer
c            Hamiltonian(m,m-48*1) = Ht ! down layer
c       Spin-orbital couple energy
c      ! left-right couple energy
c            Hamiltonian(m,m-len2+1) = Hso
c            Hamiltonian(m,m-len2-1) = Hso
c      ! up-down couple energy
c            Hamiltonian(m,m-len2-48*1) = Hso
c            Hamiltonian(m,m-len2+48*1) = Hso
c         superconductor piring energy
c            Hamiltonian(m,m-3*len2) = del ! In do loop has increase 1 layer
            !write(15,*)m
c      end do
c      end do
c==============================================================================
C================================Save Matrix with .dat file type===============
c      open(11,file="hamiltonian.dat")
c      !write(11,*)Hamiltonian(24,24)
c      close(11)
C===============================Disp something to verify code correct==========
      counter = 0
      open(16,file="position.dat")
      open(17,file='diag-value.dat')
      do m = 1,48*48*4
            do l=1,48*48*4
                  !write(17,*)Hamiltonian(m,m)
                  if (Hamiltonian(m,l) .ne. 0) then
                       counter = counter + 1
                       ! write(16,"('('I6','I6')')")m,l
                       !write(11,"(48(F9.4))",advance='No')Hamiltonian(m,1)
                  end if
            end do
      enddo
      close(16)
      close(17)
      open(14,file='count-nozero.dat')
      if (counter/N .ne. 9.0)then
            write(*,*)"The number of non-zero is:",counter
            write(*,*)"You are code isn't correct"
            stop
      end if
      write(14,*)counter
      close(14)
      !write(11,"(48(F9.4))",advance='No')Hamiltonian
c      do m=1,10
c            do l = 1,10 
c                  write(11,*)Hamiltonian(m,l)
c            end do
c      end do
c      close(11)
      !call dsyevd('V','U',N,Hamiltonian,lda,w,work,lwork,info)     ! Need lapack and blas 
      call ssyev('V','U',N,Hamiltonian,lda,w,work,lwork,info)       ! intel fortran mkl package
      open(12,file='eig-value.dat')
      write(12,*)w
      close(12)
      do m = 1,10
            do l = 1,10
                  open(13,file='eig-vector.dat')
                  write(13,*)Hamiltonian(m,l)
            end do
      end do
      close(13)
      stop 
      end


!==============================Function construct ==========================

!      real function delta(i,eig_value,eig_vector)
!      use param
!      real eig_value(N),eig_vector(N,N)
!      do m = 1,   
!      tanh
!      end function delta
