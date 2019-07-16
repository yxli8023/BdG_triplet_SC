!      Author :YuXuanLi
!      E-mail:yxli406@gmail.com
!      Article:Phys.Rev.B 96,184508
      module param
      implicit none
      ! two direction x and y,the index are i and j
      ! i=48  j=48  
      ! include two spin and create and elimination operator so have 4
      ! types
      integer,parameter::N=48*48*4
      real,parameter::pi = 3.14159265359
      real h ! h is external magnetic field
      real V ! V is point potential
      real a ! Lettice constant
      real kb ! Boltzmann constant
      real T ! Temperature
      real t0 ! t0 is the couple strength of Hubbrd model
      real lambda ! lambda is the spin-orbital coupling strength
      complex*8::image=(0.0,1.0) ! Image number
      real mu ! mu is the chemcial potential
      real phi0 ! superconduct flux quantum
      real B ! external Magnetic
      integer,parameter::len2 = 48*48
      integer,parameter::len1 = 48
      real,parameter::eps = 0.0001 ! measure distance of two delta
      real del(len2)  
      complex*8::Ham(N,N)=(0,0)
      real delta_init        ! The initial value of delta and used in self-consistently loop
      real del_max
      ! lapack parameter 
      integer::lda = N
      integer(kind=8),parameter::lwmax = 2*48*48*4+(48*48*4)**2
      real,allocatable::w(:)
      complex*8,allocatable::work(:)
      real,allocatable::rwork(:)
      integer,allocatable::iwork(:)
      integer lwork   ! at least 2*N+N**2
      integer lrwork    ! at least 1 + 5*N +2*N**2
      integer liwork   ! at least 3 +5*N
      integer info
      end module param
!===============================Program Start=========================================
      program bdg
      use param
      integer l,m ! Both of them are loop variation
      complex*8,external::phi ! function statement
      complex*8,external::delta ! function statemen
      external::cheev
      allocate(w(N))
      allocate(work(lwmax))
      allocate(rwork(lwmax))
      allocate(iwork(lwmax))

c===============================Parameter Setting=======================
      phi0 = pi! flux quantum
      a = 1 ! lattice constant
      h = 0 ! zeeman  field
      mu = -4  ! chemical potential
      V = 0    ! point potential
      t0 = 1   ! hopping energy
      lambda = 0   ! couple energy
      B = 2*phi0/(48*48-1)  ! magnetic field
      Kb = 0.1    !
      T = 0.1     ! temperature
      delta_init = 1   ! initial order parameter
c==============================Construct Ham Matrix================
!      wr = 1.0
!      vl  = 1.0
c--------------------------- point potential element setting complete-----------
      ! add a point potential at (24,24)
!      Ham(24,24)=V ! position at spin-up oriention
!      Ham(24+48,24+48)=V ! position at spin-down oriention
      call diag()   ! setting diag to matrix
      call hopping() ! full of elementary about hopping energy
      call couple()  ! couple energy elementary put in matrix 
      call pair(1)     ! 在求解过程中要用自恰的方法，所以需要对超导配对能进行重新赋值
      call count_number()    ! 统计哈密顿矩阵中非零元素个数
!      write(*,*)image*lambda*phi(1)
      stop
!=================================== eig-value and eig-vector solve ====================================
!      call matrix_verify()   ! 检验矩阵是否为厄米矩阵
!      stop
      lwork = -1
!      liwork = -1
!      lrwork = -1
      call cheev('Vectors','Upper',N,Ham,lda,w,work,lwork,rwork,info)
!      call cheevd('Vectors','Upper',N,Ham,lda,w,work,lwork,rwork,lrwork
!     & ,iwork,liwork,info)
      lwork = min( lwmax, int( work( 1 ) ) )
!      lrwork = min( lwmax, int( rwork( 1 ) ) )
!      liwork = min( lwmax, iwork( 1 ) )
      call cheev('Vectors','Upper',N,Ham,lda,w,work,lwork,rwork,info)
!      call cheevd('Vectors','Upper',N,Ham,lda,w,work,lwork,rwork,lrwork
!     & ,iwork,liwork,info)
      if( info .GT. 0 ) then
         write(*,*)'The algorithm failed to compute eigenvalues.'
         stop
      end if
      call err()
      call gap()
!      write(*,*)Ham(30,30)
!  通过调用dggev函数可以求解得到Hamiltonian矩阵的本征值和对应的本征矢量
!  下一步通过本征值和对应本征矢来求解超导配对能delta
!  求解得到delta后再次带入到Hamiltonian矩阵中，然后再进行本征值与本征矢的求解
!  接下的过程就是重复上述过程，知道相邻两次的delta的差小于一个确定的误差
!      call f1(w)
!      call pair(delta,0)
!      do m = 1,48*48
!            write(21,*)abs(delta(m,w,Ham))
!      end do
!      open(17,'count.txt')
!      num = 0                        
!      do while(abs(maxval(abs(del))-delta_init)>eps)
!            num = num + 1
!            write(17,*)num
!            delta_init = maxval(abs(del))
!            call matrix()
!            lwork = -1           
!            call cheev('Vectors','Upper',N,Ham,lda,w,
!     &                   work,lwork,rwork,info)
!            lwork = min( lwmax, int( work( 1 ) ) )
!            call cheev('Vectors','Upper',N,Ham,lda,w,
!     &                   work,lwork,rwork,info)
!            if( info .GT. 0 ) then
!                  write(*,*)'Algorithm failed to compute eigenvalues.'
!            stop
!            end if
!            call f1(w)
!      end do
!--------------------------- open file---------------------------------------------------------
!      open(10,file='H.dat')
      open(11,file='eig-value-single.dat')
!      open(12,file='eig-vector.dat')
!      open(14,file = 'wl.dat')
!      open(15,file="eig-value-R.dat")
!      open(16,file='del.dat')


!      write(10,*)Ham
      do m = 1,N
            write(11,*)w(m) ! eigenvalue
      end do
!      write(15,*)real(w)  ! real_part for eig-value
!      write(12,*)Ham ! eigenvector
!      write(14,*)wl
!      write(16,*)del
!------------------------- close file -----------------------------------------------------------
!      close(10)
      close(11)
!      close(12)
!      close(13)
!      close(15)
!      close(16)
!      close(17)
      stop
      end


!===================================== Function Construct ======================================================================
      complex*8 function phi(y)  ! Phase function   A = (-By,0,0)
      ! dir is charactered the direction of hopping or couple (right or left)
      ! n: right is +1   left is -1
      ! In the y direction hopping and couple isn't contribute a phase change
      use param
      integer y ! y-direction lattice index
      phi = exp(-image*pi/phi0*B*y)      
      end function phi 


!=========================================== Energy Gap function ==============================================================
      complex*8 function delta(x,eig_value,eig_vector)
      ! x,y point position
      ! eig_value eigen-value for matrix
      ! eig_vector eigen vector for matrix of eigen value
!      include 'parameter.f'
      use param
      integer x
      real eig_value(N)
      complex*8 eig_vector(N,N)
      complex*8 s
      real temp
      integer kk
      s = (0,0)
      do kk = 1,len2
            temp = eig_value(kk)/(2*Kb*T)   ! Hermitian Matrix eigvalue is real number
            s = s + eig_vector(x,kk)*
     &      conjg(eig_vector(x+len2*2,kk))*tanh(temp)
      end do
      delta = 5/2.0*s
      end function delta
!********************************************Pariing energy***********************************************
      subroutine pair(input)
      use param
      complex*8,external::delta
      integer input
      
      if (input .eq. 0)then
!----------------------------------------positive energy ---------------------------------------------
            do m = 1,48*48               !   (1,3)block
                  Ham(m,len2*2+m) = delta(m,w,Ham)  ! C-up*C-down
                  Ham(len2*2+m,m) = conjg(delta(m,w,Ham))  ! E-up*E-down
            end do
            
!-----------------------------------------negative energy ---------------------------------------------
            do m = 1,48*48               !   (2,4)block
                  Ham(len2+m,len2*3+m) = -delta(m,w,Ham)  ! C-up*C-down
                  Ham(len2*3+m,len2+m) = -conjg(delta(m,w,Ham))  ! E-up*E-down
            end do
        else
!----------------------------------------positive energy ---------------------------------------------
            do m = 1,48*48               
                  Ham(m,len2*2+m) = delta_init  ! C-up*C-down   (1,3)block
                  Ham(len2*2+m,m) = delta_init  ! E-up*E-down    (3,1) block
            end do
            
!-----------------------------------------negative energy ---------------------------------------------

            do m = 1,48*48        
                  Ham(len2+m,len2*3+m) = -delta_init  ! C-up*C-down (2,4)block
                  Ham(len2*3+m,len2+m) = -delta_init  ! E-up*E-down  (4,2) block
            end do
      end if
      return
      end subroutine            




!******************************************** hopping energy *****************************************
      subroutine hopping()
      use param
      complex*8,external::phi ! function statement
!     --------------------------------------  x-direction right hopping -------------------------------
!     ---------------------- positive energy  -----------------------------------------
!     right hopping is positive       left hopping is negatie 
      do l = 1,48                 !第一列和最后一列在这里没有考虑,为了保证他们向左或者向右跳跃时数组越界 
      do m = 2+(l-1)*48,48*l-1    ! (1,1)block   spin-up
            Ham(m,m+1) = -t0*phi(l)  ! (i--->i+1)  正向跳跃
            Ham(m,m-1) = -t0*conjg(phi(l))  ! (i------>i-1) 相当于考虑了反向跳跃
      end do
      do m = len2+2+(l-1)*48,len2+48*l-1   ! (2,2)block   spin-down
            Ham(m,m+1) = -t0*phi(l)
            Ham(m,m-1) = -t0*conjg(phi(l))
      end do
!     --------------------------negative energy ------------------------
      do m = len2*2+2+(l-1)*48,len2*2+48*l-1    ! (3,3) block
            Ham(m,m+1) = t0*phi(l)
            Ham(m,m-1) = t0*conjg(phi(l))
      end do
      do m = len2*3+2+(l-1)*48,len2*3+48*l-1     ! (4,4) block
            Ham(m,m+1) = t0*phi(l)
            Ham(m,m-1) = t0*conjg(phi(l))
      end do
      end do
!----------------------------------------- y-direction up hooping -------------------------------------
!----------------------------------------- 最后一行向上跳跃并未考虑------------------------------------
! 由朗道规范可以知道，在y方向的跳跃不包含位相信息(除了最后一行)
! 所以在此处没有考虑最后一行在y方向的跳跃情况

      do m = 48+1,48*47   ! (1,1) block
            Ham(m,m+len1) = -t0 ! up(positive) y方向上的向上跳跃
            Ham(m,m-len1) = -t0 ! up(negative) y 方向上的向下跳跃
      end do      
      
      do m = len2+48+1,len2+48*47  ! (2,2)  block
            Ham(m,m+len1) = -t0
            Ham(m,m-len1) = -t0
      end do

      do m = len2*2+48+1,len2*2+48*47    ! (3,3) block
            Ham(m,m+len1) = t0
            Ham(m,m-len1) = t0
      end do

      do m = len2*3+48+1,len2*3+48*47     ! (4,4) block
            Ham(m,m+len1) = t0
            Ham(m,m-len1) = t0
      end do
! ------------ first row down-direction hopping ----------------------------------------------------
! 最后一行向第一行跳跃的过程中，人为的加上一个位相，保证整个波函数跳跃一圈的位相是2pi
!     up hopping is positive               down hopping is negative
      do m = 1,len1   ! (1,1)block
            Ham(m,m+len1) = -t0  !  第一行电子向上跳跃
            Ham(m,48*47+m) = -t0*conjg(phi(m*48))  ! 第一行电子向下跳跃
            Ham(48*47+m,m) = -t0*phi(48*m)  ! 最后一行电子向上跳跃，人为的引入一个位相，该位相与x方向坐标有关
            Ham(48*47+m,48*47+m-len1) = -t0 !最后一行的电子向下进行跳跃
      end do      

      do m = 1,len1   ! (2,2)block
            Ham(len2+m,len2+m+len1) = -t0  !  first row hopping to up direction
            Ham(len2+m,len2+48*47+m) = -t0*conjg(phi(m*48))  ! 
            Ham(len2+48*47+m,len2+m) = -t0*phi(48*m)  ! last row hopping t first row
            Ham(len2+48*47+m,len2+48*47+m-len1) = -t0 !last row hopping to down direction
      end do      
      do m = 1,len1   ! (3,3)block
            Ham(len2*2+m,len2*2+m+len1) = t0  !  first row hopping to up direction
            Ham(len2*2+m,len2*2+48*47+m) = t0*conjg(phi(m*48))  !
            Ham(len2*2+48*47+m,len2*2+m) = t0*phi(48*m)  ! last row hopping t first row
            Ham(len2*2+48*47+m,len2*2+48*47+m-len1) = t0 !last row hopping to down direction 
      end do      
      do m = 1,len1   ! (4,4)block
            Ham(len2*3+m,len2*3+m+len1) = t0  !  first row hopping to up direction
            Ham(len2*3+m,len2*3+48*47+m) = t0*conjg(phi(m*48))  ! 
            Ham(len2*3+48*47+m,len2*3+m) = t0*phi(48*m)  ! last row hopping t first row
            Ham(len2*3+48*47+m,len2*3+48*47+m-len1) = t0 !last row hopping to down direction
      end do      
!------------------------ first coloum left-direction hopping (include last coloum hopping towards first cloum cause of matrix symmetry consider)----
      do m = 1,48   ! (1,1) block
            Ham(m*48,m*48-47) = -t0*phi(m)   ! 最后一行电子向右进行跳跃
            Ham(m*48,m*48-1) = -t0*conjg(phi(m))! 最后一行电子向左跳跃
            Ham(1+(m-1)*48,1+(m-1)*48+1) = -t0*phi(m)! 第一行电子向右跳跃
            Ham(1+(m-1)*48,m*48) = -t0*conjg(phi(m))! 第一行电子向左跳跃
      end do            
      do m = 1,48   ! (2,2) block
            Ham(len2+m*48,len2+m*48-47) = -t0*phi(m)
            Ham(len2+m*48,len2+m*48-1) = -t0*conjg(phi(m))
            Ham(len2+1+(m-1)*48,len2+1+(m-1)*48+1) = -t0*phi(m)! first point hopping to right
            Ham(len2+1+(m-1)*48,len2+m*48) = -t0*conjg(phi(m))! first point hopping to left
      end do            
      do m = 1,48   ! (3,3) block
            Ham(len2*2+m*48,len2*2+m*48-47) = t0*phi(m)
            Ham(len2*2+m*48,len2*2+m*48-1) = t0*conjg(phi(m))
            Ham(len2*2+1+(m-1)*48,len2*2+1+(m-1)*48+1) = t0*phi(m)! first point hopping to right
            Ham(len2*2+1+(m-1)*48,len2*2+m*48) = t0*conjg(phi(m))! first point hopping to left
      end do            
      do m = 1,48   ! (4,4) block
            Ham(len2*3+m*48,len2*3+m*48-47) = t0*phi(m)
            Ham(len2*3+m*48,len2*3+m*48-1) = t0*conjg(phi(m))
            Ham(len2*3+1+(m-1)*48,len2*3+1+(m-1)*48+1) = t0*phi(m)! first point hopping to right
            Ham(len2*3+1+(m-1)*48,len2*3+m*48) = t0*conjg(phi(m))! first point hopping to left
      end do            
!------- 四个顶角之间的相互跳跃构成一个完整的位相--------------------------------------------
            Ham(1,48) = -t0*phi(47)               ! (1,1)--->(1,48)
            Ham(48,48*48) = -t0*phi(48)        ! (1,48)-->(48,48)
            Ham(48*47+48,48) = -t0*conjg(phi(48)) ! (48,48)-->(48,1)           
            Ham(48*48,48*47+1) = -t0*conjg(phi(48)) !(48,1)-->(1,1)     

            Ham(48,1) = -t0*conjg(phi(47))               ! (1,1)--->(1,48)
            Ham(48*48,48) = -t0*conjg(phi(48))        ! (1,48)-->(48,48)
            Ham(48,48*47+48) = -t0*phi(48) ! (48,48)-->(48,1)           
            Ham(48*47+1,48*48) = -t0*phi(48) !(48,1)-->(1,1)     
!================================================================================================
            Ham(len2+1,len2+48) = -t0*phi(47)               ! (1,1)--->(1,48)
            Ham(len2+48,len2+48*48) = -t0*phi(48*47)        ! (1,48)-->(48,48)
            Ham(len2+48*47+48,len2+48) = -t0*conjg(phi(48*47)) ! (48,48)-->(48,1)           
            Ham(len2+48*48,len2+48*47+1) = -t0*conjg(phi(48)) !(48,1)-->(1,1)     

            Ham(len2+48,len2+1) = -t0*conjg(phi(47))               ! (1,1)--->(1,48)
            Ham(len2+48*48,len2+48) = -t0*conjg(phi(48*47))        ! (1,48)-->(48,48)
            Ham(len2+48,len2+48*47+48) = -t0*phi(48*47) ! (48,48)-->(48,1)           
            Ham(len2+48*47+1,len2+48*48) = -t0*phi(48) !(48,1)-->(1,1)     
!===============================================================================================
            Ham(len2*2+1,len2*2+48) = -t0*phi(47)               ! (1,1)--->(1,48)
            Ham(len2*2+48,len2*2+48*48) = -t0*phi(48*47)        ! (1,48)-->(48,48)
            Ham(len2*2+48*47+48,len2*2+48) = -t0*conjg(phi(48*47)) ! (48,48)-->(48,1)           
            Ham(len2*2+48*48,len2*2+48*47+1) = -t0*conjg(phi(48)) !(48,1)-->(1,1)     

            Ham(len2*2+48,len2*2+1) = -t0*conjg(phi(47))               ! (1,1)--->(1,48)
            Ham(len2*2+48*48,len2*2+48) = -t0*conjg(phi(48*47))        ! (1,48)-->(48,48)
            Ham(len2*2+48,len2*2+48*47+48) = -t0*phi(48*47) ! (48,48)-->(48,1)           
            Ham(len2*2+48*47+1,len2*2+48*48) = -t0*phi(48) !(48,1)-->(1,1)     
!==============================================================================================
            Ham(len2*3+1,len2*3+48) = -t0*phi(47)               ! (1,1)--->(1,48)
            Ham(len2*3+48,len2*3+48*48) = -t0*phi(48*47)        ! (1,48)-->(48,48)
            Ham(len2*3+48*47+48,len2*3+48) = -t0*conjg(phi(48*47)) ! (48,48)-->(48,1)           
            Ham(len2*3+48*48,len2*3+48*47+1) = -t0*conjg(phi(48)) !(48,1)-->(1,1)     

            Ham(len2*3+48,len2*3+1) = -t0*conjg(phi(47))               ! (1,1)--->(1,48)
            Ham(len2*3+48*48,len2*3+48) = -t0*conjg(phi(48*47))        ! (1,48)-->(48,48)
            Ham(len2*3+48,len2*3+48*47+48) = -t0*phi(48*47) ! (48,48)-->(48,1)           
            Ham(len2*3+48*47+1,len2*3+48*48) = -t0*phi(48) !(48,1)-->(1,1)     
      return
      end subroutine


!***************************************spin-orbital couple energy*****************************************
      subroutine couple()
      use param
      complex*8,external::phi 
!XXXX                              - x direction copule (right and left couple) ---------------------------
!----------------------------------positive energy couple (1,2)block-------------------------------------------
      do l = 1,48                  ! 用l代表层数，每一层的最右边的位置上向右边的耦合没有考虑
      do m = 2+(l-1)*48,48*l-1     !(1,2)block    
            Ham(m,len2+m+1) = image*lambda*phi(l) ! spin_up-spin_down couple right  正向耦合(1--->2)
            Ham(m,len2+m-1) = -image*lambda*conjg(phi(l))  ! left反向耦合(2--->1)
      end do
! -------------------------------- positive energy copule (2,1)block-------------------------------------
      do m = 2+48*(l-1),48*l-1   !(2,1)block
            Ham(len2+m,m+1) = image*lambda*phi(l) ! spin_up-spin_down couple  right
            Ham(len2+m,m-1) = -image*lambda*conjg(phi(l))  ! left
      end do
! ------------------------------- negative energy couple (3,4) block------------------------------------
      do m = 2+(l-1)*48,48*l-1
            Ham(len2*2+m,len2*3+m+1) = -conjg(image*lambda*phi(l))
            Ham(len2*2+m,len2*3+m-1) = 
     &                               -conjg(-image*lambda*conjg(phi(l)))
      end do
! ------------------------------- negative energy couple (4,3) block------------------------------------
      do m = 2+(l-1)*48,48*l-1
          Ham(len2*3+m,len2*2+m+1) = -conjg(image*lambda*phi(l)) ! spin_up-spin_down couple  right
          Ham(len2*3+m,len2*2+m-1) = -conjg(-image*lambda*conjg(phi(l)))
      end do
      end do
! XXXX--------------------------------------------------------------------------------------------------------------


!YYYY ------------------------------------ y direction copule (up and down couple)加一层---------------------------
! -------------- 不包括最后一行向上耦合的行为 ---------------------------------
! ----------------------  positive energy couple (1,2)block ---------------------------------------------
      do m = 48+1,48*47
            Ham(m,len2+m+48) = lambda
            Ham(m,len2+m-48) = -lambda
      end do
! ----------------------------------------- (2,1) block ------------------------------------------------
      do m = 48+1,48*47
            Ham(len2+m,m+48) = -lambda
            Ham(len2+m,m-48) = lambda
      end do
! ----------------------------------------- (3,4) block -------------------------------------------------
      do m = 48+1,48*47
            Ham(len2*2+m,len2*3+m+48) = lambda
            Ham(len2*2+m,len2*3+m-48) = -lambda
      end do
! -----------------------------------------(4,3) block -----------------------------------------------
      do m = 48+1,48*47
            Ham(len2*3+m,len2*2+m+48) = -lambda
            Ham(len2*3+m,len2*2+m-48) = lambda
      end do

! --------------------------------------- right and left boundary--------------------------------------------------
!---------------------------------------- positive energy -------------------------------------------------
      do m = 1,48   !  (1,2)block
            Ham(m*48,len2+48*m-47) = image*lambda*phi(m)
            Ham(m*48,len2+m*48-1) = -image*lambda*conjg(phi(m)) ! 增加一层会使自旋取向改变一次
            Ham(48*m-47,len2+48*m-47+1) = image*lambda*phi(m)
            Ham(48*m-47,len2+m*48-47+47) = -image*lambda*conjg(phi(m))
      end do
      do m = 1,48   !  (2,1)block
            Ham(len2+m*48,48*m-47) = image*lambda*phi(m)
            Ham(len2+m*48,48*m-1) = -image*lambda*conjg(phi(m)) ! 增加一层会使自旋取向改变一次(1------>48)
            Ham(len2+48*m-47,48*m-47+1) = image*lambda*phi(m)
            Ham(len2+48*m-47,48*m-47+47) = -image*lambda*conjg(phi(m))
      end do            
!---------------------------------------negative energy---------------------------------
      do m = 1,48   !  (3,4)block
            Ham(len2*2+m*48,len2*3+m*48-47) =
     &                        -conjg(image*lambda*phi(m))
            Ham(len2*2+m*48,len2*3+m*48-1) = 
     &                        -conjg(-image*lambda*conjg(phi(m)))! 增加一层会使自旋取向改变一次
            Ham(len2*2+m*48-47,len2*3+m*48-47+1) =
     &                        -conjg(image*lambda*phi(m))
            Ham(len2*2+m*48-47,len2*3+m*48-47+47) = 
     &                        -conjg(-image*lambda*conjg(phi(m)))
      end do 
      do m = 1,48   !  (4,3)block
            Ham(len2*3+m*48,len2*2+48*m-47) =
     &                        -conjg(image*lambda*phi(m))
            Ham(len2*3+m*48,len2*2+48*m-1) = 
     &                        -conjg(-image*lambda*conjg(phi(m))) ! 增加一层会使自旋取向改变一次
            Ham(len2*3+m*48-47,len2*2+48*m-47+1) =
     &                        -conjg(image*lambda*phi(m))
            Ham(len2*3+m*48-47,len2*2+48*m-47+47) = 
     &                        -conjg(-image*lambda*conjg(phi(m)))
      end do           

!------------------------------------- up and down boundary -----------------------------------------------
! ------------------------------------ positive energy ---------------------------------------------------
! ------------------------------------ 向下耦合，第一行会和最后一行进行周期性边界条件的转换----------------
! 
      do m = 1,48       ! (1,2)block
            Ham(m,len2+m+len1) = lambda  ! 第一行向上一行进行耦合
            Ham(m,len2+48*47+m) = -lambda*conjg(phi(m*48)) ! 第一行向下进行耦合(直接耦合到最后一行) 
            Ham(48*47+m,len2+m) = lambda*phi(m*48)! 最后一行向上耦合(直接耦合到第一行) 
            Ham(48*47+m,len2+48*47-len1+m) = -lambda  ! 最后一行电子向下进行耦合
      end do
      do m = 1,48       ! (2,1)block
            Ham(len2+m,m+len1) = -lambda
            Ham(len2+m,48*47+m) = lambda*conjg(phi(m*48))
            Ham(len2+m+48*47,m) = -lambda*phi(m*48)
            Ham(len2+m+48*47,m+48*47-len1) = lambda
      end do
! ------------------------------------ negative energy ----------------------------------------------------------
      do m = 1,48       ! (3,4)block
            Ham(len2*2+m,len2*3+m+len1) =
     &                        lambda
            Ham(len2*2+m,len2*3+48*47+m) =
     &                        -conjg(lambda*conjg(phi(m*48)))
            Ham(len2*2+m+48*47,len2*3+m) =
     &                        -conjg(-lambda*phi(m*48))
            Ham(len2*2+m+48*47,len2*3+m+48*47-len1) =
     &                        -lambda
      end do
      do m = 1,48       ! (4,3)block
            Ham(len2*3+m,len2*2+m+len1) =
     &                        -lambda
            Ham(len2*3+m,len2*2+48*47+m) =
     &                        -conjg(-lambda*conjg(phi(m*48)))
            Ham(len2*3+m+48*47,len2*2+m) =
     &                        -conjg(lambda*phi(m*48))
            Ham(len2*3+m+48*47,len2*2+m+48*47-len1) =
     &                        lambda
      end do
!YYYY--------------------------------------------------------------------------------------------------------
      return
      end subroutine couple

c--------------------------- Diagonal element setting complete-----------
      subroutine diag()
      use param
       ! spin-up positive energy
      do m = 1,48*48   !(1,1)block 
            Ham(m,m) = h-mu
      end do

      ! spin-down positive energy
      do m = len2+1,len2*2    !(2,2)block
            Ham(m,m) = -h-mu
      end do
      ! spin-down negative energy
      do m = len2*2+1,len2*3        !(3,3)block
            Ham(m,m) = -(-h-mu)
      end do
      ! spin-up negative energy
      do m = len2*3+1,len2*4       !(4,4)block
            Ham(m,m) = -(h-mu)
      end do
      return
      end subroutine diag
!--------------------------------------------------------------------------------------------------------------------



      subroutine count_number()
      use param
      integer i,j
      integer c
      c = 0
      do i = 1,len2*4
            do j = 1,len2*4
                  if (abs(Ham(i,j)) .ne. 0)then
                        c = c + 1
                        write(24,*)i,j
                        write(24,*)Ham(i,j)
                        write(24,*)"=========================="
                  end if
            end do
      end do 
      write(*,*)"Non-zero number is:",c 
      return
      end subroutine count_number



      subroutine eigvalue_verify()
      use param
      integer m
      integer::ccc=0
      do m = 1,48*48*2             
            if ((abs(w(m))-abs(w(48*48*4-m))) .ne. 0 )then
                  ccc = m
                  write(*,*)"You are eigenvalue is not symmetric!!!!"
                  write(*,*)"The index of no-symmetric is:",ccc
                  stop
            end if
      end do
      return
      end subroutine eigvalue_verify


      subroutine matrix_verify()
      use param
      integer i,j
      integer ccc
      ccc = 0
      open(16,file = 'test.dat')
      do i = 1,48*48*4
            do j = 1,48*48*4
                  if (Ham(i,j) .ne. conjg(Ham(j,i)))then
                        ccc = ccc +1
!                        write(*,*)"Your matrix is not Hermitian"
                        write(16,*)i,j
                        write(16,*)Ham(i,j)
                        write(16,*)Ham(j,i)
!                        write(16,*)j,i
!                        write(16,*)Ham(i,j)
                        write(16,*)"==================================="
!                        stop
                  end if
            end do
      end do
      write(*,*)ccc
      end subroutine matrix_verify
      


      subroutine f1()
      use param
      integer m     
      complex*8,external::delta
      do m = 1,len2
            del(m)= abs(delta(m,w,Ham))
      end do
      end subroutine f1
      
      subroutine matrix()
      use param
      complex*8,external::delta
      call diag()
      call hopping()
      call couple()
      call pair(0)  ! eig-value and eig-vector have been overwrited,use new data comput pairing energy strength
      end subroutine matrix
      
      subroutine err()
      use param
      integer m
      real d
      open(18,file='eig-err-single.dat')
      do m = 1,len2*4
            d = abs(w(m) + w(len2*4-m+1))
            write(18,*)d
      end do
      close(18)
      return
      end subroutine err


      subroutine gap()
      use param
      open(19,file='gap-single.dat')
      write(19,*)"The superconduct gap is:"
      write(19,*)abs(w(4608)-w(4609))
      close(19)
      return
      end subroutine gap
