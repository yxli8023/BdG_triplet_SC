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
      real,parameter::pi = 3.14159265379
      real h ! h is external magnetic field
      real V ! V is point potential
      real a ! Lettice constant
      real kb ! Boltzmann constant
      real T ! Temperature
      real t0 ! t0 is the couple strength of Hubbrd model
      real lambda ! lambda is the spin-orbital coupling strength
      complex*16::image=(0.0,1.0) ! Image number
      real mu ! mu is the chemcial potential
      real phi0 ! superconduct flux quantum
      real B ! external Magnetic
      integer,parameter::len2 = 48*48
      integer,parameter::len1 = 48
      real,parameter::eps = 0.0000001 ! measure distance of two delta
      real del(len2)  
      complex::Ham(N,N)=(0,0)
      real::del_init = 1        ! The initial value of delta and used in self-consistently loop
      real del_max
      real Hso
      integer counter ! loop number counter
      integer num ! test number
      end module param
!===============================Program Start=========================================
      program bdg
!      include 'parameter.f'
      use param
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
      

      integer i,j ! counter used in program and can't used in function
      integer l,m ! Both of them are loop variation
      complex*16,external::phi ! function statement
      complex*8,external::delta ! function statemen
!      external cheevd
      allocate(w(N))
      allocate(work(lwmax))
      allocate(rwork(lwmax))
      allocate(iwork(lwmax))

c===============================Parameter Setting=======================
      phi0 = 2*pi! flux quantum
      a = 1 ! lattice constant
      h = 0.6
      mu = -4
      V = 10.0
      t0 = 1.0
      lambda = 1.0
      B = 2*phi0/(48*48)
      Kb = 1.0
      T = 1.0
      Hso = 1.0
c==============================Construct Ham Matrix================
!      wr = 1.0
!      vl  = 1.0
c--------------------------- point potential element setting complete-----------
      ! add a point potential at (24,24)
      Ham(24,24)=V ! position at spin-up oriention
      Ham(24+48,24+48)=V ! position at spin-down oriention

      call diag()   ! setting diag to matrix
      call hopping() ! full of elementary about hopping energy
      call couple()  ! couple energy elementary put in matrix 
      call pair(delta,1)     ! 在求解过程中要用自恰的方法，所以需要对超导配对能进行重新赋值
!      call count_number()    ! 统计哈密顿矩阵中非零元素个数
!=================================== eig-value and eig-vector solve ====================================
!      call matrix_verify()   ! 检验矩阵是否为厄米矩阵
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
!      do while(abs(maxval(abs(del))-del_init)>eps)
!            num = num + 1
!            write(17,*)num
!            del_init = maxval(abs(del))
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
      open(11,file='eig-value.dat')
!      open(12,file='eig-vector.dat')
!      open(14,file = 'wl.dat')
!      open(15,file="eig-value-R.dat")
!      open(16,file='del.dat')


!      write(10,*)Ham
      write(11,*)w  ! eigenvalue
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
      complex*16 function phi(y,dir)  ! Phase function   A = (-By,0,0)
      ! dir is charactered the direction of hopping or couple (right or left)
      ! n: right is +1   left is -1
      ! In the y direction hopping and couple isn't contribute a phase change
!      include 'parameter.f'
      use param
      integer y,dir ! position of lattice
      phi = exp(image*pi/phi0*B*y*(-1.0)**dir)
!      phi = cexp(image*pi/phi0*B*y*(-1.0)**dir)
!      phi = exp(image*pi/phi0)        
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
     &      conjg(eig_vector(x+len2*2,kk))*((exp(temp)-exp(-temp))/
     &       (exp(temp)+exp(-temp)))
      end do
      delta = 1/2.0*s
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
                  Ham(len2*2+m,m) = delta(m,w,Ham)  ! E-up*E-down
            end do

            do m = 1,48*48               !   (2,4)block
                  Ham(len2*2+m,len2+len2*2+m) = -delta(m,w,Ham)  ! C-up*C-down
                  Ham(len2+len2*2+m,len2*2+m) = -delta(m,w,Ham)  ! E-up*E-down
            end do
!-----------------------------------------negative energy ---------------------------------------------
            do m = 1,48*48               !   (3,1)block
                  Ham(len2*2+m,m) = conjg(delta(m,w,Ham))  ! C-up*C-down
                  Ham(m,len2*2+m) = conjg(delta(m,w,Ham))  ! E-up*E-down
            end do

            do m = 1,48*48               !   (4,2)block
                  Ham(len2*3+m,len2+m) = -conjg(delta(m,w,Ham))  ! C-up*C-down
                  Ham(len2+m,len2*3+m) = -conjg(delta(m,w,Ham))  ! E-up*E-down
            end do
        else
!----------------------------------------positive energy ---------------------------------------------
            do m = 1,48*48               !   (1,3)block
                  Ham(m,len2*2+m) = del_init  ! C-up*C-down
                  Ham(len2*2+m,m) = del_init  ! E-up*E-down
            end do

            do m = 1,48*48               !   (2,4)block
                  Ham(len2*2+m,len2+len2*2+m) = -del_init  ! C-up*C-down
                  Ham(len2+len2*2+m,len2*2+m) = -del_init  ! E-up*E-down
            end do
!-----------------------------------------negative energy ---------------------------------------------
            do m = 1,48*48               !   (3,1)block
                  Ham(len2*2+m,m) = del_init  ! C-up*C-down
                  Ham(m,len2*2+m) = del_init  ! E-up*E-down
            end do

            do m = 1,48*48               !   (4,2)block
                  Ham(len2*3+m,len2+m) = -del_init  ! C-up*C-down
                  Ham(len2+m,len2*3+m) = -del_init  ! E-up*E-down
            end do
      end if
      return
      end subroutine            




!******************************************** hopping energy *****************************************
      subroutine hopping()
      use param
      complex*16,external::phi ! function statement
!     --------------------------------------  x-direction right hopping -------------------------------
!     ---------------------- positive energy  -----------------------------------------
      do l = 1,48
      do m = 1+(l-1)*48,48*l-1    ! (1,1)block   spin-up
            Ham(m,m+1) = -t0*phi(l,2)  ! (i--->i+1)  正向跳跃
            Ham(m+1,m) = -t0*phi(l,1)  ! (i+1--->i) 相当于考虑了反向跳跃
      end do
      do m = len2+1+(l-1)*48,len2+48*l-1   ! (2,2)block   spin-down
            Ham(m,m+1) = -t0*phi(l,2)
            Ham(m+1,m) = -t0*phi(l,1)
      end do
!     --------------------------negative energy ------------------------
      do m = len2*2+1+(l-1)*48,len2*2+48*l-1    ! (3,3) block
            Ham(m,m+1) = t0*phi(l,2)
            Ham(m+1,m) = t0*phi(l,1)
      end do
      do m = len2*3+1+(l-1)*48,len2*3+48*l-1     ! (4,4) block
            Ham(m,m+1) = t0*phi(l,2)
            Ham(m+1,m) = t0*phi(l,1)
      end do
      end do
!     -- 上面的循环并未考虑到最后一列的点向右跳跃应该回到最左端，下面进行修改------------------------
      do m = 1,48   !  (1,1) block
            Ham(48*m,48*m+1) = 0
            Ham(48*m,1+48*(m-1)) = -t0*phi(m,2)
      end do

      do m = 1,48   !  (2,2) block
            Ham(len2+48*m,len2+48*m+1) = 0
            Ham(len2+48*m,1+len2+48*(m-1)) = -t0*phi(m,2)
      end do

      do m = 1,48   !   (3,3) block
            Ham(len2*2+48*m,len2*2+48*m+1) = 0
            Ham(len2*2+48*m,1+len2*2+48*(m-1)) = t0*phi(m,2)
      end do
 
      do m = 1,47   !   (4,4)  block
            Ham(len2*3+48*m,len2*3+48*m+1) = 0
            Ham(len2*3+48*m,1+len2*3+48*(m-1)) = t0*phi(m,2)
      end do
!----------------------------------------- y-direction up hooping -------------------------------------
!----------------------------------------- 最后一行向上跳跃并未考虑------------------------------------
      do m = 1,48*47   ! (1,1) block
            Ham(m,len1+m) = -t0 ! up(positive) y方向上的向上跳跃
            Ham(len1+m,m) = -t0 ! up(negative) y 方向上的向下跳跃
      end do      
      
      do m = len2+1,len2+48*47  ! (2,2)  block
            Ham(m,len1+m) = -t0
            Ham(len1+m,m) = -t0
      end do

      do m = len2*2+1,len2*2+48*47    ! (3,3) block
            Ham(m,len1+m) = t0
            Ham(len1+m,m) = t0
      end do

      do m = len2*3+1,len2*3+48*47     ! (4,4) block
            Ham(m,len1+m) = t0
            Ham(len1+m,m) = t0
      end do
! ------------ first row down-direction hopping (考虑对称项后正好包括了最后一行向第一行跳跃)----------------------------------------------------
! ------------- 暂时没有考虑第一行向下跳跃时的位相问题------------------------------------------------------------------------------------------
      do m = 1,len1   ! (1,1)block
            Ham(m,48*47+m) = -t0*phi(m,1)  ! 第一行下最后一行跳跃
            Ham(48*47+m,m) = -t0*phi(m,2)  ! 该项为上一项的对称相,也就代表着从最后一行向第一行跳跃的问题
      end do      
      
      do m = len2+1,len2+len1  ! (2,2)  block
            Ham(m,48*47+m) = -t0*phi(m,1)
            Ham(48*47+m,m) = -t0*phi(m,2)
      end do

      do m = len2*2+1,len2*2+len1    ! (3,3) block
            Ham(m,48*47+m) = t0*phi(m,1)
            Ham(48*47+m,m) = t0*phi(m,2)
      end do

      do m = len2*3+1,len2*3+len1     ! (4,4) block
            Ham(m,48*47+m) = t0*phi(m,1)
            Ham(48*47+m,m) = t0*phi(m,2)
      end do
!------------------------ first coloum left-direction hopping (include last coloum hopping towards first cloum cause of matrix symmetry consider)----
      do m = 1,48   ! (1,1) block
            Ham(1+(m-1)*48,m*48) = -t0*phi(m,1)
            Ham(m*48,1+(m-1)*48) = -t0*phi(m,2)
      end do            
      do m = 1,48   ! (2,2) block
            Ham(len2+1+(m-1)*48,len2+m*48) = -t0*phi(m,1)
            Ham(len2+m*48,1+len2+(m-1)*48) = -t0*phi(m,2)
      end do            
      do m = 1,48   ! (3,3) block
            Ham(len2*2+1+(m-1)*48,len2*2+m*48) = t0*phi(m,1)
            Ham(len2*2+m*48,len2*2+1+(m-1)*48) = t0*phi(m,2)
      end do            
      do m = 1,48   ! (4,4) block
            Ham(1+len2*3+(m-1)*48,len2*3+m*48) = t0*phi(m,1)
            Ham(len2*3+m*48,1+len2*3+(m-1)*48) = t0*phi(m,2)
      end do            

      return
      end subroutine


      subroutine couple()
      use param
      complex*16,external::phi 
!***************************************spin-orbital couple energy*****************************************
!XXXX --------------------------------- x direction copule (right and left couple) ---------------------------
!----------------------------------positive energy couple (1,2)block-------------------------------------------
      do l = 1,48
      do m = 1+(l-1)*48,48*l-1             ! 在这种情况下，最右边的点在耦合时应该和最左边的点进行耦合，在这种循环的时候并没有考虑到
            Ham(m,len2+m+1) = image*lambda*Hso*phi(l,2) ! spin_up-spin_down couple right  正向耦合(1--->2)
            Ham(len2+m+1,m) = -image*lambda*Hso*phi(l,1)  ! left                          反向耦合(2--->1)
!            Ham(len2+m,m+1) = -Hso  ! left
      end do
! -------------------------------- positive energy copule (2,1)block-------------------------------------
      do m = 1+48*(l-1),48*l-1
            Ham(len2+m,m+1) = image*lambda*Hso*phi(l,2) ! spin_up-spin_down couple  right
            Ham(m+1,m+len2) = -image*lambda*Hso*phi(l,1)  ! left
!            Ham(m,m+len2+1) = -Hso  ! left
      end do
! ------------------------------- negative energy couple (3,4) block------------------------------------
      do m = 1+(l-1)*48,48*l-1
            Ham(len2*2+m,len2*3+m+1) = -conjg(image*lambda*Hso*phi(l,2)) ! spin_up-spin_down couple  right
            Ham(len2*3+m+1,m+len2*2) = conjg(image*lambda*Hso*phi(l,1))  ! left
!            Ham(len2*3+m,m+len2*2+1) = Hsoc  ! left
      end do
! ------------------------------- negative energy couple (4,3) block------------------------------------
      do m = 1+(l-1)*48,48*l-1
          Ham(len2*3+m,len2*2+m+1) = conjg(image*lambda*Hso*phi(l,2)) ! spin_up-spin_down couple  right
          Ham(len2*2+m+1,m+len2*3) = -conjg(image*lambda*Hso*phi(l,1))  ! left
!            Ham(len2*2+m,m+len2*3+1) = Hsoc  ! left
      end do
      end do
!---- 以上程序将所有格点向耦合的过程进行了赋值------------------------------------------------------------------

      
! XXXX--------------------------------------------------------------------------------------------------------------

!  第一列向左方向的耦合并没有考虑到！！！！！！！！！！


!YYYY ------------------------------------ y direction copule (up and down couple)加一层---------------------------
! -------------- 不包括最后一行向上耦合的行为 ---------------------------------
! ----------------------  positive energy couple (1,2)block ---------------------------------------------
      do m = 1,48*47
            Ham(m,len2+m+48) = lambda*Hso
            Ham(len2+m+48,m) = lambda*Hso
      end do
! ----------------------------------------- (2,1) block ------------------------------------------------
      do m = 1,48*47
            Ham(len2+m,m+48) = -lambda*Hso
            Ham(m+48,len2+m) = -lambda*Hso
      end do
! ----------------------------------------- (3,4) block -------------------------------------------------
      do m = 1,48*47
            Ham(len2*2+m,len2*3+m+48) = lambda*Hso
            Ham(len2*3+m+48,len2*2+m) = lambda*Hso
      end do
! -----------------------------------------(4,3) block -----------------------------------------------
      do m = 1,48*47
            Ham(len2*3+m,len2*2+m+48) = -lambda*Hso
            Ham(len2*2+48+m,len2*3+m) = -lambda*Hso
      end do

! --------------------------------------- right and left boundary--------------------------------------------------
!---------------------------------------- positive energy -------------------------------------------------
      do m = 1,48   !  (1,2)block
            Ham(1+(m-1)*48,len2+m*48) = image*lambda*Hso*phi(m,2) ! 增加一层会使自旋取向改变一次
            Ham(len2+m*48,1+(m-1)*48) = -image*lambda*Hso*phi(m,1)
      end do      


      do m = 1,48   !  (2,1)block
            Ham(len2+1+(m-1)*48,m*48) = image*lambda*Hso*phi(m,2) ! 增加一层会使自旋取向改变一次
            Ham(m*48,len2+1+(m-1)*48) = -image*lambda*Hso*phi(m,1)
      end do      
!---------------------------------------negative energy---------------------------------
      do m = 1,48   !  (3,4)block
            Ham(len2*2+1+(m-1)*48,len2*3+m*48) = 
     &                        conjg(image*lambda*Hso*phi(m,2)) ! 增加一层会使自旋取向改变一次
            Ham(len2*3+m*48,len2*2+1+(m-1)*48) =
     &                        -conjg(image*lambda*Hso*phi(m,1))
      end do      

      do m = 1,48   !  (4,3)block
            Ham(len2*3+1+(m-1)*48,len2*2+m*48) = 
     &                        conjg(image*lambda*Hso*phi(m,2)) ! 增加一层会使自旋取向改变一次
            Ham(len2*2+m*48,len2*3+1+(m-1)*48) =
     &                        -conjg(image*lambda*Hso*phi(m,1))
      end do      




!------------------------------------- up and down boundary -----------------------------------------------
! ------------------------------------ positive energy ---------------------------------------------------
! ------------------------------------ 向下耦合，第一行会和最后一行进行周期性边界条件的转换----------------
      do m = 1,48       ! (1,2)block
            Ham(m,len2+48*47+m) = lambda*Hso*phi(m,2)
            Ham(len2+48*47+m,m) = lambda*Hso*phi(m,1)
      end do
      do m = 1,48       ! (2,1)block
            Ham(len2+m,48*47+m) = lambda*Hso*phi(m,2)
            Ham(48*47+m,len2+m) = lambda*Hso*phi(m,1)
      end do
! ------------------------------------ negative energy ----------------------------------------------------------
      do m = 1,48       ! (3,4)block
            Ham(len2*2+m,len2*2+len2+48*47+m) = 
     &       conjg(lambda*Hso*phi(m,2))
            Ham(len2*2+len2+48*47+m,len2*2+m) = 
     &       conjg(lambda*Hso*phi(m,1))
      end do
      do m = 1,48       ! (4,3)block
            Ham(len2*3+m,len2*2+len2+48*47+m) = 
     &       conjg(lambda*Hso*phi(m,2))
            Ham(len2*2+len2+48*47+m,len2*3+m) =
     &       conjg(lambda*Hso*phi(m,1))
      end do
!YYYY--------------------------------------------------------------------------------------------------------
      return
      end subroutine couple

c--------------------------- Diagonal element setting complete-----------
      subroutine diag()
      use param
            ! spin-up positive energy
            do m = 1,48*48    
                  Ham(m,m) = h-mu
            end do

            ! spin-down positive energy
            do m = len2+1,len2*2    
                  Ham(m,m) = -h-mu
            end do
            ! spin-up negative energy
            do m = len2*2+1,len2*3        
                  Ham(m,m) = -(h-mu)
            end do

            ! spin-down negative energy
            do m = len2*3+1,len2*4
                  Ham(m,m) = -(-h-mu)
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
                  end if
            end do
      end do 
      write(*,*)"Non-zero number is:",c 
      return
      end subroutine count_number



      subroutine test(ff,w,vl)
      integer m
      complex*8,external::ff
      do m = 1,48*48
            write(*,*) ff(m,w,vl)
      end do
      return
      end subroutine test


      subroutine matrix_verify()
      use param
      integer i,j
      integer ccc
      ccc = 0
      open(16,file = 'test.txt')
      do i = 1,48*48*4
            do j = 1,48*48*4
                  if (Ham(i,j) .ne. conjg(Ham(j,i)))then
                        ccc = ccc +1
!                        write(*,*)"Your matrix is not Hermitian"
                        write(16,*)i,j
                        write(16,*)Ham(i,j)
!                        write(16,*)j,i
!                        write(16,*)Ham(i,j)
                        write(16,*)"==================================="
!                        stop
                  end if
            end do
      end do
      write(*,*)ccc
      end subroutine matrix_verify
      


      subroutine f1(w)
      use param
      integer m     
      real w(N)
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
      call pair(delta,0)  ! eig-value and eig-vector have been overwrited,use new data comput pairing energy strength
      end subroutine matrix
