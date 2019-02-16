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
      complex*16,parameter::pi=(3.14159265379,0)
      complex*16 h ! h is external magnetic field
      complex*16 V ! V is point potential
      complex*16 a ! Lettice constant
      complex*16 kb ! Boltzmann constant
      complex*16 T ! Temperature
      complex*16 t0 ! t0 is the couple strength of Hubbrd model
      real lambda ! lambda is the spin-orbital coupling strength
      complex*16::image=(0.0,1.0) ! Image number
      real mu ! mu is the chemcial potential
      complex*16 phi0 ! superconduct flux quantum
      complex*16 B ! external Magnetic
      complex*16 Ham(N,N)
      integer,parameter::len2 = 48*48
      integer,parameter::len1 = 48
      real,parameter::eps = 0.0000001 ! measure distance of two delta
      complex*8 del(len2,len2)  ! stor data of pairing energy
      real del_max
      real Hsoc  ! Hso conjugate
      real Hso
      real delc       ! conjugate of pairing energy
      integer counter ! loop number counter
      integer num ! test number
      ! lapack parameter 
      integer::lda = N
      integer::ldvl = 2*N
      integer,parameter::lwmax = 50*N
      integer::ldvr = 2*N
      real wr(N),wi(N)
      complex*16,allocatable::w(:)
      complex*16,allocatable::vl(:,:)
      complex*16,allocatable::vr(:,:)
      complex*16,allocatable::work(:)
      double precision rwork(2*N)
      integer lwork, info
      end module param
!===============================Program Start=========================================
      program bdg
!      include 'parameter.f'
      use param
      integer i,j ! counter used in program and can't used in function
      integer l,m ! Both of them are loop variation
      complex*16,external::phi ! function statement
      real,external::delta ! function statemen
      external zgeev
      allocate(w(N))
      allocate(vl(2*N,N))
      allocate(vr(2*N,N))
      allocate(work(lwmax))

c===============================Parameter Setting=======================
      phi0 = (1.3,0) ! flux quantum
      a = (1,0) ! lattice constant
      h = (0.6,0)
      mu = 1.0
      V = (3.0,0)
      t0 = (1.0,0)
      lambda = 1.0
!      B = 2*phi0/(N*a**2)
      B = 2*phi0
      Kb = (1.0,0)
      T = (1.0,0)
      lambda = 0.2
      Hso = 1.0
      Hsoc = 1.6
c==============================Construct Ham Matrix================
      Ham = (0,0)   ! Initinal Matrix
      wr = 1.0
      vl  = 1.0
c--------------------------- point potential element setting complete-----------
      ! add a point potential at (24,24)
      Ham(24,24)=V ! position at spin-up oriention
      Ham(24+48,24+48)=V ! position at spin-down oriention

      call diag()   ! setting diag to matrix
      call hopping() ! full of elementary about hopping energy
      call couple()  ! couple energy elementary put in matrix 
      call pair(delta,1)     ! 在求解过程中要用自恰的方法，所以需要对超导配对能进行重新赋值
!=================================== eig-value and eig-vector solve ====================================

!      lowrk = 100
!      call sgeev('Vectors','Vectors',N,Ham,lda,wr,wi,vl,ldvl,vr,ldvr
!     & ,work,lwork,info)     ! 在这里可以得到矩阵对角化之后的本征值(w)和本征矢量(work)
!      lwork = min(lwmax,int(work( 1 )))
!      call sgeev( 'Vectors', 'Vectors', N, Ham, lda, wr, wi, vl, ldvl,
!     $            vr, ldvr, work, lwork, info )
!      open(15,file='info.txt')
!      write(15,*)info
!      close(15)       
      lwork = -1
      call zgeev('V','V',N,Ham,lda,w,vl,ldvl,vr,ldvr
     & ,work,lwork,rwork,info)     ! 在这里可以得到矩阵对角化之后的本征值(w)和本征矢量(work)
      lwork = min(lwmax,int(work(1)))
      open(13,file='info.txt')
      write(13,*)'lwork:',lwork
      write(13,*)'lwmax:',lwmax
      write(13,*)'work(1):',work(1)
!      lwork = N
      call zgeev('V','V',N,Ham,lda,w,vl,ldvl,vr,ldvr
     & ,work,lwork,rwork,info)     ! 在这里可以得到矩阵对角化之后的本征值(w)和本征矢量(work)
!      write(*,*)Ham(1,2)
!  通过调用dggev函数可以求解得到Hamiltonian矩阵的本征值和对应的本征矢量
!  下一步通过本征值和对应本征矢来求解超导配对能delta
!  求解得到delta后再次带入到Hamiltonian矩阵中，然后再进行本征值与本征矢的求解
!  接下的过程就是重复上述过程，知道相邻两次的delta的差小于一个确定的误差
!      do m = 1,len2
!            do l = 1,len2
!                  del(m,l) = Ham(m,len2*2+l)   ! take out data of pairing energy   
!            end do
!      end do
!      del_max = maxval(abs(del))
!      del_max = 1 + del_max
!      write(*,*)"wi(10):",wi(10)

!      do while (abs(del_max) > eps)
!            call pair(delta,0)  ! eig-value and eig-vector have been overwrited,use new data comput pairing energy strength
!            write(*,*)"del_max:",del_max
!            lowrk = -1
!            call sgeev('Vectors','Vectors',N,Ham,lda,wr,wi,vl,ldvl,vr,
!     &             ldvr,work,lwork,info)     ! 在这里可以得到矩阵对角化之后的本征值(w)和本征矢量(work)
!            lwork = min(lwmax,int(work( 1 )))
!            call sgeev("vectors",'Vectors',N,Ham,lda,wr,wi,vl,ldvl,vr,
!     &             ldvr,work,lwork,info)
!      end do
      



!      open(10,file='H.dat')
      open(11,file='wr.dat')
!      open(12,file='eig-vector.dat')
!      open(14,file = 'wl.dat')
!------------------------------------output Ham matrix----------------------------
!      do m = 1,48
!            do l = 1,48
!                  write(10,'(48F9.5/)')Ham(m,l)
!            end do
!      end do
      
!                  write(10,'(48F9.5/)')Ham
!------------------------------------- output eigenvalue and eigenvector of Ham matrix-------------------------------
!      write(10,*)Ham 
      write(11,*)w  ! eigenvalue
!      write(12,*)vl ! eigenvector
!      write(14,*)wl
!      close(10)
      close(11)
!      close(12)
!   ------------------------------------       program check           -----------
!      counter = 0
!      do m = 1,N
!            do l = 1,N
!                  if (Ham(l,m) .ne. 0.0) then
!                        counter = counter + 1
!                  end if      
!            end do
!      end do
      
!      do i = 1,48*48
!            write(13,*)delta(i,wr,vl)
!      end do
!      write(13,*)"Non-zero numbers:",counter       
      close(13)

!      do m = 1,10
!            write(*,*)phi(m,m)
!      end do

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
      complex*8 eig_value(N)
      complex*8 eig_vector(2*N,N)
      complex*8 s
      s = 0.0
      do m = 1,len2
!            s = s + eig_vector(x,m)*
!     &      eig_vector(x+len2,m)*tanh(eig_value(m)/(2*Kb*T))   
      end do
      delta = V/2.0*s
      end function delta
!********************************************Pariing energy***********************************************
      subroutine pair(dd,input)
!      include 'parameter.f'
      use param
      real dd
      integer input
      if (input == 0)then
!----------------------------------------positive energy ---------------------------------------------
            do m = 1,48*48               !   (1,3)block
                  Ham(m,len2*2+m) = dd(m,wr,vl)  ! C-up*C-down
                  Ham(len2*2+m,m) = dd(m,wr,vl)  ! E-up*E-down
            end do

            do m = 1,48*48               !   (2,4)block
                  Ham(len2*2+m,len2+len2*2+m) = -dd(m,wr,vl)  ! C-up*C-down
                  Ham(len2+len2*2+m,len2*2+m) = -dd(m,wr,vl)  ! E-up*E-down
            end do
!-----------------------------------------negative energy ---------------------------------------------
            do m = 1,48*48               !   (3,1)block
                  Ham(len2*2+m,m) = dd(m,wr,vl)  ! C-up*C-down
                  Ham(m,len2*2+m) = dd(m,wr,vl)  ! E-up*E-down
            end do

            do m = 1,48*48               !   (4,2)block
                  Ham(len2*3+m,len2+m) = -dd(m,wr,vl)  ! C-up*C-down
                  Ham(len2+m,len2*3+m) = -dd(m,wr,vl)  ! E-up*E-down
            end do
            
            do m = 1,len2
                  do l = 1,len2
                        del(m,l) = Ham(m,len2*2+l)   ! take out data of pairing energy   
                  end do
            end do
            del_max = maxval(abs(del))
      end if
!            write(*,*)maxval(del)
      end subroutine            




!******************************************** hopping energy *****************************************
      subroutine hopping()
      use param
      complex*16,external::phi ! function statement
!     --------------------------------------  x-direction right hopping -------------------------------
!     ---------------------- positive energy  -----------------------------------------
      do l = 1,48
      do m = 1,48*l    ! (1,1)block   spin-up
            Ham(m,m+1) = t0*phi(l,1)  ! (i--->i+1)  正向跳跃
            Ham(m+1,m) = t0*phi(l,-1)  ! (i+1--->i) 相当于考虑了反向跳跃
      end do
      do m = len2+1,len2+48*l   ! (2,2)block   spin-down
            Ham(m,m+1) = t0*phi(l,1)
            Ham(m+1,m) = t0*phi(l,-1)
      end do
!     --------------------------negative energy ------------------------
      do m = len2*2+1,len2*2+48*l    ! (3,3) block
            Ham(m,m+1) = -t0*phi(l,1)
            Ham(m+1,m) = -t0*phi(l,-1)
      end do
      do m = len2*3+1,len2*3+48*l-1     ! (4,4) block
            Ham(m,m+1) = -t0*phi(l,1)
            Ham(m+1,m) = -t0*phi(l,-1)
      end do
      end do
!     -- 上面的循环并未考虑到最后一列的点向右跳跃应该回到最左端，下面进行修改------------------------
      do m = 1,48   !  (1,1) block
            Ham(48*m,48*m+1) = 0
            Ham(48*m,1+48*(m-1)) = t0*phi(m,1)
      end do

      do m = 1,48   !  (2,2) block
            Ham(len2+48*m,len2+48*m+1) = 0
            Ham(len2+48*m,1+len2+48*(m-1)) = t0*phi(m,1)
      end do

      do m = 1,48   !   (3,3) block
            Ham(len2*2+48*m,len2*2+48*m+1) = 0
            Ham(len2*2+48*m,1+len2*2+48*(m-1)) = -t0*phi(l,1)
      end do
 
      do m = 1,47   !   (4,4)  block
            Ham(len2*3+48*m,len2*3+48*m+1) = 0
            Ham(len2*3+48*m,1+len2*3+48*(m-1)) = -t0*phi(l,1)
      end do
!----------------------------------------- y-direction up hooping -------------------------------------
!----------------------------------------- 最后一行向上跳跃并未考虑------------------------------------
      do m = 1,48*47   ! (1,1) block
            Ham(m,len1+m) = t0 ! up(positive) y方向上的向上跳跃
            Ham(len1+m,m) = t0 ! up(negative) y 方向上的向下跳跃
      end do      
      
      do m = len2+1,len2+48*47  ! (2,2)  block
            Ham(m,len1+m) = t0
            Ham(len1+m,m) = t0
      end do

      do m = len2*2+1,len2*2+48*47    ! (3,3) block
            Ham(m,len1+m) = -t0
            Ham(len1+m,m) = -t0
      end do

      do m = len2*3+1,len2*3+48*47     ! (4,4) block
            Ham(m,len1+m) = -t0
            Ham(len1+m,m) = -t0
      end do
! ------------ first row down-direction hopping (考虑对称项后正好包括了最后一行向第一行跳跃)----------------------------------------------------
! ------------- 暂时没有考虑第一行向下跳跃时的位相问题------------------------------------------------------------------------------------------
      do m = 1,len1   ! (1,1)block
            Ham(m,48*47+m) = t0*phi(m,-1)  ! 第一行下最后一行跳跃
            Ham(48*47+m,m) = t0*phi(m,1)  ! 该项为上一项的对称相,也就代表着从最后一行向第一行跳跃的问题
      end do      
      
      do m = len2+1,len2+len1  ! (2,2)  block
            Ham(m,48*47+m) = t0*phi(m,-1)
            Ham(48*47+m,m) = t0*phi(m,1)
      end do

      do m = len2*2+1,len2*2+len1    ! (3,3) block
            Ham(m,48*47+m) = -t0*phi(m,-1)
            Ham(48*47+m,m) = -t0*phi(m,1)
      end do

      do m = len2*3+1,len2*3+len1     ! (4,4) block
            Ham(m,48*47+m) = -t0*phi(m,-1)
            Ham(48*47+m,m) = -t0*phi(m,1)
      end do
!------------------------ first coloum left-direction hopping (include last coloum hopping towards first cloum cause of matrix symmetry consider)----
      do m = 1,48   ! (1,1) block
            Ham(1+(m-1)*48,m*48) = t0*phi(m*48,-1)
            Ham(m*48,1+(m-1)*48) = t0*phi(m*48,1)
      end do            
      do m = 1,48   ! (2,2) block
            Ham(len2+1+len2+(m-1)*48,m*48) = t0*phi(len2+m*48,-1)
            Ham(len2+m*48,1+len2+(m-1)*48) = t0*phi(len2+48*m,1)
      end do            
      do m = 1,48   ! (3,3) block
            Ham(len2*2+1+len2*2+(m-1)*48,m*48) = -t0*phi(len2*2+48*m,-1)
            Ham(len2*2+m*48,len2*2+1+(m-1)*48) = -t0*phi(len2*2+48*m,1)
      end do            
      do m = 1,48   ! (4,4) block
            Ham(1+len2*3+(m-1)*48,m*48) = -t0*phi(len2*3+48*m,-1)
            Ham(len2*3+m*48,1+len2*3+(m-1)*48) = -t0*phi(len2*3+48*m,1)
      end do            


      end subroutine


      subroutine couple()
      use param
      complex*16,external::phi 
!***************************************spin-orbital couple energy*****************************************
!XXXX --------------------------------- x direction copule (right and left couple) ---------------------------
!----------------------------------positive energy couple (1,2)block-------------------------------------------
      do l = 1,48
      do m = 1,48*l             ! 在这种情况下，最右边的点在耦合时应该和最左边的点进行耦合，在这种循环的时候并没有考虑到
            Ham(m,len2+m+1) = Hso*phi(l,1) ! spin_up-spin_down couple right  正向耦合(1--->2)
            Ham(len2+m+1,m) = -Hso*phi(l,-1)  ! left                          反向耦合(2--->1)
!            Ham(len2+m,m+1) = -Hso  ! left
      end do
! -------------------------------- positive energy copule (2,1)block-------------------------------------
      do m = 1,48*l
            Ham(len2+m,m+1) = Hso*phi(l,1) ! spin_up-spin_down couple  right
            Ham(m+1,m+len2) = -Hso*phi(l,-1)  ! left
!            Ham(m,m+len2+1) = -Hso  ! left
      end do
! ------------------------------- negative energy couple (3,4) block------------------------------------
      do m = 1,48*l-1
            Ham(len2*2+m,len2*3+m+1) = -conjg(Hso*phi(l,1)) ! spin_up-spin_down couple  right
            Ham(len2*3+m+1,m+len2*2) = conjg(Hso*phi(l,-1))  ! left
!            Ham(len2*3+m,m+len2*2+1) = Hsoc  ! left
      end do
! ------------------------------- negative energy couple (4,3) block------------------------------------
      do m = 1,48*l
            Ham(len2*3+m,len2*2+m+1) = conjg(Hso*phi(l,1)) ! spin_up-spin_down couple  right
            Ham(len2*2+m+1,m+len2*3) = -conjg(Hso*phi(l,-1))  ! left
!            Ham(len2*2+m,m+len2*3+1) = Hsoc  ! left
      end do
      end do
!---- 以上程序将所有格点向耦合的过程进行了赋值------------------------------------------------------------------


!  考虑到上述情况中并没有对最右边的点利用周期性边界条件
!  下面对这种情况进行修正
!  在修正的过程中同事设定了最右端格点上的耦合情况
      do m = 1,48                      !   (1,2) block   
            Ham(m*48,len2+m*48+1) = 0  !  
            Ham(m*48,len2+1+(m-1)*48) = Hso  ! 修正了错误的边界条件并对其进行正确的赋值(周期性边界条件保证了最右边点的耦合情况有些许不同)
      end do

      do m = 1,48                      !   (2,1)  block     
            
            Ham(len2+m*48,m*48+1) = 0
            Ham(len2+m*48,1+(m-1)*48) = Hso
      end do

      do m = 1,48-1                        !    (3,4) block
            Ham(len2*2+m*48,len2*3+m*48+1) = 0
            Ham(len2*2+m*48,len2*3+1+(m-1)*48) = -Hso
      end do

      do m = 1,48                        !      (4,3) block
            Ham(len2*3+m*48,len2*2+m*48+1) = 0
            Ham(len2*3+m*48,len2*3+1+(m-1)*48) = -Hso
      end do
! XXXX--------------------------------------------------------------------------------------------------------------

!  第一列向左方向的耦合并没有考虑到！！！！！！！！！！


!YYYY ------------------------------------ y direction copule (up and down couple)加一层---------------------------
! -------------- 不包括最后一行向上耦合的行为 ---------------------------------
! ----------------------  positive energy couple (1,2)block ---------------------------------------------
      do m = 1,48*47
            Ham(m,len2+m+48*1) = Hso
            Ham(len2+m+48*1,m) = Hso
      end do
! ----------------------------------------- (2,1) block ------------------------------------------------
      do m = 1,48*47
            Ham(len2+m,m+48*1) = Hso
            Ham(m+48*1,len2+m) = Hso
      end do
! ----------------------------------------- (3,4) block -------------------------------------------------
      do m = 1,48*47
            Ham(len2*2+m,len2*2+len2+m+48*1) = -Hso
            Ham(len2*2+len2+m+48*1,len2*2+m) = -Hso
      end do
! -----------------------------------------(4,3) block -----------------------------------------------
      do m = 1,48*47
            Ham(len2*3+m,len2+len2+m) = -Hso
            Ham(len2+len2+m,len2*3+m) = -Hso
      end do

! --------------------------------------- right and left boundary--------------------------------------------------
!---------------------------------------- positive energy -------------------------------------------------
      do m = 1,48   !  (1,2)block
            Ham(1+(m-1)*48,len2+m*48) = Hso ! 增加一层会使自旋取向改变一次
            Ham(len2+m*48,1+(m-1)*48) = Hso
      end do      


      do m = 1,48   !  (2,1)block
            Ham(len2+1+(m-1)*48,m*48) = Hso ! 增加一层会使自旋取向改变一次
            Ham(m*48,len2+1+(m-1)*48) = Hso
      end do      
!---------------------------------------negative energy---------------------------------
      do m = 1,48   !  (3,4)block
            Ham(len2*2+1+(m-1)*48,len2*2+len2+m*48) = -Hsoc ! 增加一层会使自旋取向改变一次
            Ham(len2*2+len2+m*48,len2*2+1+(m-1)*48) = -Hsoc
      end do      

      do m = 1,48   !  (4,3)block
            Ham(len2*3+1+(m-1)*48,len2+len2+m*48) = -Hsoc ! 增加一层会使自旋取向改变一次
            Ham(len2+len2+m*48,len2*3+1+(m-1)*48) = -Hsoc
      end do      




!------------------------------------- up and down boundary -----------------------------------------------
! ------------------------------------ positive energy ---------------------------------------------------
! ------------------------------------ 向下耦合，第一行会和最后一行进行周期性边界条件的转换----------------
      do m = 1,48       ! (1,2)block
            Ham(m,len2+48*47+m) = Hso
            Ham(len2+48*47+m,m) = Hso
      end do
      do m = 1,48       ! (2,1)block
            Ham(len2+m,48*47+m) = Hso
            Ham(48*47+m,len2+m) = Hso
      end do
! ------------------------------------ negative energy ----------------------------------------------------------
      do m = 1,48       ! (3,4)block
            Ham(len2*2+m,len2*2+len2+48*47+m) = Hso
            Ham(len2*2+len2+48*47+m,len2*2+m) = Hso
      end do
      do m = 1,48       ! (4,3)block
            Ham(len2*3+m,len2*2+len2+48*47+m) = Hso
            Ham(len2*2+len2+48*47+m,len2*3+m) = Hso
      end do
!YYYY--------------------------------------------------------------------------------------------------------
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
      end subroutine diag

