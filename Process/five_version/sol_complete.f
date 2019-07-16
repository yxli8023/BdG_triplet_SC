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
      end module param
      

      program bdg
      ! Author :YuXuanLi
      ! E-mail:yxli406@gmail.com
      ! Article:Phys.Rev.B 96,184508
      use param
c===============================Parameter Setting=======================
      h = 0.6
      mu = -4
      V = 5
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
c=============================== start down boundary ============================================================                        
c                 hopping energy
            ! left-right hopping
            Hamiltonian(1,1+1) = Ht  ! right hopping
            Hamiltonian(1,48) = Ht  ! left hopping
            ! up-down hopping
            Hamiltonian(1,48*1+1) = Ht
            Hamiltonian(1,48*47+1) = Ht
c                 spin-orbital couple
            ! left-right couple
            Hamiltonian(1,len2+2) = Hso ! right
            Hamiltonian(1,len2) = Hso ! left
            ! up-down couple
            Hamiltonian(1,len2+48*1+1) = Hso !up
            Hamiltonian(1,len2+48*47+1) = Hso ! down !cause of periodic first row and last row are nearest-neighbor
c                superconductor pairing energy
            Hamiltonian(1,len2*2+1) = del ! conjg is complex conjugate function

      do m=2,47
c                 hopping energy
            ! left-right hopping
            Hamiltonian(m,m+1) = Ht  ! right hopping
            Hamiltonian(m,m-1) = Ht  ! left hopping
            ! up-down hopping
            Hamiltonian(m,48*1+m) = Ht
            Hamiltonian(m,48*47+m) = Ht
c                 spin-orbital couple
            ! left-right couple
            Hamiltonian(m,len2+m+1) = Hso ! right
            Hamiltonian(m,len2+m) = Hso ! left
            ! up-down couple
            Hamiltonian(1,len2+48*1+m) = Hso !up
            Hamiltonian(1,len2+48*47+m) = Hso ! down !cause of periodic first row and last row are nearest-neighbor
c                superconductor pairing energy
            Hamiltonian(1,len2*2+m) = del ! conjg is complex conjugate function
      end do
c                 hopping energy
            ! left-right hopping
            Hamiltonian(48,1) = Ht  ! right hopping
            Hamiltonian(48,47) = Ht  ! left hopping
            ! up-down hopping
            Hamiltonian(48,48*1+48) = Ht
            Hamiltonian(48,48*47+48) = Ht
c                 spin-orbital couple
            ! left-right couple
            Hamiltonian(48,len2+1) = Hso ! right
            Hamiltonian(48,len2*2-1) = Hso ! left
            ! up-down couple
            Hamiltonian(48,len2+1+48*1) = Hso !up
            Hamiltonian(48,len2*2-1+48*47) = Hso ! down !cause of periodic first row and last row are nearest-neighbor
c                superconductor pairing energy
            Hamiltonian(48,len2*2+1) = del ! conjg is complex conjugate function
c=============================  down boundary complete=============================================================



c============================  start left boundary=================================================================
      do m = 2,47
c                 hopping energy
            ! left-right hopping
            Hamiltonian(1+(m-1)*48,(m-1)*48+1) = Ht  ! right hopping
            Hamiltonian(1+(m-1)*48,48+(m-1)*48) = Ht  ! left hopping
            ! up-down hopping
            Hamiltonian(1+(m-1)*48,48*1+1+(m-1)*48) = Ht
            Hamiltonian(1+(m-1)*48,48*47+1+(m-1)*48) = Ht
c                 spin-orbital couple
            ! left-right couple
            Hamiltonian(1+(m-1)*48,len2+2+(m-1)*48) = Hso ! right
            Hamiltonian(1+(m-1)*48,len2+(m-1)*48) = Hso ! left
            ! up-down couple
            Hamiltonian(1+(m-1)*48,len2+48*1+1+(m-1)*48) = Hso !up
            Hamiltonian(1+(m-1)*48,len2+48*47+1+(m-1)*48) = Hso ! down !cause of periodic first row and last row are nearest-neighbor
c                superconductor pairing energy
            Hamiltonian(1+(m-1)*48,len2*2+1+(m-1)*48) = del ! conjg is complex conjugate function
      end do

c                 hopping energy
            ! left-right hopping
            Hamiltonian(48*47+1,48*47+2) = Ht  ! right hopping
            Hamiltonian(48*47+1,48*48) = Ht  ! left hopping
            ! up-down hopping
            Hamiltonian(48*47+1,1) = Ht ! up hopping
            Hamiltonian(48*47+1,48*47+1-48) = Ht ! down hopping 
c                 spin-orbital couple
            ! left-right couple
            Hamiltonian(48*47+1,len2+48*47+1) = Hso ! right
            Hamiltonian(48*47+1,len2+48*48) = Hso ! left
            ! up-down couple
            Hamiltonian(48*47+1,len2+48*47+1) = Hso !up
            Hamiltonian(48*47+1,len2-48*47+1) = Hso ! down !cause of periodic first row and last row are nearest-neighbor
c                superconductor pairing energy
            Hamiltonian(48*47+1,len2*2+48*47+1) = del ! conjg is complex conjugate function
c==========================  left boundary complete ================================================================
                        



c=========================   up boundary start ==================================================================
      do m = 48*47+2,48*48-1
c                 hopping energy
            ! left-right hopping
            Hamiltonian(m,m+1) = Ht  ! right hopping
            Hamiltonian(m,m-1) = Ht  ! left hopping
            ! up-down hopping
            Hamiltonian(m,m-48*47) = Ht ! up hopping
            Hamiltonian(m,m-48*1) = Ht ! down hopping
c                 spin-orbital couple
            ! left-right couple
            Hamiltonian(m,len2+m+1) = Hso ! right
            Hamiltonian(m,len2+m-1) = Hso ! left
            ! up-down couple
            Hamiltonian(m,len2-48*47+m) = Hso !up
            Hamiltonian(m,len2-48*1+m) = Hso ! down !cause of periodic first row and last row are nearest-neighbor
c                superconductor pairing energy
            Hamiltonian(m,len2*2+m) = del ! conjg is complex conjugate function
      end do

c                 hopping energy
            ! left-right hopping
            Hamiltonian(48*48,48*47+1) = Ht  ! right hopping
            Hamiltonian(48*48,48*48-1) = Ht  ! left hopping
            ! up-down hopping
            Hamiltonian(48*48,len2-48*47) = Ht ! up hopping
            Hamiltonian(48*48,len2-48*1) = Ht ! down hopping 
c                 spin-orbital couple
            ! left-right couple
            Hamiltonian(48*48,len2+48*47+1) = Hso ! right
            Hamiltonian(48*48,len2+48*48-1) = Hso ! left
            ! up-down couple
            Hamiltonian(48*48,len2-48*47) = Hso !up
            Hamiltonian(48*48,len2-48*1) = Hso ! down !cause of periodic first row and last row are nearest-neighbor
c                superconductor pairing energy
            Hamiltonian(48*48,len2*2+48*48) = del ! conjg is complex conjugate function

c=========================== up boundary complete =========================================================


c=========================== right boundry start ===========================================================

      do m = 2,47
c                 hopping energy
            ! left-right hopping
            Hamiltonian(48*m,(m-1)*48+1) = Ht  ! right hopping
            Hamiltonian(48*m,m*48-1) = Ht  ! left hopping
            ! up-down hopping
            Hamiltonian(48*m,48*(m+1)) = Ht
            Hamiltonian(48*m,48*(m-1)) = Ht
c                 spin-orbital couple    //  While hopping to the next block
            ! left-right couple
            Hamiltonian(48*m,len2+(m-2)*48) = Hso ! right
            Hamiltonian(48*m,len2+48*m-1) = Hso ! left
            ! up-down couple
            Hamiltonian(48*m,len2+48*(m-1)) = Hso !up
            Hamiltonian(48*m,len2-(m-1)*48) = Hso ! down !cause of periodic first row and last row are nearest-neighbor
c                superconductor pairing energy
            Hamiltonian(48*47+1,len2*2+48*m) = del ! conjg is complex conjugate function
      end do

c========================== right boundary complete  =======================================================


      
c==========================  up boundary complete  =========================================================
! According periodic condiction only boundary point have something trouble,we can handle point without boundary
! In the eare without boundary hopping and couple condition is identical
c================================Don't consider boundary points first=========================
! spin-up eliminate operator
c       Hopping energy
      ! left-right hop energy (x-direction)
      do m=48+2,48*47-1
            Hamiltonian(m,m-1) = Ht
            Hamiltonian(m,m+1) = Ht 
c     ! up-down hop energy(y-driection)
            Hamiltonian(m,48*1+m) = Ht ! up layer have same x coordinate
            Hamiltonian(m,m-48*1) = Ht ! down layer  last: 48+48*47-1 
c       Spin-orbital couple energy
c      ! left-right couple energy
            Hamiltonian(m,len2+m-1) = Hso ! next block
            Hamiltonian(m,len2+m+1) = Hso ! next block
c      ! up-down couple energy
            Hamiltonian(m,len2-48*1+m) = Hso
            Hamiltonian(m,len2+48*1+m) = Hso
c         superconductor piring energy
            Hamiltonian(m,len2*2+m) = del ! last:48*48*2+48*47-1
      end do
c=========================================================================
! spin-down elinate operator
c       Hopping energy
      ! left-right hop energy (x-direction)
      do m=48*48+48+2,2*48*48-48-1 ! second block // decrease 1 is ensure avoid boundary point first
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
      end do
c============================================================================
! spin-down create operator
c       Hopping energy
      ! left-right hop energy (x-direction)
      do m=2*48*48+48+2,3*48*48-48-1
            Hamiltonian(m,m-1) = Ht
            Hamiltonian(m,m+1) = Ht 
c     ! up-down hop energy(y-driection)
            Hamiltonian(m,48*1+m) = Ht ! up layer
            Hamiltonian(m,m-48*1) = Ht ! down layer
c       Spin-orbital couple energy
c      ! left-right couple energy
            Hamiltonian(m,m+len2+1) = Hso
            Hamiltonian(m,m+len2-1) = Hso
c      ! up-down couple energy
            Hamiltonian(m,m+len2-48*1) = Hso
            Hamiltonian(m,m+len2+48*1) = Hso
c         superconductor piring energy
            Hamiltonian(m,m-2*len2) = del ! In do loop has increase 1 layer
      end do
c=============================================================================
! spin-up create operator
c       Hopping energy
      ! left-right hop energy (x-direction)
      do m=3*48*48+48+2,4*48*48-48-1
            Hamiltonian(m,m-1) = Ht
            Hamiltonian(m,m+1) = Ht 
c     ! up-down hop energy(y-driection)
            Hamiltonian(m,48*1+m) = Ht ! up layer
            Hamiltonian(m,m-48*1) = Ht ! down layer
c       Spin-orbital couple energy
c      ! left-right couple energy
            Hamiltonian(m,m-len2+1) = Hso
            Hamiltonian(m,m-len2-1) = Hso
c      ! up-down couple energy
            Hamiltonian(m,m-len2-48*1) = Hso
            Hamiltonian(m,m-len2+48*1) = Hso
c         superconductor piring energy
            Hamiltonian(m,m-2*len2) = del ! In do loop has increase 1 layer
      end do
c==============================================================================
! spin-down create operator
!   ??????????????????????????????????????????????????????????????????
c       Hopping energy
      ! left-right hop energy (x-direction)
      do m=3*48*48+1,4*48*48-5 !?????????
            Hamiltonian(m,m-1) = Ht
            Hamiltonian(m,m+1) = Ht 
c     ! up-down hop energy(y-driection)
            Hamiltonian(m,(m+1)) = Ht ! up layer
            Hamiltonian(m,(m-1)) = Ht ! down layer
c       Spin-orbital couple energy
c      ! left-right couple energy
            Hamiltonian(m,m-1-48*48) = Hso
            Hamiltonian(m,m+1-48*48) = Hso
c      ! up-down couple energy
            Hamiltonian(m,(m+1)-48*49) = Hso
            Hamiltonian(m,(m-1)-48*49) = Hso
c         superconductor piring energy
            Hamiltonian(m,m-3*len2) = del ! In do loop has increase 1 layer
      end do
C================================Save Matrix with .dat file type===============
c      open(11,file="hamiltonian.dat")
c      !write(11,*)Hamiltonian(24,24)
c      close(11)
C===============================Disp something to verify code correct==========
      counter = 0
c      do m = 1,48*48*4
c            do l=1,48*48*4
c                  if (Hamiltonian(m,l) .ne. 0) then
c                       ! counter = counter + 1
c                       write(11,"(48(F9.4))",advance='No')Hamiltonian(m,1)
c                  end if
c            end do
c      enddo
c      write(*,*)counter
c      write(11,"(48(F9.4))",advance='No')Hamiltonian
c      close(11)
c===============================TEST END====================================
      stop 
      end

C      real function delta(i)
C      use param
c     delta = V/2*sum(u-spin_up*v^(\dagger)-spin_down*tanh(E_n/(2*kb*t)))
c     At this,we don't known the eigen-value and egenfunction(u,v)c      
C      end function delta

C      real function phi(x)
c      ! This function is a phase factor calculator,should be a integral
C      use param
C      real x
C      end function phi
