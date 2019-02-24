!      Author :YuXuanLi
!      E-mail:yxli406@gmail.com
!      Article:Phys.Rev.B 96,184508
      module param
      implicit none
      ! include two spin and create and elimination operator so have 4
      ! types
      integer xn,yn
      ! xn : Numbers of lattice points in x-direction
      ! yn : Numbers of lattice points in y-direction
      parameter(xn = 24,yn = 24)
      integer,parameter::N=xn*yn*4
      real,parameter::pi = 3.14159265359
      real h ! h is external magnetic field
      real V ! V is point potential
      real a ! Lettice constant
      real kb ! Boltzmann constant
      real T ! Temperature
      real t0 ! t0 is the couple strength of Hubbrd model
      real soc ! soc is the spin-orbital coupling strength
      complex::image=(0.0,1.0) ! Image unit
      real mu ! mu is the chemcial potential
      real phi0 ! superconduct flux quantum
      real B ! Zeeman Magnetic
      real beta
      integer xp,yp
      integer::ZeroPoint = xn*yn*2
      integer,parameter::len2 = xn*yn
      integer,parameter::len1 = xn
      real,parameter::eps = 1e-5 ! measure distance of two delta
      complex del(xn,yn)
      complex Ham(N,N),Ham_diag(N,N)
      integer nbr(4,len2) ! 记录格点的周边信息
      ! lapack parameter
      integer::lda = N
      integer,parameter::lwmax=2*N+N**2
      real,allocatable::w(:)
      complex*8,allocatable::work(:)
      real,allocatable::rwork(:)
      integer,allocatable::iwork(:)
      integer lwork   ! at least 2*N+N**2
      integer lrwork    ! at least 1 + 5*N +2*N**2
      integer liwork   ! at least 3 +5*N
      integer info
      end module param
!===============================Program Start=============================
      program BdG
      use param
      integer l,m,i,j ! Both of them are loop variation
      complex,external::phi ! function statement
      complex,external::delta ! function statemen
      real,external::sqrt
      external::cheevd
      allocate(w(N))
      allocate(work(lwmax))
      allocate(rwork(1+5*N+2*N**2))
      allocate(iwork(3+5*N))
!===============================Parameter Setting=======================
      phi0 = pi! flux quantum
      a = 1.0 ! lattice constant
      h = 0.6 ! zeeman  field
      mu = -4  ! chemical potential
      V = 0    ! point potential
      t0 = 1   ! hopping energy
      soc = 0.5*t0   ! couple energy
      B = 2*phi0/(xn*yn)  ! magnetic field
      Kb = 1    !
      T = 1e-4     ! temperature
      beta = 1/(Kb*T)
      Ham = (0,0)
      xp = 1
      yp = 1
      open(11,file='param.dat')
      write(11,*)"h=",h
      write(11,*)"V=",V
      write(11,*)"T=",T
      write(11,*)"xp=",xp
      write(11,*)"yp=",yp
      close(11)
!=======================================================
      do m=1,xn
         do l=1,yn
            i=(l-1)*xn+m
!            write(*,*)m,l
            nbr(1,i)=i+1  ! 所有向右跳跃的点
            if(m.eq.xn)nbr(1,i)=nbr(1,i)-xn  ! 最右边的点向右跳跃应该回到最左端
            nbr(2,i)=i-1  ! 所有向左跳跃的点
            if(m.eq.1)nbr(2,i)=nbr(2,i)+xn    ! 最左边的点向左跳跃应该回到最右端
            nbr(3,i)=i+xn ! 所有向上跳跃的点
            if(l.eq.yn)nbr(3,i)=nbr(3,i)-len2 ! 最上端一行跳跃到最下端   
            nbr(4,i)=i-xn ! 所有向下跳跃的点
            if(l.eq.1)nbr(4,i)=nbr(4,i)+len2  ! 最下端一行跳跃到最上端
         enddo
      enddo
!===================Construct Ham Matrix================
      call diag()
      call hopping()
      call couple()
      call pair_Init(0)
      call potential(xp,yp)

!================eig-value and eig-vector solve ==========
      call eig()
      call Majorana()
!==================self-consistently==================
      call loop()
!--------------------------- open file------------------------
      open(11,file='eigenvalue.dat')
      do m = 1,N
            write(11,*)w(m) ! eigenvalue
      end do
!------------------------- close file ------------------------
      close(11)
!--------------  IO format control ---------------------------
100   format(48f9.6) !8位宽度输出数据,小数部分6位
      stop
      end
!========================= PROGRAM END ==================================


!========================  del Read ========================
! 在计算过程中，可能需要通过之前已经计算得到的delta来进行进一步的
! 的计算，该子过程可以将对应位置的delta赋值到矩阵中的对应位置
      subroutine delRead(str)
      use param
      character*20 str
      complex temp(xn*yn)
      integer m
      open(11,file=str)
      do m = 1,xn*yn
            read(11,*)temp(m)
            Ham(m,len2*2+m) = temp(m)
            Ham(len2*2+m,m) = conjg(temp(m))
!-------negative energy --------
            Ham(len2+m,len2*3+m) = -temp(m)  ! C-up*C-down   (2,4)block
            Ham(len2*3+m,len2+m) = -conjg(temp(m))  ! E-up*E-down    (4,2)block
      end do
      return
      end subroutine
!============================= Function Construct ========================
      complex function phi(y)  ! Phase function   A = (-By,0,0)
      use param
      integer y ! position of lattice
      phi = cexp(-image*pi/phi0*B*y)
      end function phi
!================================= Energy Gap function =======
      complex function delta(x,y,eig_val,eig_vec)
      use param
      integer x,y,kk
      real eig_val(N)
      complex s,eig_vec(N,N)
      s = (0,0)
      do kk = 1,N  ! Hermitian Matrix eigvalue is real number
            s = s + eig_vec(x,kk)*
     &      conjg(eig_vec(y+len2*2,kk))*tanh(eig_val(kk)/2*beta)
      end do
      delta = 5.0/2.0*s
      end function delta
!-------------------------- point potential -----------------
!------------------------ point potential element setting complete-----------
! add a point potential at xy
! 对于矩形的点阵模型，假设x方向格点数为xn,y方向格点数为yn
! 在(a,b)点增加point potential时,xy=(a-1)*xn+b
      subroutine potential(x,y)
      use param
      integer xy,x,y
      xy = (y-1)*xn+x
!     Positive energy (正能位置的杂质拖动正能对应的态)
      Ham(xy,xy)=Ham(xy,xy)+V ! position at spin-up oriention   (1,1)block
      Ham(xy+len2,xy+len2)=Ham(xy+len2,xy+len2)+V ! position at spin-down oriention   (2,2) block
!     Negative energy (负能对应的杂质拖动负能对应的态)
      Ham(xy+len2*2,xy+len2*2)=Ham(xy+len2*2,xy+len2*2)-V    ! (3,3)block
      Ham(xy+len2*3,xy+len2*3)=Ham(xy+len2*3,xy+len2*3)-V       !(4,4)block
      return
      end subroutine potential
!***************************** Pariing energy **********************
      subroutine pair_Init(input)
      use param
      integer m,l,input,i
      complex,external::delta
      real re
      !re = 0.5
      if (input .eq. 0) then
!--------------positive energy -----
      do m = 1,xn*yn
      !当input=0时，采用随机的delta来进行自恰过程
            call random_seed()
            call random_number(re)  ! C-up*C-down    (1,3)block
            Ham(len2*2+m,m) = re
            !call random_number(re)  ! C-down*C-up    (3,1)block
            Ham(m,len2*2+m) = re
!-------negative energy --------
            Ham(len2+m,len2*3+m) = -re  ! C-up*C-down   (2,4)block
            Ham(len2*3+m,len2+m) = -re  ! E-up*E-down    (4,2)block
      end do

      else if (input .eq. 100) then
      ! 此步骤用来在已经自恰收敛的基础上,通过求解好的delta数据来进行进一步的计算
            !通过已经求解得到的数据进行下一步的计算,其一是为了避免程序的不正常中断
            !仍然可以从原来的位置起开始计算,其二是为了增加程序的可拓展性
            call delRead("order.dat")
      else
      do i = 1,xn*yn
      ! 该过程才是进行自恰求解的关键步骤,通过求解得到的本征值和本征矢来计算出
      ! 新的delta,然后代入到Hamiltonian Matrix中进行迭代过程
            Ham(i,len2*2+i) = delta(i,i,w,Ham_diag)  ! C-up*C-down    (1,3)block
            Ham(len2*2+i,i) = conjg(delta(i,i,w,Ham_diag))  ! E-up*E-down   (3,1)block
!-----------------------------------------negative energy ---------------------------------------------
            Ham(len2+i,len2*3+i) = -delta(i,i,w,Ham_diag)  ! C-up*C-down     (2,4)block
            Ham(len2*3+i,len2+i) = -conjg(delta(i,i,w,Ham_diag))  ! E-up*E-down   (4,2)block
      end do
      end if
      return
      end subroutine pair_Init

!--------------------------- Diagonal element setting complete-----------
      subroutine diag()
      use param
            ! spin-up positive energy
      do m = 1,xn*yn   !(1,1)block  ! spin-down positive energy
            Ham(m,m) = h-mu
      !(2,2)block
            Ham(len2+m,len2+m) = -h-mu
      !(3,3)block  ! spin-down negative energy
            Ham(len2*2+m,len2*2+m) = -(-h-mu)
      !(4,4)block   ! spin-up negative energy
            Ham(len2*3+m,len2*3+m) = -(h-mu)
      end do
      return
      end subroutine diag
!=============================================================
      subroutine del_c(num)
      use param
      integer m,l,i
      integer num
      complex,external::delta
      open(12,file="order.dat")
      open(13,file="energy.dat")
      open(14,file='order_module.dat')
      do m = 1,xn
            do l = 1,yn
                  i = (m-1)*yn+l
                  del(m,l)= delta(i,i,w,Ham_diag)
                  write(12,*)del(m,l)
                  write(14,*)abs(del(m,l))
            end do
      end do
      do m = 1,N
            write(13,*)w(m) ! eigenvalue
      end do
      close(12)
      close(13)
      close(14)
      return
      end subroutine del_c
!============================================================
      subroutine matrix()
      use param
      complex,external::delta
      Ham = (0,0)
      call diag()
      call hopping()
      call couple()
! In pair_Init subroutine 1 as a parameter to chose new order as element
      call pair_Init(1)
      call potential(xp,yp)
      end subroutine matrix
!===========================================================================
      subroutine loop()
      use param
      integer num,ref
      !store value of del and use for self-consistently
      complex del_loop(xn,yn),del_err(xn,yn)
      del_loop = 0
      del_err = 0
      call del_c(0) ! 计算出一组delta
      num = 0
      ref = 3
      do while(ref .gt. 1)
      open(24,file="count.dat")
      open(25,file="check.dat")
            ref = 0
      !do while(minval(abs(del_err))>eps)
            do m=1,xn
            do l=1,yn
                  !del_err用来存储两次计算得到的delta的差值,用来进行自恰条件的比较
                  del_err(m,l)=del(m,l)-del_loop(m,l)
                  del_loop(m,l) = del(m,l)  ! 重新得到一组del_loop的值，用于下一次的比较
                  if (abs(real(del_err(m,l)))>eps) ref=ref+1
                  if (abs(imag(del_err(m,l)))>eps) ref=ref+1
            end do
            end do
            write(24,*)num
            write(25,*)ref
            num = num + 1
            call matrix()  ! Reconstruction Matrix
            call eig()
            call del_c(num) ! Get a new value for del use new w and Ham
      !end do
      end do
      close(24)
      close(25)
      return
      end subroutine

!===========================================================
      subroutine eig()
      use param
      lwork = -1
      liwork = -1
      lrwork = -1
      call cheevd('V','Upper',N,Ham,lda,w,work,lwork
     &            ,rwork,lrwork,iwork,liwork,info)
      lwork = min(2*N+N**2, int( work( 1 ) ) )
      lrwork = min(1+5*N+2*N**2, int( rwork( 1 ) ) )
      liwork = min(3+5*N, iwork( 1 ) )
      call cheevd('V','Upper',N,Ham,lda,w,work,lwork
     &            ,rwork,lrwork,iwork,liwork,info)
      Ham_diag = Ham
      if( info .GT. 0 ) then
            open(11,file="mes.txt")
            write(11,*)'The algorithm failed to compute eigenvalues.'
            close(11)
      end if

      return
      end subroutine eig
!==========================================================
      subroutine Majorana()
      use param
      complex z1(xn*yn),z2(xn*yn)
      integer m,l
      z1 = 0
      z2 = 0
      open(11,file="gamma1.dat")
      open(12,file="gamma2.dat")
      do m = 1,xn*yn
        z1(m)=(Ham(m,ZeroPoint)+Ham(m,ZeroPoint+1))/sqrt(2.0)
        z2(m)=image*(Ham(m,ZeroPoint)-Ham(m,ZeroPoint+1))/sqrt(2.0)
        write(11,*)abs(z1(m))
        write(12,*)abs(z2(m))
      end do
      close(11)
      close(12)
      return
      end subroutine Majorana
!================================= Hopping Term ====================================
      subroutine hopping()
      use param
      integer m,l,i
      complex,external::phi
      do m=1, xn
            do l=1, yn
            i=(l-1)*xn+m
            ! hopping term
            ! (1,1) block
            ham(i,nbr(1,i))=-t0*conjg(phi(l))
            ham(i,nbr(2,i))=-t0*phi(l)
            ham(i,nbr(3,i))=-t0
            ham(i,nbr(4,i))=-t0
            ! 上边界人为引入位相
            if (l.eq.yn) ham(i,nbr(3,i))=-t0*phi(yn*m)
            if (l.eq.1) ham(i,nbr(4,i))=-t0*conjg(phi(yn*m))
            ! (2,2)block
            ham(i+len2,nbr(1,i)+len2)=-t0*conjg(phi(l))
            ham(i+len2,nbr(2,i)+len2)=-t0*phi(l)
            ham(i+len2,nbr(3,i)+len2)=-t0
            ham(i+len2,nbr(4,i)+len2)=-t0
            if (l.eq.yn) ham(i+len2,nbr(3,i)+len2)=
     &         -t0*phi(yn*m)
            if (l.eq.1) ham(i+len2,nbr(4,i)+len2)=
     &         -t0*conjg(phi(yn*m))
            ! (3,3)block
            ham(i+len2*2,nbr(1,i)+len2*2)=t0*phi(l)
            ham(i+len2*2,nbr(2,i)+len2*2)=t0*conjg(phi(l))
            ham(i+len2*2,nbr(3,i)+len2*2)=t0
            ham(i+len2*2,nbr(4,i)+len2*2)=t0
            if (l.eq.yn) ham(i+len2*2,nbr(3,i)+len2*2)=
     &       t0*conjg(phi(yn*m))
            if (l.eq.1) ham(i+len2*2,nbr(4,i)+len2*2)=
     &       t0*phi(yn*m)
            ! (4,4)block
            ham(i+len2*3,nbr(1,i)+len2*3)=t0*phi(l)
            ham(i+len2*3,nbr(2,i)+len2*3)=t0*conjg(phi(l))
            ham(i+len2*3,nbr(3,i)+len2*3)=t0
            ham(i+len2*3,nbr(4,i)+len2*3)=t0 
            if (l.eq.yn) ham(i+len2*3,nbr(3,i)+len2*3)=
     &       t0*conjg(phi(yn*m))
            if (l.eq.1) ham(i+len2*3,nbr(4,i)+len2*3)=
     &       t0*phi(yn*m)
            enddo
      enddo
      return
      end subroutine
!=================================== Couple Term =================================
      subroutine couple()
      use param
      integer m,l,i
      complex,external::phi
      do m=1, xn
            do l=1, yn
            i=(l-1)*xn+m
            ! (1,2)block
            ham(i,nbr(1,i)+len2)=t0*image*0.5*conjg(phi(l))
            ham(i,nbr(2,i)+len2)=-t0*image*0.5*phi(l)
            ! y-direction
            ham(i,nbr(3,i)+len2)=soc
            ham(i,nbr(4,i)+len2)=-soc
            if (l.eq.yn) ham(i,nbr(3,i)+len2)=
     &         soc*phi(yn*m)
            if (l.eq.1) ham(i,nbr(4,i)+len2)=
     &         -soc*conjg(phi(yn*m))
            ! (2,1)block
            ham(i+len2,nbr(1,i))=t0*image*0.5*conjg(phi(l))
            ham(i+len2,nbr(2,i))=-t0*image*0.5*phi(l)
            ! y-direction
            ham(i+len2,nbr(3,i))=-soc
            ham(i+len2,nbr(4,i))=soc
            if (l.eq.yn) ham(i+len2,nbr(3,i))=
     &       -soc*phi(yn*m)
            if (l.eq.1) ham(i+len2,nbr(4,i))=
     &       soc*conjg(phi(yn*m))
            ! (3,4)block
            ham(i+len2*2,nbr(1,i)+len2+len2*2)=
     &       t0*image*0.5*phi(l)
            ham(i+len2*2,nbr(2,i)+len2+len2*2)=
     &       -t0*image*0.5*conjg(phi(l))
            !y-direction
            ham(i+len2*2,nbr(3,i)+len2+len2*2)=soc
            ham(i+len2*2,nbr(4,i)+len2+len2*2)=-soc 
            if (l.eq.yn) ham(i+len2*2,nbr(3,i)+len2*3)=
     &           soc*conjg(phi(yn*m))
            if (l.eq.1) ham(i+len2*2,nbr(4,i)+len2*3)=
     &          -soc*phi(yn*m)    
            ! (4,3) block
            ham(i+len2*3,nbr(1,i)+len2*2)=
     &       t0*image*0.5*phi(l)
            ham(i+len2*3,nbr(2,i)+len2*2)=
     &       -t0*image*0.5*conjg(phi(l))
            ! y-direction
            ham(i+len2*3,nbr(3,i)+len2*2)=-soc
            ham(i+len2*3,nbr(4,i)+len2*2)=soc 
            if (l.eq.yn) ham(i+len2*3,nbr(3,i)+len2*2)=
     &         -soc*conjg(phi(yn*m))
            if (l.eq.1) ham(i+len2*3,nbr(4,i)+len2*2)=
     &         soc*phi(yn*m)
            enddo
      enddo
      return
      end subroutine

