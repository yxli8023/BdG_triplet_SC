!      Author :YuXuanLi
!      E-mail:yxli406@gmail.com
!      Article:Phys.Rev.B 96,184508
    module param
    implicit none 
    integer xn,yn
    ! xn : Numbers of lattice points in x-direction
    ! yn : Numbers of lattice points in y-direction
    parameter(xn = 24,yn = 24)
    integer,parameter::N=xn*yn*4
    real,parameter::pi = 3.14159265359
    real h ! h is external magnetic field
    real V ! V is point potential
    real kb ! Boltzmann constant
    real T ! Temperature
    real t0 ! t0 is the couple strength of Hubbrd model
    real lambda ! lambda is the spin-orbital coupling strength
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
    h = 0.6 ! zeeman  field
    mu = -4  ! chemical potential
    V = 0    ! point potential
    t0 = 1   ! hopping energy
    lambda = 0.5*t0   ! couple energy
    B = 2*phi0/(xn*yn)  ! magnetic field
    Kb = 1    !
    T = 1e-4     ! temperature
    beta = 1/(Kb*T)
    Ham = (0,0) 
!===================Construct Ham Matrix================
    call diag()   ! setting diag to matrix
    call hopping() ! full of elementary about hopping energy
    call couple()  ! couple energy elementary put in matrix 
    call pair_Init(0)     ! 在求解过程中要用自恰的方法，所以需要对超导配对能进行重新赋值
    !call potential(xp,yp)
    call eig()
    call loop()
    call Majorana() 
    stop
    end
!=============================================================
    subroutine delc()
    use param
    integer m,l,i
    integer num     
    complex,external::delta
    open(12,file="op.dat")
    do m = 1,xn
        do l = 1,yn
            i = (m-1)*yn+l
            del(m,l)= delta(i,i,w,Ham_diag)
            write(12,*)m,l,del(m,l)
        end do
    end do
    close(12)
    
    open(14,file='opmodule.dat')
    do m = 1,xn
        do l = 1,yn
            write(14,*)m,l,abs(del(m,l))**2
        end do
    end do
    close(14)

    return
    end subroutine delc
!===========================================================
    subroutine loop()   
    use param
    integer num,ref
    !store value of del and use for self-consistently
    complex del_loop(xn,yn),del_err(xn,yn)
    del_loop = 0
    del_err = 0
    call delc() ! 计算出一组delta
    num = 0    
    ref = 3                    
    do while(ref .gt. 1)
        open(24,file="count.dat")
        open(25,file="check.dat")
            ref = 0
            do m=1,xn
                do l=1,yn
                    !del_err用来存储两次计算得到的delta的差值,用来进行自恰条件的比较
                    del_err(m,l)=del(m,l)-del_loop(m,l)
                    del_loop(m,l) = del(m,l)
                    if (abs(real(del_err(m,l)))>eps) ref=ref+1
                    if (abs(imag(del_err(m,l)))>eps) ref=ref+1
                end do
            end do
            write(24,*)num
            write(25,*)ref
            num = num + 1
            call matrix()  ! Reconstruction Matrix
            call eig()
            call delc() ! Get a new value for del use new w and Ham
    end do
    close(24)
    close(25)
    return
    end subroutine loop
!==========================================================
    subroutine Majorana()
    use param
    complex z1(xn*yn),z2(xn*yn)
    integer m,l,i
    z1 = 0
    z2 = 0
    open(11,file="gamma1.dat")
    open(12,file="gamma2.dat")
    do m = 1,xn
        do l = yn
            i = (m-1)*xn + l
            z1(i)=(Ham(ZeroPoint,i)+Ham(ZeroPoint+1,len2+i))/sqrt(2.0)
            z2(i)=image*(Ham(ZeroPoint,i)-Ham(ZeroPoint+1,len2+i))/sqrt(2.0)  
            write(11,*)m,l,abs(z1(i))
            write(12,*)m,l.abs(z2(i)) 
    end do 
    close(11)
    close(12)
    return
    end subroutine Majorana
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
        s = s + eig_vec(x,kk)*conjg(eig_vec(y+len2*2,kk))*tanh(eig_val(kk)/2*beta)
    end do
    delta = 5.0/2.0*s
    end function delta
!-------------------------- point potential -----------------
    subroutine potential(x,y)
    use param
    integer xy,x,y
    xy = (y-1)*xn+x
    Ham(xy,xy)=Ham(xy,xy)+V ! position at spin-up oriention   (1,1)block
    Ham(xy+len2,xy+len2)=Ham(xy+len2,xy+len2)+V ! position at spin-down oriention   (2,2) block
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
        !do l = 1,yn
        !i = (m-1)*x_n+l
        Ham(i,len2*2+i) = delta(i,i,w,Ham_diag)  ! C-up*C-down    (1,3)block
        Ham(len2*2+i,i) = conjg(delta(i,i,w,Ham_diag))  ! E-up*E-down   (3,1)block
!-----------------------------------------negative energy ---------------------------------------------
        Ham(len2+i,len2*3+i) = -delta(i,i,w,Ham_diag)  ! C-up*C-down     (2,4)block
        Ham(len2*3+i,len2+i) = -conjg(delta(i,i,w,Ham_diag))  ! E-up*E-down   (4,2)block
        !end do
    end do
    end if
    return
    end subroutine pair_Init         
!******************************** hopping energy ***********************
    subroutine hopping()
    use param
    complex,external::phi ! function statement
!     ------------------------  x-direction right hopping -------------------
!     ---------------------- positive energy  -----------------------------------------
!     right hopping is positive       left hopping is negatie 
    do l = 1,yn                 ! 最后一列向右跳跃没有考虑到，但是向左跳跃考虑率到了
        do m = 2+(l-1)*xn,xn*l-1    ! (1,1)block   spin-up
            Ham(m,m+1) = -t0*phi(l)  ! (i--->i+1)  正向跳跃
            Ham(m,m-1) = -t0*conjg(phi(l))  ! (i+1--->i) 相当于考虑了反向跳跃
        ! (2,2)block   spin-down
            Ham(len2+m,len2+m+1) = -t0*phi(l)
            Ham(len2+m,len2+m-1) = -t0*conjg(phi(l))
    !     --------------------------negative energy ------------------------
        ! (3,3) block
            Ham(len2*2+m,len2*2+m+1) = t0*conjg(phi(l))
            Ham(len2*2+m,len2*2+m-1) = t0*(phi(l))
        ! (4,4) block
            Ham(len2*3+m,len2*3+m+1) = t0*conjg(phi(l))
            Ham(len2*3+m,len2*3+m-1) = t0*(phi(l))
        end do
    end do
!----------------------------------------- y-direction up hooping -------------------------------------
!----------------------------------------- 最后一行向上跳跃并未考虑------------------------------------
! 由朗道规范可以知道，在y方向的跳跃不包含位相信息(除了最后一行)
! 所以在此处没有考虑最后一行在y方向的跳跃情况
    do m = xn+1,xn*(yn-1)   ! (1,1) block
        Ham(m,m+len1) = -t0 ! up(positive) y方向上的向上跳跃
        Ham(m,m-len1) = -t0 ! up(negative) y 方向上的向下跳跃
    ! (2,2)  block
        Ham(len2+m,len2+m+len1) = -t0
        Ham(len2+m,len2+m-len1) = -t0
    ! (3,3) block
        Ham(len2*2+m,len2*2+m+len1) = t0
        Ham(len2*2+m,len2*2+m-len1) = t0
    ! (4,4) block
        Ham(len2*3+m,len2*3+m+len1) = t0
        Ham(len2*3+m,len2*3+m-len1) = t0
    end do
! ------------ first row down-direction hopping ----------------------------------------------------
! 最后一行向第一行跳跃的过程中，人为的加上一个位相，保证整个波函数跳跃一圈的位相是2pi
!     up hopping is positive               down hopping is negative
    do m = 1,xn   ! (1,1)block
        Ham(m,m+len1) = -t0  !  first row hopping to up direction
        Ham(m,xn*(yn-1)+m) = -t0*(phi(m*yn))  !
        Ham(xn*(yn-1)+m,m) = -t0*conjg(phi(m*yn))  ! last row hopping t first row
        Ham(xn*(yn-1)+m,xn*(yn-1)+m-len1) = -t0 !last row hopping to down direction
! (2,2)block
        Ham(len2+m,len2+m+len1) = -t0  !  first row hopping to up direction
        Ham(len2+m,len2+xn*(yn-1)+m) = -t0*(phi(m*yn))  ! 
        Ham(len2+xn*(yn-1)+m,len2+m) = -t0*conjg(phi(m*yn))  ! last row hopping t first row
        Ham(len2+xn*(yn-1)+m,len2+xn*(yn-1)+m-len1) = -t0 !last row hopping to down direction
! (3,3)block
        Ham(len2*2+m,len2*2+m+len1) = t0  !  first row hopping to up direction
        Ham(len2*2+m,len2*2+xn*(yn-1)+m) = t0*conjg(phi(m*yn))  !
        Ham(len2*2+xn*(yn-1)+m,len2*2+m) = t0*(phi(m*yn))  ! last row hopping t first row
        Ham(len2*2+xn*(yn-1)+m,len2*2+xn*(yn-1)+m-len1) = t0 !last row hopping to down direction 
! (4,4)block
        Ham(len2*3+m,len2*3+m+len1) = t0  !  first row hopping to up direction
        Ham(len2*3+m,len2*3+xn*(yn-1)+m) = t0*conjg(phi(m*yn))  ! 
        Ham(len2*3+xn*(yn-1)+m,len2*3+m) = t0*(phi(m*yn))  ! last row hopping t first row
        Ham(len2*3+xn*(yn-1)+m,len2*3+xn*(yn-1)+m-len1) = t0 !last row hopping to down direction
    end do      
!------------------------ first coloum left-direction hopping (include last coloum hopping towards first cloum cause of matrix symmetry consider)----
    do m = 1,yn   ! (1,1) block
        Ham(m*xn,m*xn-(xn-1)) = -t0*phi(m)   ! last point hopping to right
        Ham(m*xn,m*xn-1) = -t0*conjg(phi(m))! last point hopping to left
        Ham(1+(m-1)*xn,1+(m-1)*xn+1) = -t0*phi(m)! first point hopping to right
        Ham(1+(m-1)*xn,m*xn) = -t0*conjg(phi(m))! first point hopping to left
    ! (2,2) block
        Ham(len2+m*xn,len2+m*xn-(xn-1)) = -t0*phi(m)
        Ham(len2+m*xn,len2+m*xn-1) = -t0*conjg(phi(m))
        Ham(len2+1+(m-1)*xn,len2+1+(m-1)*xn+1) = -t0*phi(m)  ! first point hopping to right
        Ham(len2+1+(m-1)*xn,len2+m*xn) = -t0*conjg(phi(m))  ! first point hopping to left
    ! (3,3) block
        Ham(len2*2+m*xn,len2*2+m*xn-(xn-1)) = t0*conjg(phi(m))
        Ham(len2*2+m*xn,len2*2+m*xn-1) = t0*(phi(m))
        Ham(len2*2+1+(m-1)*xn,len2*2+1+(m-1)*xn+1) = t0*conjg(phi(m))! first point hopping to right
        Ham(len2*2+1+(m-1)*xn,len2*2+m*xn) = t0*(phi(m))! first point hopping to left
    ! (4,4) block
        Ham(len2*3+m*xn,len2*3+m*xn-(xn-1)) = t0*conjg(phi(m))
        Ham(len2*3+m*xn,len2*3+m*xn-1) = t0*(phi(m))
        Ham(len2*3+1+(m-1)*xn,len2*3+1+(m-1)*xn+1) = t0*conjg(phi(m))! first point hopping to right
        Ham(len2*3+1+(m-1)*xn,len2*3+m*xn) = t0*(phi(m))! first point hopping to left
    end do            
    return
    end subroutine
!***************************************spin-orbital couple energy*****************************************
    subroutine couple()
    use param
    complex,external::phi 
!XXXX                              - x direction copule (right and left couple) --------
!----------------------------------positive energy couple (1,2)block----------------------
    do l = 1,yn                  ! 用l代表层数，每一层的最右边的位置上向右边的耦合没有考虑
    do m = 2+(l-1)*xn,xn*l-1     !(1,2)block    
        Ham(m,len2+m+1) = image*lambda*phi(l) ! spin_up-spin_down couple right  正向耦合(1--->2)
        Ham(m,len2+m-1) = -image*lambda*conjg(phi(l))  ! left反向耦合(2--->1)
! -------------------------------- positive energy copule (2,1)block----------------
    !(2,1)block
        Ham(len2+m,m+1) = image*lambda*phi(l) ! spin_up-spin_down couple  right
        Ham(len2+m,m-1) = -image*lambda*conjg(phi(l))  ! left
! ------------------------------- negative energy couple (3,4) block-----------------
        Ham(len2*2+m,len2*3+m+1) = image*lambda*conjg(phi(l))
        Ham(len2*2+m,len2*3+m-1) = -image*lambda*(phi(l))
! ------------------------------- negative energy couple (4,3) block-----------------------
        Ham(len2*3+m,len2*2+m+1) = image*lambda*conjg(phi(l)) ! spin_up-spin_down couple  right
        Ham(len2*3+m,len2*2+m-1) = -image*lambda*(phi(l))  ! left
    end do
    end do
! XXXX--------------------------------------------------------------------------------------------------------------
!YYYY ------------------------------------ y direction copule (up and down couple)加一层---------------------------
    do m = xn+1,xn*(yn-1)
        Ham(m,len2+m+len1) = lambda
        Ham(m,len2+m-len1) = -lambda
! ----------------------------------------- (2,1) block ------------------------------------------------
        Ham(len2+m,m+len1) = -lambda
        Ham(len2+m,m-len1) = lambda
! ----------------------------------------- (3,4) block -------------------------------------------------
        Ham(len2*2+m,len2*3+m+len1) = lambda
        Ham(len2*2+m,len2*3+m-len1) = -lambda
! -----------------------------------------(4,3) block -----------------------------------------------
        Ham(len2*3+m,len2*2+m+len1) = -lambda
        Ham(len2*3+m,len2*2+m-len1) = lambda
    end do
! ------------------------ right and left boundary---------------------------------
!------------------------ positive energy ---------------------------------------
    do m = 1,yn   !  (1,2)block
        ! 
        Ham(m*xn,len2+m*xn-(xn-1)) = image*lambda*phi(m)
        Ham(m*xn,len2+m*xn-1) = -image*lambda*conjg(phi(m)) ! 增加一层会使自旋取向改变一次
        ! 
        Ham(m*xn-(xn-1),len2+m*xn-(xn-1)+1) = image*lambda*phi(m)
        Ham(m*xn-(xn-1),len2+m*xn) = -image*lambda*conjg(phi(m))
    !  (2,1)block
        Ham(len2+m*xn,m*xn-(xn-1)) = image*lambda*phi(m)
        Ham(len2+m*xn,m*xn-1) = -image*lambda*conjg(phi(m)) ! 增加一层会使自旋取向改变一次(1------>48)
        Ham(len2+m*xn-(xn-1),m*xn-(xn-1)+1) = image*lambda*phi(m)
        Ham(len2+m*xn-(xn-1),m*xn) = -image*lambda*conjg(phi(m))
!---------------------------------------negative energy---------------------------------
    !  (3,4)block
        Ham(len2*2+m*xn,len2*3+m*xn-(xn-1)) = image*lambda*conjg(phi(m))
        Ham(len2*2+m*xn,len2*3+m*xn-1) = -image*lambda*(phi(m))! 增加一层会使自旋取向改变一次
        Ham(len2*2+m*xn-(xn-1),len2*3+m*xn-(xn-1)+1) = image*lambda*conjg(phi(m))
        Ham(len2*2+m*xn-(xn-1),len2*3+m*xn) = -image*lambda*(phi(m))
    !  (4,3)block
        Ham(len2*3+m*xn,len2*2+m*xn-(xn-1)) = image*lambda*conjg(phi(m))
        Ham(len2*3+m*xn,len2*2+m*xn-1) = -image*lambda*(phi(m)) ! 增加一层会使自旋取向改变一次
        Ham(len2*3+m*xn-(xn-1),len2*2+m*xn-(xn-1)+1) = image*lambda*conjg(phi(m))
        Ham(len2*3+m*xn-(xn-1),len2*2+m*xn) = -image*lambda*(phi(m))
    end do           
!----------------------- up and down boundary -------------------------------
! 
    do m = 1,xn       ! (1,2)block
        Ham(m,len2+m+len1) = lambda   ! up
        Ham(m,len2+xn*(yn-1)+m) = -lambda*(phi(m*yn)) ! down  
        Ham(xn*(yn-1)+m,len2+m) = lambda*conjg(phi(m*yn))!up  
        Ham(xn*(yn-1)+m,len2+xn*(yn-1)-len1+m) = -lambda  !down
    ! (2,1)block
        Ham(len2+m,m+len1) = -lambda   !up
        Ham(len2+m,xn*(yn-1)+m) = lambda*(phi(m*yn))  ! down
        Ham(len2+m+xn*(yn-1),m) = -lambda*conjg(phi(m*yn))
        Ham(len2+m+xn*(yn-1),m+xn*(yn-1)-len1) = lambda
! ----------------------------------- negative energy ----------------------------
    ! (3,4)block
        Ham(len2*2+m,len2*3+m+len1) = lambda
        Ham(len2*2+m,len2*3+xn*(yn-1)+m) = -lambda*conjg(phi(m*yn))
        Ham(len2*2+m+xn*(yn-1),len2*3+m) = lambda*(phi(m*yn))
        Ham(len2*2+m+xn*(yn-1),len2*3+m+xn*(yn-1)-len1) = -lambda
    ! (4,3)block
        Ham(len2*3+m,len2*2+m+len1) = -lambda
        Ham(len2*3+m,len2*2+xn*(yn-1)+m) = lambda*conjg(phi(m*yn))
        Ham(len2*3+m+xn*(yn-1),len2*2+m) = -lambda*(phi(m*yn))
        Ham(len2*3+m+xn*(yn-1),len2*2+m+xn*(yn-1)-len1) = lambda
    end do
!YYYY--------------------------------------------------------------------
    return
    end subroutine couple
!--------------------------- Diagonal element setting complete-----------
    subroutine diag()
    use param
        ! spin-up positive energy
    do m = 1,xn*yn
    !(1,1)block  ! spin-down positive energy
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
!===========================================================
    subroutine eig()
    use param
    integer m
    lwork = -1
    liwork = -1
    lrwork = -1
    call cheevd('V','Upper',N,Ham,lda,w,work,lwork,rwork,lrwork,iwork,liwork,info)
    lwork = min(2*N+N**2, int( work( 1 ) ) )
    lrwork = min(1+5*N+2*N**2, int( rwork( 1 ) ) )
    liwork = min(3+5*N, iwork( 1 ) )
    call cheevd('V','Upper',N,Ham,lda,w,work,lwork,rwork,lrwork,iwork,liwork,info)
    Ham_diag = Ham
    if( info .GT. 0 ) then
        open(11,file="mes.txt")
        write(11,*)'The algorithm failed to compute eigenvalues.'
        close(11)
    end if
    open(12,file="eigval.dat",status="unknown")
	do m = 1,N
		write(12,*)m,w(m)
	end do
	close(12)
    return
    end subroutine eig
!===========================================================
    subroutine matrix_verify()
    use param
    integer i,j
    integer ccc
    ccc = 0
    open(16,file = 'MatrixVerify.dat')
    do i = 1,N
        do j = 1,N
            if (Ham(i,j) .ne. conjg(Ham(j,i)))then
                ccc = ccc +1
                write(*,*)"Your matrix is not Hermitian"
                write(16,*)i,j
                write(16,*)Ham(i,j)
                write(16,*)"==================================="
                !stop
            end if
        end do
    end do
    write(16,*)ccc
    end subroutine matrix_verify
!--------------------------------------------------------------------------------------------
    subroutine count_number()
    use param
    integer i,j
    integer c
    c = 0
    open(13,file='MatrixNumber.dat')
    do i = 1,N
        do j = 1,N
            if (abs(Ham(i,j)) .ne. 0)then
                c = c + 1
            end if
        end do
    end do 
    write(13,*)"Non-zero number is:",c
    !close(13) 
    return
    end subroutine count_number

    
