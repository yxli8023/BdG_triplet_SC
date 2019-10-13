    module pub
    implicit none
    integer xn,yn,N,len2,up
    parameter(xn=20,yn=20,N=xn*yn*4,len2=xn*yn,up=xn)
    complex Ham(N,N),Ham_diag(N,N)
    real t,lam,h,mu,B,phi0
    complex::im = (0.,1.0)
    real,parameter::pi = 3.14159265359
    integer::ZeroPoint = xn*yn*2 
    real::eps = 1e-5
    integer bcx,bcy
    complex del(xn,yn)  
    !-----------------LAPACK PACKAGE PARAM
    integer::lda = N
    integer,parameter::lwmax=2*N+N**2
    real,allocatable::w(:)
    complex,allocatable::work(:)
    real,allocatable::rwork(:)
    integer,allocatable::iwork(:)
    integer lwork   ! at least 2*N+N**2
    integer lrwork    ! at least 1 + 5*N +2*N**2
    integer liwork   ! at least 3 +5*N
    integer info
    end module pub
!==================================================
    program sol
    use pub
    !================ Physics memory allocate =================
    allocate(w(N))
    allocate(work(lwmax))
    allocate(rwork(1+5*N+2*N**2))
    allocate(iwork(3+5*N))
    !---------------------
    bcx = 1
    bcy = 1
    !----------------------
    t = 1.0
    lam = 0.5
    h = 0.6
    mu = -4
    del = 0.5
    phi0 = pi
    B = 2.0*phi0/len2  
    call main()
    call eigsol()
    call loop()
    stop
    end program sol
!===========================================================================
    subroutine loop()   
    use pub
    integer num,ref
    !store value of del and use for self-consistently
    complex del_loop(xn,yn),del_err(xn,yn)
    del_loop = 0
    del_err = 0
    call delc() !计算出一组delta
    num = 0    
    ref = 3                    
    do while(ref .gt. 1)
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
        write(25,*)ref
        num = num + 1
        call main()  ! Reconstruction Matrix
        call eigsol()
        call delc() ! Get a new value for del use new w and Ham
    end do
    close(25)
    return
    end subroutine
!==================================================
    subroutine delc()
    use pub
    integer m,l,i
    complex,external::delta
    open(12,file="order.dat")
    open(14,file='order_module.dat')
    do m = 1,xn
        do l = 1,yn
            i = (m-1)*yn + l
            del(m,l)= delta(i,i,w,Ham_diag)
            write(14,*)m,l,abs(del(m,l))
            write(12,*)del(m,l)
        end do
    end do
    close(12)
    close(14)
    return
    end subroutine delc
!================================= Energy Gap function =======
    complex function delta(x,y,eig_val,eig_vec)
    use pub
    integer x,y,kk
    real eig_val(N),beta
    complex s,eig_vec(N,N)
    beta = 1e5
    s = (0,0)
    do kk = 1,N  
          s = s + eig_vec(x,kk)*conjg(eig_vec(y+len2*3,kk))*tanh(eig_val(kk)/2*beta)
    end do
    delta = 5.0/2.0*s
    end function delta
!==========================================
    subroutine Majorana()
    use pub
    complex z1(xn*yn),z2(xn*yn)
    integer m,l,i
    z1 = 0
    z2 = 0
    open(11,file="gamma1.dat")
    open(12,file="gamma2.dat")
    do m = 1,yn
        do l = 1,xn
            i = (m-1)*xn + l
            z1(i)=(Ham(ZeroPoint,i)+Ham(ZeroPoint+1,len2+i))/sqrt(2.0)
            z2(i)=image*(Ham(ZeroPoint,i)-Ham(ZeroPoint+1,len2+i))/sqrt(2.0)  
            write(11,*)m,l,abs(z1(i))
            write(12,*)m,l,abs(z2(i)) 
        end do
    end do 
    close(11)
    close(12)
    return
    end subroutine Majorana
!========================== Function Construct ========================
    complex function phi(y)  ! Landau guage   A = (-By,0,0)
    use pub
    integer y ! position of lattice
    phi = cexp(-im*B*y)       
    end function phi
!=======================================================
    subroutine main()
    use pub
    complex,external::phi
    integer m,l,i
    Ham = 0
    do i = 1,len2
        Ham(i,i) = h - mu
        Ham(len2+i,len2+i) = -h - mu
    end do
    !------------- X ---------------------------
    do l = 1,yn
        do m = 2+(l-1)*xn,xn*l-1 
            ! (1,1)  
            Ham(m,m+1) = -t*phi(l)
            Ham(m,m-1) = -t*conjg(phi(l))
            ! (2,2)   
            Ham(len2+m,len2+m+1) = -t*phi(l)
            Ham(len2+m,len2+m-1) = -t*conjg(phi(l))
        end do
    end do
    !========== X boundry ==========================
    do m = 1,yn
    ! (1,1)   
        Ham(m*xn,m*xn-(xn-1)) = -t*phi(m)*bcx  
        Ham(m*xn,m*xn-1) = -t*conjg(phi(m))        
        Ham(1+(m-1)*xn,1+(m-1)*xn+1) = -t*phi(m) 
        Ham(1+(m-1)*xn,m*xn) = -t*conjg(phi(m))*bcx  
    ! (2,2)		
        Ham(len2+m*xn,m*xn-(xn-1)+len2) = -t*phi(m)*bcx
        Ham(len2+m*xn,m*xn-1+len2) = -t*conjg(phi(m))
        Ham(len2+1+(m-1)*xn,1+(m-1)*xn+1+len2) = -t*phi(m)
        Ham(len2+1+(m-1)*xn,m*xn+len2) = -t*conjg(phi(m))*bcx
    end do
    ! ---------------------------------------------
    do m = xn+1,xn*(yn-1)
    !(1,1) 
        Ham(m,m+up) = -t
        Ham(m,m-up) = -t
    !(2,2)		
        Ham(len2+m,len2+m+up) = -t
        Ham(len2+m,len2+m-up) = -t
    end do
    !---------------------------------------------------
    do m = 1,xn
    ! (1,1)
        Ham(m,m+up) = -t
        Ham(m,xn*(yn-1)+m) = -t*phi(m*yn)*bcy 
        Ham(xn*(yn-1)+m,m) = -t*conjg(phi(m*yn))*bcy 
        Ham(xn*(yn-1)+m,xn*(yn-1)+m-up) = -t
    ! (2,2)
        Ham(len2+m,m+up+len2) = -t
        Ham(len2+m,xn*(yn-1)+m+len2) = -t*phi(m*yn)*bcy
        Ham(len2+xn*(yn-1)+m,m+len2) = -t*conjg(phi(m*yn))*bcy
        Ham(len2+xn*(yn-1)+m,xn*(yn-1)+m-up+len2) = -t
    end do
    !=============== SOC =====================================
    !------------- X ---------------------------
    do l = 1,yn
        do m = 2+(l-1)*xn,xn*l-1 
            ! (1,2)  
            Ham(m,m+1+len2) = im*lam*phi(l)
            Ham(m,m-1+len2) = -im*lam*conjg(phi(l))
            ! (2,1)   
            Ham(len2+m,m+1) = im*lam*phi(l)
            Ham(len2+m,m-1) = -im*lam*conjg(phi(l))
        end do
    end do
    !========== X boundry ==========================
    do m = 1,yn
    ! (1,2)   
        Ham(m*xn,m*xn-(xn-1)+len2) = im*lam*phi(m)*bcx  
        Ham(m*xn,m*xn-1+len2) = -im*lam*conjg(phi(m))        
        Ham(1+(m-1)*xn,1+(m-1)*xn+1+len2) = im*lam*phi(m) 
        Ham(1+(m-1)*xn,m*xn+len2) = -im*lam*conjg(phi(m))*bcx  
    ! (2,1)		
        Ham(len2+m*xn,m*xn-(xn-1)) = im*lam*phi(m)*bcx
        Ham(len2+m*xn,m*xn-1) = -im*lam*conjg(phi(m))
        Ham(len2+1+(m-1)*xn,1+(m-1)*xn+1) = im*lam*phi(m)
        Ham(len2+1+(m-1)*xn,m*xn) = -im*lam*conjg(phi(m))*bcx
    end do
    ! ---------------------------------------------
    do m = xn+1,xn*(yn-1)
    !(1,2) 
        Ham(m,m+up+len2) = lam
        Ham(m,m-up+len2) = -lam
    !(2,1)		
        Ham(len2+m,m+up) = -lam
        Ham(len2+m,m-up) = lam
    end do
    !---------------------------------------------------
    do m = 1,xn
    ! (1,2)
        Ham(m,m+up+len2) = lam
        Ham(m,xn*(yn-1)+m+len2) = -lam*phi(m*yn)*bcy 
        Ham(xn*(yn-1)+m,m+len2) = lam*conjg(phi(m*yn))*bcy 
        Ham(xn*(yn-1)+m,xn*(yn-1)+m-up+len2) = -lam
    ! (2,1)
        Ham(len2+m,m+up) = -lam
        Ham(len2+m,xn*(yn-1)+m) = lam*phi(m*yn)*bcy
        Ham(len2+xn*(yn-1)+m,m) = -lam*conjg(phi(m*yn))*bcy
        Ham(len2+xn*(yn-1)+m,xn*(yn-1)+m-up) = lam
    end do    
    !=================Pairing====================================
    do m = 1,yn
        do l = 1,xn
            i = (m-1)*xn + l
            !   (1,4)
            Ham(i,len2*3+i) = del(m,l)
            !   (2,3)
            Ham(len2+i,len2*2+i) = -del(m,l)
            !   (4,1)
            Ham(len2*3+i,i) = conjg(del(m,l))
            !   (3,2)
            Ham(len2*2+i,len2+i) = -conjg(del(m,l))
        end do
    enddo
    !--------------------------------------------------------------
    do m = 1,len2
        do l =1,len2
            !   (3,3)
            Ham(len2*2+m,len2*2+l) = -conjg(Ham(m,l))
            !   (4,4)
            Ham(len2*3+m,len2*3+l) = -conjg(Ham(len2+m,len2+l))
            !   (3,4)
            Ham(len2*2+m,len2*3+l) = -conjg(Ham(m,len2+l))
            !   (4,3)
            Ham(len2*3+m,len2*2+l) = -conjg(Ham(len2+m,l))
        enddo
    enddo
    call matrixVerify()
    return
    end subroutine main
!============================================================
    subroutine matrixVerify()
    use pub
    integer i,j
    do i = 1,N
        do j = 1,N
            if (Ham(i,j) .ne. conjg(Ham(j,i)))then
                open(16,file = 'hermitian.dat')
                write(16,*)i,j
                write(16,*)Ham(i,j)
                write(16,*)Ham(j,i)
                write(16,*)"===================="
                close(16)
                stop
            end if
        end do
    end do
    return
    end subroutine matrixVerify
    !================================================================
    subroutine eigsol()
    use pub
    integer m
    lwork = -1
    liwork = -1
    lrwork = -1
    call cheevd('V','U',N,Ham,lda,w,work,lwork,rwork,lrwork,iwork,liwork,info)
    lwork = min(2*N+N**2, int( work( 1 ) ) )
    lrwork = min(1+5*N+2*N**2, int( rwork( 1 ) ) )
    liwork = min(3+5*N, iwork( 1 ) )
    call cheevd('V','U',N,Ham,lda,w,work,lwork,rwork,lrwork,iwork,liwork,info)
    Ham_diag = Ham
    if( info .GT. 0 ) then
        open(11,file="mes.tt",status="unknown")
        write(11,*)'The algorithm failed to compute eigenvalues.'
        close(11)
    end if
    open(12,file="eigval.dat",status="unknown")
    do m = 1,N
        write(12,*)m,w(m)
    end do
    close(12)
    return
    end subroutine eigsol
