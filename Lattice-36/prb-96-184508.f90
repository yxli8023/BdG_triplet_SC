!======================================================================================================================================================
    module code_param
    implicit none
    real,parameter::pi = acos(-1.0)
    complex::im = (0.,1.0)
    integer xn,yn,hn, len2,up,ZeroPoint,bcx,bcy
    real t0,lam,h0,mu,B0,phi0,U0,eps
    parameter(xn = 36,yn = xn,hn = xn * yn * 4 , len2 = xn * yn,up = xn,ZeroPoint = xn * yn * 2, eps = 1e-5 , U0 = 5.0)
    complex Ham(hn,hn),Ham_diag(hn,hn)
    complex del(xn,yn)  
    !-----------------LAPACK PACKAGE PARAM
    integer::lda = hn
    integer,parameter::lwmax = 2*hn + hn**2
    real,allocatable::matval(:)
    complex,allocatable::work(:)
    real,allocatable::rwork(:)
    integer,allocatable::iwork(:)
    integer lwork   ! at least 2*N+N**2
    integer lrwork    ! at least 1 + 5*N +2*N**2
    integer liwork   ! at least 3 +5*N
    integer info
    end module code_param
!======================================================================================================================================================
    program sol
    use code_param
    !================ Physics memory allocate =================
    allocate(matval(hn))
    allocate(work(lwmax))
    allocate(rwork(1 + 5*hn + 2*hn**2))
    allocate(iwork(3 + 5*hn))
    !================ Physics memory allocate =================
    bcx = 1
    bcy = 1
    !----------------------
    t0 = 1.0
    lam = 0.5
    h0 = 0.6
    mu = -4
    del = 0.5
    phi0 = pi
    B0 = 2.0*phi0/len2  
    call Mat_Set()
    call eigsol()
    call loop()
    call Majorana()
    stop
    end program sol
!======================================================================================================================================================
    subroutine loop()   
    use code_param
    integer num,ref
    !store value of del and use for self-consistently
    complex del_loop(xn,yn),del_err(xn,yn)
    del_loop = 0
    del_err = 0
    call delc() !   计算出一组 delta
    num = 0    
    ref = 3                    
    do while(ref .gt. 1)
        open(25,file = "check.dat")
        !  迭代次数减少,说明程序正在收敛
        ref = 0
        do m= 1,xn
            do l= 1,yn
                !del_err用来存储两次计算得到的delta的差值,用来进行自恰条件的比较
                del_err(m,l) = del(m,l) - del_loop(m,l)
                del_loop(m,l) = del(m,l)
                if (abs(real(del_err(m,l))) > eps) ref = ref + 1
                if (abs(imag(del_err(m,l))) > eps) ref = ref + 1
            end do
        end do
        write(25,"(10I15.6)")ref
        num = num + 1
        call Mat_Set()  ! Reconstruction Matrix
        call eigsol()
        call delc() ! Get a new value for del use new w and Ham
    end do
    close(25)
    return
    end subroutine
!======================================================================================================================================================
    subroutine delc()
    use code_param
    integer m,l,i
    complex,external::delta
    open(12,file = "order.dat")
    open(14,file = "order_module.dat")
    do m = 1,xn
        do l = 1,yn
            i = (m - 1) * yn + l
            del(m,l) = delta(i,i,matval,Ham_diag)
            write(14,"(10F15.8)")abs(del(m,l))
            write(12,"(10F15.8)")real(del(m,l)),aimag(del(m,l))
        end do
    end do
    close(12)
    close(14)
    return
    end subroutine delc
!======================================================================================================================================================
!  能隙方程自洽
    complex function delta(x,y,eig_val,eig_vec)
    use code_param
    integer x,y,kk
    real eig_val(hn),beta
    complex s,eig_vec(hn,hn)
    beta = 1e5
    s = (0,0)
    do kk = 1,hn 
          s = s + eig_vec(x,kk) * conjg(eig_vec(y + len2 * 3, kk)) * tanh(eig_val(kk)/2 * beta)
    end do
    delta = U0/2.0 * s
    end function delta
!======================================================================================================================================================
subroutine Majorana()
    !  获取零能模波函数绝对值
    use code_param
    complex z1(xn * yn),z2(xn * yn)
    integer m
    z1 = 0
    z2 = 0
    open(11,file = "gamma1.dat")
    open(12,file = "gamma2.dat")
    do m = 1,xn*yn
        z1(m) = (Ham(ZeroPoint,m) + Ham(ZeroPoint + 1 , len2 + m))/sqrt(2.0)
        z2(m) = im * (Ham(ZeroPoint,m) - Ham(ZeroPoint + 1 , len2 + m))/sqrt(2.0)  
        write(11,"(10F15.8)")abs(z1(m))
        write(12,"(10F15.8)")abs(z2(m)) 
    end do 
    close(11)
    close(12)
    return
end subroutine Majorana
!======================================================================================================================================================
complex function phi(iy)  ! Landau guage   A = (-By,0,0)
    use code_param
    integer iy ! position of lattice
    phi = cexp(-im * B0 * iy)       
end function phi
!======================================================================================================================================================
subroutine Mat_Set()
    use code_param
    complex,external::phi
    integer m,l,i
    Ham = 0
    do i = 1 , len2
        Ham(i,i) = h0 - mu
        Ham(len2 + i , len2 + i) = -h0 - mu
    end do
    !------------- X ---------------------------
    do l = 1,yn
        do m = 2+(l-1)*xn,xn*l-1 
            ! (1,1)  
            Ham(m,m + 1) = -t0 *phi(l)
            Ham(m,m-1) = -t0 *conjg(phi(l))
            ! (2,2)   
            Ham(len2 + m , len2 + m + 1) = -t0 *phi(l)
            Ham(len2 + m , len2 + m-1) = -t0 *conjg(phi(l))
        end do
    end do
    !========== X boundry ==========================
    do m = 1,yn
    ! (1,1)   
        Ham(m*xn,m*xn-(xn-1)) = -t0 *phi(m) * bcx  
        Ham(m*xn,m*xn-1) = -t0 *conjg(phi(m))        
        Ham(1+(m-1)*xn,1+(m-1)*xn + 1) = -t0 *phi(m) 
        Ham(1+(m-1)*xn,m*xn) = -t0 *conjg(phi(m)) * bcx  
    ! (2,2)		
        Ham(len2 + m*xn,m*xn-(xn-1) + len2) = -t0 *phi(m) * bcx
        Ham(len2 + m*xn,m*xn-1 + len2) = -t0 *conjg(phi(m))
        Ham(len2 + 1+(m-1)*xn,1+(m-1)*xn + 1 + len2) = -t0 *phi(m)
        Ham(len2 + 1+(m-1)*xn,m*xn + len2) = -t0 *conjg(phi(m)) * bcx
    end do
    ! ---------------------------------------------
    do m = xn + 1,xn * (yn-1)
    !(1,1) 
        Ham(m,m + up) = -t0
        Ham(m,m - up) = -t0
    !(2,2)		
        Ham(len2 + m , len2 + m + up) = -t0
        Ham(len2 + m , len2 + m - up) = -t0
    end do
    !---------------------------------------------------
    do m = 1,xn
    ! (1,1)
        Ham(m,m + up) = -t0
        Ham(m,xn * (yn-1) + m) = -t0 *phi(m*yn) * bcy 
        Ham(xn * (yn-1) + m,m) = -t0 *conjg(phi(m*yn)) * bcy 
        Ham(xn * (yn-1) + m,xn * (yn-1) + m - up) = -t0
    ! (2,2)
        Ham(len2 + m,m + up + len2) = -t0
        Ham(len2 + m,xn * (yn-1) + m + len2) = -t0 *phi(m*yn) * bcy
        Ham(len2+xn * (yn-1) + m,m + len2) = -t0 *conjg(phi(m*yn)) * bcy
        Ham(len2+xn * (yn-1) + m,xn * (yn-1) + m - up + len2) = -t0
    end do
    !=============== SOC =====================================
    !------------- X ---------------------------
    do l = 1,yn
        do m = 2+(l-1)*xn,xn*l-1 
            ! (1,2)  
            Ham(m,m + 1 + len2) = im*lam * phi(l)
            Ham(m,m-1 + len2) = -im*lam * conjg(phi(l))
            ! (2,1)   
            Ham(len2 + m,m + 1) = im*lam * phi(l)
            Ham(len2 + m,m-1) = -im*lam * conjg(phi(l))
        end do
    end do
    !========== X boundry ==========================
    do m = 1,yn
    ! (1,2)   
        Ham(m*xn,m*xn-(xn-1) + len2) = im*lam * phi(m) * bcx  
        Ham(m*xn,m*xn-1 + len2) = -im*lam * conjg(phi(m))        
        Ham(1+(m-1)*xn,1+(m-1)*xn + 1 + len2) = im*lam * phi(m) 
        Ham(1+(m-1)*xn,m*xn + len2) = -im*lam * conjg(phi(m)) * bcx  
    ! (2,1)		
        Ham(len2 + m*xn,m*xn-(xn-1)) = im*lam * phi(m) * bcx
        Ham(len2 + m*xn,m*xn-1) = -im*lam * conjg(phi(m))
        Ham(len2 + 1+(m-1)*xn,1+(m-1)*xn + 1) = im*lam * phi(m)
        Ham(len2 + 1+(m-1)*xn,m*xn) = -im*lam * conjg(phi(m)) * bcx
    end do
    ! ---------------------------------------------
    do m = xn + 1,xn * (yn-1)
    !(1,2) 
        Ham(m,m + up + len2) = lam
        Ham(m,m - up + len2) = -lam
    !(2,1)		
        Ham(len2 + m,m + up) = -lam
        Ham(len2 + m,m - up) = lam
    end do
    !---------------------------------------------------
    do m = 1,xn
    ! (1,2)
        Ham(m,m + up + len2) = lam
        Ham(m,xn * (yn-1) + m + len2) = -lam * phi(m*yn) * bcy 
        Ham(xn * (yn-1) + m,m + len2) = lam * conjg(phi(m*yn)) * bcy 
        Ham(xn * (yn-1) + m,xn * (yn-1) + m - up + len2) = -lam
    ! (2,1)
        Ham(len2 + m,m + up) = -lam
        Ham(len2 + m,xn * (yn-1) + m) = lam * phi(m*yn) * bcy
        Ham(len2+xn * (yn-1) + m,m) = -lam * conjg(phi(m*yn)) * bcy
        Ham(len2+xn * (yn-1) + m,xn * (yn-1) + m - up) = lam
    end do    
    !=================Pairing====================================
    do m = 1,yn
        do l = 1,xn
            i = (m-1)*xn + l
            !   (1,4)
            Ham(i , len2 * 3 + i) = del(m,l)
            !   (2,3)
            Ham(len2 + i , len2 * 2 + i) = -del(m,l)
            !   (4,1)
            Ham(len2 * 3 + i,i) = conjg(del(m,l))
            !   (3,2)
            Ham(len2 * 2 + i , len2 + i) = -conjg(del(m,l))
        end do
    enddo
    !--------------------------------------------------------------
    do m = 1 , len2
        do l = 1 , len2
            !   (3,3)
            Ham(len2 * 2 + m , len2 * 2 + l) = -conjg(Ham(m,l))
            !   (4,4)
            Ham(len2 * 3 + m , len2 * 3 + l) = -conjg(Ham(len2 + m , len2 + l))
            !   (3,4)
            Ham(len2 * 2 + m , len2 * 3 + l) = -conjg(Ham(m , len2 + l))
            !   (4,3)
            Ham(len2 * 3 + m , len2 * 2 + l) = -conjg(Ham(len2 + m,l))
        enddo
    enddo
    !--------------------------------------------------------------
    call Is_Hermatian()  
    return
end subroutine Mat_Set
!======================================================================================================================================================
subroutine Is_Hermatian()
    use code_param
    integer i,j
    do i = 1,hn
        do j = 1,hn
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
end subroutine Is_Hermatian
!======================================================================================================================================================
subroutine eigsol()
    use code_param
    integer m
    lwork = -1
    liwork = -1
    lrwork = -1
    call cheevd('V','U',hn,Ham,lda,matval,work,lwork,rwork,lrwork,iwork,liwork,info)
    lwork = min(2*hn+hn**2, int( work( 1 ) ) )
    lrwork = min(1+5*hn+2*hn**2, int( rwork( 1 ) ) )
    liwork = min(3+5*hn, iwork( 1 ) )
    call cheevd('V','U',hn,Ham,lda,matval,work,lwork,rwork,lrwork,iwork,liwork,info)
    Ham_diag = Ham
    if( info .GT. 0 ) then
        open(11,file="mes.tt",status="unknown")
        write(11,*)'The algorithm failed to compute eigenvalues.'
        close(11)
    end if
    open(12,file="eigval.dat",status="unknown")
    do m = 1,hn
        write(12,"(I10,10F15.8)")m,matval(m)
    end do
    close(12)
    return
end subroutine eigsol
