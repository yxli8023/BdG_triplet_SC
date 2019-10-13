program ex01
    implicit none
    integer xn,yn
    integer i,x,y
    parameter(xn=48,yn=48)
    real a(xn*yn)
    open(11,file="order_module.dat")
    do i = 1,xn*yn
        read(11,*)a(i)
        if (a(i) .lt. 0.1) then
        x = i/xn 
        y = i - x*xn
        write(*,*)"The energy zero mode in:"
        write(*,*)x+1,y
        end if
    end do
    close(11)
    stop 
end program ex01
