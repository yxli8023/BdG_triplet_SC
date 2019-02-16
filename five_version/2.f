      program ex01
      implicit none
      real a
      integer i,j
      do i = 1,48
            do j = 1,48*i
                  write(*,*)48*48+j+1
            end do
      end do
      stop
      end
