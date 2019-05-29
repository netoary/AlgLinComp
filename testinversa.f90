include 'methods.f90'
include 'functions.f90'


program testinversa
    use functions
    IMPLICIT NONE
    real :: x1, x2, x3, deltaX, deltaX1, deltaX2, p
    real, allocatable, dimension(:,:) :: a, a1, a2
    integer :: i,j

    allocate(a(2,2))
    allocate(a1(2,2))
    allocate(a2(2,2))

    a(1,1) = 1
    a(2,1) = 3
    a(1,2) = 2
    a(2,2) = 2


    call inversa(a, a1, 2)

    a(1,1) = 1
    a(2,1) = 3
    a(1,2) = 2
    a(2,2) = 2
    call inversaALT(a, a2, 2)
    do i=1, 2
        do j=1, 2
            write (*,*) a(i,j)
        end do
    end do
    write (*,*) "a1"
    write (*,*) a1
    write (*,*) "a2"
    write (*,*) a2

end program testinversa