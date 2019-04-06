include 'methods.f90'

program eigenvalues
implicit none
    integer :: n, m, i, j
    real, allocatable, dimension(:) :: lambda
    real, allocatable, dimension(:,:) :: a, x
    real :: tolerance

    open (1, file='sistema_2_2.txt', status='old', action='read')

    read(1,*) n, m
    allocate(a(n,m))

    do i=1, m
        read(1,*) ( a(i,j) , j= 1 , n)
    end do

    read(1, *) tolerance

    close(1)


    allocate(x(n, n))
    allocate(lambda(n))

    x = 0
    lambda = 0
    !call powerMethod(a, x, n, tolerance, lambda(1))
    call jocobi(a, x, n, tolerance, lambda)

    open(2, file='RESUL.txt', status='replace')
    do i=1, n
        write(2,*) (a(i, j), j=1, n)
    end do
    do i=1, n
        write(2,*) (x(i, j), j=1, n)
    end do
    write(2, *) (lambda(i), i=1, n)
    write(2, *) tolerance
    close(2)

    call sleep(10000)
    
end program eigenvalues