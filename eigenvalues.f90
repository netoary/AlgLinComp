include 'methods.f90'

program eigenvalues
implicit none
    integer :: n, m, i, j
    real, allocatable, dimension(:) :: x, lambda
    real, allocatable, dimension(:,:) :: a
    real :: tolerance

    open (1, file='sistema_2_1.txt', status='old', action='read')

    read(1,*) n, m
    allocate(a(n,m))

    do i=1, m
        read(1,*) ( a(i,j) , j= 1 , n)
    end do

    read(1, *) tolerance

    close(1)


    allocate(x(n))
    allocate(lambda(n))
    lambda = 0
    call powerMethod(a, x, n, tolerance, lambda(1))

    open(2, file='RESUL.txt', status='replace')
    do i=1, n
        write(2,*) (a(i, j), j=1, n)
    end do
    write(2, *) (x(i), i=1, n)
    write(2, *) (lambda(i), i=1, n)
    write(2, *) tolerance
    close(2)

    call sleep(10000)
    
end program eigenvalues