include 'methods.f90'

program eigenvalues
implicit none
    integer :: n, m, i, j
    real, allocatable, dimension(:) :: lambda
    real, allocatable, dimension(:,:) :: a, x
    real :: tolerance
    Character(len = 100) :: fileName


    print*,"Digite o nome do arquivo com o sistema a ser resolvido:"
    read(*,*) fileName
    open (1, file=fileName, status='old', action='read')

    read(1,*) n, m
    allocate(a(n,m))

    do i=1, m
        read(1,*) ( a(j,i) , j= 1 , n)
    end do

    read(1, *) tolerance

    close(1)


    allocate(x(n, n))
    allocate(lambda(n))

    x = 0
    lambda = 0
    call powerMethod(a, x, n, tolerance, lambda(1))
    !call jacobi(a, x, n, tolerance, lambda)

    print*, "lambda = ", lambda
    print*, "x = ", x

    open(2, file='RESUL_' // fileName //'.txt', status='replace')
    write(2, *) "A"
    do i=1, n
        write(2,*) (a(i, j), j=1, n)
    end do
    write(2, *) "X"
    do i=1, n
        write(2,*) (x(i, j), j=1, n)
    end do
    write(2, *) "lambda"
    write(2, *) (lambda(i), i=1, n)
    write(2, *) "tolerância"
    write(2, *) tolerance
    close(2)

    call sleep(10000)

end program eigenvalues