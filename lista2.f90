include 'methods.f90'

program lista2
    integer n, m, i, j
    real, allocatable, dimension(:) :: b, x
    real, allocatable, dimension(:,:) :: a, d, c
    real lambda

    ! abertura do arquivo sistema.txt e leitura
    open (1, file='sistema.txt', status='old', action='read')
    read(1,*) n, m
    allocate(a(n,m))
    do i=1, m
        read(1,*) ( a(i,j) , j= 1 , n)
    end do
    allocate(b(n))
    read(1,*) ( b(i) , i= 1 , n)
    close(1)


    allocate(d(n,n))
    d(1,1) = 1.0
    d(1,2) = 0.0
    d(1,3) = 0.0
    d(2,1) = 0.0
    d(2,2) = 1.0
    d(2,3) = 0.0
    d(3,1) = 0.0
    d(3,2) = 0.0
    d(3,3) = 1.0
    allocate(c(n,n))
    allocate(x(n))

    call powerMethod(a, x, n, 0.001, lambda)

    open(2, file='RESUL.txt', status='replace')
    do i=1, n
        write(2,*) (a(i, j), j=1, n)
    end do
    write(2, *) (x(i), i=1, n)
    write(2, *) (lambda)
    close(2)

end program lista2