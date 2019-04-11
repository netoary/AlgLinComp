include 'methods.f90'

program questao2
    integer n, i, j
    real, allocatable, dimension(:) :: b, x, y
    real, allocatable, dimension(:,:) :: a, d, c
    real lambda, beta

    beta = 0.6
    ! abertura do arquivo sistema.txt e leitura
    open (1, file='sistema_2.txt', status='old', action='read')
    read(1,*) n
    allocate(x(n))
    read(1,*) ( x(j) , j= 1 , n)

    allocate(y(n))
    read(1,*) ( y(i) , i= 1 , n)
    close(1)

    y = y * beta
    allocate(b(2))
    call minimosQuadrados(b, x, y, n)

    open(2, file='Q2.txt', status='replace')
    write(2, *) "X"
    write(2,*) (x(j), j=1, n)
    write(2, *) "Y"
    write(2, *) (y(i), i=1, n)
    write(2, *) "Coeficiente a"
    write(2, *) (b(1))
    write(2, *) "Coeficiente b"
    write(2, *) (b(2))
    close(2)

end program questao2