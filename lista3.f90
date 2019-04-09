include 'methods.f90'

program lista3
    integer n, i, j
    real, allocatable, dimension(:) :: b, x, y
    real, allocatable, dimension(:,:) :: a, d, c
    real lambda

    ! abertura do arquivo sistema.txt e leitura
    open (1, file='sistemaM.txt', status='old', action='read')
    read(1,*) n
    allocate(x(n))
    read(1,*) ( x(j) , j= 1 , n)

    allocate(y(n))
    read(1,*) ( y(i) , i= 1 , n)
    close(1)

    allocate(b(2))
    call minimosQuadrados(b, x, y, n)

    open(2, file='RESUL.txt', status='replace')
    write(2,*) (x(j), j=1, n)
    write(2, *) (y(i), i=1, n)
    write(2, *) (b(i), i=1, 2)
    close(2)

end program lista3