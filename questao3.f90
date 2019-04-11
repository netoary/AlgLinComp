include 'methods.f90'

program questao3
    integer :: n, m, i, j, k
    real :: det, tolerancia
    real, allocatable, dimension(:) :: b, x
    real, allocatable, dimension(:,:) :: a, aOriginal, at

    Character(len = 100) :: fileName


    print*,"Digite o nome do arquivo com o sistema a ser resolvido:"
    read(*,*) fileName
    open (1, file=fileName, status='old', action='read')
    read(1,*) n, m
    allocate(a(n,m))
    do i=1, m
        read(1,*) ( a(i, j) , j= 1 , n)
    end do
    allocate(b(n))
    read(1,*) ( b(i) , i= 1 , n)
    close(1)
    allocate(x(n))
    allocate(aOriginal(n, n))
    x = 0.0
    
    tolerancia = 0.00001
    aOriginal = a
    call interativoJacobi(a, b, x, n, tolerancia, k)

    open(2, file='Q3.txt', status='replace')
    write(2, *) "A"
    do i=1, n
        write(2,*) (aOriginal(i, j), j=1, n)
    end do
    write(2, *) "Jacobi interativo"
    write(2, *) "X"
    write(2, *) (x(i), i=1, n)
    write(2, *) "Tolerância"
    write(2, *) tolerancia
    write(2, *) "Número de iterações"
    write(2, *) k
    close(2)

    print*, "x = ", x

    call sleep(20000)
end program questao3


