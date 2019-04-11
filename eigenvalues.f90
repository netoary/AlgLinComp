include 'methods.f90'

program eigenvalues
implicit none
    integer :: n, m, i, j, k
    real, allocatable, dimension(:) :: lambda
    real, allocatable, dimension(:,:) :: a, x
    real :: tolerance, alfa
    Character(len = 100) :: fileName

    alfa = 0.18
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

    print*, a
    a = a * alfa

    print*, a
    allocate(x(n, n))
    allocate(lambda(n))

    x = 0
    lambda = 0
    !call powerMethod(a, x, n, tolerance, lambda(1), k)
    call jacobi(a, x, n, tolerance, lambda, k)

    print*, "lambda = ", lambda
    print*, "x = ", x

    open(2, file='RESUL_' // fileName //'.txt', status='old')
    write(2, *)
    write(2, *) "Questão 1 - b"
    write(2, *) "A"
    do i=1, n
        write(2,*) (a(i, j), j=1, n)
    end do
    write(2, *) "Auto-vetor do maior auto-valor"
    do i=1, n
         write(2,*) (x(i, j), j=1, n)
    end do
    write(2,*) (x(j, 1), j=1, n)
    write(2, *) "Todos os auto-valores"
    write(2, *) (lambda(i), i=1, n)
    !write(2, *) (lambda(1))
    write(2, *) "Tolerância"
    write(2, *) tolerance
    write(2, *) "Numero de iterações"
    write(2, *) k
    close(2)

    call sleep(10000)

end program eigenvalues