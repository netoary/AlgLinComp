include 'methods.f90'

program algcomp
    integer :: n, m, i, j, k
    real :: det, tolerancia, alfa
    real, allocatable, dimension(:) :: b, x
    real, allocatable, dimension(:,:) :: a, d, c, at, aInversa, aOriginal

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

    allocate(d(n,n))
    d(1,1) = 1.0
    d(1,2) = 0.0
    d(1,3) = 0.0
    d(2,1) = 5.0
    d(2,2) = 1.0
    d(2,3) = 0.0
    d(3,1) = 0.0
    d(3,2) = 0.0
    d(3,3) = 1.0
    allocate(c(n,n))
    allocate(x(n))
    x = 0.0
    allocate(aInversa(n,n))
    allocate(aOriginal(n,n))
    aOriginal = a

    alfa = 0.18
    a = a * alfa
    tolerancia = 0.00001
    !call decomposicaoLU(a, b, x, n)
    !call decomposicaoCholesky(a, b, x, n)
    !call GaussSeidel(a, b, x, n, tolerancia, k)
    call determinante(a, n, det)

    !call eliminacaoGaussJordan(a, b, x, n)
    !call interativoJacobi(a, b, x, n, tolerancia, k)
    !call inversa(a, aInversa, n)

    allocate(at(n,n))
    call transposta(a, at, n)

    open(2, file='RESUL_' // fileName //'.txt', status='replace')
    ! write(2, *) "A"
    ! do i=1, n
    !     write(2,*) (a(i, j), j=1, n)
    ! end do
    ! write(2, *) "A^-1"
    ! do i=1, n
    !     write(2,*) (aInversa(j,i), j=1, n)
    ! end do
    ! write(2, *) "Gauss Seidel"
    ! write(2, *) "X"
    ! write(2, *) (x(i), i=1, n)
    ! write(2, *) "Tolerância"
    ! write(2, *) tolerancia
    ! write(2, *) "Número de iterações"
    ! write(2, *) k
    write(2, *) "Determinante"
    write(2, *) det
    close(2)

    print*, "x = ", x

    call sleep(20000)
end program algcomp


