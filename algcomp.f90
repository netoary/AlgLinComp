! FALTA
! prepare as mesmas para indicar casos onde as decomposições não são possíveis
! na decomposicaoLU e decomposicaoCholesky
include 'methods.f90'

program algcomp
    integer :: n, m, i, j
    real :: det
    real, allocatable, dimension(:) :: b, x
    real, allocatable, dimension(:,:) :: a, d, c, at, aInversa

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
    x = 0.0
    allocate(aInversa(n,n))

    !call eliminacaoGauss(a, b, x, n)
    !call eliminacaoGaussJordan(a, b, x, n)
    !call decomposicaoLU(a, b, x, n)
    !call decomposicaoCholesky(a, b, x, n)
    !call interativoJacobi(a, b, x, n)
    call GaussSeidel(a, b, x, n)
    !call inversa(a, aInversa, n)
    !call determinante(a, n, det)

    allocate(at(n,n))
    call transposta(a, at, n)

    open(2, file='RESUL.txt', status='replace')
    do i=1, n
        write(2,*) (a(i, j), j=1, n)
    end do
    write(2, *) (x(i), i=1, n)
    write(2, *) det
    close(2)

end program algcomp


