include 'methods.f90'

program questao1
    integer :: n, m, i, j, k
    real :: det, tolerancia, alfa
    real, allocatable, dimension(:) :: b, x, lambda
    real, allocatable, dimension(:,:) :: a, eigenvalues, at, aInversa, aOriginal

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
    !read(1, *) tolerance
    close(1)

    allocate(eigenvalues(n,m))

    allocate(x(n))
    allocate(aOriginal(n,n))
    allocate(lambda(n))

    x = 0.0
    aOriginal = a

    alfa = 0.18
    a = a * alfa
    aOriginal = a
    tolerancia = 0.00001
    
    open(2, file='Q1.txt', status='replace')
    write(2, *) "A"
    do i=1, n
        write(2,*) (a(i, j), j=1, n)
    end do

    call powerMethod(a, eigenvalues, n, tolerancia, lambda(1), k)
    write(2, *) "Auto-vetor transposto do maior auto-valor"
    write(2,*) (eigenvalues(j, 1), j=1, n)
    write(2, *) "Maior auto-valor"
    write(2, *) (lambda(1))
    write(2, *) "Tolerância"
    write(2, *) tolerancia
    write(2, *) "Numero de iterações"
    write(2, *) k
    write(2, *) "----------------------------------------------"

    a = aOriginal
    call jacobi(a, eigenvalues, n, tolerancia, lambda, k)
    write(2, *) "Auto-vetor do maior auto-valor"
    do i=1, n
         write(2,*) (eigenvalues(i, j), j=1, n)
    end do
    write(2, *) "Todos os auto-valores"
    write(2, *) (lambda(i), i=1, n)
    write(2, *) "Tolerância"
    write(2, *) tolerancia
    write(2, *) "Numero de iterações"
    write(2, *) k
    write(2, *) "----------------------------------------------"

    a = aOriginal
    call decomposicaoLU(a, b, x, n)
    write(2, *) "Decomposição LU"
    write(2, *) "X"
    write(2, *) (x(i), i=1, n)
    write(2, *) "----------------------------------------------"

    a = aOriginal
    call decomposicaoCholesky(a, b, x, n)
    write(2, *) "Decomposição Cholesky"
    write(2, *) "X"
    write(2, *) (x(i), i=1, n)
    write(2, *) "----------------------------------------------"

    a = aOriginal
    call GaussSeidel(a, b, x, n, tolerancia, k)
    write(2, *) "Decomposição Gauss-Seidel"
    write(2, *) "X"
    write(2, *) (x(i), i=1, n)
    write(2, *) "Tolerância"
    write(2, *) tolerancia
    write(2, *) "Numero de iterações"
    write(2, *) k
    write(2, *) "----------------------------------------------"

    a = aOriginal
    call determinante(a, n, det)
    write(2, *) "Determinante"
    write(2, *) det
    write(2, *) "----------------------------------------------"

    allocate(at(n,n))
    call transposta(a, at, n)

    open(2, file='RESUL_' // fileName //'.txt', status='replace')
    write(2, *) "A"
    do i=1, n
        write(2,*) (a(i, j), j=1, n)
    end do
    write(2, *) "A^-1"
    do i=1, n
        write(2,*) (aInversa(j,i), j=1, n)
    end do
    write(2, *) "Gauss Seidel"
    write(2, *) "X"
    write(2, *) (x(i), i=1, n)
    write(2, *) "Tolerância"
    write(2, *) tolerancia
    write(2, *) "Número de iterações"
    write(2, *) k
    write(2, *) "Determinante"
    write(2, *) det
    close(2)
    call sleep(20000)
end program questao1


