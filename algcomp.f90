program algcomp
    integer n, i, j
    real, allocatable, dimension(:) :: b, x, vamos
    real, allocatable, dimension(:,:) :: a, d, c, at

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

    call teste(a, d, c, n)

    !call eliminacaoGauss(a, b, x, n)
    !call eliminacaoGaussJordan(a, b, x, n)
    call decomposicaoCholesky(a, b, x, n)

    allocate(at(n,n))
    call transposta(a, at, n)


    open(2, file='RESUL.txt', status='replace')
    do i=1, n
        write(2,*) (a(i, j), j=1, n)
    end do
    write(2, *) (x(i), i=1, n)
    close(2)

end program algcomp

subroutine retroSubs(a, b, x, n)
    integer :: n, i, j
    real :: a(n,n), b(n), x(n), soma

    x = 0.0

    do i = n, 1, -1
        soma = 0.0
        do j = i+1, n
            soma = soma + a(i,j) * x(j)
        end do
        x(i) = (b(i) - soma)/a(i,i)
    end do
end subroutine

subroutine eliminacaoGauss(a, b, x, n)
    integer :: n, i, j, k, pivo, r, aux
    integer :: p(n)
    real :: a(n,n), b(n), x(n), soma, mult

    x = 0.0

    do i=1, n
        p(i) = i
    end do

    do k = 1, (n-1)
        pivo = a(k,k)
        r = k
        if (pivo == 0) then
            do i= k+1, n
                if(a(i,k)==0) then
                    if (i == n) then
                        open(2, file='RESUL.txt', status='replace')
                            write(2,*) 'A matriz A é singular.'
                        close(2)
                        stop 'A matriz A é singular.'
                    end if
                else
                    r = i
                end if
            end do
        end if

        if (r /= k) then
            aux = p(k)
            p(k) = p(r)
            p(r) = aux
            do j=1, n !troca linha r por linha k
                aux = a(k,j)
                a(k,j) = a(r,j)
                a(r,j) = aux
            end do
        end if

        do i=(k+1), n !eliminação ate chegar numa matriz triangular superior
            mult = a(i,k)/a(k,k)
            a(i,k) = a(i,k) - mult * a(k,k)
            do j=(k+1), n
                a(i,j) = a(i,j) - mult*a(k,j)
            end do
            b(i) = b(i) - mult * b(k)
        end do
    end do

    call retroSubs(a,b, x,n)

end subroutine

subroutine eliminacaoGaussJordan(a, b, x, n)
    integer :: n, i, j, k
    real :: a(n,n), b(n), x(n), mult

    x = 0.0

    call eliminacaoGauss(a, b, x, n)

    do k = n, 2, -1

        do i=(k-1), 1, -1 !eliminação ate chegar numa matriz diagonal
            mult = a(i,k)/a(k,k)
            a(i,k) = a(i,k) - mult * a(k,k)
            do j=(k-1), 1, -1
                a(i,j) = a(i,j) - mult*a(k,j)
            end do
            b(i) = b(i) - mult * b(k)
        end do
    end do

    do i=1, n
        if (a(i,i) /= 0.0) then
            x(i) = b(i)/a(i,i)
        end if
    end do

end subroutine

subroutine decomposicaoCholesky(a, b,x, n)
    integer :: n, i, j, k
    real :: a(n,n), b(n), x(n), l(n,n), soma, y(n), u(n,n)

    x = 0.0
    y = 0.0

    do i = 1, n
        soma = 0.0
        do k = 1, i-1
            soma = soma + l(i,k)**2
        end do
        l(i,i) = sqrt(a(i,i) - soma)
        do j = (i+1), n
            soma = 0.0
            do k = 1, i-1
                soma = soma + l(i,k)*l(j,k)
            end do
            l(j,i) = (a(i,j) - soma)/l(i,i)
        end do
    end do

    do i=1, n !Resolução do sistema Ly=b
        soma = 0.0
        do j=1, (i-1)
            soma = soma + l(i,j)*y(j)
        end do
        y(i) = (b(i) - soma)/l(i,i)
    end do

    call transposta(l, u, n)

    do i=n, 1, -1 !Resolução do sistema Ux=y
        soma = 0.0
        do j=(i+1), n
            soma = soma + u(i,j)*x(j)
        end do
        x(i) = (y(i)-soma)/u(i,i)
    end do

end subroutine

subroutine multiplicaMatrizes(a, b, c, n)
    !multiplicação de matrizes quadradas, precisa generalizar?
    integer :: n, i, j, k
    real :: a(n,n), b(n,n), c(n,n)

    c = 0.0

    do i=1, n
        do j=1, n
            do k=1, n
                c(i,j) = c(i,j) + a(i,k)*b(k,j)
            end do
        end do
    end do
end subroutine

subroutine multiplicaVetor(a, b, c, n)
    !multiplicação de vetores de tamanho n
    integer :: n, i, j, k, m
    real :: a(n), b(n), c

    c = 0.0
    do j=1, n
        c = c + a(j)*b(j)
    end do

end subroutine

subroutine transposta(a, at, n)
    integer :: n
    real :: a(n,n), at(n,n)

    do i = 1, n
        at(i,i) = a(i,i)
        do j = (i+1), n
            at(i,j) = a(j,i)
            at(j,i) = a(i,j)
        end do
    end do
end subroutine

subroutine teste(a, b, c, n)
    !multiplicação de matrizes quadradas, precisa generalizar?
    integer :: n, i, j, k
    real :: a(n,n), b(n,n), c(n,n)

    call multiplicaMatrizes(a, b, c, n)
end subroutine