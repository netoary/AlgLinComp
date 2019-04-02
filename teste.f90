program teste
    implicit none
    integer :: n, m, i, j
    real :: det
    real, allocatable, dimension(:) :: b, x, x1, x2, x3, x4, x5, x6
    real, allocatable, dimension(:,:) :: a, aOriginal, aInversa

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

    allocate(aOriginal(n,n))
    aOriginal = a

    allocate(x(n))
    x = 0.0
    allocate(x1(n))
    x1 = 0.0
    allocate(x2(n))
    x2 = 0.0
    allocate(x3(n))
    x3 = 0.0
    allocate(x4(n))
    x4 = 0.0
    allocate(x5(n))
    x5 = 0.0
    allocate(x6(n))
    x6 = 0.0
    allocate(aInversa(n,n))

    call eliminacaoGauss(a, b, x1, n)
    call eliminacaoGaussJordan(a, b, x2, n)
    call decomposicaoLU(a, b, x3, n)
    call decomposicaoCholesky(a, b, x4, n)
    call jacobi(a, b, x5, n)
    call GaussSeidel(a, b, x6, n)
    call inversa(a, aInversa, n)
    call determinante(a, n, det)

    open(2, file='RESUL.txt', status='replace')
    do i=1, n
        write(2,*) (aOriginal(i, j), j=1, n)
    end do
    write(2, *)
    write(2, *) "eliminacaoGauss"
    write(2, *) (x1(i), i=1, n)
    write(2, *)
    write(2, *) "eliminacaoGaussJordan"
    write(2, *) (x2(i), i=1, n)
    write(2, *)
    write(2, *) "decomposicaoLU"
    write(2, *) (x3(i), i=1, n)
    write(2, *)
    write(2, *) "decomposicaoCholesky"
    write(2, *) (x4(i), i=1, n)
    write(2, *)
    write(2, *) "jacobi"
    write(2, *) (x5(i), i=1, n)
    write(2, *)
    write(2, *) "GaussSeidel"
    write(2, *) (x6(i), i=1, n)
    write(2, *)

    write(2, *) det
    close(2)


end program teste

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

subroutine decomposicaoLU(a, b, x, n)
    integer :: n, i, j, k, pivo, r, aux
    integer :: p(n)
    real :: a(n,n), b(n), x(n), soma, mult, y(n)


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

        do i=(k+1), n !eliminação
            mult = a(i,k)/a(k,k)
            a(i,k) = mult
            do j=(k+1), n
                a(i,j) = a(i,j) - mult*a(k,j)
            end do
        end do
    end do

     do i=1, n !Resolução do sistema Ly=c
        soma = 0
        do j=1, (i-1)
            soma = soma + a(i,j)*y(j)
        end do
        y(i) = b(i) - soma
    end do

    do i=n, 1, -1 !Resolução do sistema Ux=y
        soma = 0
        do j=(i+1), n
            soma = soma + a(i,j)*x(j)
        end do
        x(i) = (y(i)-soma)/a(i,i)
    end do

end subroutine

subroutine determinante(a, n, det)
    integer :: n, i
    real :: a(n,n), b(n), x(n), det

    b = 0.0
    x = 0.0
    det = 1.0
    call decomposicaoLU(a, b, x, n)
    do i=1, n
        det = det * a(i,i)
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

subroutine jacobi(a, b, x, n)
    integer :: n, i, j, k
    real :: a(n,n), b(n), xZero(n), tolerancia, r, xNovo(n), xNovoT, x(n), xT, soma

    k = 1
    r = 100.0
    tolerancia = 10**(-5)
    xZero = 1.0
    !xZero(1) = 1.0

    do while (r > tolerancia .OR. k <= 1000)
        do i=1, n
            soma = 0.0
            do j=1, n
                if (j /= i) then
                    soma = soma + a(i,j)*xZero(j)
                end if
            end do
            xNovo(i) = (b(i) - soma) / a(i,i)
        end do

        xNovoT = 0.0
        call multiplicaVetor(xNovo, xNovo, xNovoT, n)
        xNovoT = sqrt(xNovoT)

        x = xNovo - xZero
        call multiplicaVetor(x, x, xT, n)
        xT = sqrt(xT)

        r = xT / xNovoT
        xZero = xNovo
        k = k + 1
    end do

    x = xZero

end subroutine

subroutine GaussSeidel(a, b, x, n)
    integer :: n, i, j, k
    real :: a(n,n), b(n), xZero(n), tolerancia, r, xNovo(n), xNovoT, x(n), xT, soma

    k = 1
    r = 100.0
    tolerancia = 10**(-5)
    xZero = 1.0
    !xZero(1) = 1.0

    do while (r > tolerancia .OR. k <= 1000)
        do i=1, n
            soma = 0.0
            do j=1, (i-1)
                soma = soma + a(i,j)*xNovo(j)
            end do
            do j=i+1, n
                soma = soma + a(i,j)*xZero(j)
            end do
            xNovo(i) = (b(i) - soma) / a(i,i)
        end do

        xNovoT = 0.0
        call multiplicaVetor(xNovo, xNovo, xNovoT, n)
        xNovoT = sqrt(xNovoT)

        x = xNovo - xZero
        call multiplicaVetor(x, x, xT, n)
        xT = sqrt(xT)

        r = xT / xNovoT
        xZero = xNovo
        k = k + 1
    end do
    x = xZero
end subroutine

subroutine inversa(a, aInversa, n)
    integer :: n, i, j, k, l
    real :: a(n,n), b(n,n), x(n), mult, aInversa(n,n)

    x = 0.0
    aInversa = 0.0

    do i=1, n
        do j=1, n
            if (i==j) then
                b(i,j) = 1.0
            else
                b(i,j) = 0.0
            end if
        end do
    end do

    do k = 1, (n-1)
        do i=(k+1), n !eliminação ate chegar numa matriz triangular superior
            mult = a(i,k)/a(k,k)
            a(i,k) = a(i,k) - mult * a(k,k)
            do j=(k+1), n
                a(i,j) = a(i,j) - mult*a(k,j)
            end do
            do l=1, n
                b(l,i) = b(l,i) - mult * b(l,k)
            end do
        end do
    end do
    do k = n, 2, -1
        do i=(k-1), 1, -1 !eliminação ate chegar numa matriz diagonal
            mult = a(i,k)/a(k,k)
            a(i,k) = a(i,k) - mult * a(k,k)
            do j=(k-1), 1, -1
                a(i,j) = a(i,j) - mult*a(k,j)
            end do
            do l=n, 1, -1
                b(l,i) = b(l,i) - mult * b(l,k)
            end do
        end do
    end do
    do i=1, n
        if (a(i,i) /= 0.0) then
            do j=1, n
                b(j,i) = b(j,i)/a(i,i)
            end do
        end if
    end do

    aInversa = b

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