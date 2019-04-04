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

subroutine decomposicaoCholesky(a, b,x, n)
    integer :: n, i, j, k
    real :: a(n,n), a_transposta(n,n), b(n), x(n), l(n,n), soma, y(n), u(n,n), a_soma

    x = 0.0
    y = 0.0

    call transposta(a, a_transposta, n)

    if (.not.all(a.EQ.a_transposta)) then
        open(2, file='RESUL.txt', status='replace')
        write(2,*) 'A matriz A não é simétrica.'
        close(2)
        stop 'A matriz A não é simétrica.'
    end if

    do i = 1, n
        soma = 0.0
        do k = 1, i-1
            soma = soma + l(i,k)**2
        end do
        a_soma = a(i,i) - soma
        if (a_soma < 0) then
            open(2, file='RESUL.txt', status='replace')
            write(2,*) 'A matriz A não é positiva definida.'
            close(2)
            stop 'A matriz A não é positiva definida.'
        end if
        l(i,i) = sqrt(a_soma)
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

subroutine multiplicaVetorMatriz(a, b, c, n)
    !multiplicação de vetores de tamanho n
    integer :: n, i, j
    real :: a(n, n), b(n), c(n)
    
    c = 0
    do i=1, n
        do j=1, n
            c(i) = c(i) + (a(i,j) * b(j))
        end do
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

subroutine powerMethod(a, x, n, tol, lambda)
    integer :: n, i, k
    real :: a(n,n), x(n), y(n), lambda, lastLambda, r, tol

    r = tol + 1
    x = 1.0
    lastLambda = 1
    k = 1

    do while(r >= tol)
        call multiplicaVetorMatriz(a, x, y, n)
        !y = MATMUL(a, x)
        
        lambda = y(1)
        do i = 2, n
            x(i) = y(i) / lambda
        end do
        r = abs(lambda - lastLambda) / abs(lambda)
        lastLambda = lambda

        print*, "x(", k, ") = ", x
        print*, "lambda = ", lambda
        print*, "r = ", r
        print*,""

        k = k + 1
        call sleep(1)
    end do
end subroutine


subroutine jocobi
end subroutine

subroutine teste(a, b, c, n)
    !multiplicação de matrizes quadradas, precisa generalizar?
    integer :: n, i, j, k
    real :: a(n,n), b(n,n), c(n,n)

    call multiplicaMatrizes(a, b, c, n)
end subroutine