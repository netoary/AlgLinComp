subroutine retroSubs(a, b, x, n)
    integer :: n, i, j
    real :: a(n,n), b(n), x(n), soma

    x(n) = b(n) / a(n, n)
    do i = n-1, 1, -1
        soma = 0.0
        do j = i + 1, n
            soma = soma + a(j, i) * x(j)
        end do
        x(i) = (b(i) - soma) / a(i, i)
    end do
end subroutine

subroutine eliminacaoGauss(a, b, x, n)
    integer :: n, i, j, k, pivo
    real :: p(n, n), m(n, n)
    real :: aux(n)
    real :: a(n,n), b(n), x(n), soma, mult

    p = 0
    call identidade(p, n)

    do k = 1, (n-1)

        print*, "A=", a
        print*, "B=", b
        print*,""
        pivo = a(k,k)
        if (pivo == 0) then
            call identidade(p, n)
            do i= k+1, n
                if(a(i,k)==0) then
                    if (i == n) then
                        open(2, file='RESUL.txt', status='replace')
                            write(2,*) 'A matriz A é singular.'
                        close(2)
                        stop 'A matriz A é singular.'
                    end if
                else
                    aux = p(k,:)
                    p(k,:) = p(i,:)
                    p(i,:) = aux
                    exit
                end if
            end do
            print*, "P=", p
            a = matmul(a, p)
            b = matmul(b, p)
            print*, "PA=", a
            print*, "PB=", b
            pivo = a(k,k)
        end if

        call identidade(m, n)

        do i=(k+1), n !eliminação ate chegar numa matriz triangular superior
            m(k, i) = -a(k, i)/pivo
        end do
        a = matmul(a, m)
        b = matmul(b, m)
        print*, "AM=", a
        print*, "BM=", b
        print*,""
        print*,""
        print*,""
    end do
    call retroSubs(a, b, x, n)
end subroutine

subroutine eliminacaoGaussALT(a, b, x, n)
    integer :: n, i, j, k, pivo, r
    integer :: p(n)
    real :: a(n,n), b(n), x(n), soma, mult, aux, m(n,n), at(n,n)

    x = 0.0

    do i=1, n
        p(i) = i
    end do

    call transposta(a, at, n)

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

        if (at(k,k) == 0.0) then
            STOP 'falta pivoteamento'
            aux = b(k)
            b(k) = b(k+1)
            b(k+1) = aux
            do j=1, n !troca linha r por linha k
                aux = a(j,k)
                a(j,k) = a(j,k+1)
                a(j,k+1) = aux
            end do
        end if


        call identidade(m, n)
        do i=(k+1), n !eliminação ate chegar numa matriz triangular superior
            m(i,k) = -a(i,k)/a(k,k)
            !m(k,i) = -a(k,i)/a(k,k)
        end do
        write(*,*) m
        !write(*,*) a
        !write(*,*) b

        a = MATMUL(m, a)
        b = MATMUL(m, b)


    end do
    call transposta(a, at, n)

    call retroSubs(at, b, x, n)

end subroutine

subroutine multiplicaVetorMatriz(a, b, c, n)
    !multiplicação de vetores de tamanho n
    integer :: n, i, j, k, m
    real :: a(n), b(n), c

    c = 0.0
    do j=1, n
        c = c + a(j)*b(j)
    end do

end subroutine

subroutine eliminacaoGaussJordan(a, b, x, n)
    integer :: n, i, j, k
    real :: a(n,n), b(n), x(n), mult

    x = 0.0

    call eliminacaoGaussALT(a, b, x, n)

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

subroutine decomposicaoCholesky(a, b, x, n)
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

subroutine interativoJacobi(a, b, x, n, tolerancia, k)
    integer :: n, i, j, k
    real :: a(n,n), b(n), xZero(n), tolerancia, r, xNovo(n), xNovoT, x(n), xT, soma, linha, coluna

    k = 1
    r = tolerancia + 1
    xZero = 1.0
    !xZero(1) = 1.0

    do i=1, n
        linha = 0
        coluna = 0
        do j=1, n
            if (i /= j) then
                linha = linha + a(i,j)
                coluna = coluna + a(j,i)
            end if
        end do
        if (a(i,i) < linha .OR. a(i,i) < coluna) then
            open(2, file='Q3.txt', status='replace')
            write(2,*) 'A matriz A não é diagonal dominante.'
            close(2)
            stop 'A matriz A não é diagonal dominante.'
        end if
    end do

    do while (r > tolerancia)
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

        xt = 0.0
        x = xNovo - xZero
        call multiplicaVetor(x, x, xT, n)
        xT = sqrt(xT)

        r = xT / xNovoT
        xZero = xNovo
        k = k + 1
    end do

    x = xZero

end subroutine

subroutine GaussSeidel(a, b, x, n, tolerancia, k)
    integer :: n, i, j, k
    real :: a(n,n), b(n), xZero(n), tolerancia, r, xNovo(n), xNovoT, x(n), xT, soma, linha, coluna

    k = 1
    r = 100.0
    xZero = 1.0
    !xZero(1) = 1.0

    do i=1, n
        linha = 0
        coluna = 0
        do j=1, n
            if (i /= j) then
                linha = linha + a(i,j)
                coluna = coluna + a(j,i)
            end if
        end do
        if (a(i,i) < linha .OR. a(i,i) < coluna) then
            open(2, file='RESUL.txt', status='replace')
            write(2,*) 'A matriz A não é diagonal dominante.'
            close(2)
            stop 'A matriz A não é diagonal dominante.'
        end if
    end do

    do while (r > tolerancia)
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
    !write(*,*) b

    do k = n, 2, -1
        do i=(k-1), 1, -1 !eliminação ate chegar numa matriz diagonal
            mult = a(i,k)/a(k,k)
            a(i,k) = a(i,k) - mult * a(k,k)
            do j=(k-1), 1, -1
                a(i,j) = a(i,j) - mult*a(k,j)
            end do
            do l=n, 1, -1
                b(l,i) = b(l,i) - mult * b(l,k)
                !b(i,l) = b(i,l) - mult * b(k,l)
            end do
        end do
    end do
    !write(*,*) a

    do i=1, n
        if (a(i,i) /= 0.0) then
            do j=1, n
                b(j,i) = b(j,i)/a(i,i)
            end do
        end if
    end do
    !write(*,*) b

    aInversa = b

end subroutine

subroutine inversaALT(a, aInversa, n)
    integer :: n, i, j, k, l
    real :: a(n,n), b(n,n), x(n), mult, aInversa(n,n), m(n,n), at(n,n)

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

    call transposta(a, at, n)
    write(*,*) b
    write(*,*)

    do k = 1, (n-1)
        call identidade(m, n)
        do i=(k+1), n !eliminação ate chegar numa matriz triangular superior
            m(k,i) = -at(k,i)/at(k,k)
        end do
        write(*,*)
        write(*,*) m
        write(*,*) b
        at = MATMUL(m, at)
        b = MATMUL(m, b)
        write(*,*)
        write(*,*) b
    end do
    !write(*,*) m
    write(*,*)
    !write(*,*) a
    write(*,*)
    !write(*,*) b

    do k = n, 2, -1
        call identidade(m, n)
        do i=(k-1), 1, -1 !eliminação ate chegar numa matriz diagonal

            m(k,i) = -at(k,i)/at(k,k)
            !mult = a(i,k)/a(k,k)
            !a(i,k) = a(i,k) - mult * a(k,k)
            !do j=(k-1), 1, -1
            !    a(i,j) = a(i,j) - mult*a(k,j)
            !end do
            !do l=n, 1, -1
            !    b(l,i) = b(l,i) - mult * b(l,k)
            !    !b(i,l) = b(i,l) - mult * b(k,l)
            !end do
        end do
        at = MATMUL(m, at)
        b = matmul(m, b)
    end do
    !write(*,*) a

    do i=1, n
        if (a(i,i) /= 0.0) then
            do j=1, n
                b(i,j) = b(i,j)/a(i,i)
            end do
        end if
    end do
    !write(*,*) b

    aInversa = b

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

subroutine powerMethod(a, x, n, tol, lambda, k)
    integer :: n, i, k
    real :: a(n,n), x(n), y(n), lambda, lastLambda, r, tol

    r = tol + 1
    x = 1.0
    lastLambda = 1
    k = 0
    
    print*, "x(", k, ") = ", x
    do while(r > tol)
        y = MATMUL(a, x)
        k = k + 1
        print*, "y(", k, ") = ", y

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

        !call sleep(1)
    end do
end subroutine

subroutine identidade(a, n)
    integer :: n, i
    real :: a(n,n)

    a = 0
    do i = 1, n
        a(i, i) = 1
    end do
end subroutine

subroutine jacobi(a, x, n, tol, lambda, k)
    integer :: n, i, j, k, iMax, jMax
    real :: a(n,n), a_transposta(n,n), p(n, n), pT(n,n), x(n, n), lambda(n), tol, phi
    logical :: done

    call transposta(a, a_transposta, n)

    if (.not.all(a.EQ.a_transposta)) then
        open(2, file='RESUL.txt', status='replace')
        write(2,*) 'A matriz A não é simétrica.'
        close(2)
        stop 'A matriz A não é simétrica.'
    end if

    k = 0
    PI=4.D0*DATAN(1.D0)

    call identidade(x, n)

    done = .false.
    do while (.true.)
        done = .true.
        aMax = 0
        p = 0
        do i = 1, n
            do j = 1, n
                if (i == j) then
                    p(i, j) = 1
                    cycle
                end if

                if (done .and. abs(a(i, j)) > tol) then
                    done = .false.
                end if

                p(i, j) = 0
                if (abs(a(i, j)) > abs(aMax)) then
                    aMax = a(i, j)
                    iMax = i
                    jMax = j
                end if
            end do
        end do

        if(done) then
            exit
        end if

        if (a(iMax, iMax) /= a(jMax, jMax)) then
            phi = atan((2 * a(iMax, jMax)) / (a(iMax, iMax) - a(jMax, jMax))) / 2
        else
            phi = PI/4
        end if

        p(iMax, jMax) = -sin(phi)
        p(jMax, iMax) = sin(phi)

        p(iMax, iMax) = cos(phi)
        p(jMax, jMax) = cos(phi)


        print*, "iMax = ", iMax
        print*, "jMax = ", jMax
        print*, "aMax = ", aMax
        print*, "phi(", k, ") = ", phi
        print*, "a(", k, ") = ", a
        print*, "x(", k, ") = ", x
        print*, "p(", k, ") = ", p

        call transposta(p, pT, n)

        a = matmul(a, p)
        a = matmul(pT, a)
        x = matmul(x, p)

        print*, "a(", k+1, ") = ", a
        print*, "x(", k+1, ") = ", x
        print*,""
        !call sleep(1)


        k = k + 1
    end do

    do i = 1, n
        lambda(i) = a(i, i)
    end do
end subroutine

subroutine minimosQuadrados(b, x, y, n)
    integer n, i, j
    real b(2), x(n), y(n), a(2,2), soma, c(2), aI(2,2), det

    b = 0.0
    soma = 0.0

    do i=1, n
        soma = soma + 1
    end do
    a(1,1) = soma

    soma = 0.0
    do i=1, n
        soma = soma + x(i)
    end do
    a(2,1) = soma
    a(1,2) = soma

    soma = 0.0
    do i=1, n
        soma = soma + (x(i))**2
    end do
    a(2,2) = soma

    soma = 0.0
    do i=1, n
        soma = soma + y(i)
    end do
    c(1) = soma

    soma = 0.0
    do i=1, n
        soma = soma + x(i)*y(i)
    end do
    c(2) = soma

    det = a(1,1)*a(2,2) - (a(1,2)*a(2,1))
    if (det == 0) then
        STOP
    end if
    aI(1,1) = a(2,2)/det
    aI(1,2) = -a(1,2)/det
    aI(2,1) = -a(2,1)/det
    aI(2,2) = a(1,1)/det

    !call multiplicaVetorMatriz(aI, c, b, 2)
    b = matmul(aI, c)
    write(*,*) a
    write(*,*) c
    write(*,*) b


    !do i = 1, n
    !    lambda(i) = a(i, i)
    !end do
end subroutine