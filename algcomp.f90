program algcomp
    integer n, i, j
    real, allocatable, dimension(:) :: b, x
    real, allocatable, dimension(:,:) :: a, d, c

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

    call eliminacaoGaussJordan(a, b, x, n)


    open(2, file='RESUL.txt', status='replace')
    do i=1, n
        write(2,*) (a(i, j), j=1, n)
    end do
    write(2, *) (x(i), i=1, n)
    close(2)

end program algcomp

subroutine retroSubs(a, b, x, n)
    integer :: n, i, j
    real :: a(n,n), b(n), x(n), sum

    x = 0.0

    do i = n, 1, -1
        sum = 0.0
        do j = i+1, n
            sum = sum + a(i,j) * x(j)
        end do
        x(i) = (b(i) - sum)/a(i,i)
    end do
end subroutine

subroutine eliminacaoGauss(a, b, x, n)
    integer :: n, i, j, k, pivo, r, aux
    integer :: p(n)
    real :: a(n,n), b(n), x(n), sum, mult

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
    integer :: p(n)
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


subroutine teste(a, b, c, n)
    !multiplicação de matrizes quadradas, precisa generalizar?
    integer :: n, i, j, k
    real :: a(n,n), b(n,n), c(n,n)

    call multiplicaMatrizes(a, b, c, n)
end subroutine