module functions
implicit none
interface
function IFunction(x) result (y)
    REAL(8), intent(in) :: x
    REAL(8)             :: y
end function IFunction

function IFunctionD(x, t) result (y)
    REAL(8), intent(in) :: x, t
    REAL(8)             :: y
end function IFunctionD

function IFunctionDD(x, t, xL) result (y)
    REAL(8), intent(in) :: x, t, xL
    REAL(8)             :: y
end function IFunctionDD

function IFunctionMD(x) result (y)
    REAL(8), dimension(:), intent(in) :: x
    REAL(8), allocatable, dimension(:) :: y
end function IFunctionMD

function IJacobian(x) result (y)
    REAL(8), dimension(:), intent(in) :: x
    REAL(8), allocatable, dimension(:,:) :: y
end function IJacobian
end interface

contains

subroutine bissecao(f, tol, a, b)
    procedure(IFunction) ::  f
    REAL(8) :: meio, a, b, tol, delta, y
    integer :: i

    i = 0
    do while (abs(b - a) > tol)
        write (*,*) ""
        write (*,*) "k=", i
        meio = (a + b) / 2.0
        write (*,*) "meio=", meio
        y = f(meio)
        write (*,*) "y=", y
        if (y >= 0) then
            b = meio
        else
            a = meio
        end if
        i = i + 1
    end do
    write (*,*) i
    a = meio
    b = meio
    write (*,*) a

end subroutine

subroutine metodoDeNewtonOriginal(func, dfunc, tol, x0, nIter, a)
    procedure(IFunction) :: func, dfunc
    REAL(8) :: tol, x0, xK, xK_1, f, df, tolK, a
    integer :: nIter, k

    xK_1 = x0
    do k = 1, nIter
        write (*,*) "k=", k

        f = func(xK_1)
        write (*,*) "f=", f

        df = dfunc(xK_1) !call diferencaCentral(func, xK_1, 0.25, df)
        write (*,*) "df=", df

        xK = xK_1 - f/df
        write (*,*) "xK=", xK

        tolK = abs(xK - xK_1)
        if (tolK.lt.tol) then
            a = xK
            return
        endif
        xK_1 = xK
    end do
    write (*,*) 'Convergencia não atingida.'
end subroutine

subroutine metodoDeNewtonSecante(func, tol, x0, nIter, deltaX, a)
    procedure(IFunction) :: func
    REAL(8) :: tol, x0, xK, xK_1, xK_p1, tolK, a, fA, fI, deltaX
    integer :: nIter, k

    fA = func(x0)
    xK_1 = x0
    xK = xK_1 + deltaX

    do k = 1, nIter
        write (*,*) ""
        write (*,*) "k=", k

        fI = func(xK)
        !write (*,*) "f=", f


        xK_p1 = xK - ((fI * (xK - xK_1))/(fI-fA))
        write (*,*) "xK-1=", xK_1
        write (*,*) "xK=", xK
        write (*,*) "xK+1=", xK_p1

        tolK = abs(xK_p1 - xK)
        write (*,*) "tolK=", tolK

        if (tolK.lt.tol) then
            a = xK
            return
        else
            fA = fI
        endif
        xK_1 = xK
        xK = xK_p1
    end do
    write (*,*) 'Convergencia não atingida.'
end subroutine

subroutine interpolacaoInversa(func, tol, nIter, a, x)
    procedure(IFunction) :: func
    REAL(8), intent(in) :: tol
    integer, intent(in) :: nIter

    REAL(8), dimension(3) :: x, y, auxX
    REAL(8) :: xK, xK_1, tolK, aux
    integer :: k, i, j

    REAL(8), intent(out) :: a


    xK_1 = 10E+36
    auxX = x
    do k = 1, nIter
        write (*,*) ""
        write (*,*) "k=", k

        write (*,*) "x1=", auxX(1)
        write (*,*) "x2=", auxX(2)
        write (*,*) "x3=", auxX(3)
        y(1) = func(auxX(1))
        y(2) = func(auxX(2))
        y(3) = func(auxX(3))
        write (*,*) "y1=", y(1)
        write (*,*) "y2=", y(2)
        write (*,*) "y3=", y(3)

        xK = (y(2) * y(3) * auxX(1)) / ((y(1) - y(2)) * (y(1) - y(3))) + &
             (y(1) * y(3) * auxX(2)) / ((y(2) - y(1)) * (y(2) - y(3))) + &
             (y(1) * y(2) * auxX(3)) / ((y(3) - y(1)) * (y(3) - y(2)))
        write (*,*) "x*=", xK

        tolK = abs(xK - xK_1)
        write (*,*) "tol=", tolK

        if (tolK.lt.tol) then
            a = xK
            return
        else
            i = MAXLOC(abs(y), DIM = 1)
            auxX(i) = xK
            y(i) = func(xK)

            do i = 1, 3
                aux = auxX(i)
                j = MINLOC(auxX(i:), DIM = 1) + i - 1

                auxX(i) = auxX(j)
                auxX(j) = aux
            end do
        endif

        xK_1 = xK
    end do
    write (*,*) 'Convergencia não atingida.'
end subroutine

REAL(4) function euclidianModule(A)
    REAL(8), dimension(:) :: A
    integer :: i, j

    do i = 1, size(A)
        euclidianModule = euclidianModule + A(i) ** 2
    end do
    euclidianModule = sqrt(euclidianModule)
end function

function invert(a) result(a_1)
    REAL(8), dimension(:,:), intent(in) :: a
    REAL(8), dimension(:,:), allocatable :: a_1
    integer :: n

    allocate(a_1, mold = a)
    n = size(a(1,:))
    call matrixinv(a, a_1, n)
    !a_1 = TRANSPOSE(a_1)
end function

subroutine matrixinv(a,b,n)
    ! subroutine to calculate the inverse of a matrix using Gauss-Jordan elimination
    ! the inverse of matrix a(n,n) is calculated and stored in the matrix b(n,n)
    integer :: i,j,k,l,m,n,irow

    real(8):: big,a(n,n),b(n,n),dum



    !build the identity matrix
    do i = 1,n
        do j = 1,n
            b(i,j) = 0.0
        end do

        b(i,i) = 1.0
        end do
    do i = 1,n ! this is the big loop over all the columns of a(n,n)

    ! in case the entry a(i,i) is zero, we need to find a good pivot; this pivot

        ! is chosen as the largest value on the column i from a(j,i) with j = 1,n

        big = a(i,i)
    do j = i,n

        if (a(j,i).gt.big) then

        big = a(j,i)

    irow = j
        end if
        end do

    ! interchange lines i with irow for both a() and b() matrices

    if (big.gt.a(i,i)) then

        do k = 1,n

        dum = a(i,k) ! matrix a()

    a(i,k) = a(irow,k)

    a(irow,k) = dum

    dum = b(i,k) ! matrix b()

    b(i,k) = b(irow,k)

    b(irow,k) = dum
        end do
        end if

    ! divide all entries in line i from a(i,j) by the value a(i,i);

    ! same operation for the identity matrix

    dum = a(i,i)

    do j = 1,n

        a(i,j) = a(i,j)/dum

    b(i,j) = b(i,j)/dum
        end do

    ! make zero all entries in the column a(j,i); same operation for indent()

    do j = i+1,n

        dum = a(j,i)

    do k = 1,n

        a(j,k) = a(j,k) - dum*a(i,k)

    b(j,k) = b(j,k) - dum*b(i,k)
        end do

    end do

    end do



    ! substract appropiate multiple of row j from row j-1

    do i = 1,n-1

        do j = i+1,n

        dum = a(i,j)

    do l = 1,n

        a(i,l) = a(i,l)-dum*a(j,l)

    b(i,l) = b(i,l)-dum*b(j,l)

    end do

    end do

end do
end subroutine




subroutine metodoDeNewtonMD(func, Jac, tol, x, nIter)
    procedure (IFunctionMD) :: func
    procedure (IJacobian) :: Jac
    REAL(4) :: tol, tolK, a, b
    REAL(8), dimension(:) :: x
    REAL(8), allocatable, dimension(:) :: F, deltaX
    REAL(8), allocatable, dimension(:,:) :: J, J_1
    integer :: nIter, k

    tolK = 0
    allocate(deltaX, mold = x)
    allocate(F, source = func(x))
    allocate(J, source = Jac(x))
    allocate(J_1, mold = J)

    do k = 1, nIter
        write (*,*) ""
        write (*,*) 'k=', k
        write (*,*) 'F=', F
        write (*,*) 'J=', J
        J_1 = invert(J)
        write (*,*) 'J_1=', J_1
        deltaX = - matmul(J_1, F)
        write (*,*) 'x=', x
        write (*,*) 'deltaX=', deltaX
        x = x + deltaX

        b = euclidianModule(deltaX)
        a = euclidianModule(x)
        tolK = b / a

        write (*,*) 'tolK=', tolK, ' - mod(deltaX) = ', b, ' - mod(x) = ', a
        if (tolK.lt.tol) then
            return
        end if

        F = func(x)
        J = Jac(x)
    end do
    write (*,*) 'Convergencia não atingida.'
end subroutine

function transposeVector(v) result(vT)
    REAL(8), dimension(:), intent(in) :: v
    REAL(8), allocatable, dimension(:,:) :: vT
    integer :: i
    allocate(vT(1, size(v)))

    do i = 1, size(v)
        vT(1, i) = v(i)
    end do
end function

subroutine metodoDeBroydenMD(func, B, tol, x, nIter)
    procedure (IFunctionMD) :: func
    REAL(4) :: tol, tolK
    REAL(8), dimension(:) :: x
    REAL(8), dimension(:,:) :: B
    REAL(8), allocatable, dimension(:) :: F, deltaX, Y
    REAL(8), allocatable, dimension(:,:) :: J, B_1, J_1, deltaX_T
    integer :: nIter, k

    allocate(deltaX, mold = x)
    allocate(F, source = func(x))
    allocate(J, source = B)
    allocate(B_1, source = B)
    allocate(Y, mold = x)
    allocate(J_1, mold = J)
    allocate(deltaX_T(1, size(deltaX)))

    write (*,*) 'b=', B
    do k = 1, nIter
        write (*,*) ""
        write (*,*) 'k=', k


        write (*,*) 'F=', F
        write (*,*) 'J=', J
        J_1 = invert(J)
        write (*,*) 'J_1=', J_1

        deltaX = - matmul(J_1, F)
        write (*,*) 'X=', x
        write (*,*) 'deltaX=', deltaX

        x = x + deltaX
        write (*,*) 'newX=', x
        Y = func(x) - F
        write (*,*) 'Y=', Y!, '- func = ', func(x), ' - F = ',F

        tolK = euclidianModule(deltaX) / euclidianModule(x)

        write (*,*) 'tolK=', tolK!, ",", euclidianModule(deltaX), ",", euclidianModule(x)
        if (tolK.lt.tol) then
            return
        else
            deltaX_T = transposeVector(deltaX)
            !write (*,*) 'deltaX_T=',  MATMUL(reshape(Y - MATMUL(B_1, deltaX), shape=(/ 2, 1 /)), deltaX_T)
            B = B_1 + MATMUL(reshape(Y - MATMUL(B_1, deltaX), shape=(/ 2, 1 /)), deltaX_T) / dot_product(deltaX, deltaX)

            write (*,*) 'B=', B
        end if


        F = func(x)
        J = B_1
        B_1 = B
        call inversa(J, J_1, size(J_1(1,:)))

        deltaX = - matmul(J_1, F)
    end do
    write (*,*) 'Convergencia não atingida.'
end subroutine

subroutine minimosQuadrados(func, Jac, tol, x, nIter)
    procedure (IFunctionMD) :: func
    procedure (IJacobian) :: Jac
    REAL(8) :: tol, tolK
    REAL(8), dimension(:) :: x
    REAL(8), allocatable, dimension(:) :: F, deltaX
    REAL(8), allocatable, dimension(:,:) :: J, J_1, J_t, J_t1
    integer :: nIter, k

    allocate(deltaX, mold = x)
    allocate(F, source = func(x))
    write (*,*) 'F=', F

    allocate(J, source = Jac(x))
    write (*,*) 'J=', J

    allocate(J_t, mold = J)
    allocate(J_t1, mold = J)
    call transposta(J, J_t, size(J_t(1,:)))
    write (*,*) 'J_t=', J_t

    J_t1 = matmul(J_t, J)

    allocate(J_1, mold = J)
    J_1 = invert(J)
    write (*,*) 'J_1=', J_1

    J = matmul(J_1, J_t)
    deltaX = - matmul(J, F)

    do k = 2, nIter
        write (*,*) ""
        write (*,*) 'k=', k
        x = x + deltaX
        tolK = euclidianModule(deltaX) / euclidianModule(x)

        write (*,*) 'x=', x
        write (*,*) 'deltaX=', deltaX
        if (tolK.lt.tol) then
            return
        end if


        F = func(x)
        J = Jac(x)
        call transposta(J, J_t, size(J_t(1,:)))
        J_t1 = matmul(J_t, J)
        call inversa(J, J_1, size(J_1(1,:)))
        J = matmul(J_1, J_t)
        deltaX = - matmul(J, F)
    end do
    write (*,*) 'Convergencia não atingida.'
end subroutine

function readTable(fileName, n) result(table)
    Character(len=*) :: fileName
    REAL(8), dimension(:,:), allocatable :: table
    integer :: i, j
    integer, intent(out) :: n

    open (1, file=fileName, status='old', action='read')
    read(1,*) n
    allocate(table(n, 2))
    do i=1, n
        read(1,*) j, table(i, 1), table(i, 2)
    end do
    close (1)
end function

function integracaoPesos(func, n, a, b, table) result(area)
    procedure(IFunction) :: func
    integer :: i, n
    REAL(8), dimension(:, :) :: table
    REAL(8) :: area, a, b, L

    area = 0
    L = (b - a)
    do i = 1, n
        write(*,*) ""
        write(*,*) "i = ", i
        write(*,*) "W = ", table(i, 1)
        write(*,*) "X = ", table(i, 2)
        area = area + (table(i, 1) * func(table(i, 2)))
    end do
end function

subroutine integracaoPolinomial(func, points, a, b, area)
    procedure(IFunction) :: func
    integer :: n, i
    REAL(8) :: a, b, area, L
    REAL(8), dimension(:,:), allocatable :: table
    character(len=*) :: points

    L = (b-a)
    table = readTable("tabelas/polinomial/" // points //".txt", n)
    do i = 1, n
        table(i, 1) = table(i, 1) * L
        table(i, 2) = (a + table(i, 2) * L)
    end do
    area = integracaoPesos(func, n, a, b, table)
end subroutine

subroutine integracaoQuadratura(func, points, a, b, area)
    procedure(IFunction) :: func
    integer :: i, n
    REAL(8) :: a, b, area, L
    REAL(8), dimension(:,:), allocatable :: table
    character(len=*) :: points

    L = (b-a)
    table = readTable("tabelas/quadratura/" // points //".txt", n)
    do i = 1, n
        table(i, 2) = 1.0/2.0 * (a + b + table(i, 2)*L)
    end do
    area = integracaoPesos(func, n, a, b, table)
    area = L/2 * area
end subroutine

subroutine diferencaCentral(f, x, deltaX, fLinha)
    procedure(IFunction) :: f
    REAL(8) :: x, deltaX, f1, f2, fLinha

    f1 = f(x+deltaX)
    f2 = f(x-deltaX)

    fLinha = (f1-f2)/(2*deltaX)

end subroutine

subroutine passoFrente(f, x, deltaX, fLinha)
    procedure(IFunction) :: f
    REAL(8) :: x, deltaX, f1, f2, fLinha

    f1 = f(x+deltaX)
    f2 = f(x)

    fLinha = (f1-f2)/(deltaX)

end subroutine

subroutine passoTras(f, x, deltaX, fLinha)
    procedure(IFunction) :: f
    REAL(8) :: x, deltaX, f1, f2, fLinha

    f1 = f(x)
    f2 = f(x-deltaX)

    fLinha = (f1-f2)/(deltaX)

end subroutine

subroutine interpolacaoRichard(f, x, deltaX1, deltaX2, fLinha, p)
    procedure(IFunction) :: f
    REAL(8) :: x, deltaX1, deltaX2, d1, d2, q, fLinha, p

    q = deltaX1/deltaX2

    call passoFrente(f, x, deltaX1, d1)
    call passoFrente(f, x, deltaX2, d2)

    fLinha = d1 + (d1 - d2)/((q**(-p))-1)

end subroutine

subroutine euler(f, xK, xZero, tZero, k, deltaT)
    procedure(IFunctionD) :: f
    REAL(8) :: xK, xZero, tZero, deltaT, y
    integer :: i, k

    !do i=1, k
        !call fL(xZero, y, tZero)
    !    y = f(xZero, tZero)
    !    xK = xZero + y * deltaT
    !    tZero = i*deltaT
    !    xZero = xK
        !write(*,*) tZero
    !end do

    i = 1
    do while (tZero<k)
        !call fL(xZero, y, tZero)
        y = f(xZero, tZero)
        xK = xZero + y * deltaT
        tZero = i*deltaT
        xZero = xK
        i = i + 1
    end do
end subroutine

subroutine rungeKutta2(f, xK, xZero, tZero, k, deltaT)
    procedure(IFunctionD) :: f
    REAL(8) :: xK, xZero, tZero, deltaT, y1, y2
    integer :: i, k

    !do i=1, k
        !call fL(xZero, y1, tZero)
        !call fL(xZero+deltaT*y1, y2, tZero+deltaT)
    !    y1 = f(xZero, tZero)
    !    y2 = f(xZero+deltaT*y1, tZero+deltaT)
    !    xK = xZero + (y1+y2) * deltaT/2.0
    !    tZero = i*deltaT
    !    xZero = xK
    !    write(*,*) tZero
    !end do
    i = 1
    do while (tZero<k)
        !call fL(xZero, y1, tZero)
        !call fL(xZero+deltaT*y1, y2, tZero+deltaT)
        y1 = f(xZero, tZero)
        y2 = f(xZero+deltaT*y1, tZero+deltaT)
        xK = xZero + (y1+y2) * deltaT/2.0
        tZero = i*deltaT
        xZero = xK
        i = i + 1
    end do

end subroutine

subroutine rungeKutta4(f, xK, xZero, tZero, k, deltaT)
    procedure(IFunctionD) :: f
    REAL(8) :: xK, xZero, tZero, deltaT, y1, y2, y3, y4
    integer :: i, k

    !do i=1, k
        !call fL(xZero, y1, tZero)
        !call fL(xZero+deltaT/2.0*y1, y2, tZero+deltaT/2.0)
        !call fL(xZero+deltaT/2.0*y2, y3, tZero+deltaT/2.0)
        !call fL(xZero+deltaT*y3, y4, tZero+deltaT)
    !    y1 = f(xZero, tZero)
    !    y2 = f(xZero+deltaT/2.0*y1, tZero+deltaT/2.0)
    !    y3 = f(xZero+deltaT/2.0*y2, tZero+deltaT/2.0)
    !    y4 = f(xZero+deltaT*y3, tZero+deltaT)
    !    xK = xZero + (y1+2*y2+2*y3+y4) * deltaT/6.0
    !    tZero = i*deltaT
    !    xZero = xK
        !write(*,*) tZero
    !end do


    i = 1
    do while (tZero<k)
        y1 = f(xZero, tZero)
        y2 = f(xZero+deltaT/2.0*y1, tZero+deltaT/2.0)
        y3 = f(xZero+deltaT/2.0*y2, tZero+deltaT/2.0)
        y4 = f(xZero+deltaT*y3, tZero+deltaT)
        xK = xZero + (y1+2*y2+2*y3+y4) * deltaT/6.0
        tZero = i*deltaT
        xZero = xK
        i = i + 1
    end do
end subroutine

subroutine taylor2(f, xK, xZero, xLZero, tZero, k, deltaT)
    procedure(IFunctionDD) :: f
    REAL(8) :: xK, xZero, tZero, deltaT, xLZero, x, y
    integer :: i, k

    !do i=1, k
        !call fLL(xZero, y, tZero, xLZero)
    !    y = f(xZero, tZero, xLZero)
    !    xK = xZero + xLZero*deltaT + y * (deltaT**2)/2.0
    !    xLZero = xLZero + y * deltaT
    !    tZero = i*deltaT
    !    xZero = xK
    !end do

    i = 1
    do while (tZero<k)
        y = f(xZero, tZero, xLZero)
        xK = xZero + xLZero*deltaT + y * (deltaT**2)/2.0
        xLZero = xLZero + y * deltaT
        tZero = i*deltaT
        xZero = xK
        i = i + 1
    end do
end subroutine


subroutine rungeKuttaNystrom(f, xK, xZero, xLZero, tZero, k, deltaT)
    procedure(IFunctionDD) :: f
    REAL(8) :: xK, xZero, tZero, deltaT, xLZero, y1, y2, y3, y4
    integer :: i, k

    !do i=1, k
        !call fLL(xZero, y1, tZero, xLZero)
    !    y1 = f(xZero, tZero, xLZero)
    !    y1 = deltaT/2.0*y1
        !call fLL(xZero+xLZero+y1/2, y2, tZero+deltaT/2.0, xLZero+y1)
    !    y2 = f(xZero+xLZero+y1/2, tZero+deltaT/2.0, xLZero+y1)
    !    y2 = deltaT/2.0*y2
        !call fLL(xZero+xLZero+y1/2, y3, tZero+deltaT/2.0, xLZero+y2)
    !    y3 = f(xZero+xLZero+y1/2, tZero+deltaT/2.0, xLZero+y2)
    !    y3 = deltaT/2.0*y3
        !call fLL(xZero+deltaT*(xLZero+y3), y4, tZero+deltaT, xLZero+2*y3)
    !    y4 = f(xZero+deltaT*(xLZero+y3), tZero+deltaT, xLZero+2*y3)
    !    y4 = deltaT/2.0*y4

    !    xK = xZero + (xLZero + (y1+y2+y3)/3) * deltaT
    !    xLZero = xLZero + (y1 + 2*y2 + 2*y3 + y4)/3
    !    tZero = i*deltaT
    !    xZero = xK
    !end do

    i = 1
    do while (tZero<k)
        y1 = f(xZero, tZero, xLZero)
        y1 = deltaT/2.0*y1
        !call fLL(xZero+xLZero+y1/2, y2, tZero+deltaT/2.0, xLZero+y1)
        y2 = f(xZero+xLZero+y1/2, tZero+deltaT/2.0, xLZero+y1)
        y2 = deltaT/2.0*y2
        !call fLL(xZero+xLZero+y1/2, y3, tZero+deltaT/2.0, xLZero+y2)
        y3 = f(xZero+xLZero+y1/2, tZero+deltaT/2.0, xLZero+y2)
        y3 = deltaT/2.0*y3
        !call fLL(xZero+deltaT*(xLZero+y3), y4, tZero+deltaT, xLZero+2*y3)
        y4 = f(xZero+deltaT*(xLZero+y3), tZero+deltaT, xLZero+2*y3)
        y4 = deltaT/2.0*y4

        xK = xZero + (xLZero + (y1+y2+y3)/3) * deltaT
        xLZero = xLZero + (y1 + 2*y2 + 2*y3 + y4)/3
        tZero = i*deltaT
        xZero = xK
        i = i + 1
    end do
end subroutine

end module functions