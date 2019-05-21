module functions
implicit none
interface
function IFunction(x) result (y)
    real, intent(in) :: x
    real             :: y
end function IFunction


function IFunctionMD(x) result (y)
    real, dimension(:), intent(in) :: x
    real, allocatable, dimension(:) :: y
end function IFunctionMD

function IJacobian(x) result (y)
    real, dimension(:), intent(in) :: x
    real, allocatable, dimension(:,:) :: y
end function IJacobian
end interface

contains

subroutine fL(x, y, t)
    real x, y, t
    y = -2*t*(x**2)

end subroutine

subroutine fLL(x, y, t, xL)
    real :: x, y, t, xL
    y = - 9.81 - (1 * xL * abs(xL))

end subroutine

subroutine bissecao(f, tol, a, b)
    real :: f, meio, a, b, tol, delta, y
    integer :: i

    i = 0
    do while (abs(b-a) > tol)
        meio = (a+b) / 2.0
        y = f(meio)
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

end subroutine

subroutine metodoDeNewtonOriginal(func, dfunc, tol, x0, nIter, a)
    procedure(IFunction) :: func, dfunc
    real :: tol, x0, xK, xK_1, f, df, tolK, a
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
    real :: tol, x0, xK, xK_1, xK_p1, tolK, a, fA, fI, deltaX
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
    real, intent(in) :: tol
    integer, intent(in) :: nIter

    real, dimension(3) :: x, y, auxX
    real :: xK, xK_1, tolK, aux
    integer :: k, i, j

    real, intent(out) :: a


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

real function euclidianModule(A)
    real, dimension(:) :: A
    integer :: i, j

    do i = 1, size(A)
        euclidianModule = euclidianModule + A(i) ** 2
    end do
    euclidianModule = sqrt(euclidianModule)
end function

subroutine metodoDeNewtonMD(func, Jac, tol, x, nIter)
    procedure (IFunctionMD) :: func
    procedure (IJacobian) :: Jac
    real :: tol, tolK
    real, dimension(:) :: x
    real, allocatable, dimension(:) :: F, deltaX
    real, allocatable, dimension(:,:) :: J, J_1
    integer :: nIter, k

    allocate(deltaX, mold = x)
    allocate(F, source = func(x))
    write (*,*) 'F=', F

    allocate(J, source = Jac(x))
    write (*,*) 'J=', J

    allocate(J_1, mold = J)
    call inversa(J, J_1, size(J_1(1,:)))
    write (*,*) 'J_1=', J_1

    deltaX = - matmul(J_1, F)

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
        call inversa(J, J_1, size(J_1(1,:)))

        deltaX = - matmul(J_1, F)
    end do
    write (*,*) 'Convergencia não atingida.'
end subroutine

function transposeVector(v) result(vT)
    real, dimension(:), intent(in) :: v
    real, allocatable, dimension(:,:) :: vT
    integer :: i
    allocate(vT(1, size(v)))

    do i = 1, size(v)
        vT(1, i) = v(i)
    end do
end function

subroutine metodoDeBroydenMD(func, B, tol, x, nIter)
    procedure (IFunctionMD) :: func
    real :: tol, tolK
    real, dimension(:) :: x
    real, dimension(:,:) :: B
    real, allocatable, dimension(:) :: F, deltaX, Y
    real, allocatable, dimension(:,:) :: J, B_1, J_1, deltaX_T
    integer :: nIter, k

    allocate(deltaX, mold = x)
    allocate(F, source = func(x))
    allocate(J, source = B)
    allocate(B_1, source = B)
    allocate(Y, mold = x)
    allocate(J_1, mold = J)
    allocate(deltaX_T(1, size(deltaX)))

    x = x + deltaX
    write (*,*) 'b=', B
    do k = 2, nIter
        write (*,*) ""
        write (*,*) 'k=', k


        write (*,*) 'F=', F
        write (*,*) 'J=', J
        call inversa(J, J_1, size(J_1(1,:)))
        write (*,*) 'J_1=', J_1

        deltaX = - matmul(J_1, F)
        write (*,*) 'deltaX=', deltaX

        x = x + deltaX
        write (*,*) 'x=', x
        Y = func(x) - F
        write (*,*) 'Y=', Y

        tolK = euclidianModule(deltaX) / euclidianModule(x)

        write (*,*) 'tolK=', tolK, ",", euclidianModule(deltaX), ",", euclidianModule(x)
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

function readTable(fileName, n) result(table)
    Character(len=*) :: fileName
    real, dimension(:,:), allocatable :: table
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
    real, dimension(:, :) :: table
    real :: area, a, b, L

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
    real :: a, b, area, L
    real, dimension(:,:), allocatable :: table
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
    real :: a, b, area, L
    real, dimension(:,:), allocatable :: table
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
    real :: x, deltaX, f1, f2, fLinha

    f1 = f(x+deltaX)
    f2 = f(x-deltaX)

    fLinha = (f1-f2)/(2*deltaX)

end subroutine

subroutine passoFrente(f, x, deltaX, fLinha)
    procedure(IFunction) :: f
    real :: x, deltaX, f1, f2, fLinha

    f1 = f(x+deltaX)
    f2 = f(x)

    fLinha = (f1-f2)/(deltaX)

end subroutine

subroutine passoTras(f, x, deltaX, fLinha)
    procedure(IFunction) :: f
    real :: x, deltaX, f1, f2, fLinha

    f1 = f(x)
    f2 = f(x-deltaX)

    fLinha = (f1-f2)/(deltaX)

end subroutine

subroutine interpolacaoRichard(f, x, deltaX1, deltaX2, fLinha, p)
    procedure(IFunction) :: f
    real :: x, deltaX1, deltaX2, d1, d2, q, fLinha, p

    q = deltaX1/deltaX2

    call passoFrente(f, x, deltaX1, d1)
    call passoFrente(f, x, deltaX2, d2)

    fLinha = d1 + (d1 - d2)/((q**(-p))-1)

end subroutine

subroutine euler(xK, xZero, tZero, k, deltaT)
    real :: xK, xZero, tZero, deltaT, y
    integer :: i, k

    do i=1, k
        call fL(xZero, y, tZero)
        xK = xZero + y * deltaT
        tZero = i*deltaT
        xZero = xK
    end do

end subroutine

subroutine rungeKutta2(xK, xZero, tZero, k, deltaT)
    real :: xK, xZero, tZero, deltaT, y1, y2
    integer :: i, k

    do i=1, k
        call fL(xZero, y1, tZero)
        call fL(xZero+deltaT*y1, y2, tZero+deltaT)
        xK = xZero + (y1+y2) * deltaT/2.0
        tZero = i*deltaT
        xZero = xK
    end do

end subroutine

subroutine rungeKutta4(xK, xZero, tZero, k, deltaT)
    real :: xK, xZero, tZero, deltaT, y1, y2, y3, y4
    integer :: i, k

    do i=1, k
        call fL(xZero, y1, tZero)
        call fL(xZero+deltaT/2.0*y1, y2, tZero+deltaT/2.0)
        call fL(xZero+deltaT/2.0*y2, y3, tZero+deltaT/2.0)
        call fL(xZero+deltaT*y3, y4, tZero+deltaT)
        xK = xZero + (y1+2*y2+2*y3+y4) * deltaT/6.0
        tZero = i*deltaT
        xZero = xK
    end do

end subroutine

subroutine taylor2(xK, xZero, xLZero, tZero, k, deltaT)
    real :: xK, xZero, tZero, deltaT, xLZero, x, y
    integer :: i, k

    do i=1, k
        call fLL(x, y, tZero, xLZero)
        xK = xZero + xLZero*deltaT + y * (deltaT**2)/2.0
        xLZero = xLZero + y * deltaT
        tZero = i*deltaT
        xZero = xK
    end do

end subroutine


subroutine rungeKuttaNystrom(xK, xZero, xLZero, tZero, k, deltaT)
    real :: xK, xZero, tZero, deltaT, xLZero, y1, y2, y3, y4
    integer :: i, k

    do i=1, k
        call fLL(xZero, y1, tZero, xLZero)
        y1 = deltaT/2.0*y1
        call fLL(xZero+xLZero+y1/2, y2, tZero+deltaT/2.0, xLZero+y1)
        y2 = deltaT/2.0*y2
        call fLL(xZero+xLZero+y1/2, y3, tZero+deltaT/2.0, xLZero+y2)
        y3 = deltaT/2.0*y3
        call fLL(xZero+deltaT*(xLZero+y3), y4, tZero+deltaT, xLZero+2*y3)
        y4 = deltaT/2.0*y4

        xK = xZero + (xLZero + (y1+y2+y3)/3) * deltaT
        xLZero = xLZero + (y1 + 2*y2 + 2*y3 + y4)/3
        tZero = i*deltaT
        xZero = xK
    end do

end subroutine

end module functions