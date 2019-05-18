module functions
implicit none
interface
function IFunction(x) result (y)
    real, intent(in) :: x
    real             :: y
end function IFunction
end interface

contains


subroutine fL(x, y, t)
    real x, y, t
    y = t + x

end subroutine

subroutine fLL(x, y, t, xL)
    real x, y, t, xL
    y = t + x

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

subroutine metodoDeNewtonOriginal(func, tol, x0, nIter, a)
    procedure(IFunction) :: func
    real :: tol, x0, xK, xK_1, f, df, tolK, a
    integer :: nIter, k

    xK_1 = x0
    do k = 1, nIter
        write (*,*) "k=", k

        f = func(xK_1)
        write (*,*) "f=", f

        call diferencaCentral(func, xK_1, 0.25, df)
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