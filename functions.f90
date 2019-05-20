subroutine f(x, y)
    real x, y
    !y = x*x - 4 * cos(x)
    !y = x*x
    y = sin(x)
    !y = x**3 + 1/exp(x)
    !y = x**(1/3) + log(x)
    !y = 1 - exp(-(x/5)**2)

end subroutine

subroutine fL(x, y, t)
    real x, y, t
    y = -2*t*(x**2)

end subroutine

subroutine fLL(x, y, t, xL)
    real :: x, y, t, xL
    y = - 9.81 - (1 * xL * abs(xL))

end subroutine

subroutine bissecao(tol, a, b)
    real :: meio, a, b, tol, delta, y
    integer :: i

    i = 0
    do while (abs(b-a) > tol)
        meio = (a+b) / 2.0
        call f(meio, y)
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

subroutine diferencaCentral(x, deltaX, fLinha)
    real :: x, deltaX, f1, f2, fLinha

    call f(x+deltaX, f1)
    call f(x-deltaX, f2)

    fLinha = (f1-f2)/(2*deltaX)

end subroutine

subroutine passoFrente(x, deltaX, fLinha)
    real :: x, deltaX, f1, f2, fLinha

    call f(x+deltaX, f1)
    call f(x, f2)

    fLinha = (f1-f2)/(deltaX)

end subroutine

subroutine passoTras(x, deltaX, fLinha)
    real :: x, deltaX, f1, f2, fLinha

    call f(x, f1)
    call f(x-deltaX, f2)

    fLinha = (f1-f2)/(deltaX)

end subroutine

subroutine interpolacaoRichard(x, deltaX1, deltaX2, fLinha, p)
    real :: x, deltaX1, deltaX2, d1, d2, q, fLinha

    q = deltaX1/deltaX2

    call passoFrente(x, deltaX1, d1)
    call passoFrente(x, deltaX2, d2)

    fLinha = d1 + (d1 - d2)/((q**(-p))-1)

end subroutine

subroutine euler(xK, xZero, tZero, k, deltaT)
    real :: xK, xZero, tZero, deltaT
    integer :: i, k

    do i=1, k
        call fL(xZero, y, tZero)
        xK = xZero + y * deltaT
        tZero = i*deltaT
        xZero = xK
    end do

end subroutine

subroutine rungeKutta2(xK, xZero, tZero, k, deltaT)
    real :: xK, xZero, tZero, deltaT
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
    real :: xK, xZero, tZero, deltaT
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
    real :: xK, xZero, tZero, deltaT
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
    real :: xK, xZero, tZero, deltaT
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