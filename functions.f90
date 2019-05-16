subroutine f(x, y)
    real x, y
    !y = x*x - 4 * cos(x)
    !y = x*x
    y = sin(x)

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