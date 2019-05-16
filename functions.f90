subroutine f(x, y)
    real x, y
    y = x*x - 4 * cos(x)

end subroutine

subroutine bissecao(tol, a, b)
    real :: meio, a, b, tol, delta
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

