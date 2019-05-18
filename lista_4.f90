include 'methods.f90'
include 'functions.f90'

program lista_4
    use functions
    IMPLICIT NONE
    integer n, i, j
    real :: a, b, tol, y

    tol = 0.0001
    a = 0
    b = 10
    !call bissecao(f1, tol, a, b)
    call metodoDeNewtonOriginal(f1, tol, 10.0, 100, a)
    !call metodoDeNewtonSecante(f1, tol, 10.0, 100, 0.001, a)
    !call f(b, y)

    write (*,*) a
    call sleep(10000)

contains
    
    real function f1(x) result(y)
        real, intent(in) :: x
        !y = x*x - 4 * cos(x)
        !y = x*x
        !y = sin(x)
        y = x**2 - 4 * cos(x)
    end function

end program lista_4