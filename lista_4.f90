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
    !call metodoDeNewtonOriginal(f1, df1, tol, 10.0, 100, a)
    !call metodoDeNewtonSecante(f1, tol, 10.0, 100, 0.001, a)
    call interpolacaoInversa(f1, 5E-4, 100, a, [3.0, 5.0, 10.0])
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

    real function df1(x) result(y)
        real, intent(in) :: x
        !y = x*x - 4 * cos(x)
        !y = x*x
        !y = sin(x)
        y = 2*x + 4 * sin(x)
    end function

end program lista_4