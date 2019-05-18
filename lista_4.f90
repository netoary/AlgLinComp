include 'methods.f90'
include 'functions.f90'

program lista_4
    use functions
    IMPLICIT NONE
    integer n, i, j
    real :: a, b, tol, y
    real, dimension(2) :: x


    tol = 0.0001
    a = 0
    b = 10
    x = [2, 3]
    !call bissecao(f1, tol, a, b)
    !call metodoDeNewtonOriginal(f1, df1, tol, 10.0, 100, a)
    !call metodoDeNewtonSecante(f1, tol, 10.0, 100, 0.001, a)
    !call interpolacaoInversa(f1, 5E-4, 100, a, [3.0, 5.0, 10.0])
    !call metodoDeNewtonMD(f2, J2, 1E-4, x, 100)
    call metodoDeBroydenMD(f2, B2(), 1E-10, x, 100)
    !call f(b, y)

    write (*,*) x
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

    
    function f2(x) result(y)
        real, dimension(:), intent(in) :: x
        real, allocatable, dimension(:) :: y
        
        allocate(y(2))
        y(1) = x(1) + 2 * x(2) - 2
        y(2) = (x(1) ** 2) + (4 * (x(2) ** 2)) - 4
    end function

    function J2(x) result(J)
        real, dimension(:), intent(in) :: x
        real, allocatable, dimension(:, :) :: J
        
        allocate(J(2,2))
        J(1,:) = [1.0, 2 * x(1)]
        J(2,:) = [2.0, 8 * x(2)]
    end function
    
    function B2() result(J)
        real, allocatable, dimension(:, :) :: J
        allocate(J(2,2))
        J(1,:) = [1.0, 4.0]
        J(2,:) = [2.0, 24.0]
    end function


end program lista_4