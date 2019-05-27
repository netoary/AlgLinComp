include 'methods.f90'
include 'functions.f90'

program lista_4
    use functions
    IMPLICIT NONE
    REAL(8) :: a, b, tol
    REAL(8), dimension(2) :: x


    !TESTE SLIDE
    tol = 0.0001
    a = 0
    b = 10
    !call bissecao(ft, tol, a, b)
    !call metodoDeNewtonOriginal(ft, dft, tol, 10.0, 100, a)
    !call metodoDeNewtonSecante(ft, tol, 10.0, 100, 0.001, a)
    !call interpolacaoInversa(ft, 5E-3, 1000, a, [1.0, 5.0, 10.0])
    !write (*,*) a


    !Questão 1
    tol = 1E-5
    a = 500
    b = 1000
    !call bissecao(f1, tol, a, b)
    !call metodoDeNewtonOriginal(f1, df1, tol, 500.0d0, 100, a)
    !call metodoDeNewtonSecante(f1, tol, 10.0d0, 100, 0.001d0, a)
    !call interpolacaoInversa(f1, tol, 1000, a, [100.0d0, 600.0d0, 1000.0d0])
    write (*,*) a

    !Questão 2
    tol = 1E-5
    a = 1
    b = 0
    !call bissecao(f2, tol, a, b)
    !call metodoDeNewtonOriginal(f2, df2, tol, 0.0, 100, a)
    !call metodoDeNewtonSecante(f2, tol, 10.0, 100, 0.001, a)
    !call interpolacaoInversa(f2, tol, 100, a, [1.0, 5.0, 10.0])
    !write (*,*) a

    !Questão 3
    tol = 1E-5
    a = -10
    b = 10
    x = [0, 1]
    !call metodoDeNewtonMD(f5, J2, tol, x, 100)
    !call metodoDeBroydenMD(f5, B2(), tol, x, 100)
    !call minimosQuadrados(fQ, JQ, tol, x, 100)
    !write (*,*) a

    call sleep(10000)

contains
    REAL(8) function ft(x) result(y)
        REAL(8), intent(in) :: x
        y = x**2 - 4*cos(x)
    end function

    REAL(8) function dft(x) result(y)
        REAL(8), intent(in) :: x
        y = 2*x + 4*sin(x)
    end function


    REAL(8) function f1(x) result(y)
        REAL(8), intent(in) :: x
        REAL(8) :: g, k
        
        g = 9.806
        k = 0.00341
        y = log10(cosh(x * sqrt(g * k))) - 50.0
    end function

    REAL(8) function df1(x) result(y)
        REAL(8), intent(in) :: x
        y = (sqrt(9.806*0.00341)*sinh(x*sqrt(9.806*0.00341)))/cosh(x*sqrt(9.806*0.00341))
    end function

    REAL(8) function f2(x) result(y)
        REAL(8), intent(in) :: x
        y = 4*cos(x)-exp(2*x)
    end function

    REAL(8) function df2(x) result(y)
        REAL(8), intent(in) :: x
        y = 4*(-sin(x))-2*exp(2*x)
    end function


    function f5(x) result(y)
        REAL(8), dimension(:), intent(in) :: x
        REAL(8), allocatable, dimension(:) :: y

        allocate(y(2))
        y(1) = x(1) + 2 * x(2) - 2
        y(2) = (x(1) ** 2) + (4 * (x(2) ** 2)) - 4
    end function

    function J2(x) result(J)
        REAL(8), dimension(:), intent(in) :: x
        REAL(8), allocatable, dimension(:, :) :: J

        allocate(J(2,2))
        J(1,:) = [1d0, 2 * x(1)]
        J(2,:) = [2d0, 8 * x(2)]
    end function

    function B2() result(J)
        REAL(8), allocatable, dimension(:, :) :: J
        allocate(J(2,2))
        J(1,:) = [1.0, 4.0]
        J(2,:) = [2.0, 24.0]
    end function

    function fQ(x) result(y)
        REAL(8), dimension(:), intent(in) :: x
        REAL(8), allocatable, dimension(:) :: y

        allocate(y(3))
        y(1) = x(1) + 2 * x(2) - 2
        y(2) = (x(1) ** 2) + (4 * (x(2) ** 2)) - 47
        y(3) = (x(1) ** 2) + (4 * (x(2) ** 2)) - 4
    end function

    function JQ(x) result(J)
        REAL(8), dimension(:), intent(in) :: x
        REAL(8), allocatable, dimension(:, :) :: J
        allocate(J(3,2))
        J(1,:) = [1.0, 4.0]
        J(2,:) = [2.0, 24.0]
        J(3,:) = [2.0, 24.0]
    end function

end program lista_4