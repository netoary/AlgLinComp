include 'methods.f90'
include 'functions.f90'

program lista_4
    use functions
    IMPLICIT NONE
    REAL(8) :: a, b
    REAL(4) :: tol
    REAL(8), dimension(2) :: x
    REAL(8) :: c(3,3), ci(3,3)


    !TESTE SLIDE
    !tol = 0.0001
    !a = 0
    !b = 10
    !call bissecao(ft, tol, a, b)
    !call metodoDeNewtonOriginal(ft, dft, tol, 10.0, 100, a)
    !call metodoDeNewtonSecante(ft, tol, 10.0, 100, 0.001, a)
    !call interpolacaoInversa(ft, 5E-3, 1000, a, [1.0, 5.0, 10.0])
    !write (*,*) a


    !Questão 1
    ! tol = 1E-5
    ! a = 500
    ! b = 1000
    !call bissecao(f1, tol, a, b)
    !call metodoDeNewtonOriginal(f1, df1, tol, 500.0d0, 100, a)
    !call metodoDeNewtonSecante(f1, tol, 10.0d0, 100, 0.001d0, a)
    !call interpolacaoInversa(f1, tol, 1000, a, [100.0d0, 600.0d0, 1000.0d0])
    ! write (*,*) a

    !Questão 2
    ! tol = 1E-5
    ! a = 1
    ! b = 0
    !call bissecao(f2, tol, a, b)
    !call metodoDeNewtonOriginal(f2, df2, tol, 0.0, 100, a)
    !call metodoDeNewtonSecante(f2, tol, 10.0, 100, 0.001, a)
    !call interpolacaoInversa(f2, tol, 100, a, [1.0, 5.0, 10.0])
    !write (*,*) a

    !Questão 3
    tol = 1E-4
    x = [2,5]
    call metodoDeNewtonMD(f3, J3, tol, x, 100)

    ! a(1,:) = [0.008333333333333333d0, -0.14166666666666666d0, 0.25d0]
    ! a(2,:) = [0.008333333333333333d0, 0.10833333333333334d0, -0.25d0]
    ! a(3,:) = [-0.016666666666666666d0, 0.5333333333333333d0, 0d0]
    ! b(:) = [17, 0, -2]
    ! write (*,*) a
    ! x = matmul(a, b)

    !call metodoDeBroydenMD(f3_, J2(x), tol, x, 20)
    !call minimosQuadrados(fQ, JQ, tol, x, 100)

    write (*,*) x

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


    function f3_(inputs) result(funcs)
        REAL(8), dimension(:), intent(in) :: inputs
        REAL(8), allocatable, dimension(:) :: funcs

        REAL(8) :: x, y, z

        x = inputs(1)
        y = inputs(2)

        allocate(funcs(2))
        funcs(1) = x + 2*y -2
        funcs(2) = x**2 + 4*y**2 -4
    end function

    function f3(inputs) result(funcs)
        REAL(8), dimension(:), intent(in) :: inputs
        REAL(8), allocatable, dimension(:) :: funcs

        REAL(8) :: x, y, z

        x = inputs(1)
        y = inputs(2)
        z = inputs(3)

        allocate(funcs(3))
        funcs(1) = 16.0d0 * x**4 + 16 * y**4 + z**4 - 16
        funcs(2) = x**2 + y**2 + z**2 - 3.0d0
        funcs(3) = x**3 - y - z - 1.0d0
    end function

    function J3(inputs) result(J)
        REAL(8), dimension(:), intent(in) :: inputs
        REAL(8), allocatable, dimension(:, :) :: J

        REAL(8) :: x, y, z

        x = inputs(1)
        y = inputs(2)
        z = inputs(3)

        allocate(J(3,3))
        J(1,1) = 64 * (x**3)
        J(2,1) = 64 * (y**3)
        J(3,1) = 4 * (z**3)

        J(1,2) = 2 * x
        J(2,2) = 2 * y
        J(3,2) = 2 * z

        J(1,3) = 3 * x**2
        J(2,3) = -1
        J(3,3) = 1

        ! J(1,1) = 64 * (x**3)
        ! J(1,2) = 64 * (y**3)
        ! J(1,3) = 4 * (z**3)

        ! J(2,1) = 2 * x
        ! J(2,2) = 2 * y
        ! J(2,3) = 2 * z

        ! J(3,1) = 3 * x**2
        ! J(3,2) = -1
        ! J(3,3) = 1


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
        J(1,:) = [1d0, 2d0]
        J(2,:) = [2 * x(1), 8 * x(2)]
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