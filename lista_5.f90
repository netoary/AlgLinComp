include 'methods.f90'
include 'functions.f90'

program lista_5
    use functions
    IMPLICIT NONE
    REAL(8) :: a, b
    REAL(8) :: area

    a = 0
    b = 5
    
    call integracaoPolinomial(f1, "7", a, b, area)

    write (*,*) area
    call sleep(10000)

contains
REAL(8) function f1(x) result(y)
    REAL(8), intent(in) :: x
    REAL, PARAMETER :: Pi = 3.1415927
    y = 1/sqrt(2*Pi) * exp(-1.0d0/2.0d0 * x**2)
end function

end program lista_5