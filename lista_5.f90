include 'methods.f90'
include 'functions.f90'

program lista_5
    use functions
    IMPLICIT NONE
    real :: a, b
    real :: area

    a = 0
    b = 10
    
    call integracaoQuadratura(f1, "10", a, b, area)

    write (*,*) area
    call sleep(10000)

contains
    real function f1(x) result(y)
        real, intent(in) :: x
        y = x
    end function
end program lista_5