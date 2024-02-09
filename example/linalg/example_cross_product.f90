program demo_cross_product
    use stdlib_linalg, only: cross_product
    implicit none
    real :: a(3), b(3), c(3)
    a = [1., 0., 0.]
    b = [0., 1., 0.]
    c = cross_product(a, b)
    print *, "c = ", c
    print *, "sum(c) = ", sum(c)
    if (abs(sum(c) - 1.0) > 1e-8) error stop
    !c = [0., 0., 1.]
end program demo_cross_product
