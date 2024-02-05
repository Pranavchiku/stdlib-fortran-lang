# #!/bin/bash

set -ex

# git clean -dfx

# FC=gfortran cmake .
# make fortran_stdlib
# make
# ./example/math/example_clip_integer
# ./example/math/example_clip_real
# ./example/math/example_diff
# ./example/math/example_gcd
# ./example/math/example_linspace_complex
# ./example/math/example_linspace_int16
# ./example/math/example_logspace_complex
# ./example/math/example_logspace_int
# ./example/math/example_logspace_rstart_cbase
# ./example/math/example_math_all_close
# ./example/math/example_math_arange
# ./example/math/example_math_arg
# ./example/math/example_math_argd
# ./example/math/example_math_argpi
# ./example/math/example_math_is_close


# git clean -dfx

FC=lfortran cmake .
make fortran_stdlib
cp src/*.mod example/math
make VERBOSE=1
./example/math/example_clip_integer
./example/math/example_clip_real
./example/math/example_diff
./example/math/example_gcd
./example/math/example_linspace_complex
./example/math/example_linspace_int16
./example/math/example_logspace_complex
./example/math/example_logspace_int
./example/math/example_logspace_rstart_cbase
./example/math/example_math_all_close
./example/math/example_math_arange
./example/math/example_math_arg
./example/math/example_math_argd
./example/math/example_math_argpi
./example/math/example_math_is_close

# git clean -dfx
