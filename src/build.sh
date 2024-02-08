set -ex
FC="lfortran --realloc-lhs" cmake .
make fortran_stdlib
cp src/*.mod example/linalg
make VERBOSE=1

./example/linalg/example_eye1		
./example/linalg/example_is_symmetric
# ./example/linalg/example_cross_product	
./example/linalg/example_eye2		
./example/linalg/example_is_triangular
./example/linalg/example_diag1		
./example/linalg/example_is_diagonal		
# ./example/linalg/example_kronecker_product
./example/linalg/example_diag2		
./example/linalg/example_is_hermitian	
./example/linalg/example_outer_product
./example/linalg/example_diag3		
./example/linalg/example_is_hessenberg	
./example/linalg/example_trace
./example/linalg/example_diag4		
./example/linalg/example_is_skew_symmetric

git clean -dfx
