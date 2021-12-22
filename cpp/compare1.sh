echo "Compare B-spline curve with generation points"
rm *.out
rm *.in

g++ num_bspline.cpp -o num_bspline.out -O2
g++ prec_bspline.cpp -o prec_bspline.out -O2
g++ bspline.cpp  -o bspline.out -O2

python3 generate.py $2
./num_bspline.out $1
./prec_bspline.out $1
./bspline.out $1

python3 compare.py