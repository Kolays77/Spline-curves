echo "Compare B-spline curve without generation points. Files needed cv_x.in, cv_y.in"
g++ num_bspline.cpp -o num_bspline.out 
g++ prec_bspline.cpp -o prec_bspline.out
g++ bspline.cpp  -o bspline.out 

./num_bspline.out $1
./prec_bspline.out $1
./bspline.out $1

python3 compare.py