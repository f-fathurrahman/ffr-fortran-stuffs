FC="gfortran"
#FC="flang"
$FC -c -O3 random.f90
$FC -c -O3 -Wall position.f90
$FC -O3 -Wall main.f90 position.o random.o -o main.x
