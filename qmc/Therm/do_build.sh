gfortran -c -O3 -Wall m_random.f90
gfortran -c -O3 -Wall m_position.f90
gfortran -c -O3 -Wall m_output.f90
gfortran -C -O3 -Wall main.f90 m_output.o m_position.o m_random.o -o main.x
