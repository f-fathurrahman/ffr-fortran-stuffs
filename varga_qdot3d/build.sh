gfortran -O3 -c m_qd3d.f90
gfortran -O3 -c xc.f90
gfortran -O3 3dqdot.f90 xc.o m_qd3d.o