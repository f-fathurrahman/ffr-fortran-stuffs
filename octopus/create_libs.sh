LIB_DIRS="basic   hamiltonian  ions  math  poisson  species  sternheimer  td
basis_set  grid  opt_control  scf   states   system       utils"
for dir in $LIB_DIRS
do
    echo $dir
    cd $dir
    echo "liboct_$dir.a"
    rm -fv "liboct_$dir.a"
    ar rcs "liboct_$dir.a" *.o
    cd ../
done

cd main
rm main.o
ar rcs liboct_main *.o
cd ../
