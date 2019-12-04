gfortran -c *.FOR
gfortran -c DW4UNIX.F

cd culib4

gfortran -c *.FOR

cd ../culib8

gfortran -c *.FOR

cd ..

gfortran *.o culib8/*.o -o DWUCK4.exe

./DWUCK4.exe < DW4TST.DAT

#diff -y output.txt DW4TST.LIS 
