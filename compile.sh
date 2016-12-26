#Delete old output data
rm -r output
#Prepare output directory
mkdir output
cd output
mkdir datafile distribution
cd ..
#Delete old nohup.out and ifort files
rm nohup.out
rm -r ifortfile
mkdir ifortfile
#Compile fortran files using OpenMP
ifort -fast -openmp -openmp-report2 -parallel -traceback -CB -gen_interfaces -warn all -o RPL.exe \
program/parameter.f90 \
program/global.f90 \
program/misc.f90 \
program/diffusion.f90 \
program/energy.f90 \
program/radtrans.f90 \
program/main.f90
#Move fortran files
mv *.mod ifortfile/
mv *__genmod.f90 ifortfile/

