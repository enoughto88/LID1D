module parameters_mod

   implicit none
   !Program parameter
   character(len=80)             :: TOPDIR  = '/home/kawashima/RPL1D'
   integer,parameter             :: ITMAX   = 1000000                   !Maximum number of time steps
   integer,parameter             :: ISMP    = 1000                      !Sampling interval
   double precision,parameter    :: UPS     = 1.0d-20                   ![-] Truncation error
   !Grid parameter
   integer,parameter             :: NX      = 192                       ![-]     Number of cells
   double precision,parameter    :: XL      = 2.0d-3                    ![m]     Length of calculation field in x-direction
   double precision,parameter    :: DXL     = XL/dble(NX)               ![m]     Length of cell in x-direction
   double precision,parameter    :: DTNOD   = 1.0d-7                    ![s]     Time step interval
   !Physical parameter
   double precision,parameter    :: PI      = 3.1415926535893d0
   double precision,parameter    :: ECH     = 1.60217733d-19            ![C]     Elementary charge
   double precision,parameter    :: EION    = 15.76d0                   ![eV]    First ionization energy
   double precision,parameter    :: NELE0   = 0.0d0                     ![m-3]   Inflow electron number density
   double precision,parameter    :: TELE0   = 2.0d0                     ![eV]    Inflow electron temperature
   double precision,parameter    :: ILAS0   = 1.0d11                    ![Wm-2]  Inflow laser intensity
   double precision,parameter    :: NNEU0   = 2.45d25                   ![m-3]   Inflow neutral number density
   double precision,parameter    :: UGRD0   = 2.0d3                     ![ms-1]  Initial grid speed
   double precision,parameter    :: DELE0   = 1.83d-1                   ![m2s-1] Electron diffusion coefficient
   double precision,parameter    :: DAMB0   = 3.88d-4                   ![m2s-1] Ambipolar diffusion coefficient
   double precision,parameter    :: KABS0   = 6.44d-22                  ![m2]    Absorption cross section!6.44d-22
   double precision,parameter    :: KABS1   = 1.07d-44                  ![m5]    E-I absorption coefficient!1.07d-44
   double precision,parameter    :: ALPHA0  = 4.0d8                     !Acceleration coefficient for wave speed
   double precision,parameter    :: ALPHA1  = 6.0d4                     !Acceleration coefficient for wave speed
   double precision,parameter    :: ERATIO  = 0.8d0                     !ERATIO
   double precision,parameter    :: TeUV    = 2.0d0                     ![eV] Electron temperature assumed in Pabs

end module parameters_mod


