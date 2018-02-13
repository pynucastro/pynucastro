# executable
#FC := pgf95
#FC := pgf95_local
FC := gfortran
#FC := Cray

ifdef ACC
   acc_suf := .acc
endif

ifdef OMP
   omp_suf := .omp
endif

suf=$(FC)$(omp_suf)$(acc_suf)