# alpha 3 states


************************************************************************
*****  LANDAU WANG ALGORITHM FOR THE ALPHA-3-STATES MODEL **************
************************************************************************

There are two Wang-Landau programs and both present bugs/problems descibed in 'issues'. There is also a python program ('dos_mean.py'), which calculates the density of states using combinatorics.

The first Wang-Landau code is not able to visit all the energy bins of the histogram for a number of agents N > 14. Nevertheless it delivers results that are in agreement with the exact solution (calculated by 'dos_mean.py') for N=14. The comparison can be observed in 'Images'. The second Wang-Landau code output (which works for any N) can also be seen in 'Images' for N=14 and for N=100.

CONTENTS of the 'alpha3states_wl' folder: --------------------------------

 *  'Code' (subfolder):
      It contains two versions of the Wang-Landau algorithm the first one uses an N-dimensional array with the explicit state of each spin, the second one uses a 3 integers tuple with the number of agents in each state. It also contains a third program ('dos_mean.py') in python, which calculates the density of states using combinatorics. 
      'main.f90'
      'module.f90'
      'input.dat' (output name: nom, size: N (must be even), random seed: iseed,
                   initial modifier: ln_f, alpha value: alpha, number of MC attempts: MCtot,
                   flatness tolerance: tol, modifier threshold: eps, numerical tolerance: eps2))
      'r1279.o', 'ran2.o', secs.o' (random generator files)
      
      'main_2.f90'
      'module_2.f90'
      'input_2.dat'  (identical to input.dat)
      
      'dos_mean.py'

 *  'Images' (subfolder):
      'ln_g_N10.png', 'ln_g_N12.png', 'ln_g_N14.png' (Plots of the ln(states density(E)) .vs. E)
      'wl_N14.png' (A comparison of the results obtained for N=14 with the three programs)
      'wl_N100.png' (A comparison of the results obtained for N=100 with the 'dos_mean.py' and the second Wang-Landau program)
      'exe_times.png' (Execution times for the different programs for diffrerent values of N and alpha = 0.2)
    
 *  'Scripts'
      'ln_g.gnu' (For the output visualization of the first Wang-Landau code) Change the filename and 'N = 10' 
                 for the desired size
      'ln_g_N14_a0.2.gnu' and 'ln_g_N100_a0.2.gnu' (For the output visualization of the second Wang-Landau code)
                                                    Each one is specific for a size N and it selects the
                                                    appropiate yrange so the figure can be observed.
                                                    
      'dos_ln_g_N14_a0.2.gnu' and 'dos_ln_g_N100_a0.2.gnu' (For the output visualization of the python code)


COMPILING INSTRUCTIONS (from local folder 'Code'): ------------------------------------------------

gfortran -g -fcheck=all -Wall -fbacktrace module.f90 main.f90 r1279.o ran2.o secs.o -o main.exe
./main.exe (execution)  (Same for main_2.f90)

    
The outputs are  .csv files with 2 columns. The first one are the energy bins, 
and the second one correspond to the logarithm of the states density corresponding to 
the given energy.

