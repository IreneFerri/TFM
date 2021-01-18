!
!    
! -------------------------------------------------------------------------
!                              landau_wang.f90
! --------------------------------------------------------------------------
!
! Monte Carlo simulator of a system of spins that can take 3 values 
! corresponding to the vectors (-1,0),(0,alpha) and (1,0) 
! respectively, where alpha is a real parameter. Energy is defined as the
! negative sum of dot products over all pairs (fully connected graph).
! The program proposes sequencial changes and uses the histogram entropic
! Landau-Wang algorithm to visit all states and find the first order transition
! of the system for alpha above the tricritical point.
! It stores the energies and the logarithm of the states density array in an 
! external file.
!
! -----------------------------------------------------------------------------
! -----------------------------------------------------------------------------
!      MAIN
! -----------------------------------------------------------------------------
      program landau_wang
      use my_func   ! Functions and subroutines in module.f90
!
      implicit none
!
      integer*16, dimension(:), allocatable :: H  ! Histogram of visited configurations
      real*16, dimension (:), allocatable :: log_g      ! Logarithm of the DOS for every possible combination of agents
      real*16, dimension (:), allocatable :: final_log_g   ! Logarithm of the DOS for every different energy
      real*16, dimension (:), allocatable :: E_array   ! Energy corresponding to every possible combination of agents
      real*16, dimension (:), allocatable :: final_E_array, unique_E    !  Arrays for every different energy of the system (sorted)
      real :: TIME1,TIME2   ! Efficiency control
      real*16 :: min_E_val, max_E_val, deg  ! Auxiliar variables 
      real*16 :: E0, Ef  ! Maximum and minimum possible values of energy for the system
      integer*16 :: N, N2   ! Number of agents and its square
      real*16 :: alpha, alpha2, ln_f   ! Neutrality value, its square and the modifier for the log(DOS) 
      integer*16 :: i, j, k, istep, counter, iseed     ! Loop indexes, random gen seed  
      integer*16 :: E_bin, new_E_bin   ! Indexes of the energy in the E_array 
      integer*8 :: MCtot, check_H  ! Number of total steps and integer*16 :: i, j, k, N, istep, counter, iseed, check_H   ! Loop indexes, random gen seed and 
      real*16 :: energy, new_energy, DE    ! Current state energy, proposed state energy and difference of energy between them
      integer*16 :: num0, num1, num_1, new_num0, new_num1, new_num_1, n1, n_1, n0  ! Number of agents in each state (microstate)
      integer*16 :: g_dim , final_g_dim   ! Number of combinations of agents and number of non-degenerated values of the energy 
      real*16 ::  tol, eps, eps2 ! Tolerance forhistogram flatness, tolerance for the final value of the modifier, tolerance for the maximum and minimum values of energy
      real*16 :: prob_flip  ! Probability of accepting a given spin-flip proposal
      character*4:: nom    ! Output file name
      character*30:: date, x1, x2   ! Auxiliar string variables
! -----------------------------------------------------------------------------
!   READ DATA
! ------------------------------------------------------------------------------
      NAMELIST/DADES/nom,N,iseed,ln_f,alpha,&
                      MCtot, check_H, tol, eps, eps2! Initial data
      OPEN(UNIT=10,FILE='input.dat')        ! Input
      READ(10,DADES)
      CLOSE(10)
!
      write(x1,'(I4.4)') N  ! converting integer to string 
      write(x2,'(F4.2)') alpha  ! converting float to string 
!
!     Calculate DOS preliminar dimension ---------------------------------------
! -------------------------------------------------------------------
      g_dim = ((N+2)*(N+1)/2) ! Number of possible combinations of agents 
      allocate(E_array(g_dim))  ! It stores the values of the energy corresponding all possible combinations
      allocate(unique_E(g_dim))  ! It will copy the non-degenerated energy values
!     
      N2 = N*N
      alpha2 = alpha*alpha
! -----------------------------------------------------------------------------
! -----------------------------------------------------------------------------
!     Print input data
      write(*,*) 'Number of spins = ', N
      write(*,*) 'MCS = ', MCtot
      write(13,*) '# N = ', N, 'TOL = ', tol, 'Steps = ', MCtot,&
                   'epsilon = ', eps
! -----------------------------------------------------------------------------
! -------------------------------------------------------------------------
      call CPU_TIME(TIME1)          ! Time control
! -------------------------------------------------------------------
!     INITIAL CONFIGURATION 
! -------------------------------------------------------------------
      call setr1279(iseed)   
! --- Comment or uncomment next lines to start from T = 0 or T --> inf ---------
!
! --- Initial Configuration. Fully ordered +1 ----------------------- 
      num0 = 0 ; num1 = 0 ; num_1 = N
! --- Initial Configuration. Random spins -1, 0 or +1 ----------------------- 
!      num1 = int(N/3) ; num0 = num1 ; num_1 = N-num0-num1
      energy=-0.5*((float(num1-num_1))**2+alpha2*float(num0**2)-float(num1)-float(num_1)-alpha2*float(num0)) ! Mean-field Hamiltonian
! -- Energy_array and limits ----------------------------------------------
      i = 0
      do n0 = 0, N
        do n1 = 0, int(N-n0)
          i = i+1
          n_1 = N - n1 - n0
          E_array(i) = -0.5*((float(n1-n_1))**2+alpha2*float(n0**2)-float(n1)-float(n_1)-alpha2*float(n0))
        enddo
      enddo
      E0 = -float(N*(N-1))/2.0 - eps2 ! Ground State. eps2 is to include numerical errors from the computer 
      Ef = maxval(E_array) + eps2 ! Maximum possible value of energy
!                                                       
! - Corresponding energy index (E_bin) ---------------------------------------------
!    
      E_bin = find_index(N, num0, num1)
      print*, 'ALPHA = ', alpha, 'E0 = ', E0, 'Ef = ', Ef
      print*, 'Initial energy = ', energy, 'E_bin = ', E_bin
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! --- Initialize energy histogram and states density arrays -------------
      allocate(H(g_dim))    ! Histogram
      allocate(log_g(g_dim))  ! Logarithm of the density of states (DOS)
!    
      do i = 1, g_dim
        H(i) = 0
        log_g(i) = 0.0
      end do
! ------------------------------------------------------------------
! --  WANG-LANDAU MARKOV CHAIN (non-stationary) ---------------------------------------------- 
! ------------------------------------------------------------------
      istep = 0
      counter = 0
      do while ((ln_f > eps).and.(istep < MCtot))
        do k = 1, check_H ! ------------- HISTOGRAM LOOP
          istep = istep + 1
          do i = 1,N         ! Each MC has N flip attempts -------------
            counter = counter + 1
            call system_evolves(N, alpha2, num_1, num0, num1, &   ! Single Spin flip proposal 
                                DE, new_num_1, new_num0, new_num1) 
!
!  ACCEPT OR REJECT  ------------------------------------------------------------------------------
            new_energy = energy + DE
            new_E_bin = find_index(N, new_num0, new_num1)  ! Find the bin corresponding to the new energy
            prob_flip = log_g(E_bin) - log_g(new_E_bin)   ! The algorithm accepts the proposed flip just if the new configuration has been less visited than the current one
            call acceptance(prob_flip,Ef,E0,&
                            num1, num_1, num0, new_num1, new_num_1, new_num0, &
                            E_bin, new_E_bin, energy, DE, & ! Input
!
                             num1, num_1, num0, new_energy, new_E_bin)  ! Output
!
            energy = new_energy; E_bin = new_E_bin  ! Actualize energy and energy bin     
            log_g(E_bin) = log_g(E_bin) + ln_f
            H(E_bin) = H(E_bin) + 1
          end do  ! Ends steps (1 to N) loop ------------------------------
        end do !  Ends histogram check loop ------------------------------
!
!   CHECK FLATNESS   ----------------------------------------------------------------------------
        if (minval(H) > tol*(counter/g_dim)) then   ! Just considering the histogram's bottom
          print*, 'Flat!'
          ln_f = 0.5*ln_f ! f -> sqrt(f) (actualize f)
          do j = 1, g_dim 
            H(j) = 0     ! Reinitialize histogram
          end do
          counter = 0
        endif
      end do  ! ends imc loop --------------------------------------------------
!
      print*,'Number of steps needed:', istep, 'ln_f =', ln_f
!
!  ------------------------------------------------------------------ 
! -- ELIMINATE DEGENERACY AND SORT ENERGY VALUES  ---------------------------------
!  ------------------------------------------------------------------
      min_E_val = minval(E_array)-1
      max_E_val = maxval(E_array)
      i = 0
      do while (min_E_val<max_E_val)
        i = i+1
        min_E_val = minval(E_array, mask=E_array>min_E_val)
        unique_E(i) = min_E_val
      enddo
      allocate(final_E_array(i), source=unique_E(1:i))  ! It stores all possible values of energy, non degenerated
      final_g_dim = size(final_E_array)                     ! and sorted from smallest to largest
      allocate(final_log_g(final_g_dim))
!
!  -- Sum the log of states corresponding to the same energy -------------------------   
      do i = 1, final_g_dim   ! Run over the non-degenerated energy array
        deg = 0     
        do j = 1, g_dim       ! Run over the degenerated energy array
          if (E_array(j) == final_E_array(i)) then
            deg = deg + log_g(j)
          endif
        enddo
        final_log_g(i) = deg
      enddo
!
!----------------------------------------------------------------------
! WRITE the ln(g(E))
!----------------------------------------------------------------------
      OPEN(UNIT=11,FILE=nom//"_N"//trim(x1)//"_a"//trim(x2)//"log.csv")  ! ln(g(E)) file
      do i = 1, final_g_dim
        write(11,*) final_E_array(i), final_log_g(i)
      enddo
      close(11)
!----------------------------------------------------------------------
!
      CALL CPU_TIME(TIME2)            ! Time control
      CALL FDATE(DATE)
      WRITE(*,*) DATE
      WRITE(*,*) 'CPUTIME=',TIME2-TIME1
      END program landau_wang 

