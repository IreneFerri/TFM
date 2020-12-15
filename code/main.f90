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
! After obtaining the states density array it uses it to calculate the
! canonical probability, the partition function and the average energy for 
! different temperatures and stores them in a file
! 
! -----------------------------------------------------------------------------
! -----------------------------------------------------------------------------
!      MAIN
! -----------------------------------------------------------------------------
      program landau_wang
      use my_func
!
      implicit none
!
      integer*16, dimension(:), allocatable :: spins , H
      real*16, dimension (:), allocatable :: log_g, g, E_array
      real*16, dimension (:), allocatable :: final_E_array, unique
      real :: TIME1,TIME2
      real*16 :: min_E_val, max_E_val
      real*16 :: E0, Ef,  alpha2, ln_f 
      integer*16 :: chosen_agent, i, j, k, N, istep, counter, iseed, check_H
      integer*16 :: h_zero, b
      integer*16 :: E_bin, new_E_bin 
      integer*8 :: MCtot
      real*16 :: DE,energy, new_energy, p_new_energy
      integer*16 :: num0, num1, num_1, new_num0, new_num1, new_num_1, n1, n_1, n0
      integer*16 :: spin, new_spin, p_g_dim , g_dim, chosen_agent_spin
      integer*16 :: N2
      real*16 ::  alpha, tol, eps, eps2
      character*2:: nom
      character*30:: date, x1, x2
! -----------------------------------------------------------------------------
!   READ DATA
! ------------------------------------------------------------------------------
      NAMELIST/DADES/nom,N,iseed,ln_f,alpha,&
                      MCtot, check_H, tol, eps,eps2 ! Initial data
      OPEN(UNIT=10,FILE='input.dat')        ! Input
      READ(10,DADES)
      CLOSE(10)
!
      write(x1,'(I4.4)') N  ! converting integer to string 
      write(x2,'(F4.2)') alpha  ! converting float to string 
!
!     Calculate DOS preliminar dimension ---------------------------------------
! -------------------------------------------------------------------
      p_g_dim = (((N+2)*(N+1)/2) - (N/2+1))/2 + (N/2+1) ! Not valid if N is odd; it considers 
!                                                     half triangle and the elements along the height
      allocate(spins(N))  ! It stores the state of every agent
      allocate(E_array(p_g_dim))  ! It stores the values of the energy corresponding to half triangle
      allocate(unique(p_g_dim))
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
!     LANDAU - WANG 
! -------------------------------------------------------------------
        call setr1279(iseed)   
! --- Comment or uncomment next lines to start from T = 0 or T --> inf ---------
!
! --- Initial Configuration. Fully ordered -1 ----------------------- 
      call initial_T0(N, spins, num1, num_1, num0)
! --- Initial Configuration. Random spins -1, 0 or +1 ----------------------- 
!      call initial_inf(N, spins, num1, num_1, num0)
! -------------------------------------------------------------------------
      energy=-0.5*((float(num1-num_1))**2+alpha2*float(num0**2)-float(num1)-float(num_1)-alpha2*float(num0)) ! Mean-field Hamiltonian
! -- Energy_array and limits ----------------------------------------------
      i = 0
      do n0 = 0, N
        do n1 = 0, int(N-n0)/2
          i = i+1
          n_1 = N - n1 - n0
          E_array(i) = -0.5*((float(n1-n_1))**2+alpha2*float(n0**2)-float(n1)-float(n_1)-alpha2*float(n0))
        enddo
      enddo
      E0 = -float(N*(N-1))/2.0 - eps2  ! This is to include numerical errors from the computer
      Ef = maxval(E_array) + eps2
! -- Eliminate degeneration and sort values ---------------------------------
      min_E_val = minval(E_array)-1
      max_E_val = maxval(E_array)
      i = 0
      do while (min_E_val<max_E_val)
        i = i+1
        min_E_val = minval(E_array, mask=E_array>min_E_val)
        unique(i) = min_E_val
      enddo
      allocate(final_E_array(i), source=unique(1:i))  ! It stores all possible values of energy, non degenerated
      g_dim = size(final_E_array)

      do i = 1, g_dim
        print*, final_E_array(i)
      enddo

!                                                       and sorted from smallest to largest
! - Corresponding energy index (E_bin) ---------------------------------------------
!    
      E_bin = find_index(g_dim,final_E_array, energy, 1, g_dim)  ! Find the bin corresponding to the current energy

      print*, 'ALPHA = ', alpha, 'E0 = ', E0, 'Ef = ', Ef
      print*, 'Initial energy = ', energy, 'E_bin = ', E_bin
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! --- Initialize energy histogram and states density arrays -------------
      allocate(H(g_dim))
      allocate(log_g(g_dim))
      allocate(g(g_dim))
!    
      do i = 1, g_dim
        H(i) = 0
        log_g(i) = 0.0
      end do
! ------------------------------------------------------------------
! --  MARKOV CHAIN (non-stationary) ---------------------------------------------------------   
! ------------------------------------------------------------------
      istep = 0
      counter = 0
      b=0
      DO while ((ln_f > eps).and.(istep < MCtot))
        do k = 1, check_H ! ------------- HISTOGRAM LOOP
          istep = istep + 1
          DO i = 1,N         ! Each MC has N flip attempts -------------
!            print*, spins, energy, E_bin, '**',chosen_agent, b
            counter = counter + 1
            call system_evolves(N, spins, alpha2, num1, num_1, num0, & ! input
                                new_num1, new_num_1, new_num0, DE,chosen_agent, new_spin, b) ! Output
!  Accept or reject ------------------------------------------------------------------------------
            chosen_agent_spin = spins(chosen_agent)
            p_new_energy = energy + DE
            new_E_bin =find_index(g_dim, final_E_array,p_new_energy,1,g_dim)  !Find the bin corres    ponding to the new energy
            call acceptance(g_dim,log_g,Ef,E0,&
                            num1, num_1, num0, new_num1, new_num_1, new_num0, &
                            chosen_agent_spin, new_spin, E_bin, new_E_bin, energy, DE, & ! Input
!
                             num1, num_1, num0, new_energy, spin, new_E_bin)  ! Output
!
            energy = new_energy; E_bin = new_E_bin  ! Actualize energy          
            spins(chosen_agent) = spin
            log_g(E_bin) = log_g(E_bin) + ln_f
            H(E_bin) = H(E_bin) + 1
          end do  ! Ends steps (1 to N) loop ------------------------------
        end do !  Ends histogram check loop ------------------------------
        print*, H
!   Check flatness ----------------------------------------------------------------------------
        if (minval(H) > tol*(counter/g_dim)) then   ! Just considering the histogram's bottom
          print*, 'Flat!'
          ln_f = 0.5*ln_f ! f -> sqrt(f) (actualize f)
          do j = 1, g_dim ! Reinitialize histogram
            H(j) = 0
          end do
          counter = 0
        endif
      end do  ! ends imc loop --------------------------------------------------
! Check non-visited configurations           
      print*, 'Non-visited configurations ****************'
      do h_zero = 1, g_dim
        if (H(h_zero) == 0) then
          print*, 'index: ',h_zero, 'energy =',final_E_array(h_zero)
        endif
      enddo
!
      print*,'Number of steps needed:', istep, 'ln_f =', ln_f
!----------------------------------------------------------------------
! Write the ln(g(E))
!----------------------------------------------------------------------
      OPEN(UNIT=11,FILE=nom//"_N"//trim(x1)//"_a"//trim(x2)//"log.csv")  ! ln(g(E)) file
      do i = 1, g_dim
        write(11,*) final_E_array(i), log_g(i)
      enddo
      close(11)
!----------------------------------------------------------------------
!  Normalization 
!------------------------------------------------------------------
!      max_log_g = maxval(log_g)
!      log_g = log_g - max_log_g
!      g = exp((log_g))/sum(log_g)
!
!      C = 2.0/g(1)  ! g(E_min = 2) ! For alpha<1 
!      g = C*g
!
      CALL CPU_TIME(TIME2)            ! Time control
      CALL FDATE(DATE)
      WRITE(*,*) DATE
      WRITE(*,*) 'CPUTIME=',TIME2-TIME1
      END program landau_wang 

