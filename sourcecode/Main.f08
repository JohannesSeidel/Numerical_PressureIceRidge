program Numerical_Ridge
  ! change log:
  ! 29.12. -> added termination if condition for option 1 in the main loop
  !
  use iso_fortran_env
  use omp_lib
  use mod_Constants
  use mod_Global
  use mod_ImportExport
  use mod_Particle
  use mod_PredCorr
  use mod_OverlapComputation
  use mod_Force
  use mod_Output
  use mod_Neighbour
  use mod_Functions

  !===========!
  ! variables !
  !===========!
  implicit none
  double precision   :: length, width, thickness, density, youngs
  double precision   :: coh_coeff, gamma_n, mu
  integer            :: n_ice, n_elem, n_struct
  integer, parameter :: n_wall = 3
  double precision   :: dt, frac, m_min, Y_max
  integer            :: wcount, winterval !write count, write interval
  integer            :: fps
  integer            :: option
  integer            :: i
  integer            :: t1, t2, t3, clock_rate, clock_max
  double precision   :: calc_time
  logical            :: termination
  character (200)    :: dummy, path


  call getcwd (path)
  if (path(1:1) == '/') then !path starts with /
    call execute_command_line ('clear')
    call execute_command_line ('(cd output && rm -r vtk/*.vtk *.txt csv/*.csv)')
    write (*,*) "Running on Linux!"
    write (*,*)
  else
    call execute_command_line ('del output\vtk\*.vtk output\*.txt output\csv\*.csv')
    write (*,*) "Running on Windows!"
    write (*,*)
  end if

  write (*, 1001) "!=================!"
  write (*, 1001) "! Numerical Ridge !"
  write (*, 1001) "!=================!"
  write (*, 1001)
  write (*, '(1x, 2a)') "compiler version: ", compiler_version()
  write (*, '(1x, 2a)') "compiler options: ", compiler_options()
  write (*, 1001)
  write (*, 1001) "What do you want to do?"
  write (*, 1001) "1) New ridge"
  write (*, 1001) "2) Punch test"
  write (*, 1001) "3) Model test"
  write (*, 1001)
  1001 format (1x, a)
  read (*, *) option


  select case (option)
  case (1) ! New Ridge
    !================!
    ! read from file !
    !================!
    open (20, file = './Simulation_Settings.txt')
    read (20, *) dummy, rdge%width, dummy, rdge%keel_width, &
      dummy, rdge%keel_height, dummy, rdge%length, dummy, rdge%porosity, &
      dummy, length, dummy, width, dummy, thickness, dummy, density, dummy, youngs, &
      dummy, coh_coeff, dummy, gamma_n, dummy, mu, dummy, frac
    close (20)

    ! Calculate number of needed ice elements
    n_ice = nint (rdge%porosity * (0.5 * (rdge%width + rdge%keel_width) * rdge%keel_height * rdge%length) & 
      & / ( length * width * thickness ) ) !length, width, thickness are mean values

    !n_ice = 1
    write (*,'(1x,a,i6)') "Die Anzahl der Eisbloecke ist ", n_ice

    !================!
    ! Initialization !
    !================!
    n_struct = 2
    n_elem   = n_ice + n_wall + n_struct

    !==========================!
    ! Initial element position !
    !==========================!
    call InitializeParticles (density, length, mu, n_ice, n_struct, n_wall, thickness, width, youngs)
    call MakeStructure (n_ice, n_struct, n_wall, option)
    call UpdateBoundingBoxes (n_elem)
    call InitialSort
  case (2) ! punch test
    n_struct = 1
    call ImportDomain (coh_coeff, frac, gamma_n, mu, n_ice, n_struct, n_elem)
    call MakeStructure (n_ice, n_struct, n_wall, option)
    call UpdateBoundingBoxes (n_elem)
    call InitialSort
  case (3) ! model test
    n_struct = 1
    call ImportDomain (coh_coeff, frac, gamma_n, mu, n_ice, n_struct, n_elem)
    call MakeStructure (n_ice, n_struct, n_wall, option)
    call UpdateBoundingBoxes (n_elem)
    call InitialSort
  end select

  !===============!
  ! Vorrechnungen !
  !===============!

  m_min = minval ( [(ice(i)%mass, i = 1, n_ice)] )
  Y_max = maxval ( [(ice(i)%youngs, i = 1, n_ice)] )
  dt    = frac * 2.0 * sqrt (m_min / (2.0 * Y_max))
  write (*, '(1x,a,d8.2,a)') "dt = ", dt, " [s]"

  fps       = 5 !graphical output files per second
  winterval = nint (1.0 / (dt * fps))
  wcount    = 1

  !==================!
  ! graphical output !
  !==================!

  !initial configuration output at time step 0
  step = 0
  call Graphics (n_elem, n_ice, n_struct, n_wall)

  open (25, file = './output/txt_debugging.txt', status = 'replace')
  open (26, file = './output/csv_debugging.csv', status = 'replace')
  if (option /= 1) then
    open (30, file = "./output/csv/xForce.csv", status = "replace")
    open (31, file = "./output/csv/yForce.csv", status = "replace")
    open (32, file = "./output/csv/zForce.csv", status = "replace")
  end if
  open (33, file = "./output/csv/kinetic_energy.csv", status = "replace")

  termination = .false.
  call system_clock (t1, clock_rate, clock_max)

  do while (.not. termination) 
    step = step + 1

    call Predictor (dt, n_ice)
    call UpdateParticles (dt, n_ice, n_elem, n_struct, option)
    call ForceComputation (coh_coeff, gamma_n, dt, n_elem, n_ice, n_struct, n_wall, option)
    call Corrector (dt, n_ice, n_struct, option)

    if (step == wcount * winterval) then 
      call Graphics (n_elem, n_ice, n_struct, n_wall)
      wcount = wcount + 1
      write (*, '(1x,a,f7.4,a)') "Calculated real time: ", step*dt, " s"
    end if

    call TerminationCheck (dt, n_elem, option, termination)
  end do

  call system_clock (t2, clock_rate, clock_max)
  calc_time = real(t2-t1)/real(clock_rate)

  close (25); close (26)
  if (option /= 1) then
    close (30); close (31); close (32)
  end if
  close (33)

  call Output (option, frac, dt, step*dt, n_ice, step, calc_time, m_min, Y_max)

  if (option == 1) then
    call ExportDomain (coh_coeff, frac, gamma_n, mu, n_elem, n_ice, density, youngs, &
      length, width, thickness)
  end if
end program Numerical_Ridge
