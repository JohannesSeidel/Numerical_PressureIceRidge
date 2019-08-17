module mod_Output
  implicit none
contains
  subroutine Graphics (n_elem, n_ice, n_struct, n_wall)
    use mod_Global

    implicit none

    integer :: i, j, k
    integer, intent (in) :: n_ice, n_struct, n_wall, n_elem
    integer, save :: gcount
    character (80) :: filename = ''
    integer :: points, faces

    ! initialization of the graphical counter (step is a global variable)
    if (step == 0) gcount = 1

    !============!
    ! vtk-output !
    !============!

    !==== POINTS ====!
    points = 0
    !points of all ice elements + consolidated layer
    do i = 1, n_ice + 1
      points = points + elements(i)%vertices
    end do
    !points of structure
    do i = n_ice + n_wall + 1, n_elem
      points = points + elements(i)%vertices
    end do

    write (filename,'(a,i5.5,a)') './output/vtk/output',gcount,'.vtk'
    gcount = gcount + 1
    open (unit = 21, file = filename, status = 'replace')

    write (21, '(a)') '# vtk DataFile Version 2.0'
    write (21, '(a)') 'Startpositionen'
    write (21, '(a)') 'ASCII'
    write (21, '(a)') 'DATASET POLYDATA'
    write (21, *)     'POINTS ', points, ' float'
    write (21, *)

    !ice elements + consolidated layer
    do i = 1, n_ice + 1
      do j = 1, elements(i)%vertices
        write (21, 1000) elements(i)%vert_coord_space(1,j), elements(i)%vert_coord_space(2,j), elements(i)%vert_coord_space(3,j)
      end do
    end do
    !structure
    do i = n_ice + n_wall + 1, n_elem
      do j = 1, elements(i)%vertices
        write (21, 1000) elements(i)%vert_coord_space(1,j), elements(i)%vert_coord_space(2,j), elements(i)%vert_coord_space(3,j)
      end do
    end do

    write (21, *)

    !==== FACES ====!
    faces = 0
    !faces of all ice elements + consolidated layer
    do i = 1, n_ice + 1
      faces = faces + elements(i)%faces
    end do
    !faces of structure
    do i = n_ice + n_wall + 1, n_elem
      faces = faces + elements(i)%faces
    end do

    write (21, *) 'POLYGONS', faces, faces * 5
    k = - 1
    do i = 1, n_ice + 1
      do j = 1, elements(i)%faces
        ! for vtk format, the first vertex index has to be zero
        write (21, 1001) 4, elements(i)%face_vertex_table(1, j) + k, elements(i)%face_vertex_table(2, j) + k, &
          elements(i)%face_vertex_table(3, j) + k, elements(i)%face_vertex_table(4, j) + k
      end do
      k = k + elements(i)%vertices
    end do

    do i = n_ice + n_wall + 1, n_elem
      do j = 1, elements(i)%faces
        ! for vtk format, the first vertex index has to be zero
        write (21, 1001) 4, elements(i)%face_vertex_table(1, j) + k, elements(i)%face_vertex_table(2, j) + k, &
          elements(i)%face_vertex_table(3, j) + k, elements(i)%face_vertex_table(4, j) + k
      end do
      k = k + elements(i)%vertices
    end do

    write (21, *)
    
    !====== POINT DATA ===========!
    write (21,'(a,i6)') "POINT_DATA", points
    write (21,'(a)') "VECTORS velocity float"

    do i = 1, n_ice + 1 ! +1, because the consolidated layer is elements(n_ice+1)
      do j = 1, elements(i)%vertices
        write (21,1002) elements(i)%rtd1
      end do
    end do
    do i = 1, n_struct
      do j = 1, structure(i)%vertices
        write (21,1002) structure(i)%rtd1
      end do
    end do

    close (21)

    1000 format (3(1x, f8.4))
    1002 format (3(1x, f10.6))
    1001 format (5(1x, i6))
  end subroutine Graphics

  subroutine Output (option, frac, dt, t_end, n_ice, n_steps, calc_time, m_min, Y_max)
    implicit none
    double precision, intent (in) :: frac, dt, t_end, calc_time, m_min, Y_max
    integer,          intent (in) :: n_ice, n_steps, option
    character (200)               :: path

    open (26, file = './output/results.txt', status = 'replace')
    if (option == 1) then
      write (26, *) "RIDGE CREATION"
    else if (option == 2) then
      write (26, *) "PUNCH TEST"
    else if (option == 3) then
      write (26, *) "MODEL TEST"
    end if
    write (26, *) 'Simulation successfully ended'
    write (26, '(1x,a,f10.3)')      'frac:             ', frac
    write (26, '(1x,a,d10.3,a)')    'm_min:            ', m_min, " kg"
    write (26, '(1x,a,d10.3,a)')    'Y_max:            ', Y_max, " N/m^2"
    write (26, '(1x,a,d10.3,a)')    'dt:               ', dt, " s"
    write (26, '(1x,a,f10.3,a)')    't_end:            ', t_end, " s"
    write (26, '(1x,a,i10)')        'n_ice:            ', n_ice
    write (26, '(1x,a,i10)')        'n_steps:          ', n_steps
    write (26, '(1x, a, f10.3, a)') 'calculation time: ', calc_time, " s"
    close (26)

    !compress all output files in ./output and move the compressed file to ./output/save
    call getcwd ( path )
    if (path(1:1) == '/') then !Operating system is Linux
      if (option==1) &
        call execute_command_line ('(cd output && tar -cf Ridge_Creation-$(date +%F_%H%M%S).tar.gz ridge csv vtk *.txt)')
      if (option==2) call execute_command_line ('(cd output && tar -cf Punch_Test-$(date +%F_%H%M%S).tar.gz csv vtk *.txt)')
      if (option==3) call execute_command_line ('(cd output && tar -cf Model_Test-$(date +%F_%H%M%S).tar.gz csv vtk *.txt)')
      call execute_command_line ('mv ./output/*.tar.gz ./output/save')
    else !OS is Windows / Mac
      write (*,*) "compress ridge csv vtk and *.txt and move folder to save"
      !call execute_command_line ('cd output')
      !call execute_command_line ('set backupFilename=%DATE:~6,4%%DATE:~3,2%%DATE:~0,2%')
      !call execute_command_line ('compact Save%backupFilename%.zip *.backup')
    end if
  end subroutine Output
end module mod_Output
