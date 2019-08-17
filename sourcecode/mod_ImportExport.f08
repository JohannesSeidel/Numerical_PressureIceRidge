module mod_ImportExport
contains
  subroutine ImportDomain (coh_coeff, frac, gamma_n, mu, n_ice, n_struct, n_elem)
    use mod_Global
    use mod_Particle
    implicit none
    integer, parameter :: n_wall = 3
    double precision :: length, width, thickness, density, youngs
    double precision :: coh_coeff, frac, gamma_n, mu
    integer :: i, j, n_elem, n_ice, n_struct
    character :: dummy

    open (20, file = './output/ridge/Test_Settings.txt')
    read (20,*) dummy, n_ice
    read (20,*) dummy, length
    read (20,*) dummy, width
    read (20,*) dummy, thickness
    read (20,*) dummy, density
    read (20,*) dummy, youngs
    read (20,*) dummy, rdge%width
    read (20,*) dummy, rdge%keel_width
    read (20,*) dummy, rdge%keel_height
    read (20,*) dummy, rdge%length
    read (20,*) dummy, rdge%porosity
    close (20)

    n_elem   = n_ice + n_wall + n_struct

    open (50, file = './output/ridge/simulation.txt')
    read (50,*) dummy, coh_coeff
    read (50,*) dummy, gamma_n
    read (50,*) dummy, mu
    read (50,*) dummy, nPairs
    read (50,*) dummy, frac
    close (50)

    call InitializeParticles (density, length, mu, n_ice, n_struct, n_wall, thickness, width, youngs)

    open (50, file = './output/ridge/rubble%r.txt')
    do i = 1, n_ice
      read (50, *) ice(i)%r
      read (50, *) ice(i)%rtd1
      read (50, *) ice(i)%rtd2
    end do
    close (50)

    open (50, file = './output/ridge/rubble%q.txt')
    do i = 1, n_ice
      read (50, *) ice(i)%q%w,    ice(i)%q%x,    ice(i)%q%y,    ice(i)%q%z
      read (50, *) ice(i)%qtd1%w, ice(i)%qtd1%x, ice(i)%qtd1%y, ice(i)%qtd1%z
      read (50, *) ice(i)%qtd2%w, ice(i)%qtd2%x, ice(i)%qtd2%y, ice(i)%qtd2%z
    end do
    close (50)

    open (50, file = './output/ridge/rubble%wltA.txt')
    do i = 1, n_ice
      read (50, *) ice(i)%length, ice(i)%width, ice(i)%thickness, ice(i)%A
      ice(i)%vert_coord_body = reshape (&
        [ -0.5*ice(i)%width, -0.5*ice(i)%length,  0.5*ice(i)%thickness, &
        & -0.5*ice(i)%width,  0.5*ice(i)%length,  0.5*ice(i)%thickness, &
        &  0.5*ice(i)%width,  0.5*ice(i)%length,  0.5*ice(i)%thickness, &
        &  0.5*ice(i)%width, -0.5*ice(i)%length,  0.5*ice(i)%thickness, &
        & -0.5*ice(i)%width, -0.5*ice(i)%length, -0.5*ice(i)%thickness, &
        & -0.5*ice(i)%width,  0.5*ice(i)%length, -0.5*ice(i)%thickness, &
        &  0.5*ice(i)%width,  0.5*ice(i)%length, -0.5*ice(i)%thickness, &
        &  0.5*ice(i)%width, -0.5*ice(i)%length, -0.5*ice(i)%thickness], &
        shape (ice(i)%vert_coord_body) )
    end do
    close (50)

    do i = 1, n_ice
      !Update vertex coordinates (the function 'Rotate_q' is declared in mod_Quaternion)
      do j = 1, ice(i)%vertices
        ice(i)%vert_coord_space(:, j) = ice(i)%r + Rotate_q (ice(i)%vert_coord_body(:, j), ice(i)%q)
      end do
    end do

    do i = 1, n_ice
      ice(i)%mass = ice(i)%length * ice(i)%width * ice(i)%thickness * ice(i)%density
    end do

    open (50, file = './output/ridge/overlapPairs.txt')
    do i = 1, nPairs
      read (50, *) overlapPairs(i)%owner
      read (50, *) overlapPairs(i)%f_t_old
      read (50, *) overlapPairs(i)%overlap_volume
      read (50, *) overlapPairs(i)%overlap_volume_old
      read (50, *) overlapPairs(i)%overlap_area
      read (50, *) overlapPairs(i)%overlap_centroid
      read (50, *) overlapPairs(i)%force_direction
    end do
    close (50)

    open (50, file = './output/ridge/endPoints.txt')
    do i = 1, 2*(n_ice + n_wall)
      read (50, *) xEndPoints(i)%mData, yEndPoints(i)%mData, zEndPoints(i)%mData

      !die Pointer in den EndPoints-Arrays werden wieder der richtigen BoundingBox zugeordert
      if (xEndPoints(i)%mData > 0) then
        xEndPoints(i)%mValue => boundingBoxes(xEndPoints(i)%mData)%mMax(1)
      else
        xEndPoints(i)%mValue => boundingBoxes(abs(xEndPoints(i)%mData))%mMin(1)
      end if

      if (yEndPoints(i)%mData > 0) then
        yEndPoints(i)%mValue => boundingBoxes(yEndPoints(i)%mData)%mMax(2)
      else
        yEndPoints(i)%mValue => boundingBoxes(abs(yEndPoints(i)%mData))%mMin(2)
      end if

      if (zEndPoints(i)%mData > 0) then
        zEndPoints(i)%mValue => boundingBoxes(zEndPoints(i)%mData)%mMax(3)
      else
        zEndPoints(i)%mValue => boundingBoxes(abs(zEndPoints(i)%mData))%mMin(3)
      end if
    end do
    close (50)
  end subroutine ImportDomain

  subroutine ExportDomain (coh_coeff, frac, gamma_n, mu, n_elem, n_ice, density, youngs, &
      length, width, thickness)
    use mod_Global
    implicit none
    double precision, intent (in) :: density, youngs, length, width, thickness
    double precision, intent (in) :: coh_coeff, gamma_n, mu, frac
    integer,          intent (in) :: n_elem, n_ice
    integer                       :: i, counter
    integer, parameter            :: n_wall = 3 

    open (50, file = './output/ridge/rubble%r.txt', status = 'replace')
    do i = 1, n_ice
      write (50, *) ice(i)%r
      write (50, *) ice(i)%rtd1
      write (50, *) ice(i)%rtd2
    end do
    close (50)

    open (50, file = './output/ridge/rubble%q.txt', status = 'replace')
    do i = 1, n_ice
      write (50, *) ice(i)%q%w,    ice(i)%q%x,    ice(i)%q%y,    ice(i)%q%z
      write (50, *) ice(i)%qtd1%w, ice(i)%qtd1%x, ice(i)%qtd1%y, ice(i)%qtd1%z
      write (50, *) ice(i)%qtd2%w, ice(i)%qtd2%x, ice(i)%qtd2%y, ice(i)%qtd2%z
    end do
    close (50)

    open (50, file = './output/ridge/rubble%wltA.txt', status = 'replace')
    do i = 1, n_ice
      write (50, *) ice(i)%length, ice(i)%width, ice(i)%thickness, ice(i)%A
    end do
    close (50)

    !remove the push bars from the pair manager (overlapPairs) and
    !from the EndPoints
    !look if overlapPairs(i) has owner structure. If true then
    !copy the the last overlapPair into the i-th place and decrease
    !i by one to check in the next loop if the last Pair has a
    !structure itself
    counter = 0
    i = 0
    do while (i < nPairs)
      i = i + 1
      if (any (overlapPairs(i)%owner > n_ice+n_wall)) then
        overlapPairs(i) = overlapPairs(nPairs-counter)
        overlapPairs(nPairs-counter)%owner = [0, 0]
        counter = counter + 1
        i = i - 1
      end if
    end do
    nPairs = nPairs - counter

    open (50, file = './output/ridge/overlapPairs.txt', status = 'replace')
    do i = 1, nPairs
      write (50, *) overlapPairs(i)%owner
      write (50, *) overlapPairs(i)%f_t_old
      write (50, *) overlapPairs(i)%overlap_volume
      write (50, *) overlapPairs(i)%overlap_volume_old
      write (50, *) overlapPairs(i)%overlap_area
      write (50, *) overlapPairs(i)%overlap_centroid
      write (50, *) overlapPairs(i)%force_direction
    end do
    close (50)

    counter = 0
    i = 0
    do while (i < n_elem*2)
      i = i + 1
      if (abs (xEndPoints(i)%mData) > n_ice+n_wall) then
        counter = counter + 1
        xEndPoints( i : (n_elem*2-counter) ) = xEndPoints( i+1 : (n_elem*2-counter+1))
        xEndPoints(n_elem*2-counter+1)%mData = 0
        i = i - 1
      end if
    end do

    counter = 0
    i = 0
    do while (i < n_elem*2)
      i = i + 1
      if (abs (yEndPoints(i)%mData) > n_ice+n_wall) then
        counter = counter + 1
        yEndPoints( i : (n_elem*2-counter) ) = yEndPoints( i+1 : (n_elem*2-counter+1))
        yEndPoints(n_elem*2-counter+1)%mData = 0
        i = i - 1
      end if
    end do

    counter = 0
    i = 0
    do while (i < n_elem*2)
      i = i + 1
      if (abs (zEndPoints(i)%mData) > n_ice+n_wall) then
        counter = counter + 1
        zEndPoints( i : (n_elem*2-counter) ) = zEndPoints( i+1 : (n_elem*2-counter+1))
        zEndPoints(n_elem*2-counter+1)%mData = 0
        i = i - 1
      end if
    end do

    open (50, file = './output/ridge/endPoints.txt', status = 'replace')
    do i = 1, 2*(n_ice+n_wall)
      write (50, *) xEndPoints(i)%mData, yEndPoints(i)%mData, zEndPoints(i)%mData
    end do
    close (50)

    open (50, file = './output/ridge/simulation.txt', status = 'replace')
    write (50,*) "coh_coeff    ", coh_coeff
    write (50,*) "gamma_n      ", gamma_n
    write (50,*) "mu           ", mu
    write (50,*) "nPairs       ", nPairs
    write (50,*) "frac         ", frac
    close (50)

    open (50, file = './output/ridge/Test_Settings.txt', status = 'replace')
    write (50,*) "n_ice ", n_ice
    write (50,*)
    write (50,*) "ice_length    ", length
    write (50,*) "ice_width     ", width
    write (50,*) "ice_thickness ", thickness
    write (50,*) "ice_density   ", density
    write (50,*) "ice_youngs    ", youngs
    write (50,*)
    write (50,*) "rdge%width         ", rdge%width
    write (50,*) "ridge%keel_width   ", rdge%keel_width
    write (50,*) "ridge%keel_height  ", rdge%keel_height
    write (50,*) "rdge%length        ", rdge%length
    write (50,*) "rdge%porosity      ", rdge%porosity
    close (50)
  end subroutine ExportDomain
end module mod_ImportExport
