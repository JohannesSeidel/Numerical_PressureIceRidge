module mod_Neighbour

contains
  subroutine CalculateNeighbours
    use mod_Global
    implicit none
    call SweepAndPrune (xEndPoints)
    call SweepAndPrune (yEndPoints)
    call SweepAndPrune (zEndPoints)
  end subroutine CalculateNeighbours

  subroutine SweepAndPrune (endPoints)
    !basically an "insertionsort"-algorithm
    use mod_Global
    implicit none
    type (endPoint), intent (in out) :: endPoints(:)
    type (endPoint) :: tmp, swapper
    integer         :: i, j, k, iID, jID, dummy
    logical         :: flag

    do i = 2, size(endPoints)
      j = i - 1
      tmp = endPoints(i)

      do while (j >= 1 .and. endPoints(j)%mValue > tmp%mValue)

        !when two entries swap conduct a full 3D-overlap test
        swapper = endPoints(j)
        iID = abs(tmp%mData)
        jID = abs(swapper%mData)

        if (tmp%mData < 0.0 .and. swapper%mData > 0.0) then

          call AABBOverlap (boundingBoxes(iID), boundingBoxes(jID), flag)     

          if (flag .eqv. .true.) then
            !Even if they overlap, it is possible that this pair is already registered during the sorting of
            !another axis. Therefore, one has to test first if this pair is not already known. Each 
            !array index is checked for [i,j] and [j,i]. If every check for not being equal returns .true.
            !a new pair is added

            if (any ([(all ([iID, jID] == overlapPairs(k)%owner), k = 1, nPairs)]) .or. &
              & any ([(all ([jID, iID] == overlapPairs(k)%owner), k = 1, nPairs)])) then
              
              !is known
            
            else

              !is not known

              nPairs = nPairs + 1
              
              ! THE ERROR I WAS WORKING ON FOR TWO WEEKS !!!
              ! iID is the element index which is passed to the subroutine
              ! Compute_contact_area (in mod_OverlapComputation). If iID is the
              ! ID of a wall or structure, then this could result into an error
              ! in the computation of the contact normals.
              ! The walls and structures have big length-thickness-ratio,
              ! therfore elements(iID)%r's transversal displacement relative to
              ! the overlap centroid could be so big, that the program flips the
              ! normal vector of a face to the wrong direction
              if (iID > jID) then
                dummy = iID
                iID   = jID
                jID   = dummy
              end if

              overlapPairs(nPairs)%owner              = [iID, jID]
              overlapPairs(nPairs)%f_t_old            = 0.0d0
              overlapPairs(nPairs)%force_direction    = 0.0d0
              overlapPairs(nPairs)%overlap_volume     = 0.0d0
              overlapPairs(nPairs)%overlap_volume_old = 0.0d+00
              overlapPairs(nPairs)%overlap_centroid   = 0.0d0
              overlapPairs(nPairs)%overlap_area       = 0.0d0

            end if
          end if
        end if

        if (tmp%mData > 0 .and. swapper%mData < 0) then
          !remove i and j from pairList
          do k = 1, nPairs
            if (all(overlapPairs(k)%owner == [iID, jID]) .or. all(overlapPairs(k)%owner == [jID, iID])) then
              overlapPairs(k) = overlapPairs(nPairs)
              nPairs = nPairs - 1
              exit
            end if
          end do
        end if

        endPoints(j + 1) = endPoints(j)
        j = j - 1

        if (j == 0) exit
      end do

      endPoints(j + 1) = tmp
    end do
  end subroutine SweepAndPrune

  subroutine AABBOverlap (a, b, flag)
    use mod_Global
    implicit none
    type (aabb), intent (in)  :: a, b
    logical,     intent (out) :: flag
    if (a%mMax(1) < b%mMin(1) .or. b%mMax(1) < a%mMin(1)) then 
      flag = .false.
      return
    end if
    if (a%mMax(2) < b%mMin(2) .or. b%mMax(2) < a%mMin(2)) then 
      flag = .false.
      return
    end if
    if (a%mMax(3) < b%mMin(3) .or. b%mMax(3) < a%mMin(3)) then 
      flag = .false.
      return
    end if
    ! only if the three if conditions above are not fulfilled the aabbs overlap and flag ist set .true.
    flag = .true.
  end subroutine AABBOverlap

  subroutine InitialSort
    use mod_Global
    implicit none
    call InsertionSort (xEndPoints)
    call InsertionSort (yEndPoints)
    call InsertionSort (zEndPoints)
  end subroutine InitialSort

  pure subroutine InsertionSort (a)
    ! This algorithm is called during initialization. It just sorts the 
    ! endPoint-arrays without calling the aabb-overlap, since in the
    ! initial positions no overlaps are expected
    use mod_Global
    implicit none
    type (endPoint), intent (in out) :: a(:)
    type (endPoint) :: temp
    integer         :: i, j
    do i = 2, size (a)
      j = i - 1
      temp = a(i)
      do while (j >= 1 .and. a(j)%mValue > temp%mValue)
        a(j+1) = a(j)
        j = j - 1
        if (j == 0) exit
      end do
      a(j+1) = temp
    end do
  end subroutine InsertionSort

end module mod_Neighbour
