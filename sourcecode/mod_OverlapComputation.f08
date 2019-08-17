module mod_OverlapComputation
  !
  ! this module contains all the procedures for the overlap computation:
  ! a) contact detection, e.g. "is there a contact?"
  ! b) inherited vertices
  ! c) generated vertices
  ! d) determination of the faces of the overlap polyhedron
  ! e) determination of the contact area and normal
  !
  ! The conducted steps are:
  ! 1) Find the inherited vertices and compute the generated vertices from the intersections
  !    of the (triangular/square) faces of the two polyhedra
  ! 2) Determine the generated faces and inherited faces based on the face intersection
  !    information
  ! 3) With the vertex coordinates and the vertex index list for the faces, compute the volume
  !    of the overlap polyhedron and its centroid
  ! 4) Join the intersection line segments to determine the contact line, and then determine
  !    the contact area by constructing the contact triangles from the centroid of the overlap
  !    polygon and each pair of successive vertices on the contact line
  !
  !    Record of revisions:
  !
  !      Date         Programmer          Description of change
  !     ======        ==========          =====================
  !     21.1.15       Johannes            Original code
  !     23.1.15                           a) is working more or less
  !     27.1.                             b) coded, but not tested yet
  !     5.2.                              c) coded
  !     12.2.                             added "compute overlap polyhedron faces"
  !     16.2.                             started compute overlap volume
  !     19.2.                             completed overlap volume and force direction
  !
  implicit none
  private
  public  :: OverlapComputation
  double precision, parameter :: EPS = 10.0d-14
  integer :: currentPair

contains
  subroutine OverlapComputation (elem1, elem2, k, force_direction, overlap_area, overlap_centroid, overlap_volume)
    use type_Iceblock
    use mod_Global

    implicit none

    integer,          intent (in)  :: elem1, elem2, k
    double precision, intent (out) :: force_direction(:)
    double precision, intent (out) :: overlap_centroid(:)
    double precision, intent (out) :: overlap_volume
    double precision, intent (out) :: overlap_area
    
    integer          :: intersect_point_pair_idx(2, 50)
    double precision :: vert_gen(3, 50)
    double precision :: overlap_face_equation(4, 20)
    double precision :: overlap_face_verts(10, 10, 3)
    integer          :: vert_idx
    double precision :: inherited_vert(2, 20, 3)
    integer :: inherited_count(2)
    integer :: contact_face_pair(2, 50)
    integer :: nverts(20)
    integer :: current_face
    integer :: num_int_pair

    currentPair = k

    call ComputeGeneratedVertices (elem1, elem2, contact_face_pair, intersect_point_pair_idx, &
      & num_int_pair, vert_idx, vert_gen)

    if (num_int_pair < 3) then !not enough generated contact pairs -> no overlap and no need for further computation
      overlap_volume = 0.0
      return
    end if

    call Compute_inherited_vertices (elem1, elem2, inherited_count, inherited_vert)

    call Compute_overlap_polyhedron_faces (contact_face_pair, elem1, elem2, inherited_count, intersect_point_pair_idx, &
      inherited_vert, num_int_pair, vert_gen, current_face, nverts, overlap_face_equation, overlap_face_verts)

    call Compute_overlap_volume (current_face, nverts, overlap_centroid, overlap_face_verts, overlap_volume)

    call Compute_contact_area (elem1, intersect_point_pair_idx, num_int_pair, &
      overlap_area, overlap_centroid, vert_gen, force_direction)

    if (overlap_area == 0.0) then
      overlap_volume = 0.0
      return
    end if

    overlapPairs(k)%force_direction  = force_direction
    overlapPairs(k)%overlap_area     = overlap_area
    overlapPairs(k)%overlap_centroid = overlap_centroid
    overlapPairs(k)%overlap_volume   = overlap_volume
  end subroutine OverlapComputation

  !============================!
  ! COMPUTE GENERATED VERTICES !
  !============================!
  subroutine ComputeGeneratedVertices (elem1, elem2, contact_face_pair, intersect_point_pair_idx, &
      & num_int_pair, vert_idx, vert_gen)
    !
    ! compute all the generated vertices by 'brute force', i.e. by computing
    ! the intersections of all the faces of one element with all the faces of the other
    ! element
    !
    ! INPUT:
    ! 'elem1' -> ID of element i
    ! 'elem2' -> ID of element j
    ! OUTPUT:
    ! 'vert_gen' -> vector containing all generated vertices
    ! 'intersect_point_pair_idx' -> array (dimension(2,n_pairs)) numbering the generated vertices

    use type_Iceblock
    use mod_Global

    implicit none
    integer,          intent (in)  :: elem1, elem2
    integer,          intent (out) :: contact_face_pair(:, :)
    integer,          intent (out) :: intersect_point_pair_idx(:, :)
    integer,          intent (out) :: num_int_pair
    integer,          intent (out) :: vert_idx
    double precision, intent (out) :: vert_gen(:, :)
    integer          :: iface_idx, jface_idx
    integer          :: inter_type
    double precision :: v_int1(3), v_int2(3), v_check(3)
    integer          :: i, j, k, dummy
    double precision :: intersect_point_pair(2, 50, 3)
    double precision :: distance, minimum
    logical          :: new


    contact_face_pair    = 0
    intersect_point_pair = 0.0
    num_int_pair         = 0

    do iface_idx = 1, elements(elem1)%faces
      do jface_idx = 1, elements(elem2)%faces
        call FindFaceIntersection (elem1, elem2, iface_idx, jface_idx, inter_type, v_int1, v_int2)
        if (inter_type == 1) then  
          !if two intersection points v_int1 and v_int2 exist
          num_int_pair = num_int_pair + 1
          !record [iface_idx,jface_idx] in a list of contact face pairs
          contact_face_pair(:, num_int_pair) = [iface_idx, jface_idx]
          !record the two points in a list of pairs of intersection points
          intersect_point_pair(1, num_int_pair, :) = v_int1
          intersect_point_pair(2, num_int_pair, :) = v_int2
        end if
      end do
    end do
    
    if (num_int_pair == 3 .and. &
      any ( [( norm2 (intersect_point_pair(1,i,:)-intersect_point_pair(2,i,:)) < EPS, i = 1, 3 )] ) ) then
      ! the distance between contact pair is so small, that the overlap volume
      ! would be very small as well
      write (25,*) "norm2 (v1-v2) < EPS"
      num_int_pair = 2
      return
    end if

    if (num_int_pair < 3) then
      return
    end if

    !index the intersection points as generated vertices
    !assign the first two intersection points as the first two generated vertices
    vert_gen = 0.0
    intersect_point_pair_idx = 0
    vert_gen(:, 1) = intersect_point_pair(1, 1, :)
    vert_gen(:, 2) = intersect_point_pair(2, 1, :)
    !assign the indices of the generated vertices for the intersection point pairs
    intersect_point_pair_idx(:, 1) = [1, 2]
    vert_idx = 2
    do i = 2, num_int_pair
      do j = 1, 2
        v_check = intersect_point_pair(j, i, :) 
        minimum = 10.0
        new     = .true.
        dummy   = vert_idx
        !the following loop looks checks intersect_point_pair for double entries
        !of v_check. If v_check is already known, k, which is the index number
        !of the already known vertex is saved into intersect_point_pair_idx. If
        !v_check isn't already known, it's added to vert_gen and its new index
        !is added to intersect_point_pair_idx 
        do k = 1, dummy
          distance = norm2 (v_check - vert_gen(:,k))
          if (distance < EPS .and. distance < minimum) then
            intersect_point_pair_idx(j,i) = k
            minimum = distance
            new     = .false.
          end if
          if (k == dummy .and. new) then
            vert_idx = vert_idx + 1
            vert_gen(:,vert_idx) = v_check
            intersect_point_pair_idx(j,i) = vert_idx
          end if
        end do
      end do
    end do
  end subroutine ComputeGeneratedVertices

  subroutine FindFaceIntersection (elem1, elem2, iface_idx, jface_idx, inter_type, v_int1, v_int2)
    use mod_Functions
    !
    ! this subroutine calculates the two intersection points (if they exist) which result from a
    ! face-face intersection.
    ! It returns as well the value 1 or 0 for the intersection type. As long there are less than
    ! two intersection points inter_type is 0.
    ! 
    ! INPUT:
    ! 1) elem1, elem2          -> Element IDs of element 1 and element 2
    ! 2) iface_idx, jface_idx  -> The face IDs of both elements
    !    Example: face 1 of element 1, and face 1 of element 2
    !
    ! OUTPUT:
    ! 1) inter_type     -> The type of intersection. Initially it is 0. When two intersection points
    !                      are found it is set to 1.
    ! 2) v_int, v_int2  -> The two intersection points, which form a line segment for the overlap volume
    !
    ! possible cases:
    ! 1) No intersection point after first call   -> possible face-vertex intersection or no intersection at all
    ! 2) One intersection point after first call  -> possible vertex-face intersection or edge-edge intersection
    ! 3) Two intersection points after first call -> possible overlap
    ! 3.1) ... + zero intersection points after second call
    ! 3.2) ... + one intersection point after second call
    ! 3.3) ... + two intersection points after second call -> the only case with an actual overlap
    !
    ! After four intersection points are found they are ordered in ascending order. Then the vertices 2 and 3
    ! form the intersection segment

    implicit none
    integer,          intent (in)  :: elem1, elem2, iface_idx, jface_idx
    integer,          intent (out) :: inter_type
    double precision, intent (out) :: v_int1(:), v_int2(:)
    integer :: num_inter_pts
    integer :: axis
    double precision :: v_i1(3), v_i2(3), v_i3(3), v_i4(3)
    double precision :: intersect_points(3, 4)

    v_i1 = 0.0
    v_i2 = 0.0
    v_i3 = 0.0
    v_i4 = 0.0
    v_int1 = 0.0
    v_int2 = 0.0
    intersect_points = 0.0
    inter_type = 0

    call FacePlaneIntersection (elem1, elem2, iface_idx, jface_idx, num_inter_pts, v_i1, v_i2)
    if (num_inter_pts == 2) then
      ! if two intersection points v_i1 and v_i2 exist there is a chance that the face of element 2
      ! is in contact with plane of element 1
      call FacePlaneIntersection (elem2, elem1, jface_idx, iface_idx, num_inter_pts, v_i3, v_i4)
      if (num_inter_pts == 2) then
        ! if two intersection points v_i3 and v_i4 exist
        ! - Choose a coordinate axis to align v_i1-v_i4
        ! - Arrange (v_i1, v_i2) and (v_i3, v_i4) in ascending direction along that axis
        axis = 1 ! x-axis
        if (abs (v_i1(axis) - v_i2(axis)) < EPS) axis = 2 !y-axis
        if (abs (v_i1(axis) - v_i2(axis)) < EPS) axis = 3 !z-axis
        if (v_i1(axis) > v_i2(axis)) call Swap_vertices(v_i1, v_i2)
        if (v_i3(axis) > v_i4(axis)) call Swap_vertices(v_i3, v_i4)
        if (maxval ([v_i1(axis), v_i3(axis)]) < minval ([v_i2(axis), v_i4(axis)])) then
          !the line segments (v_i1, v_i2) and (v_i3, v_i4) overlap
          !return intersection segment between max(V_i1,V_i3) and min(V_i2,V_i4)
          intersect_points = reshape ([v_i1, v_i2, v_i3, v_i4], shape (intersect_points))
          call Bubble_sort (intersect_points, axis)
          v_int1 = intersect_points(:, 2)
          v_int2 = intersect_points(:, 3)
          inter_type = 1
          !check for edge-edge intersection -> an coordinate difference of the order of 10e-14 can be 
          !regard as "acceptable" tolerance for the same point
          if (norm2 (v_int1 - v_int2) < EPS) inter_type = 0
        end if
      end if
    end if
  end subroutine FindFaceIntersection

  subroutine FacePlaneIntersection (elem1, elem2, iface_idx, jface_idx, num_inter_pts, v_ix, v_iy)
    !
    ! This subroutine calculates the intersection points for a given face
    ! of one element with the plane which is given by a face of an other
    ! element.
    !
    ! Input: Element 1, Element 2, Face_ID and Plane_ID
    ! Output: Intersection points V_ix and V_iy
    !
    use type_Iceblock
    use mod_Global
    implicit none
    integer,          intent (in)  :: elem1, elem2, iface_idx, jface_idx
    double precision, intent (out) :: v_ix(:), v_iy(:)
    integer,          intent (out) :: num_inter_pts
    integer          :: i, max_vertices, vert_idx
    double precision :: d
    double precision :: lambda
    double precision :: n(3)
    double precision :: v_1(3), v_2(3)

    num_inter_pts = 0
    v_ix   = 0.0
    v_iy   = 0.0
    max_vertices = size (elements(elem1)%face_vertex_table(:, iface_idx), 1) !max vertices of face_idx of element1
    vert_idx = 1

    d = elements(elem2)%face_equation(4, jface_idx)
    n = elements(elem2)%face_equation(1:3, jface_idx)

    do i = 1, max_vertices
      v_1 = elements(elem1)%vert_coord_space(:, elements(elem1)%face_vertex_table(vert_idx, iface_idx))
      if (vert_idx == max_vertices) vert_idx = 0 !=> vertex 5 is vertex 1
      v_2 = elements(elem1)%vert_coord_space(:, elements(elem1)%face_vertex_table(vert_idx + 1, iface_idx))
      vert_idx = vert_idx + 1

      lambda = (d - dot_product (n, v_1)) / (dot_product (n, (v_2 - v_1))) 

      if (0.0 <= lambda .and. lambda <= 1.0) then
        if (num_inter_pts == 0) then
          num_inter_pts = 1
          v_ix = v_1 + lambda * (v_2 - v_1)
          ! intersection point 1 is found
        else
          num_inter_pts = 2
          v_iy = v_1 + lambda * (v_2 - v_1)
          ! the maximum of possible intersection points is two. After the second point is found
          ! it is not necessary to jump back to the start of the inter_points loop. The "return"
          ! statement returns the conrol back to the procedure that invoced this procedure
          return
        end if
      !else if (abs (lambda) == 0.0 .or. lambda == 1.0) then ! degenerate cases: 
        !write (25,*) "degenerate case: lambda = 1 or 0"
        ! 1) one of two vertices lays on plane
        ! -> go through all edges and check for norm2 (v_ix-v_iy) < EPS and merge them
        ! 2) point on plane
        ! -> go through all edges and check for norm2 (v_ix-v_iy) < EPS and merge them
        ! 3) edge parallel to plane
        ! -> if any lambda = infinity (and .or. NaN?!) then check if there is an inherited vertex
        ! -> if no inherited vertex, then this face is outside and doesn't contribute to OP
      end if
    end do
  end subroutine FacePlaneIntersection

  !============================!
  ! COMPUTE INHERITED VERTICES !
  !============================!
  subroutine Compute_inherited_vertices (elem1, elem2, inherited_count, inherited_vert)
    implicit none
    integer,          intent (in)  :: elem1, elem2
    integer,          intent (out) :: inherited_count(:) !index 1 for element 1, index 2 for element 2
    double precision, intent (out) :: inherited_vert(:, :, :) !row 1 for element 1, row 2 for element 2
    inherited_count = 0
    inherited_vert  = 0.0d0
    call FindInheritedVertices (elem1, elem2, 1, inherited_count, inherited_vert)
    call FindInheritedVertices (elem2, elem1, 2, inherited_count, inherited_vert)
  end subroutine Compute_inherited_vertices

  subroutine FindInheritedVertices (elem1, elem2, row, inherited_count, inherited_vert)
    !
    ! This subroutine computes the inherited vertices
    ! of elem1 in elem2
    !
    use type_Iceblock
    use mod_Global

    implicit none
    integer,          intent (in)     :: elem1, elem2, row
    integer,          intent (in out) :: inherited_count(:)
    double precision, intent (in out) :: inherited_vert(:, :, :)
    integer          :: i, j
    double precision :: d, distance
    double precision :: n(3), v(3)

    outer_loop: do i = 1, elements(elem1)%vertices
      v = elements(elem1)%vert_coord_space(:, i)
      inner_loop: do j = 1, elements(elem2)%faces
        n = elements(elem2)%face_equation(1:3, j)
        d = elements(elem2)%face_equation(4, j)
        distance = dot_product (n, v) - d
        if (distance > -EPS) cycle outer_loop
        ! if this loop doesn't exit for every face the vertex i is inside element j
        ! => therefore, it is an inherited vertex
        if (j == elements(elem2)%faces) then
          ! Record v in the list fo vertices inherited from P_1
          inherited_count(row) = inherited_count(row) + 1
          inherited_vert(row, inherited_count(row), :) = v
        end if
      end do inner_loop
    end do outer_loop
  end subroutine FindInheritedVertices

  !=====================================!
  ! COMPUTE FACES OF OVERLAP POLYHEDRON !
  !=====================================!
  subroutine Compute_overlap_polyhedron_faces (contact_face_pair, elem1, elem2, inherited_count, intersect_point_pair_idx, &
      inherited_vert, num_int_pair, vert_gen, current_face, nverts, overlap_face_equation, overlap_face_verts)
    !   this subroutine calculates the faces of the overlap polyhedron
    !   INPUT: elem1, elem2, contact_face_pair, intersect_point_pair_idx
    !   OUTPUT: overlap_face_equation, overlap_face_verts
    !   
    !   overlap_face_equation: in the first row the element index is saved, in the second row the overlaping face index
    !   overlap_face_verts : each column stores the vertex indices for one face. The column index corresponds to the index of
    !                         overlap_face_equation.

    implicit none

    integer,          intent (in)  :: contact_face_pair(:, :)
    integer,          intent (in)  :: elem1, elem2
    integer,          intent (in)  :: inherited_count(:)
    integer,          intent (in)  :: intersect_point_pair_idx(:, :)
    double precision, intent (in)  :: inherited_vert(:, :, :)
    integer,          intent (in)  :: num_int_pair
    double precision, intent (in)  :: vert_gen(:, :)
    integer,          intent (out) :: current_face
    integer,          intent (out) :: nverts(:)
    double precision, intent (out) :: overlap_face_equation(:, :)
    double precision, intent (out) :: overlap_face_verts(:, :, :)
    integer :: found_face_idx(20)

    nverts = 0
    current_face = 0
    overlap_face_equation = 0
    overlap_face_verts = 0

    call FindFaces (nverts, current_face, found_face_idx, overlap_face_equation, overlap_face_verts, elem1, &
      num_int_pair, 1, contact_face_pair, intersect_point_pair_idx, vert_gen, inherited_count, inherited_vert)

    call FindFaces (nverts, current_face, found_face_idx, overlap_face_equation, overlap_face_verts, elem2, &
      num_int_pair, 2, contact_face_pair, intersect_point_pair_idx, vert_gen, inherited_count, inherited_vert)
  end subroutine Compute_overlap_polyhedron_faces

  subroutine FindFaces (nverts, current_face, found_face_idx, overlap_face_equation, overlap_face_verts, elem_id, &
      num_int_pair, row, contact_face_pair, intersect_point_pair_idx, vert_gen, inherited_count, inherited_vert)
    !
    ! INPUT:
    ! 1) Element id "element" -> for reference to face_equation
    ! 2) Number of intersection pairs "num_int_pair" -> used a the upper limit in the do-loop cycling trough contact_face_pair
    ! 3) "row" -> can be 1 or 2. The first row in contact_face_pair is occupied by the face indices of element 1. The second row ...
    ! 4) Array with the face pair indices "contact_face_pair" -> stores the face indices which are in contact in the same column
    ! 5) Array with the point indices (coordinates can be looked up in vert_gen) -> two points which are connected as a segment are
    !    stored in the same column
    ! 6) List with generated vertices "vert_gen" -> the column number corresponds to the indices in "intersect_point_pair_idx"
    !
    ! IN-OUTPUT:
    ! remark: this subroutine is called two times. Each time the following three variables are updated
    ! 1) List with the face equations of each face of the overlap polyhedron "overlap_face_equation"
    ! 2) List in which each column the coordinates of all vertices belonging to one face are stored. The column index corresponds to
    !    the column of "overlap_face_equation"
    ! 3) The number of calculated faces "current_face"
    ! 4) An vector containing for each face the number of vertices "nverts"
    !

    use type_Iceblock
    use mod_Global

    implicit none
    integer,          intent (in out) :: nverts(:)
    integer,          intent (in out) :: current_face
    integer,          intent (in out) :: found_face_idx(:)
    double precision, intent (in out) :: overlap_face_equation(:, :)
    double precision, intent (in out) :: overlap_face_verts(:, :, :)
    integer,          intent (in)     :: elem_id, num_int_pair, row 
    integer,          intent (in)     :: inherited_count(:)
    double precision, intent (in)     :: inherited_vert(:, :, :)
    integer,          intent (in)     :: contact_face_pair(:, :)
    integer,          intent (in)     :: intersect_point_pair_idx(:, :)
    double precision, intent (in)     :: vert_gen(:, :)
    integer :: found_vert_idx(10, 20)
    integer :: face_num, face_idx, vert_idx, max_face_idx
    integer :: i, j, k

    found_face_idx = 0
    found_vert_idx = 0
    max_face_idx   = maxval (contact_face_pair(row, :)) !number of faces of element

    !=================!
    ! GENERATED FACES !
    !=================!
    do j = 1, num_int_pair
      do face_idx = 1, max_face_idx
        !has face F_i an intersection with polygon P_2?
        if (contact_face_pair(row, j) == face_idx) then
          !is face F_i already known? If not...
          if (all (face_idx /= found_face_idx) ) then
            current_face                           = current_face + 1
            found_face_idx(current_face)           = face_idx
            overlap_face_equation(:, current_face) = elements(elem_id)%face_equation(:, face_idx)
            face_num = current_face
          else !F_i is already known
            !the face index then is the first column where it appears
            known_face: do k = 1, current_face
              if (face_idx == found_face_idx(k)) then
                face_num = k
                exit known_face
              end if
            end do known_face
          end if
          ! do loop over each row of "intersect_point_pair_idx"
          do i = 1, 2
            if (all (intersect_point_pair_idx(i, j) /= found_vert_idx(:, face_num)) ) then
              nverts(face_num)                                  = nverts(face_num) + 1
              found_vert_idx(nverts(face_num), face_num)        = intersect_point_pair_idx(i, j)
              overlap_face_verts(nverts(face_num), face_num, :) = vert_gen(:, intersect_point_pair_idx(i, j))
            else
              cycle
            end if
          end do
        end if
      end do
    end do

    !=================!
    ! INHERITED FACES !
    !=================!
    do i = 1, inherited_count(row)
      do vert_idx = 1, size (elements(elem_id)%vert_coord_space, 2)
        if ( norm2 (inherited_vert(row, i, :) - elements(elem_id)%vert_coord_space(:, vert_idx)) < EPS ) then
          do j = 1, size (elements(elem_id)%vertex_face_table, 1)
            face_idx = elements(elem_id)%vertex_face_table(j, vert_idx)
            !is F_i in the list of faces of P_0
            if (any (found_face_idx == face_idx)) then
              do face_num = 1, current_face
                if (found_face_idx(face_num) == face_idx) then
                  nverts(face_num) = nverts(face_num) + 1
                  overlap_face_verts(nverts(face_num), face_num, :) = inherited_vert(row, i, :)
                end if
              end do
            else
              !register F_i as new entry in the face list of P_0
              current_face = current_face + 1
              found_face_idx(current_face) = face_idx
              overlap_face_equation(:, current_face) = elements(elem_id)%face_equation(:, face_idx)
              !register V_k as a vertex of the face F_i
              nverts(current_face) = nverts(current_face) + 1
              overlap_face_verts(nverts(current_face), current_face, :) = inherited_vert(row, i, :)
            end if
          end do
        end if
      end do
    end do
  end subroutine

  !============================!
  ! COMPUTE THE OVERLAP VOLUME !
  !============================!
  subroutine Compute_overlap_volume (current_face, nverts, overlap_centroid, overlap_face_verts, overlap_volume)
    !
    ! INPUT: 
    ! overlap_face_verts, intersect_point_pair_idx
    !
    ! OUTPUT :: overlap_volume
    ! 
    use mod_Functions
    use mod_Global
    use type_Iceblock
    implicit none
    double precision, intent (in out) :: overlap_face_verts(:, :, :)
    integer,          intent (in)     :: current_face
    integer,          intent (in)     :: nverts(:)
    double precision, intent (out)    :: overlap_centroid(:)
    double precision, intent (out)    :: overlap_volume
    integer          :: i
    integer          :: face, vert
    double precision :: volume
    double precision :: origin(3)
    double precision :: v_a(3), v_b(3), v_c(3), v_i(3), v_1(3), v_2(3)
    double precision, allocatable :: cos_theta(:, :)
    integer,          allocatable :: sorted_idx(:)

    if (any ( [(nverts(face), face = 1, current_face)] < 3 ) ) then
      write (25,*) "face has less than three vertices"
      overlap_centroid = overlapPairs(currentPair)%overlap_centroid
      overlap_volume   = overlapPairs(currentPair)%overlap_volume
      return
    end if

    overlap_volume   = 0.0
    overlap_centroid = 0.0
    origin           = 0.0

    ! compute origin -> sum of all vertices devided by number of vertices
    ! the origin is simultaneously the centroid of the overlap polyhedron
    do face = 1, current_face
      origin = origin + sum (overlap_face_verts(1:nverts(face), face, :), 1) / nverts(face)
    end do  
    origin = origin / current_face
    overlap_centroid = origin

    ! sorting of the vertices and calculating the overlap volume
    do face = 1, current_face
      v_1 = overlap_face_verts(1, face, :)
      v_2 = overlap_face_verts(2, face, :)
      v_b = (v_2 - v_1) / norm2 (v_2 - v_1)
      !
      ! cos(theta_i) has just to be of size nverts-2 because the first two vertices are used to calculate v_1 and v_2.
      ! the first row stores the vertex index and the second row the cos(theta_i) value. Both information are necessary
      ! because in the next step one has to sort this array's cos(theta_i) values in ascending order. The vertex index
      ! has to be a clear assingement. The subroutine, coded above, is able to sort this array.
      !
      allocate (cos_theta(2, nverts(face)-2)); cos_theta = 0.0
      allocate (sorted_idx(nverts(face))); sorted_idx = 0
      do vert = 3, nverts(face)
        v_i = overlap_face_verts(vert, face, :)
        cos_theta(1, vert-2) = vert
        cos_theta(2, vert-2) = dot_product ((v_i - v_1), v_b) / norm2 (v_i - v_1)
      end do
      ! sort cos_theta according to its 2nd row in ascending order
      call Bubble_sort (cos_theta, 2)

      !save sorted indices (1st row of cos_theta) to sorted_idx
      sorted_idx(1:2) = [1, 2]
      do i = 3, nverts(face) 
        sorted_idx(i) = nint (cos_theta(1,i-2))
      end do

      v_a = overlap_face_verts(1, face, :)
      ! a face containing n vertices is split up in n-2 triangulars
      do i = 2, nverts(face)-1
        v_b = overlap_face_verts(sorted_idx(i), face, :)
        v_c = overlap_face_verts(sorted_idx(i+1), face, :)

        volume   = 1.0/6.0 * abs (dot_product ( Cross (v_a - origin, v_b - origin), v_c - origin ))
        overlap_volume   = overlap_volume + volume
      end do
      deallocate (cos_theta)
      deallocate (sorted_idx)
    end do
  end subroutine Compute_overlap_volume

  !======================!
  ! COMPUTE CONTACT AREA !
  !======================!
  subroutine Compute_contact_area (elem_id, intersect_point_pair_idx, num_int_pair, &
      overlap_area, overlap_centroid, vert_gen, force_direction)
    use mod_Functions
    use mod_Global
    use type_Iceblock
    implicit none
    integer,          intent (in)     :: elem_id
    integer,          intent (in)     :: num_int_pair
    double precision, intent (in)     :: overlap_centroid(:)
    double precision, intent (in)     :: vert_gen(:, :)
    double precision, intent (out)    :: force_direction(:)
    double precision, intent (out)    :: overlap_area
    integer,          intent (in out) :: intersect_point_pair_idx(:, :)
    integer,          allocatable     :: contact_line(:)
    double precision, allocatable     :: normals(:,:)
    integer          :: i, j
    double precision :: v_1(3), v_2(3)
    integer          :: v_end

    double precision, allocatable :: dummy(:,:)
    allocate (dummy(size(intersect_point_pair_idx,1), size(intersect_point_pair_idx,2)))
    dummy = intersect_point_pair_idx

    allocate (contact_line(num_int_pair + 1))
    contact_line = 0

    ! copy the first point pair into the contact line array
    contact_line(1:2) = intersect_point_pair_idx(1:2, 1)
    intersect_point_pair_idx(:, 1) = 0 !clear the entry
    v_end = contact_line(2)
    line_idx: do i = 2, num_int_pair
      intersect_idx: do j = 2, num_int_pair
        ! find V(end) as entry j in the list of pairs intersect_point_pair_idx
        if (v_end == intersect_point_pair_idx(1, j)) then
          ! assing the other vertex of this pair to V(end)
          v_end = intersect_point_pair_idx(2, j)
          contact_line(i + 1) = v_end
          intersect_point_pair_idx(:, j) = 0
          exit intersect_idx
        else if (v_end == intersect_point_pair_idx(2, j)) then
          v_end = intersect_point_pair_idx(1, j)
          contact_line(i + 1) = v_end
          intersect_point_pair_idx(:, j) = 0
          exit intersect_idx
        else
          cycle intersect_idx
        end if
      end do intersect_idx
    end do line_idx

    if (contact_line(1) /= v_end) then
      !reasons may be:
      ! 1) intersect_point_pair_idx and num_int_pair doesn't correlate
      ! 2) each index didn't appear twice
      force_direction = overlapPairs(currentPair)%force_direction
      overlap_area    = overlapPairs(currentPair)%overlap_area
      return
    end if
    
    if (count (contact_line == 0) /= 0) then
      !reasons may be:
      ! 1) intersect_point_pair_idx and num_int_pair doesn't correlate
      ! 2) each index didn't appear twice
      force_direction = overlapPairs(currentPair)%force_direction
      overlap_area    = overlapPairs(currentPair)%overlap_area
      return
    end if

    ! compute triangle normals and resuling from these the force direction
    allocate (normals(3, num_int_pair))

    do i = 1, num_int_pair
      v_1 = vert_gen(:, contact_line(i))
      v_2 = vert_gen(:, contact_line(i + 1))
      normals(:, i) = 0.5 * Cross ((v_1 - overlap_centroid), (v_2 - overlap_centroid))

      if (dot_product (normals(:, i), (elements(elem_id)%r - overlap_centroid)) < 0.0) then
        normals(:, i) = -1 * normals(:, i)
      end if
    end do

    force_direction = sum (normals, 2) / (norm2 (sum (normals, 2)))
    overlap_area = norm2 (sum (normals, 2))
  end subroutine Compute_contact_area
end module mod_OverlapComputation
