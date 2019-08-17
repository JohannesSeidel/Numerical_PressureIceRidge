module type_ContactData
  type :: overlapPair
    integer :: owner(2)
    double precision :: f_t_old(3)
    double precision :: overlap_volume, overlap_volume_old
    double precision :: overlap_area
    double precision :: overlap_centroid(3)
    double precision :: force_direction(3)
  end type
end module type_ContactData
