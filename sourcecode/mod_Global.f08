module mod_Global
  use type_Iceblock
  use type_Ridge
  use type_Element
  use type_Wall
  use type_AABB
  use type_ContactData

  implicit none

  type (element),          allocatable :: elements(:) ! pointer array to ALL elements
  type (iceblock), target, allocatable :: ice(:) ! array for all rubble elements
  type (wall),     target, allocatable :: boundaries(:) ! array for all boundaries, e.g. consol layer, basin walls
  type (wall),     target, allocatable :: structure(:)
  type (ridge)                         :: rdge

  type (aabb), target,     allocatable :: boundingBoxes(:)
  type (endPoint),         allocatable :: xEndPoints(:), yEndpoints(:), zEndPoints(:)
  type (overlapPair),      allocatable :: overlapPairs(:)
  integer                              :: nPairs

  integer :: n_steps
  integer :: step
end module mod_Global
