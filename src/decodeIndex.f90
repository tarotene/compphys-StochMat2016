PROGRAM main

  IMPLICIT NONE

  INTEGER(kind = 4) :: len_x, len_z
  INTEGER(kind = 4), ALLOCATABLE :: spin(:, :), x(:), z(:)

  WRITE(*, *) "len_x, len_z = ?"
  READ(*, *) len_x, len_z
  ALLOCATE(spin(0:len_x + 1, 0: len_z + 1))
  ALLOCATE(x(1:len_x * len_z), z(1:len_x * len_z))
  CALL foldSpinArray(x(1:len_x * len_z), z(1:len_x * len_z))
  CALL decodeIndex(1, spin)

  WRITE(*, *) spin(1:len_x, 1:len_z)

CONTAINS
  SUBROUTINE decodeIndex(i_st, spin)
    INTEGER(kind = 4), INTENT(in) :: i_st
    INTEGER(kind = 4), INTENT(out) :: spin(0:, 0:)

    INTEGER(kind = 4) :: str_st, i_v

    str_st = i_st
    DO i_v = 1, len_x * len_z, 1
       spin(x(i_v), z(i_v)) = IAND(str_st, 1) * 2 - 1
       str_st = ISHFT(str_st, -1)
    END DO
  END SUBROUTINE decodeIndex

  SUBROUTINE foldSpinArray(x, z)
    INTEGER(kind = 4), INTENT(out) :: x(1:len_x * len_z)
    INTEGER(kind = 4), INTENT(out) :: z(1:len_x * len_z)

    INTEGER(kind = 4) :: i_v

    x(1:len_x * len_z) = (/ (MOD(i_v - 1, len_x) + 1, i_v = 1, len_x * len_z) /)
    z(1:len_x * len_z) = (/ ((i_v - 1) / len_x + 1, i_v = 1, len_x * len_z) /)
  END SUBROUTINE foldSpinArray
END PROGRAM main
