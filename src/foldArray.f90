PROGRAM main
  IMPLICIT NONE

  INTEGER(kind = 4) :: len_x, len_z, i_v, x, z

  WRITE(0, '(a)') "len_x, len_z = ?"
  READ(*, *) len_x, len_z

  DO
     WRITE(0, '(a)') "i_v = ?"
     READ(*, *) i_v
     CALL foldArray(len_x, len_z, i_v, x, z)
     WRITE(0, '(a, i3, a, i3)') "x = ", x, ", z = ", z
  END DO

CONTAINS
  SUBROUTINE foldArray(len_x, len_z, i_v, x, z)
    IMPLICIT NONE

    INTEGER(kind = 4), INTENT(in) :: len_x, len_z, i_v
    INTEGER(kind = 4), INTENT(out) :: x, z

    x = MOD(i_v, len_x)
    z = INT(i_v / len_x) + 1
    IF ( x == 0 ) THEN
       x = len_x
       z = i_v / len_x
    END IF
  END SUBROUTINE foldArray
END PROGRAM main
