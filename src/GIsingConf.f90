PROGRAM main

  IMPLICIT NONE
  INTEGER(kind = 4) :: len_x, len_z, n_st
  INTEGER(kind = 4) :: i_st, j_st, x, z
  INTEGER(kind = 4), ALLOCATABLE :: st_spin(:, :, :), spin(:, :)

  WRITE(0, '(a)') "len_x, len_z = ?"
  READ(*, *) len_x, len_z
  n_st = 2 ** (len_x * len_z)
  ALLOCATE(spin(1:len_x, 1:len_z), st_spin(1:n_st, 1:len_x, 1:len_z))

  DO i_st = 1, n_st, 1
     CALL DecodeIndex(i_st, spin(1:len_x, 1:len_z))
     st_spin(i_st, 1:len_x, 1:len_z) = spin(1:len_x, 1:len_z)
  END DO

  DO i_st = 1, n_st, 1
     CALL EncodeSpin(st_spin(i_st, 1:len_x, 1:len_z), j_st)
     IF (i_st == j_st) THEN
        WRITE(0, '(a, i3, a, i3)') "Correct: ", i_st, " ＝ ", j_st
     ELSE
        WRITE(0, '(a, i3, a, i3)') "Correct: ", i_st, " ≠ ", j_st
     END IF
  END DO

CONTAINS
  SUBROUTINE foldArray(len_x, i_v, x, z)
    IMPLICIT NONE

    INTEGER(kind = 4), INTENT(in) :: len_x, i_v
    INTEGER(kind = 4), INTENT(out) :: x, z

    x = MOD(i_v, len_x)
    z = INT(i_v / len_x) + 1
    IF ( x == 0 ) THEN
       x = len_x
       z = i_v / len_x
    END IF
  END SUBROUTINE foldArray

  SUBROUTINE DecodeIndex(i_st, spin)
    IMPLICIT NONE

    INTEGER(kind = 4), INTENT(in) :: i_st
    INTEGER(kind = 4), INTENT(inout) :: spin(1:, 1:)

    INTEGER(kind = 4) :: len_x, len_z, i_v, m, n

    len_x = SIZE(spin, 1)
    len_z = SIZE(spin, 2)

    m = n - 2 * i_st + 2
    DO i_v = 1, len_x * len_z, 1
       n = (n - m) / 2
       m = MOD(n, 2)
       CALL foldArray(len_x, i_v, x, z)
       spin(x, z) = 2 * m - 1
    END DO
  END SUBROUTINE DecodeIndex

  SUBROUTINE EncodeSpin(spin, j_st)
    IMPLICIT NONE

    INTEGER(kind = 4), INTENT(in) :: spin(1:, 1:)
    INTEGER(kind = 4), INTENT(out) :: j_st

    INTEGER(kind = 4) :: len_x, len_z, i_v

    len_x = SIZE(spin, 1)
    len_z = SIZE(spin, 2)

    j_st = 1
    DO i_v = 1, len_x * len_z, 1
       CALL foldArray(len_x, i_v, x, z)
       j_st = j_st + (2 ** (i_v - 1)) * (spin(x, z) + 1) / 2
    END DO
  END SUBROUTINE EncodeSpin
END PROGRAM main
