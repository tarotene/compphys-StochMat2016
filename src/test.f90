MODULE global
  IMPLICIT NONE

  INTEGER(kind = 4), SAVE :: len_x, len_z, n_st
  INTEGER(kind = 4), SAVE, ALLOCATABLE :: x(:), z(:)
  DOUBLE PRECISION, SAVE :: J, beta
  DOUBLE PRECISION, SAVE :: deltaE(-1:1, -1:1, -1:1, -1:1, -1:1)
  DOUBLE PRECISION, SAVE :: prob(-1:1, -1:1, -1:1, -1:1, -1:1)
END MODULE global

PROGRAM main
  USE global
  IMPLICIT NONE

  DOUBLE PRECISION, ALLOCATABLE :: M(:, :)

  WRITE(0, '(a)') "len_x, len_z = ?"
  READ(*, *) len_x, len_z
  ALLOCATE(x(1:len_x * len_z), z(1:len_x * len_z))
  CALL foldSpinArray(x(1:len_x * len_z), z(1:len_x * len_z))
  n_st = ISHFT(1, len_x * len_z)

  ALLOCATE(M(1:n_st, 1:n_st))

  WRITE(0, '(a)') "beta, J = ?"
  READ(*, *) beta, J

  CALL makeDeltaEArray
  CALL makeProbArray

  CALL constructMatrixM(prob(-1:1, -1:1, -1:1, -1:1, -1:1), M(1:n_st, 1:n_st))
  CALL indicateMatrixM(M(1:n_st, 1:n_st))

CONTAINS
  SUBROUTINE constructMatrixM(prob, M)
    DOUBLE PRECISION, INTENT(in) :: prob(-1:1, -1:1, -1:1, -1:1, -1:1)
    DOUBLE PRECISION, INTENT(out) :: M(1:, 1:)

    INTEGER(kind = 4) :: i_st, j_st, i_v, east, west, south, north
    INTEGER(kind = 4) :: spin(0:len_x + 1, 0:len_z + 1)

    M(1:n_st, 1:n_st) = 0.0d0

    DO i_st = 1, n_st, 1
       CALL decodeIndex(i_st, spin(0:len_x + 1, 0:len_z + 1))
       DO i_v = 1, len_x * len_z, 1
          CALL setDirection(  spin(0:len_x + 1, 0:len_z + 1), &
               x(i_v), z(i_v), east, west, south, north)
          M(i_st, i_st) = M(i_st, i_st) + &
               (1.0d0 - prob(spin(x(i_v), z(i_v)), east, west, south, north)) &
               / DBLE(len_x * len_z)
       END DO
    END DO

    DO i_st = 1, n_st, 1
       CALL decodeIndex(i_st, spin(0:len_x + 1, 0:len_z + 1))
       DO i_v = 1, len_x * len_z, 1
          CALL setDirection(  spin(0:len_x + 1, 0:len_z + 1), &
               x(i_v), z(i_v), east, west, south, north)
          M(IBCHNG(i_st, i_v - 1), i_st) &
               = prob(spin(x(i_v), z(i_v)), east, west, south, north)
       END DO
    END DO
  END SUBROUTINE constructMatrixM

  SUBROUTINE indicateMatrixM(M)
    DOUBLE PRECISION, INTENT(in) :: M(1:, 1:)

    INTEGER(kind = 4) :: i_st, j_st

    DO i_st = 1, n_st, 1
       DO j_st = 1, n_st, 1
          WRITE(*, '(a, i3, a, i3, a, e9.2)') "M(", j_st, ", ", i_st, ") =", M(j_st, i_st)
       END DO
    END DO
  END SUBROUTINE indicateMatrixM

  SUBROUTINE makeDeltaEArray
    INTEGER(kind = 4) :: center, east, west, south, north

    DO center = -1, 1, 1
       DO east = -1, 1, 1
          DO west = -1, 1, 1
             DO south = -1, 1, 1
                DO north = -1, 1, 1
                   deltaE(center, east, west, south, north)=&
                        2 * J * center * (east + west + south + north)
                END DO
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE makeDeltaEArray

  SUBROUTINE makeProbArray
    INTEGER(kind = 4) :: center, east, west, south, north

    WHERE ( deltaE .GT. 0 )
       prob(-1:1, -1:1, -1:1, -1:1, -1:1) = &
            EXP(- beta * deltaE(-1:1, -1:1, -1:1, -1:1, -1:1))
    END WHERE
  END SUBROUTINE makeProbArray

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

  SUBROUTINE encodeSpin(spin, j_st)
    INTEGER(kind = 4), INTENT(inout) :: spin(0:, 0:)
    INTEGER(kind = 4), INTENT(out) :: j_st

    INTEGER(kind = 4) :: i_v

    j_st = 1
    DO i_v = 1, len_x * len_z, 1
       j_st = j_st + ISHFT((spin(x(i_v), z(i_v)) + 1) / 2, i_v - 1)
    END DO
  END SUBROUTINE encodeSpin

  SUBROUTINE foldSpinArray(x, z)
    INTEGER(kind = 4), INTENT(out) :: x(1:len_x * len_z)
    INTEGER(kind = 4), INTENT(out) :: z(1:len_x * len_z)

    INTEGER(kind = 4) :: i_v

    x(1:len_x * len_z) = (/ (MOD(i_v - 1, len_x) + 1, i_v = 1, len_x * len_z) /)
    z(1:len_x * len_z) = (/ ((i_v - 1) / len_x + 1, i_v = 1, len_x * len_z) /)
  END SUBROUTINE foldSpinArray

  SUBROUTINE setDirection(spin, x, z, east, west, south, north)
    INTEGER(kind = 4), INTENT(in) :: x, z
    INTEGER(kind = 4), INTENT(in) :: spin(0:, 0:)
    INTEGER(kind = 4), INTENT(out) :: east, west, south, north

    east = INT(spin(x + 1, z))
    west = INT(spin(x - 1, z))

    IF ( x == len_x ) THEN
       east = spin(1, z)
    END IF

    IF ( x == 1 ) THEN
       west = spin(len_x, z)
    END IF

    north = spin(x, z + 1)
    south = spin(x, z - 1)

    IF (z == len_z) THEN
       north = 0
    END IF

    IF (z == 0) THEN
       south = 0
    END IF
  END SUBROUTINE setDirection
END PROGRAM main
