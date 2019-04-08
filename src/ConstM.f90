MODULE simulateSpin
  IMPLICIT NONE

CONTAINS
  SUBROUTINE zboundAntiParallel(spin)
    INTEGER(kind = 4), INTENT(inout) :: spin(0:, 0:)

    INTEGER(kind = 4) :: len_x, len_z

    len_x = SIZE(spin, dim=1) - 2
    len_z = SIZE(spin, dim=2) - 2

    spin(1:len_x, 0) = -1
    spin(1:len_x, len_z + 1) = 1
  END SUBROUTINE zboundAntiParallel

  SUBROUTINE zboundParallel(spin)
    INTEGER(kind = 4), INTENT(inout) :: spin(0:, 0:)

    INTEGER(kind = 4) :: len_x, len_z

    len_x = SIZE(spin, dim=1) - 2
    len_z = SIZE(spin, dim=2) - 2

    spin(1:len_x, 0) = -1
    spin(1:len_x, len_z + 1) = 1
  END SUBROUTINE zboundParallel

  SUBROUTINE zboundFree(spin)
    INTEGER(kind = 4), INTENT(inout) :: spin(0:, 0:)

    INTEGER(kind = 4) :: len_x, len_z

    len_x = SIZE(spin, dim=1) - 2
    len_z = SIZE(spin, dim=2) - 2

    spin(1:len_x, 0) = 0d0
    spin(1:len_x, len_z + 1) = 0d0
  END SUBROUTINE zboundFree

  SUBROUTINE makeProbArray(beta, deltaE, prob)
    DOUBLE PRECISION, INTENT(in) :: beta, deltaE(-1:1, -1:1, -1:1, -1:1, -1:1)
    DOUBLE PRECISION, INTENT(out) :: prob(-1:1, -1:1, -1:1, -1:1, -1:1)

    INTEGER(kind = 4) :: center, east, west, south, north

    prob(-1:1, -1:1, -1:1, -1:1, -1:1) = 1.0d0
    DO center = -1, 1, 1
       DO east = -1, 1, 1
          DO west = -1, 1, 1
             DO south = -1, 1, 1
                DO north = -1, 1, 1
                   IF ( deltaE(center, east, west, south, north) > 0 ) THEN
                      prob(center, east, west, south, north)=&
                           EXP(- beta * deltaE(center, east, west, south, north))
                   END IF
                END DO
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE makeProbArray

  SUBROUTINE makeDeltaEArray(J, deltaE)
    DOUBLE PRECISION, INTENT(in) :: J
    DOUBLE PRECISION, INTENT(out) :: deltaE(-1:1, -1:1, -1:1, -1:1, -1:1)

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

  SUBROUTINE setDirection(spin, x, z, east, west, south, north)
    INTEGER(kind = 4), INTENT(in) :: x, z
    INTEGER(kind = 4), INTENT(in) :: spin(0:, 0:)
    INTEGER(kind = 4), INTENT(out) :: east, west, south, north

    INTEGER(kind = 4) :: len_x, len_z

    len_x = SIZE(spin, dim=1) - 2
    len_z = SIZE(spin, dim=2) - 2

    east = INT(spin(x + 1, z))
    west = INT(spin(x - 1, z))

    IF ( x == len_x ) THEN
       east = INT(spin(1, z))
    END IF

    IF ( x == 1 ) THEN
       west = INT(spin(len_x, z))
    END IF

    north = INT(spin(x, z + 1))
    south = INT(spin(x, z - 1))
  END SUBROUTINE setDirection
END MODULE simulateSpin

MODULE calcMatrix
  USE simulateSpin

  IMPLICIT NONE

CONTAINS
  SUBROUTINE foldArray(len_x, i_v, x, z)
    INTEGER(kind = 4), INTENT(in) :: len_x, i_v
    INTEGER(kind = 4), INTENT(out) :: x, z

    x = MOD(i_v, len_x)
    z = INT(i_v / len_x) + 1
    IF ( x == 0 ) THEN
       x = len_x
       z = i_v / len_x
    END IF
  END SUBROUTINE foldArray

  SUBROUTINE decodeIndex(i_st, spin)
    INTEGER(kind = 4), INTENT(in) :: i_st
    INTEGER(kind = 4), INTENT(inout) :: spin(1:, 1:)

    INTEGER(kind = 4) :: len_x, len_z, x, z, i_v, m, n

    len_x = SIZE(spin, dim=1)
    len_z = SIZE(spin, dim=2)

    m = n - 2 * i_st + 2
    DO i_v = 1, len_x * len_z, 1
       n = (n - m) / 2
       m = MOD(n, 2)
       CALL foldArray(len_x, i_v, x, z)
       spin(x, z) = 2 * m - 1
    END DO
  END SUBROUTINE decodeIndex

  SUBROUTINE encodeSpin(spin, j_st)
    INTEGER(kind = 4), INTENT(in) :: spin(1:, 1:)
    INTEGER(kind = 4), INTENT(out) :: j_st

    INTEGER(kind = 4) :: len_x, len_z, x, z, i_v

    len_x = SIZE(spin, dim=1)
    len_z = SIZE(spin, dim=2)

    j_st = 1
    DO i_v = 1, len_x * len_z, 1
       CALL foldArray(len_x, i_v, x, z)
       j_st = j_st + (2 ** (i_v - 1)) * (spin(x, z) + 1) / 2
    END DO
  END SUBROUTINE encodeSpin

SUBROUTINE constructMatrixM(len_x, len_z, prob, M)
  INTEGER(kind = 4), INTENT(in) :: len_x, len_z
  DOUBLE PRECISION, INTENT(in) :: prob(-1:1, -1:1, -1:1, -1:1, -1:1)
  DOUBLE PRECISION, INTENT(out) :: M(1:, 1:)

  INTEGER(kind = 4) :: i_st, j_st, i_v, x, z, east, west, south, north
  INTEGER(kind = 4) :: n_st, spin(0:len_x + 1, 0:len_z + 1)

  n_st = SIZE(M, dim=1)
  M(1:n_st, 1:n_st) = 0d0

  CALL zboundFree(spin(0:len_x + 1, 0:len_z + 1))

  DO i_st = 1, n_st, 1
     CALL decodeIndex(i_st, spin(1:len_x, 1:len_z))
     DO i_v = 1, len_x * len_z, 1
        CALL foldArray(len_x, i_v, x, z)
        CALL setDirection(spin(0:len_x + 1, 0:len_z + 1), x, z, east, west, south, north)
        M(i_st, i_st) = M(i_st, i_st) + 1 - prob(spin(x, z), east, west, south, north)
     END DO
  END DO

  DO i_st = 1, n_st, 1
     CALL decodeIndex(i_st, spin(1:len_x, 1:len_z))
     DO i_v = 1, len_x * len_z, 1
        CALL foldArray(len_x, i_v, x, z)
        CALL setDirection(spin(0:len_x + 1, 0:len_z + 1), x, z, east, west, south, north)
        spin(x, z) = - spin(x, z)
        CALL encodeSpin(spin(1:len_x, 1:len_z), j_st)
        spin(x, z) = - spin(x, z)
        M(j_st, i_st) = prob(spin(x, z), east, west, south, north)
     END DO
  END DO

  M(1:n_st, 1:n_st) = M(1:n_st, 1:n_st) / DBLE(len_x * len_z)
END SUBROUTINE constructMatrixM

END MODULE calcMatrix

PROGRAM main
  USE calcMatrix
  USE simulateSpin

  IMPLICIT NONE
  INTEGER(kind = 4) :: len_x, len_z, n_st
  INTEGER(kind = 4) :: i_st, j_st, i_v, x, z, east, west, south, north
  INTEGER(kind = 4), ALLOCATABLE :: st_spin(:, :, :), spin(:, :)
  DOUBLE PRECISION :: beta, J
  DOUBLE PRECISION, ALLOCATABLE :: deltaE(:, :, :, :, :), prob(:, :, :, :, :)
  DOUBLE PRECISION, ALLOCATABLE :: M(:, :)

  DOUBLE PRECISION, ALLOCATABLE :: EIGR(:), EIGI(:), VL(:, :), VR(:, :)
  INTEGER(kind = 4), ALLOCATABLE :: WORK(:)
  INTEGER(kind = 4) :: LWORK, INFO

  WRITE(0, '(a)') "len_x, len_z = ?"
  READ(*, *) len_x, len_z
	n_st = 2 ** (len_x * len_z)
	ALLOCATE(M(1:n_st, 1:n_st))
  ALLOCATE(spin(0:len_x + 1, 0:len_z + 1))
  ALLOCATE(deltaE(-1:1,-1:1,-1:1,-1:1,-1:1), prob(-1:1,-1:1,-1:1,-1:1,-1:1))

  ALLOCATE(EIGR(1:n_st), EIGI(1:n_st), VL(1:n_st, 1:n_st), VR(1:n_st, 1:n_st))
  LWORK = 4 * n_st
  ALLOCATE(WORK(1:LWORK))

  WRITE(0, '(a)') "beta, J = ?"
  READ(*, *) beta, J

  CALL makeDeltaEArray(J, deltaE(-1:1,-1:1,-1:1,-1:1,-1:1))
  CALL makeProbArray(beta, deltaE(-1:1,-1:1,-1:1,-1:1,-1:1), prob(-1:1,-1:1,-1:1,-1:1,-1:1))

	CALL constructMatrixM(len_x, len_z, prob(-1:1, -1:1, -1:1, -1:1, -1:1), M(1:n_st, 1:n_st))

  !
  ! DO i_st = 1, n_st, 1
  !     WRITE(*, '(a, i3, a, i3, a, e)')  "M(", j_st,", ", i_st, ") = ", M(j_st, i_st)
  !   END DO
  ! END DO
END PROGRAM main
