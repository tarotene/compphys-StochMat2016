MODULE global
  IMPLICIT NONE

  INTEGER(kind = 4), SAVE :: len_x, len_z, n_st
  INTEGER(kind = 4), SAVE, ALLOCATABLE :: x(:), z(:)
  DOUBLE PRECISION, SAVE :: J, beta
  DOUBLE PRECISION, SAVE :: deltaE(-1:1, -1:1, -1:1, -1:1, -1:1)
  DOUBLE PRECISION, SAVE :: prob(-1:1, -1:1, -1:1, -1:1, -1:1)
END MODULE global

MODULE simulateSpin
  USE global
  IMPLICIT NONE

CONTAINS
  SUBROUTINE makeProbArray
    prob(-1:1, -1:1, -1:1, -1:1, -1:1) = 1.0d0
    WHERE ( deltaE .GT. 0.0d0 )
       prob(-1:1, -1:1, -1:1, -1:1, -1:1) = &
            EXP(- beta * deltaE(-1:1, -1:1, -1:1, -1:1, -1:1))
    END WHERE
  END SUBROUTINE makeProbArray

  SUBROUTINE makeDeltaEArray
    INTEGER(kind = 4) :: center, east, west, south, north

    DO center = -1, 1, 1
       DO east = -1, 1, 1
          DO west = -1, 1, 1
             DO south = -1, 1, 1
                DO north = -1, 1, 1
                   deltaE(center, east, west, south, north)=&
                        2 * J * DBLE(center * (east + west + south + north))
                END DO
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE makeDeltaEArray

  SUBROUTINE foldSpinArray(x, z)
    INTEGER(kind = 4), INTENT(out) :: x(1:len_x * len_z)
    INTEGER(kind = 4), INTENT(out) :: z(1:len_x * len_z)

    INTEGER(kind = 4) :: i_v

    x(1:len_x * len_z) = (/ (MOD(i_v - 1, len_x) + 1, i_v = 1, len_x * len_z) /)
    z(1:len_x * len_z) = (/ ((i_v - 1) / len_x + 1, i_v = 1, len_x * len_z) /)
  END SUBROUTINE foldSpinArray

  SUBROUTINE setDirection(spin, x, z, east, west, south, north)
    INTEGER(kind = 4), INTENT(in) :: x, z
    INTEGER(kind = 4), INTENT(in) :: spin(1:, 1:)
    INTEGER(kind = 4), INTENT(out) :: east, west, south, north

    IF ( x == len_x ) THEN
       east = spin(1, z)
		 ELSE
			 east = spin(x + 1, z)
    END IF

    IF ( x == 1 ) THEN
       west = spin(len_x, z)
		 ELSE
			 west = spin(x - 1, z)
    END IF

    IF (z == len_z) THEN
       north = 0
		 ELSE
			 north = spin(x, z + 1)
    END IF

    IF (z == 1) THEN
       south = 0
		 ELSE
			 south = spin(x, z - 1)
    END IF
  END SUBROUTINE setDirection

  SUBROUTINE shiftUpperHalf(spin)
    INTEGER(kind = 4), INTENT(inout) :: spin(1:, 1:)
    INTEGER(kind = 4) :: half_len_z
    INTEGER(kind = 4), ALLOCATABLE :: temp_spin(:, :)

    half_len_z = len_z / 2

    ALLOCATE(temp_spin(1:len_x, 1:half_len_z))

    temp_spin(1:len_x, 1:half_len_z) = spin(1:len_x, 1:half_len_z)
    temp_spin(1:len_x, 1:half_len_z) = CSHIFT(temp_spin(1:len_x, 1:half_len_z), 1, 1)
    spin(1:len_x, 1:half_len_z) = temp_spin(1:len_x, 1:half_len_z)
  END SUBROUTINE shiftUpperHalf

  SUBROUTINE calcTotalEnergy(spin, en_tot)
    INTEGER(kind = 4), INTENT(in) :: spin(1:, 1:)
    DOUBLE PRECISION, INTENT(out) :: en_tot

    INTEGER(kind = 4) :: x, z, east, west, south, north

    en_tot = 0.0d0

    DO z = 1, len_z, 1
       DO x = 1, len_x, 1
          CALL setDirection(spin, x, z, east, west, south, north)
          ! WRITE(*, *) spin(x, z), east, west, south, north
          en_tot = en_tot - J * DBLE(spin(x, z) * (east + west + south + north))
       END DO
    END DO

    en_tot = en_tot / 2
  END SUBROUTINE calcTotalEnergy
END MODULE simulateSpin

MODULE calcMatrix
  USE global
  USE simulateSpin

  IMPLICIT NONE
CONTAINS

  SUBROUTINE decodeIndex(i_st, spin)
    INTEGER(kind = 4), INTENT(in) :: i_st
    INTEGER(kind = 4), INTENT(out) :: spin(1:, 1:)

    INTEGER(kind = 4) :: str_st, i_v

    str_st = i_st - 1
    DO i_v = 1, len_x * len_z, 1
       spin(x(i_v), z(i_v)) = IAND(str_st, 1) * 2 - 1
       str_st = ISHFT(str_st, -1)
    END DO
  END SUBROUTINE decodeIndex

  SUBROUTINE encodeSpin(spin, j_st)
    INTEGER(kind = 4), INTENT(inout) :: spin(1:, 1:)
    INTEGER(kind = 4), INTENT(out) :: j_st

    INTEGER(kind = 4) :: i_v

    j_st = 1
    DO i_v = 1, len_x * len_z, 1
      j_st = j_st + ISHFT((spin(x(i_v), z(i_v)) + 1) / 2, i_v - 1)
    END DO
  END SUBROUTINE encodeSpin

  SUBROUTINE constructMatrixM(M)
    DOUBLE PRECISION, INTENT(out) :: M(1:, 1:)

    INTEGER(kind = 4) :: i_st, j_st, i_v, east, west, south, north
    INTEGER(kind = 4) :: spin(1:len_x, 1:len_z)

    M(1:n_st, 1:n_st) = 0.0d0

    DO i_st = 1, n_st, 1
       CALL decodeIndex(i_st, spin(1:len_x, 1:len_z))
      !  WRITE(*, *) 1
       DO i_v = 1, len_x * len_z, 1
          ! WRITE(*, *) 2
          CALL setDirection(  spin(1:len_x, 1:len_z), &
                              x(i_v), z(i_v), east, west, south, north)
          ! WRITE(*, *) 3
          ! WRITE(*, '(a, i5, a, i5, a, e9.2)') "M(", i_st, ", ", i_st, ") = ", M(i_st, i_st)
          ! WRITE(*, '(a, i5, a, i5, a, i5, a, i5, a, i5, a)') "center = ", center, &
          ! ", direction = (", east, ", ", west, ", ", south, ", ", north, ")."
          ! WRITE(*, '(a, i5, a, i5, a, i5, a, i5, a, i5, a, e9.2)') &
          !  "prob(", center, &
          !  ", ", east, ", ", west, ", ", south, ", ", north, &
            ! ") = ", prob(center, east, west, south, north)
          M(i_st, i_st) = M(i_st, i_st) &
           + (1.0d0 - prob(spin(x(i_v), z(i_v)), east, west, south, north)) &
            / DBLE(len_x * len_z)
          ! WRITE(*, *) 4
       END DO
    END DO

    DO i_st = 1, n_st, 1
       CALL decodeIndex(i_st, spin(1:len_x, 1:len_z))
       DO i_v = 1, len_x * len_z, 1
          CALL setDirection(  spin(1:len_x, 1:len_z), &
                              x(i_v), z(i_v), east, west, south, north)
          j_st = IEOR(ISHFT(1, i_v - 1), i_st - 1) + 1
          M(j_st, i_st) = prob(spin(x(i_v), z(i_v)), east, west, south, north) &
           / DBLE(len_x * len_z)
       END DO
    END DO
  END SUBROUTINE constructMatrixM

  SUBROUTINE indicateMatrixM(M)
    DOUBLE PRECISION, INTENT(in) :: M(1:, 1:)

    INTEGER(kind = 4) :: i_st, j_st

    DO i_st = 1, n_st, 1
			DO j_st = 1, n_st, 1
				WRITE(*, '(a, i5, a, i5, a, e9.2)') "M(", j_st, ", ", i_st, ") =", M(j_st, i_st)
			END DO
		END DO
	END SUBROUTINE indicateMatrixM

  SUBROUTINE constructMatrixV(V)
    DOUBLE PRECISION, INTENT(inout) :: V(1:, 1:)

    INTEGER(kind = 4) :: i_st, j_st
    INTEGER(kind = 4) :: spin(1:len_x, 1:len_z)

    V(1:n_st, 1:n_st) = 0.0d0

    DO i_st = 1, n_st, 1
       CALL decodeIndex(i_st, spin(1:len_x, 1:len_z))
       CALL shiftUpperHalf(spin(1:len_x, 1:len_z))
       CALL encodeSpin(spin(1:len_x, 1:len_z), j_st)
       V(j_st, i_st) = 1.0d0
    END DO
  END SUBROUTINE constructMatrixV

	SUBROUTINE indicateMatrixV(V)
    DOUBLE PRECISION, INTENT(in) :: V(1:, 1:)

    INTEGER(kind = 4) :: i_st, j_st

    DO i_st = 1, n_st, 1
			DO j_st = 1, n_st, 1
				WRITE(*, '(a, i5, a, i5, a, e9.2)') "V(", j_st, ", ", i_st, ") =", V(j_st, i_st)
			END DO
		END DO
	END SUBROUTINE indicateMatrixV

  SUBROUTINE constructVectorE(E)
    DOUBLE PRECISION, INTENT(out) :: E(1:)

    INTEGER(kind = 4) :: i_st, spin(1:len_x, 1:len_z)

		!TODO: ISHIFT, IANDを用いて改変．)
    E(1:n_st) = 0.0d0

    DO i_st = 1, n_st, 1
      CALL decodeIndex(i_st, spin(1:len_x, 1:len_z))
      CALL calcTotalEnergy(spin(1:len_x, 1:len_z), E(i_st))
    END DO
  END SUBROUTINE constructVectorE

  SUBROUTINE indicateVecotrE(E)
    DOUBLE PRECISION, INTENT(in) :: E(1:)

    INTEGER(kind = 4) :: i_st

    DO i_st = 1, n_st, 1
				WRITE(*, '(a, i5, a, e9.2)') "E(", i_st, ") =", E(i_st)
		END DO
	END SUBROUTINE indicateVecotrE
END MODULE calcMatrix

PROGRAM main
  USE global
  USE calcMatrix
  USE simulateSpin

  IMPLICIT NONE
  INTEGER(kind = 4) :: i_st, j_st, i_v
  INTEGER(kind = 4) :: east, west, south, north
  INTEGER(kind = 4), ALLOCATABLE :: spin(:, :)
  DOUBLE PRECISION :: EP
  DOUBLE PRECISION, ALLOCATABLE :: M(:, :), V(:, :), E(:)
  DOUBLE PRECISION, ALLOCATABLE :: T(:, :, :), temp_mat(:, :), ETP(:)

  DOUBLE PRECISION, ALLOCATABLE :: WR(:), WI(:), VL(:, :), P_eq(:, :), P_st(:, :)
  INTEGER(kind = 4), ALLOCATABLE :: WORK(:)
  INTEGER(kind = 4) :: INFO

  WRITE(0, '(a)') "len_x, len_z = ?"
  READ(*, *) len_x, len_z
  ALLOCATE(x(1:len_x * len_z), z(1:len_x * len_z))
  CALL foldSpinArray(x(1:len_x * len_z), z(1:len_x * len_z))
  !CHECK: foldSpinArray
  ! DO i_v = 1, len_x * len_z, 1
  !    WRITE(*, '(a, i5, a, i5, a, i5, a, i5)') "x(", i_v, ") = ", x(i_v), ", z(", i_v, ") = ", z(i_v)
  ! END DO

	n_st = ISHFT(1, len_x * len_z)
	ALLOCATE(   M(1:n_st, 1:n_st), V(1:n_st, 1:n_st), E(1:n_st), &
              T(0:len_x * len_z, 1:n_st, 1:n_st), temp_mat(1:n_st, 1:n_st), &
              ETP(0:len_x * len_z))

  ALLOCATE(   WR(1:n_st), WI(1:n_st), VL(1:n_st, 1:n_st), &
              P_eq(1:n_st, 1:n_st), P_st(1:n_st, 1:n_st), WORK(1:4*n_st))

  WRITE(0, '(a)') "beta, J = ?"
  READ(*, *) beta, J
  !CHECK: (beta, J)
  ! WRITE(*, '(a, e9.2, a, e9.2, a)') "(beta, J) = (", beta, ", ", J, ")"

  CALL makeDeltaEArray
  CALL makeProbArray

	CALL constructMatrixM(M(1:n_st, 1:n_st))
  CALL indicateMatrixM(M(1:n_st, 1:n_st))
	CALL constructMatrixV(V(1:n_st, 1:n_st))
	CALL indicateMatrixV(V(1:n_st, 1:n_st))
  CALL constructVectorE(E(1:n_st))
  CALL indicateVecotrE(E(1:n_st))

  temp_mat(1:n_st, 1:n_st) = M(1:n_st, 1:n_st)
  CALL dgeev("V", "V", n_st, temp_mat(1:n_st, 1:n_st), n_st, WR, WI, &
              VL, n_st, P_eq, n_st, WORK(1:4*n_st), 4*n_st, INFO)
	WRITE(*, '(a, i5)') "INFO(dgeev(M)) = ", INFO
  WRITE(*, '(a, e9.2)') "Norm(P_eq) = ", SUM(P_eq(1:n_st, 1))
  P_eq(1:n_st, 1) = P_eq(1:n_st, 1) / SUM(P_eq(1:n_st, 1))

  !CHECK: P_eq(lambda = 1)
  ! DO i_st = 1, n_st, 1
  !   WRITE(*, '(a, i5, a, e9.2)') "(M * P_eq - P_eq)(", i_st, ") = ", &
  !   dot_PRODUCT(M(i_st, 1:n_st), P_eq(1:n_st, 1)) - P_eq(i_st, 1)
  ! END DO

  temp_mat(1:n_st, 1:n_st) = V(1:n_st, 1:n_st)
  DO i_v = 1, len_x * len_z, 1
    temp_mat(1:n_st, 1:n_st) = MATMUL(M(1:n_st, 1:n_st), temp_mat(1:n_st, 1:n_st))
    T(i_v, 1:n_st, 1:n_st) = temp_mat(1:n_st, 1:n_st)
  END DO

  DEALLOCATE(M)

  temp_mat(1:n_st, 1:n_st) = T(len_x * len_z, 1:n_st, 1:n_st)
  CALL dgeev("V", "V", n_st, temp_mat(1:n_st, 1:n_st), n_st, WR, WI, &
              VL, n_st, P_st, n_st, WORK(1:4*n_st), 4*n_st, INFO)
	WRITE(*, '(a, i5)') "INFO(dgeev(T4)) = ", INFO
  WRITE(*, '(a, e9.2)') "Norm(P_st) = ", SUM(P_st(1:n_st, 1))
  P_st(1:n_st, 1) = P_st(1:n_st, 1) / SUM(P_st(1:n_st, 1))

  DEALLOCATE(temp_mat)

  !CHECK: P_st(lambda = 1)
  ! DO i_st = 1, n_st, 1
  !   WRITE(*, '(a, i5, a, e9.2)') "(T4 * P_st - P_st)(", i_st, ") = ", &
  !   dot_PRODUCT(T4(i_st, 1:n_st), P_st(1:n_st, 1)) - P_st(i_st, 1)
  ! END DO

  EP = dot_PRODUCT(E(1:n_st), P_st(1:n_st, 1))
  ETP(0) = dot_PRODUCT(E(1:n_st), MATMUL(V(1:n_st, 1:n_st), P_st(1:n_st, 1)))
  DO i_v = 1, len_x * len_z
    ETP(i_v) = &
    dot_PRODUCT(E(1:n_st), &
    MATMUL(T(i_v, 1:n_st, 1:n_st), P_st(1:n_st, 1)))
  END DO
  DEALLOCATE(E, V, T)
	!
  ! DO i_st = 1, n_st, 1
  !   WRITE(*, '(a, i5, a, e9.2)') "dP(", i_st, ") = ", P_st(i_st, 1) - P_eq(i_st, 1)
  ! END DO

  DEALLOCATE(P_eq, P_st)

  WRITE(*, '(a, e13.6)') "EP = ", EP
  DO i_v = 0, len_x * len_z, 1
    WRITE(*, '(a, i5, a, e13.6)') "ETP(", i_v, ") = ", ETP(i_v)
  END DO
END PROGRAM main
