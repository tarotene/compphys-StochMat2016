PROGRAM main
  IMPLICIT NONE
  INTEGER (kind = 4) :: i, j, N
  DOUBLE PRECISION :: beta
  DOUBLE PRECISION :: acc2, acc6, rej2, rej6
  DOUBLE PRECISION, ALLOCATABLE :: M(:, :), V(:, :)
  DOUBLE PRECISION :: DEV, DE1, DE2, DE3, DE4

  DOUBLE PRECISION, ALLOCATABLE :: EIGR(:), EIGI(:), VL(:,:), VR(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: WORK(:), PstR(:), PstL(:)
  DOUBLE PRECISION :: DE(0:5),Nom
  INTEGER (kind = 4) :: LWORK,INFO

  ALLOCATE(M(1:16, 1:16), V(1:16, 1:16))
  ALLOCATE(EIGR(1:16), EIGI(1:16),VL(1:16,1:16),VR(1:16,1:16))
  ALLOCATE(WORK(1:16*4), PstR(1:16), PstL(1:16))

  ! 各種基本行列への値の代入

  ! D(1, 2) = -6.0d0
  ! D(1, 3) = -6.0d0
  ! D(1, 4) = -4.0d0
  ! D(1, 5) = -6.0d0
  ! D(1, 6) = -8.0d0
  ! D(1, 7) = -12.0d0
  ! D(1, 8) = -6.0d0
  ! D(1, 9) = -6.0d0
  ! D(1, 10) = -12.0d0
  ! D(1, 11) = -8.0d0
  ! D(1, 12) = -6.0d0
  ! D(1, 13) = -4.0d0
  ! D(1, 14) = -6.0d0
  ! D(1, 15) = -6.0d0
  ! D(1, 16) = 0d0
  !
  ! D(2, 3) = 0d0
  ! D(2, 4) = 2.0d0
  ! D(2, 5) = 0d0
  ! D(2, 6) = -2.0d0
  ! D(2, 7) = -6.0d0
  ! D(2, 8) = 0d0
  ! D(2, 9) = 0d0
  ! D(2, 10) = -6.0d0
  ! D(2, 11) = -2.0d0
  ! D(2, 12) = 0d0
  ! D(2, 13) = 2.0d0
  ! D(2, 14) = 0d0
  ! D(2, 15) = 0d0
  ! D(2, 16) = 6.0d0
  !
  ! D(3, 4) = 2.0d0
  ! D(3, 5) = 0d0
  ! D(3, 6) = -2.0d0
  ! D(3, 7) = -6.0d0
  ! D(3, 8) = 0d0
  ! D(3, 9) = 0d0
  ! D(3, 10) = -6.0d0
  ! D(3, 11) = -2.0d0
  ! D(3, 12) = 0d0
  ! D(3, 13) = 2.0d0
  ! D(3, 14) = 0d0
  ! D(3, 15) = 0d0
  ! D(3, 16) = 6.0d0
  !
  ! D(4, 5) = -2.0d0
  ! D(4, 6) = -4.0d0
  ! D(4, 7) = -8.0d0
  ! D(4, 8) = -2.0d0
  ! D(4, 9) = -2.0d0
  ! D(4, 10) = -8.0d0
  ! D(4, 11) = -4.0d0
  ! D(4, 12) = -2.0d0
  ! D(4, 13) = 0d0
  ! D(4, 14) = -2.0d0
  ! D(4, 15) = -2.0d0
  ! D(4, 16) = 4.0d0
  !
  ! D(5, 6) = -2.0d0
  ! D(5, 7) = -6.0d0
  ! D(5, 8) = 0d0
  ! D(5, 9) = 0d0
  ! D(5, 10) = -6.0d0
  ! D(5, 11) = -2.0d0
  ! D(5, 12) = 0d0
  ! D(5, 13) = 2.0d0
  ! D(5, 14) = 0d0
  ! D(5, 15) = 0d0
  ! D(5, 16) = 6.0d0
  !
  ! D(6, 7) = -4.0d0
  ! D(6, 8) = 2.0d0
  ! D(6, 9) = 2.0d0
  ! D(6, 10) = -4.0d0
  ! D(6, 11) = 0d0
  ! D(6, 12) = 2.0d0
  ! D(6, 13) = 4.0d0
  ! D(6, 14) = 2.0d0
  ! D(6, 15) = 2.0d0
  ! D(6, 16) = 8.0d0
  !
  ! D(7, 8) = 6.0d0
  ! D(7, 9) = 6.0d0
  ! D(7, 10) = 0d0
  ! D(7, 11) = 4.0d0
  ! D(7, 12) = 6.0d0
  ! D(7, 13) = 8.0d0
  ! D(7, 14) = 6.0d0
  ! D(7, 15) = 6.0d0
  ! D(7, 16) = 12.0d0
  !
  ! D(8, 9) = 0d0
  ! D(8, 10) = -6.0d0
  ! D(8, 11) = -2.0d0
  ! D(8, 12) = 0d0
  ! D(8, 13) = 2.0d0
  ! D(8, 14) = 0d0
  ! D(8, 15) = 0d0
  ! D(8, 16) = 6.0d0
  !
  ! D(9, 10) = -6.0d0
  ! D(9, 11) = -2.0d0
  ! D(9, 12) = 0d0
  ! D(9, 13) = 2.0d0
  ! D(9, 14) = 0d0
  ! D(9, 15) = 0d0
  ! D(9, 16) = 6.0d0
  !
  ! D(10, 11) = 4.0d0
  ! D(10, 12) = 6.0d0
  ! D(10, 13) = 8.0d0
  ! D(10, 14) = 6.0d0
  ! D(10, 15) = 6.0d0
  ! D(10, 16) = 12.0d0
  !
  ! D(11, 12) = 2.0d0
  ! D(11, 13) = 4.0d0
  ! D(11, 14) = 2.0d0
  ! D(11, 15) = 2.0d0
  ! D(11, 16) = 8.0d0
  !
  ! D(12, 13) = 2.0d0
  ! D(12, 14) = 0d0
  ! D(12, 15) = 0d0
  ! D(12, 16) = 6.0d0
  !
  ! D(13, 14) = -2.0d0
  ! D(13, 15) = -2.0d0
  ! D(13, 16) = 4.0d0
  !
  ! D(14, 15) = 0d0
  ! D(14, 16) = 6.0d0
  !
  ! D(15, 16) = 6.0d0
  !
  ! DO i = 1, 16, 1
  !    DO j = 1, 16, 1
  !       IF ( i > j ) THEN
  !          D(i, j) = - D(j, i)
  !       END IF
  !    END DO
  ! END DO

  WRITE(*, *) "beta"
  READ(*, *) beta

  M(1:16, 1:16) = 0d0
  CALL set_matM(M, beta)

  ! acc2 = EXP(- 2 * beta)
  ! acc6 = EXP(- 6 * beta)
  ! rej2 = 1 - EXP(- 2 * beta)
  ! rej6 = 1 - EXP(- 6 * beta)
  !
  ! M(1,1) = rej6
  ! M(2,2) = rej2*0.25d0 + rej6*0.25d0
  ! M(3,3) = rej2*0.25d0 + rej6*0.25d0
  ! M(4,4) = rej2
  ! M(1,2:3) = 0.25d0
  ! M(2:3,1) = acc6*0.25d0
  ! M(4,2:3) = 0.25d0
  ! M(2:3,4) = acc2*0.25d0
  !
  ! M(1,5) = 0.25d0
  ! M(2,6) = 0.25d0
  ! M(3,7) = 0.25d0
  ! M(4,8) = 0.25d0
  !
  ! M(5,1) = acc6*0.25d0
  ! M(6,2) = acc2*0.25d0
  ! M(7,3) = acc6*0.25d0
  ! M(8,4) = acc2*0.25d0
  !
  ! M(5,5) = rej2*0.25d0 + rej6*0.25d0
  ! M(8,8) = rej2*0.25d0 + rej6*0.25d0
  ! M(5,6:7) = 0.25d0
  ! M(6,5) = acc2*0.25d0
  ! M(7,5) = acc6*0.25d0
  ! M(8,6:7) = 0.25d0
  ! M(6,8) = acc2*0.25d0
  ! M(7,8) = acc6*0.25d0
  !
  ! DO i = 1, 8, 1
  !    DO j  = 1, 8, 1
  !       M(i + 8, j + 8) = M(9 - i, 9 - j)
  !    END DO
  ! END DO
  !
  ! M(1,9) = 0.25d0
  ! M(2,10) = 0.25d0
  ! M(3,11) = 0.25d0
  ! M(4,12) = 0.25d0
  ! M(5,13) = acc2*0.25d0
  ! M(6,14) = acc2*0.25d0
  ! M(7,15) = acc6*0.25d0
  ! M(8,16) = acc6*0.25d0
  !
  ! DO i = 1, 8, 1
  !    DO j  = 1, 8, 1
  !       M(i + 8, j) = M(9 - i, 17 - j)
  !    END DO
  ! END DO

  V(1:16, 1:16) = 0d0
	call set_matV(V)
  ! V(1,1) = 1.0d0
  ! V(2,3) = 1.0d0
  ! V(3,2) = 1.0d0
  ! V(4,4) = 1.0d0
	!
  ! V(5,5) = 1.0d0
  ! V(6,7) = 1.0d0
  ! V(7,6) = 1.0d0
  ! V(8,8) = 1.0d0
	!
  ! V(9,9) = 1.0d0
  ! V(10,11) = 1.0d0
  ! V(11,10) = 1.0d0
  ! V(12,12) = 1.0d0
	!
  ! V(13,13) = 1.0d0
  ! V(14,15) = 1.0d0
  ! V(15,14) = 1.0d0
  ! V(16,16) = 1.0d0

  ! 非平衡定常状態の固有ベクトル（固有値1）の定義
  Ev(1:16) = 0d0
  P1(1:16, 1:16) = MATMUL(M(1:16, 1:16), V(1:16, 1:16))
  P2(1:16, 1:16) = MATMUL(M(1:16, 1:16), P1(1:16, 1:16))
  P3(1:16, 1:16) = MATMUL(M(1:16, 1:16), P2(1:16, 1:16))
  P4(1:16, 1:16) = MATMUL(M(1:16, 1:16), P3(1:16, 1:16))

  !CALL diag(N, P4, Ev)
  A=DBLE(P4)
  CALL dgeev('V','V',N,A,N,EIGR,EIGI,VL,N,VR,N,WORK,LWORK,INFO)
  WRITE(6,*) 'INFO=',INFO

  !do I = 1, 16, 1
  !  write(6,*) I,EIGR(I),EIGI(I)
  !end do

  PstR=VR(1:16,1)
  PstL=VL(1:16,1)
  DO I = 1, 16, 1
     WRITE(6,*) PstR(I)
  END DO
  Nom=SUM(PstR)
  PstR=PstR/Nom
  !  A=dble(P4)
  !  write(6,*) dot_product(PstL,MATMUL(A,PstR))

  DE=0.0d0
  DO J = 1, 16, 1
     DE(1)=DE(1)+DBLE(dot_PRODUCT(D(J,1:16),V(1:16,J)))*PstR(J)
  END DO
  DO J = 1, 16, 1
     DE(2)=DE(2)+DBLE(dot_PRODUCT(D(J,1:16),P1(1:16,J)))*PstR(J)
  END DO
  DO J = 1, 16, 1
     DE(3)=DE(3)+DBLE(dot_PRODUCT(D(J,1:16),P2(1:16,J)))*PstR(J)
  END DO
  DO J = 1, 16, 1
     DE(4)=DE(4)+DBLE(dot_PRODUCT(D(J,1:16),P3(1:16,J)))*PstR(J)
  END DO
  DO J = 1, 16, 1
     DE(5)=DE(5)+DBLE(dot_PRODUCT(D(J,1:16),P4(1:16,J)))*PstR(J)
  END DO
  WRITE(6,*) (DE(I),I=1,5)

  PR(1:16) = MATMUL(D(1:16, 1:16), P4(1:16, N))
  PLV(1:16) = MATMUL(V(1:16, 1:16), P4(1:16, N))
  PL1(1:16) = MATMUL(M(1:16, 1:16), PLV(1:16))
  PL2(1:16) = MATMUL(M(1:16, 1:16), PL1(1:16))
  PL3(1:16) = MATMUL(M(1:16, 1:16), PL2(1:16))
  PL4(1:16) = MATMUL(M(1:16, 1:16), PL3(1:16))

  DEV = dot_PRODUCT(PLV(1:16), PR(1:16))
  DE1 = dot_PRODUCT(PL1(1:16), PR(1:16))
  DE2 = dot_PRODUCT(PL2(1:16), PR(1:16))
  DE3 = dot_PRODUCT(PL3(1:16), PR(1:16))
  DE4 = dot_PRODUCT(PL4(1:16), PR(1:16))

  WRITE(*, *) DEV, DE1, DE2, DE3, DE4

  STOP
END PROGRAM main
!---------------------

SUBROUTINE set_matM(M, beta)
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(inout) :: M(1:, 1:)
  DOUBLE PRECISION, INTENT(in) :: beta

  INTEGER (kind = 4) :: i, j
  DOUBLE PRECISION :: acc2, acc6, rej2, rej6

  acc2 = EXP(- 2 * beta)
  acc6 = EXP(- 6 * beta)
  rej2 = 1 - EXP(- 2 * beta)
  rej6 = 1 - EXP(- 6 * beta)

  M(1,1) = rej6
  M(2,2) = rej2*0.25d0 + rej6*0.25d0
  M(3,3) = rej2*0.25d0 + rej6*0.25d0
  M(4,4) = rej2
  M(1,2:3) = 0.25d0
  M(2:3,1) = acc6*0.25d0
  M(4,2:3) = 0.25d0
  M(2:3,4) = acc2*0.25d0

  M(1,5) = 0.25d0
  M(2,6) = 0.25d0
  M(3,7) = 0.25d0
  M(4,8) = 0.25d0

  M(5,1) = acc6*0.25d0
  M(6,2) = acc2*0.25d0
  M(7,3) = acc6*0.25d0
  M(8,4) = acc2*0.25d0

  M(5,5) = rej2*0.25d0 + rej6*0.25d0
  M(8,8) = rej2*0.25d0 + rej6*0.25d0
  M(5,6:7) = 0.25d0
  M(6,5) = acc2*0.25d0
  M(7,5) = acc6*0.25d0
  M(8,6:7) = 0.25d0
  M(6,8) = acc2*0.25d0
  M(7,8) = acc6*0.25d0

  DO i = 1, 8, 1
     DO j  = 1, 8, 1
        M(i + 8, j + 8) = M(9 - i, 9 - j)
     END DO
  END DO

  M(1,9) = 0.25d0
  M(2,10) = 0.25d0
  M(3,11) = 0.25d0
  M(4,12) = 0.25d0
  M(5,13) = acc2*0.25d0
  M(6,14) = acc2*0.25d0
  M(7,15) = acc6*0.25d0
  M(8,16) = acc6*0.25d0

  DO i = 1, 8, 1
     DO j  = 1, 8, 1
        M(i + 8, j) = M(9 - i, 17 - j)
     END DO
  END DO

END SUBROUTINE set_matM

SUBROUTINE set_matV(V)
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(inout) :: V(1:, 1:)

  V(1,1) = 1.0d0
  V(2,3) = 1.0d0
  V(3,2) = 1.0d0
  V(4,4) = 1.0d0

  V(5,5) = 1.0d0
  V(6,7) = 1.0d0
  V(7,6) = 1.0d0
  V(8,8) = 1.0d0

  V(9,9) = 1.0d0
  V(10,11) = 1.0d0
  V(11,10) = 1.0d0
  V(12,12) = 1.0d0

  V(13,13) = 1.0d0
  V(14,15) = 1.0d0
  V(15,14) = 1.0d0
  V(16,16) = 1.0d0
END SUBROUTINE set_matV
