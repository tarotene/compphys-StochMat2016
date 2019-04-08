PROGRAM main
  IMPLICIT NONE
  INTEGER :: i, j, N
  DOUBLE PRECISION :: beta
  COMPLEX (KIND(0d0)) :: acc2, acc6, rej2, rej6
  COMPLEX (KIND(0d0)), ALLOCATABLE :: M(:, :), V(:, :), D(:, :), Ev(:)
  COMPLEX (KIND(0d0)), ALLOCATABLE :: P(:, :), PR(:), PL(:)
  COMPLEX (KIND(0d0)) :: DEV, DE1, DE2, DE3, DE4

  N=16
  ALLOCATE(M(1:N, 1:N), V(1:N, 1:N), D(1:N, 1:N))
  ALLOCATE(Ev(1:N), P(1:N, 1:N), PR(1:N), PL(1:N))

  ! 各種基本行列への値の代入
  D(1:N, 1:N) = dcmplx(0d0,0d0)
  D(1,1) = dcmplx(0, 0d0)
  D(2,2) = dcmplx(0, 0d0)
  D(3,3) = dcmplx(0, 0d0)
  D(4,4) = dcmplx(0, 0d0)
  D(1,2:3) = dcmplx(-6.0d0, 0d0)
  D(2:3,1) = dcmplx(6.0d0, 0d0)
  D(4,2:3) = dcmplx(-2.0d0, 0d0)
  D(2:3,4) = dcmplx(2.0d0, 0d0)

  D(1,5) = dcmplx(-6.0d0, 0d0)
  D(2,6) = dcmplx(-2.0d0, 0d0)
  D(3,7) = dcmplx(-6.0d0, 0d0)
  D(4,8) = dcmplx(-2.0d0, 0d0)

  D(5,1) = dcmplx(6.0d0, 0d0)
  D(6,2) = dcmplx(2.0d0, 0d0)
  D(7,3) = dcmplx(6.0d0, 0d0)
  D(8,4) = dcmplx(2.0d0, 0d0)

  D(5,5) = dcmplx(0d0, 0d0)
  D(8,8) = dcmplx(0d0, 0d0)
  D(5,6) = dcmplx(-2.0d0, 0d0)
  D(5,7) = dcmplx(-6.0d0, 0d0)
  D(6,5) = dcmplx(2.0d0, 0d0)
  D(7,5) = dcmplx(6.0d0, 0d0)
  D(8,6) = dcmplx(-2.0d0, 0d0)
  D(8,7) = dcmplx(-6.0d0, 0d0)
  D(6,8) = dcmplx(2.0d0, 0d0)
  D(7,8) = dcmplx(6.0d0, 0d0)

  DO i = 1, 8, 1
     DO j  = 1, 8, 1
        D(i + 8, j + 8) = D(9 - i, 9 - j)
     END DO
  END DO

  D(1,9) = dcmplx(-6.0d0, 0d0)
  D(2,10) = dcmplx(-6.0d0, 0d0)
  D(3,11) = dcmplx(-2.0d0, 0d0)
  D(4,12) = dcmplx(-2.0d0, 0d0)
  D(5,13) = dcmplx(2.0d0, 0d0)
  D(6,14) = dcmplx(2.0d0, 0d0)
  D(7,15) = dcmplx(6.0d0, 0d0)
  D(8,16) = dcmplx(6.0d0, 0d0)

  DO i = 1, 8, 1
     DO j  = 1, 8, 1
        D(i + 8, j) = D(9 - i, 17 - j)
     END DO
  END DO

  WRITE(*, *) "beta"
  READ(*, *) beta

  M(1:N, 1:N) = dcmplx(0d0,0d0)

  acc2 = dcmplx(EXP(- 2 * beta), 0)
  acc6 = dcmplx(EXP(- 6 * beta), 0)
  rej2 = dcmplx(1 - EXP(- 2 * beta), 0)
  rej6 = dcmplx(1 - EXP(- 6 * beta), 0)

  M(1,1) = rej6
  M(2,2) = rej2*0.25d0 + rej6*0.25d0
  M(3,3) = rej2*0.25d0 + rej6*0.25d0
  M(4,4) = rej2
  M(1,2:3) = dcmplx(0.25d0, 0d0)
  M(2:3,1) = acc6*0.25d0
  M(4,2:3) = dcmplx(0.25d0, 0d0)
  M(2:3,4) = acc2*0.25d0

  M(1,5) = dcmplx(0.25d0, 0d0)
  M(2,6) = dcmplx(0.25d0, 0d0)
  M(3,7) = dcmplx(0.25d0, 0d0)
  M(4,8) = dcmplx(0.25d0, 0d0)

  M(5,1) = acc6*0.25d0
  M(6,2) = acc2*0.25d0
  M(7,3) = acc6*0.25d0
  M(8,4) = acc2*0.25d0

  M(5,5) = rej2*0.25d0 + rej6*0.25d0
  M(8,8) = rej2*0.25d0 + rej6*0.25d0
  M(5,6:7) = dcmplx(0.25d0, 0d0)
  M(6,5) = acc2*0.25d0
  M(7,5) = acc6*0.25d0
  M(8,6:7) = dcmplx(0.25d0, 0d0)
  M(6,8) = acc2*0.25d0
  M(7,8) = acc6*0.25d0

  DO i = 1, 8, 1
     DO j  = 1, 8, 1
        M(i + 8, j + 8) = M(9 - i, 9 - j)
     END DO
  END DO

  M(1,9) = dcmplx(0.25d0, 0d0)
  M(2,10) = dcmplx(0.25d0, 0d0)
  M(3,11) = dcmplx(0.25d0, 0d0)
  M(4,12) = dcmplx(0.25d0, 0d0)
  M(5,13) = acc2*0.25d0
  M(6,14) = acc2*0.25d0
  M(7,15) = acc6*0.25d0
  M(8,16) = acc6*0.25d0

  DO i = 1, 8, 1
     DO j  = 1, 8, 1
        M(i + 8, j) = M(9 - i, 17 - j)
     END DO
  END DO

  V(1:N, 1:N) = dcmplx(0d0,0d0)
  V(1,1) = dcmplx(1.0d0, 0d0)
  V(2,3) = dcmplx(1.0d0, 0d0)
  V(3,2) = dcmplx(1.0d0, 0d0)
  V(4,4) = dcmplx(1.0d0, 0d0)

  V(5,5) = dcmplx(1.0d0, 0d0)
  V(6,7) = dcmplx(1.0d0, 0d0)
  V(7,6) = dcmplx(1.0d0, 0d0)
  V(8,8) = dcmplx(1.0d0, 0d0)

  V(9,9) = dcmplx(1.0d0, 0d0)
  V(10,11) = dcmplx(1.0d0, 0d0)
  V(11,10) = dcmplx(1.0d0, 0d0)
  V(12,12) = dcmplx(1.0d0, 0d0)

  V(13,13) = dcmplx(1.0d0, 0d0)
  V(14,15) = dcmplx(1.0d0, 0d0)
  V(15,14) = dcmplx(1.0d0, 0d0)
  V(16,16) = dcmplx(1.0d0, 0d0)

  ! 非平衡定常状態の固有ベクトル（固有値1）の定義
  Ev(1:N) = dcmplx(0d0,0d0)
  CALL diag(N, P, Ev)
  PR(1:N) = P(1:N, N)

  DEV = dot_PRODUCT(MATMUL(V(1:N, 1:N), PR(1:N)), MATMUL(D(1:N, 1:N), PR(1:N)))

  WRITE(*, *) DEV

  STOP
END PROGRAM main
!---------------------

SUBROUTINE diag(N,A,Ev)
  ! sikinote
  !date      : 2015/07/07
  !            2015/08/21
  !developer : sikino & fernandeskun
  IMPLICIT NONE
  INTEGER,INTENT(in)::N
  COMPLEX(KIND(0d0)),INTENT(inout)::A(1:N,1:N)
  COMPLEX(KIND(0d0)),INTENT(out)::Ev(1:N)

  INTEGER::ilo,ihi,info,lwork,turn(1:N),tmp,i
  DOUBLE PRECISION::SCALE(1:N),rwork(1:N)
  COMPLEX(KIND(0d0))::tau(1:N-1),w(1:N),z(1:N,1:N),Q(1:N,1:N),vr(1:N,1:N),tw(1:3)
  COMPLEX(KIND(0d0)),ALLOCATABLE::work(:)

  tau(1:N-1)=dcmplx(0d0,0d0)
  w(1:N)=dcmplx(0d0,0d0)
  z(1:N,1:N)=dcmplx(0d0,0d0)
  Q(1:N,1:N)=dcmplx(0d0,0d0)
  vr(1:N,1:N)=dcmplx(0d0,0d0)
  tw(1:3)=dcmplx(0d0,0d0)
  Ev(1:N)=dcmplx(0d0,0d0)

  !Equilibrate matrix A to equilibrated matrix A' to improve accuracy.
  !            i   i io  i   o    o     o     o
  CALL zgebal('P', N, A, N, ilo, ihi, scale, info)
  IF(info.NE.0)THEN
     WRITE(6,'(A,i0)')" At zgebal error, info --> ",info
     WRITE(6,'(A)')" Program stop"
     STOP
  ENDIF

  !Size Query
  CALL zgehrd(N, ilo, ihi, A, N, tau, tw, -1, info)
  lwork=NINT(DBLE(tw(1)))
  ALLOCATE(work(1:lwork)); work=dcmplx(0d0,0d0)

  !Degenerate matrix A to upper Hessenberg matrix H.
  !           i   i    i  io  i   o    i      i      o
  CALL zgehrd(N, ilo, ihi, A, N, tau, work, lwork, info)
  IF(info.NE.0)THEN
     WRITE(6,'(A,i0)')" At zgehrd error, info --> ",info
     WRITE(6,'(A)')" Program stop"
     STOP
  ENDIF
  DEALLOCATE(work)

  Q=a
  !Size Query
  CALL zunghr(N, ilo, ihi, Q, N, tau, tw, -1, info)
  lwork=NINT(DBLE(tw(1)))
  ALLOCATE(work(1:lwork)); work=dcmplx(0d0,0d0)

  !Make complex unitary matrix Q from upper Hessenberg matrix H.
  !           i   i    i  io  i   i    i      i      o
  CALL zunghr(N, ilo, ihi, Q, N, tau, work, lwork, info)
  IF(info.NE.0)THEN
     WRITE(6,'(A,i0)')" At zunghr error, info --> ",info
     WRITE(6,'(A)')" Program stop"
     STOP
  ENDIF
  DEALLOCATE(work)

  z=Q
  !Size Query
  CALL zhseqr('S', 'V', N, ilo, ihi, A, N, Ev, z, N, tw, -1, info)
  lwork=NINT(DBLE(tw(1)))
  ALLOCATE(work(1:lwork)); work=dcmplx(0d0,0d0)

  !Get eigenvalue of upper Hessenberg matrix H and Get Schur vector.
  !                     i   i    i  io  i   o  o  i   i      i      o
  CALL zhseqr('S', 'V', N, ilo, ihi, A, N, Ev, z, N, work, lwork, info)
  IF(info.NE.0)THEN
     WRITE(6,'(A,i0)')" At zhseqr error, info --> ",info
     WRITE(6,'(A)')" Program stop"
     STOP
  ENDIF
  DEALLOCATE(work)

  !Get right eigenvector X from upper triangular matrix T.
  ALLOCATE(work(1:2*N))
  vr=z
  !                        i  i  i         o  i  i   o   i      i      i
  CALL ztrevc('R', 'B', 0, N, A, N, 0, 1, vr, N, N, tmp, work, rwork, info)
  IF(info.NE.0)THEN
     WRITE(6,'(A,i0)')" At zhseqr error, info --> ",info
     WRITE(6,'(A)')" Program stop"
     STOP
  ENDIF
  DEALLOCATE(work)

  !Transrate right eigenvector X of Equilibrated matrix A' to right eigenvector of matrix A
  !                     i   i    i     i    i   o  i   o
  CALL zgebak('P', 'R', N, ilo, ihi, scale, N, vr, N, info)
  IF(info.NE.0)THEN
     WRITE(6,'(A,i0)')" At zhseqr error, info --> ",info
     WRITE(6,'(A)')" Program stop"
     STOP
  ENDIF

  A=vr

!swap Eigenvectol as same arrangement for Eigenvalue
  CALL sortdp2(N,Ev,turn)

  Q=A
  DO i=1,N
     tmp=turn(i)
     A(1:N,i)=Q(1:N,tmp)
  ENDDO
  RETURN

  !sort Eigenvalue of real part from small to big.
CONTAINS
  SUBROUTINE sortdp2(N,DATA,turn)
    IMPLICIT NONE
    INTEGER::i,ti,j,N,turn(1:N)
    COMPLEX(KIND(0d0))::DATA(1:N),tmp

    DO i=1,N
       turn(i)=i
    ENDDO

    DO i=1,N-1
       DO j=i+1,N
          IF(DBLE(DATA(i)) > DBLE(DATA(j)))THEN
             tmp=DATA(i)
             DATA(i)=DATA(j)
             DATA(j)=tmp

             ti=turn(i)
             turn(i)=turn(j)
             turn(j)=ti
          END IF
       END DO
    END DO

    RETURN
  END SUBROUTINE sortdp2
END SUBROUTINE diag

SUBROUTINE set_matD(D)
  IMPLICIT NONE
  COMPLEX (KIND(0d0)), INTENT(inout) :: D(1:, 1:)

  INTEGER (kind = 4) :: i, j

  D(1,1) = dcmplx(0, 0d0)
  D(2,2) = dcmplx(0, 0d0)
  D(3,3) = dcmplx(0, 0d0)
  D(4,4) = dcmplx(0, 0d0)
  D(1,2:3) = dcmplx(-6.0d0, 0d0)
  D(2:3,1) = dcmplx(6.0d0, 0d0)
  D(4,2:3) = dcmplx(-2.0d0, 0d0)
  D(2:3,4) = dcmplx(2.0d0, 0d0)

  D(1,5) = dcmplx(-6.0d0, 0d0)
  D(2,6) = dcmplx(-2.0d0, 0d0)
  D(3,7) = dcmplx(-6.0d0, 0d0)
  D(4,8) = dcmplx(-2.0d0, 0d0)

  D(5,1) = dcmplx(6.0d0, 0d0)
  D(6,2) = dcmplx(2.0d0, 0d0)
  D(7,3) = dcmplx(6.0d0, 0d0)
  D(8,4) = dcmplx(2.0d0, 0d0)

  D(5,5) = dcmplx(0d0, 0d0)
  D(8,8) = dcmplx(0d0, 0d0)
  D(5,6) = dcmplx(-2.0d0, 0d0)
  D(5,7) = dcmplx(-6.0d0, 0d0)
  D(6,5) = dcmplx(2.0d0, 0d0)
  D(7,5) = dcmplx(6.0d0, 0d0)
  D(8,6) = dcmplx(-2.0d0, 0d0)
  D(8,7) = dcmplx(-6.0d0, 0d0)
  D(6,8) = dcmplx(2.0d0, 0d0)
  D(7,8) = dcmplx(6.0d0, 0d0)

  DO i = 1, 8, 1
     DO j  = 1, 8, 1
        D(i + 8, j + 8) = D(9 - i, 9 - j)
     END DO
  END DO

  D(1,9) = dcmplx(-6.0d0, 0d0)
  D(2,10) = dcmplx(-6.0d0, 0d0)
  D(3,11) = dcmplx(-2.0d0, 0d0)
  D(4,12) = dcmplx(-2.0d0, 0d0)
  D(5,13) = dcmplx(2.0d0, 0d0)
  D(6,14) = dcmplx(2.0d0, 0d0)
  D(7,15) = dcmplx(6.0d0, 0d0)
  D(8,16) = dcmplx(6.0d0, 0d0)
END SUBROUTINE set_matD

SUBROUTINE set_matM(N, M, beta)
  IMPLICIT NONE

  INTEGER (kind = 4), INTENT(in) :: N
  COMPLEX (KIND(0d0)), INTENT(inout) :: M(1:, 1:)
  DOUBLE PRECISION, INTENT(in) :: beta

  INTEGER (kind = 4) :: i, j
  COMPLEX (KIND(0d0)) :: acc2, acc6, rej2, rej6

  acc2 = dcmplx(EXP(- 2 * beta), 0)
  acc6 = dcmplx(EXP(- 6 * beta), 0)
  rej2 = dcmplx(1 - EXP(- 2 * beta), 0)
  rej6 = dcmplx(1 - EXP(- 6 * beta), 0)

  M(1,1) = rej6
  M(2,2) = rej2*0.25d0 + rej6*0.25d0
  M(3,3) = rej2*0.25d0 + rej6*0.25d0
  M(4,4) = rej2
  M(1,2:3) = dcmplx(0.25d0, 0d0)
  M(2:3,1) = acc6*0.25d0
  M(4,2:3) = dcmplx(0.25d0, 0d0)
  M(2:3,4) = acc2*0.25d0

  M(1,5) = dcmplx(0.25d0, 0d0)
  M(2,6) = dcmplx(0.25d0, 0d0)
  M(3,7) = dcmplx(0.25d0, 0d0)
  M(4,8) = dcmplx(0.25d0, 0d0)

  M(5,1) = acc6*0.25d0
  M(6,2) = acc2*0.25d0
  M(7,3) = acc6*0.25d0
  M(8,4) = acc2*0.25d0

  M(5,5) = rej2*0.25d0 + rej6*0.25d0
  M(8,8) = rej2*0.25d0 + rej6*0.25d0
  M(5,6:7) = dcmplx(0.25d0, 0d0)
  M(6,5) = acc2*0.25d0
  M(7,5) = acc6*0.25d0
  M(8,6:7) = dcmplx(0.25d0, 0d0)
  M(6,8) = acc2*0.25d0
  M(7,8) = acc6*0.25d0

  DO i = 1, 8, 1
     DO j  = 1, 8, 1
        M(i + 8, j + 8) = M(9 - i, 9 - j)
     END DO
  END DO

  M(1,9) = dcmplx(0.25d0, 0d0)
  M(2,10) = dcmplx(0.25d0, 0d0)
  M(3,11) = dcmplx(0.25d0, 0d0)
  M(4,12) = dcmplx(0.25d0, 0d0)
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

SUBROUTINE set_matV(N, V)
  IMPLICIT NONE

  INTEGER (kind = 4), INTENT(in) :: N
  COMPLEX (KIND(0d0)), INTENT(inout) :: V(1:, 1:)

  V(1,1) = dcmplx(1.0d0, 0d0)
  V(2,3) = dcmplx(1.0d0, 0d0)
  V(3,2) = dcmplx(1.0d0, 0d0)
  V(4,4) = dcmplx(1.0d0, 0d0)

  V(5,5) = dcmplx(1.0d0, 0d0)
  V(6,7) = dcmplx(1.0d0, 0d0)
  V(7,6) = dcmplx(1.0d0, 0d0)
  V(8,8) = dcmplx(1.0d0, 0d0)

  V(9,9) = dcmplx(1.0d0, 0d0)
  V(10,11) = dcmplx(1.0d0, 0d0)
  V(11,10) = dcmplx(1.0d0, 0d0)
  V(12,12) = dcmplx(1.0d0, 0d0)

  V(13,13) = dcmplx(1.0d0, 0d0)
  V(14,15) = dcmplx(1.0d0, 0d0)
  V(15,14) = dcmplx(1.0d0, 0d0)
  V(16,16) = dcmplx(1.0d0, 0d0)
END SUBROUTINE set_matV
