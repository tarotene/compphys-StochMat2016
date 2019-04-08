PROGRAM main
  IMPLICIT NONE
  INTEGER::i,j,N
  COMPLEX(KIND(0d0)),ALLOCATABLE::A(:,:),Ev(:)
  DOUBLE PRECISION :: beta

  N=4
  ALLOCATE(A(1:N,1:N),Ev(1:N))
  A(1:N,1:N)=dcmplx(0d0,0d0)
  Ev(1:N)=dcmplx(0d0,0d0)

  READ(*, *) beta
  A(1, 1) = dcmplx(1.0d0 - EXP(- 2 * beta), 0d0)
  A(1, 2) = dcmplx(0.5d0, 0d0)
  A(1, 3) = dcmplx(0.5d0, 0d0)
  A(1, 4) = dcmplx(0d0, 0d0)
  A(2, 1) = dcmplx(0.5d0 * EXP(- 2 * beta), 0d0)
  A(2, 2) = dcmplx(0d0, 0d0)
  A(2, 3) = dcmplx(0d0, 0d0)
  A(2, 4) = dcmplx(0.5d0 * EXP(- 2 * beta), 0d0)
  A(3, 1) = dcmplx(0.5d0 * EXP(- 2 * beta), 0d0)
  A(3, 2) = dcmplx(0d0, 0d0)
  A(3, 3) = dcmplx(0d0, 0d0)
  A(3, 4) = dcmplx(0.5d0 * EXP(- 2 * beta), 0d0)
  A(4, 1) = dcmplx(0d0, 0d0)
  A(4, 2) = dcmplx(0.5d0, 0d0)
  A(4, 3) = dcmplx(0.5d0, 0d0)
  A(4, 4) = dcmplx(1.0d0 - EXP(- 2 * beta), 0d0)

  CALL diag(N,A,Ev)

  WRITE(6,'(A,e12.5e2,A,e12.5e2,A)')"(",DBLE(Ev(N)),",",dimag(Ev(N)),")"

  DO j=1,N
     WRITE(6,'(A,e12.5e2,A,e12.5e2,A,$)')"(",DBLE(A(j,N)),",",dimag(A(j,N)),")  "
  END DO

  WRITE(*, *) EXP(- 2 * beta)

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
