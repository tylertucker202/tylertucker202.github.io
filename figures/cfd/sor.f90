!*****************
DO nsor = 1,nstop
!*****************
perr = 0.0
! STEP 1: predict pressure correction
DO i = 1,nz
DO k = 1,nx
 IF(wet(i,k))THEN
 q1 = dq(i,k)
 term1 = qstar(i,k) + & 
  &      at(i,k)*dq(i-1,k) + ab(i,k)*dq(i+1,k) + & 
  &      aw(i,k)*dq(i,k-1) + ae(i,k)*dq(i,k+1)
 q2 = (1.0-omega)*q1 + omega*term1/atot(i,k) 
 dq(i,k) = q2
 perr = MAX(ABS(q2-q1),perr)
 END IF
IF(perr <= peps)THEN
  nstop = nsor
  GOTO 33
END IF
!********************
END DO
!********************
