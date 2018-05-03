! STEP 4: CG. ITERATION
! Init r and p, and q
DO i = 0,nz + 1
DO k = 0,nx + 1
rr(i,k) = 0
pp(i,k) = 0
qq(i,k) = 0
END DO
END DO

DO i = 1,nz
DO k = 1,nx
IF(wet(i,k))THEN
rr(i,k) = qstar(i,k) + ab(i,k)*dq(i+1, k) + at(i,k)*dq(i-1, k) + ae(i,k)*dq(i, k+1) + aw(i,k)*dq(i, k-1) - atot(i,k)*dq(i,k)
pp(i,k) = rr(i,k)
END IF
END DO
END DO

DO i = 1,nz
DO k = 1,nx
IF(wet(i,k))THEN
qq(i,k) = atot(i,k)*rr(i,k) - ab(i,k)*rr(i+1,k) - at(i,k)*rr(i-1,k) - ae(i,k)*rr(i,k+1) - aw(i,k)*rr(i, k-1)
END IF
END DO
END DO


r2 = ddot(nTot, RESHAPE(rr, [nTot]), 1, RESHAPE(rr, [nTot]), 1)
alpha = r2 / ddot(nTot, RESHAPE(pp, [nTot]), 1, RESHAPE(qq, [nTot]), 1)
! STEP 4.2: solve for dq using CG
DO WHILE ( r2 > eps )
nstop = nstop + 1
DO i = 1,nz
DO k = 1,nx
dq(i, k) = dq(i, k) + alpha * pp(i,k)
rr(i, k) = rr(i,k) - alpha * qq(i,k)

END DO
END DO
r2New = ddot(nTot, RESHAPE(rr, [nTot]), 1, RESHAPE(rr, [nTot]), 1)
beta = r2New/r2
DO i = 1,nz
DO k = 1,nx

pp(i,k) = rr(i,k) + beta * pp(i,k)
qq(i,k) = atot(i,k)*rr(i,k) - ab(i,k)*rr(i+1,k) - at(i,k)*rr(i-1,k) - ae(i,k)*rr(i,k+1) - aw(i,k)*rr(i, k-1) + beta*qq(i,k)

END DO
END DO
r2 = r2New
alpha = r2New / ddot(nTot, RESHAPE(pp, [nTot]), 1, RESHAPE(qq, [nTot]), 1)
END DO
