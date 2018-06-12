DO i = 1,nz
DO k = 1,nx
IF(wet(i,k))THEN
rr(i,k) = qstar(i,k) + ab(i,k)*dq(i+1, k) + at(i,k)*dq(i-1, k) + ae(i,k)*dq(i, k+1) + aw(i,k)*dq(i, k-1) - atot(i,k)*dq(i,k)
pp(i,k) = rr(i,k)
rrstar(i,k) = rr(i,k)
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

rtrstar = sdot(nTot, RESHAPE(rr, [nTot]), 1, RESHAPE(rrstar, [nTot]), 1)
r2 = sdot(nTot, RESHAPE(rr, [nTot]), 1, RESHAPE(rr, [nTot]), 1)
alphaden = dot_product(RESHAPE(qq, [nTot]), RESHAPE(rrstar, [nTot]))
alpha = rtrstar / alphaden
r2norm = snrm2(nTot, RESHAPE(rr, [nTot]),1)

! STEP 4.2: solve for dq using CG
DO WHILE ( r2norm > eps )
nstop = nstop + 1
DO i = 1,nz
DO k = 1,nx
IF(wet(i,k))THEN
ss(i, k) = rr(i, k) - alpha * qq(i,k)
END IF
END DO
END DO

DO i = 1,nz
DO k = 1,nx
IF(wet(i,k))THEN
ee(i, k) = atot(i,k)*ss(i,k) - ab(i,k)*ss(i+1,k) - at(i,k)*ss(i-1,k) - ae(i,k)*ss(i,k+1) - aw(i,k)*ss(i, k-1)
END IF
END DO
END DO

omeganum = sdot(nTot, RESHAPE(ee, [nTot]), 1, RESHAPE(ss, [nTot]), 1)
omegaden = dot_product(RESHAPE(ee, [nTot]), RESHAPE(ee, [nTot]))
omega = omeganum/omegaden

DO i = 1,nz
DO k = 1,nx
IF(wet(i,k))THEN
dq(i, k) = dq(i, k) + alpha * pp(i, k) + omega * ss(i, k)
rrnew(i, k) = ss(i, k) - omega * ee(i, k)
END IF
END DO
END DO

betanum = sdot(nTot, RESHAPE(rrnew, [nTot]), 1, RESHAPE(rrstar, [nTot]), 1)
betaden = sdot(nTot, RESHAPE(rr, [nTot]), 1, RESHAPE(rrstar, [nTot]), 1)
beta = (alpha/omega) * (betanum/ betaden)

DO i = 1,nz
DO k = 1,nx
IF(wet(i,k))THEN
pp(i,k) = rrnew(i, k) + beta * (pp(i, k) - omega * qq(i, k))
END IF
END DO
END DO

r2norm = snrm2(nTot, RESHAPE(rr, [nTot]),1)

DO i = 1,nz
DO k = 1,nx
IF(wet(i,k))THEN
qq(i,k) = atot(i,k)*pp(i,k) - ab(i,k)*pp(i+1,k) - at(i,k)*pp(i-1,k) - ae(i,k)*pp(i,k+1) - aw(i,k)*pp(i, k-1)
rr(i,k) = rrnew(i,k)
END IF
END DO
END DO

rtrstar = sdot(nTot, RESHAPE(rr, [nTot]), 1, RESHAPE(rrstar, [nTot]), 1)
r2 = sdot(nTot, RESHAPE(rr, [nTot]), 1, RESHAPE(rr, [nTot]), 1)
alphaden = sdot(nTot, RESHAPE(qq, [nTot]), 1, RESHAPE(rrstar, [nTot]), 1)
alpha = rtrstar / alphaden
r2norm = snrm2(nTot, RESHAPE(rr, [nTot]),1)

EXIT
END IF

END DO
