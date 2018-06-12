! File poisson-c-grid.f90
DO i = 1,nz
DO k = 1,nx
IF(wet(i,k))THEN
 Ax(i,k) = at(i,k)*x(i-1,k) + ab(i,k)*x(i+1,k) + aw(i,k)*x(i,k-1) + ae(i,k)*x(i,k+1)
END IF
END DO
END DO
