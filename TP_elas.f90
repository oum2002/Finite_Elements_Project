program TP_elas
    use mesh
    implicit none

    integer  :: i, r, c
    integer  :: n1, n2, n3
    integer  :: i_center
    real(8)  :: ybar, zbar
    real(8)  :: dist, dist_min
    real(8)  :: y1, y2, y3, z1, z2, z3
    real(8)  :: selem, y_c, z_c
    real(8)  :: gradN(2,3), M_mat(3,3), M_inv(3,3), B_mat(2,3), Kel(3,3), qel(3)
    real(8), allocatable :: Rsys(:,:), vecs(:,:), dep(:,:), coor(:,:)

    nord = 1
    call rw_mesh

    allocate(Rsys(nno,nno), vecs(nno,1), dep(nno,1), coor(2,nno))
    Rsys = 0.0d0
    vecs = 0.0d0
    dep  = 0.0d0

    gradN      = 0.0d0
    gradN(1,2) = 1.0d0
    gradN(2,3) = 1.0d0

    ybar = sum(ym) / nno
    zbar = sum(zm) / nno
    i_center = 1
    dist_min = (ym(1)-ybar)**2 + (zm(1)-zbar)**2
    do i = 2, nno
        dist = (ym(i)-ybar)**2 + (zm(i)-zbar)**2
        if (dist < dist_min) then
            dist_min = dist
            i_center = i
        end if
    end do
    write(*,*) 'Central node:', i_center, ' y=', ym(i_center), ' z=', zm(i_center)

    coor(1,:) = ym - ybar
    coor(2,:) = zm - zbar

    ! Assembly loop (use recentered coordinates)
    do i = 1, nelem
        n1 = nodelem(1,i) ; n2 = nodelem(2,i) ; n3 = nodelem(3,i)
        y1 = coor(1,n1) ; z1 = coor(2,n1)
        y2 = coor(1,n2) ; z2 = coor(2,n2)
        y3 = coor(1,n3) ; z3 = coor(2,n3)

        M_mat(1,:) = [1.0d0, y1, z1]
        M_mat(2,:) = [1.0d0, y2, z2]
        M_mat(3,:) = [1.0d0, y3, z3]

        call invermat(3, M_mat, M_inv)

        B_mat = matmul(gradN, M_inv)
        selem = 0.5d0 * abs((y2-y1)*(z3-z1) - (y3-y1)*(z2-z1))
        Kel   = selem * matmul(transpose(B_mat), B_mat)

        y_c = (y1+y2+y3) / 3.0d0
        z_c = (z1+z2+z3) / 3.0d0
        qel = selem * matmul(transpose(B_mat), [-z_c, y_c])

        do r = 1, 3
            do c = 1, 3
                Rsys(nodelem(r,i), nodelem(c,i)) = &
                Rsys(nodelem(r,i), nodelem(c,i)) + Kel(r,c)
            end do
            vecs(nodelem(r,i), 1) = vecs(nodelem(r,i), 1) + qel(r)
        end do
    end do

    ! Boundary condition: omega = 0 at central node
    Rsys(i_center, :)        = 0.0d0
    Rsys(i_center, i_center) = 1.0d0
    vecs(i_center, 1)        = 0.0d0

    call linsys(nno, 1, Rsys, vecs, dep)

    call ficgmshno_0(nno, nelem, nodelem, coor, dep(:,1))

    write(*,*) 'Done.  nno=', nno, '  nelem=', nelem
    write(*,*) 'omega  min=', minval(dep), '  max=', maxval(dep)

    contains

        subroutine invermat(n, raid, raidi)
            implicit none
            integer,  intent(in)  :: n
            real(8),  intent(in)  :: raid(n,n)
            real(8),  intent(out) :: raidi(n,n)
            integer              :: ipiv(n), info, nsize
            real(8)              :: workr(n)
            real(8), allocatable :: work(:)
            real(8)              :: raidt(n,n)
            raidt = raid
            call dgetrf(n, n, raidt, n, ipiv, info)
            call dgetri(n, raidt, n, ipiv, workr, -1, info)
            nsize = int(workr(1))
            allocate(work(nsize))
            call dgetri(n, raidt, n, ipiv, work, nsize, info)
            if (info /= 0) write(*,*) 'Warning: singular matrix, info=', info
            raidi = raidt
            deallocate(work)
        end subroutine invermat

        subroutine linsys(n, m, rai, force, dep)
            implicit none
            integer*4, intent(in) :: n, m
            integer               :: ipiv(2*n), info
            real*8                :: rai(n,n), rai1(n,n), dep(n,m), force(n,m)
            rai1 = rai
            dep  = force
            call dgetrf(n, n, rai1, n, ipiv, info)
            call dgetrs('N', n, m, rai1, n, ipiv, dep, n, info)
            if (info /= 0) write(*,*) 'Warning: system not solved, info=', info
        end subroutine linsys

        subroutine ficgmshno_0(nno, nelem, nodelem, coor, dep)
            implicit none
            integer, intent(in) :: nno, nelem
            integer, intent(in) :: nodelem(3,nelem)
            real(8), intent(in) :: coor(2,nno), dep(nno)
            integer :: i
            open(unit=1, file='resu.msh', form='formatted', status='unknown')
            write(1,'(a)') '$MeshFormat'
            write(1,'(a)') '2.2 0 8'
            write(1,'(a)') '$EndMeshFormat'
            write(1,'(a)') '$Nodes'
            write(1,'(I6)') nno
            do i = 1, nno
                write(1,'(I6,3(1X,F20.10))') i, 0.0d0, coor(1,i), coor(2,i)
            end do
            write(1,'(a)') '$EndNodes'
            write(1,'(a)') '$Elements'
            write(1,'(I6)') nelem
            do i = 1, nelem
                write(1,'(I6," 2 2       1      1",3I6)') i, nodelem(1,i), nodelem(2,i), nodelem(3,i)
            end do
            write(1,'(a)') '$EndElements'
            write(1,'(a)') '$NodeData'
            write(1,'(a)') '1'
            write(1,'(a)') '"omega"'
            write(1,'(a)') '1'
            write(1,'(a)') ' 0.0'
            write(1,'(a)') '3'
            write(1,'(a)') '0'
            write(1,'(a)') '1'
            write(1,'(I6)') nno
            do i = 1, nno
                write(1,'(I6,1X,F20.10)') i, dep(i)
            end do
            write(1,'(a)') '$EndNodeData'
            close(1)
        end subroutine ficgmshno_0

end program TP_elas