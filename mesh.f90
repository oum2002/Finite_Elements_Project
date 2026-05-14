module mesh
    implicit none ! Force la d�claration de chaque variable

    integer :: nord, nce, noce, nseg, nno, nelem, nbord
    integer, allocatable :: ics(:), icseg(:,:), iaseg(:), ibseg(:), ind_consec(:)
    integer, allocatable :: nodelem(:,:), ncelem(:), nodebord(:,:), ibord(:), nelbord(:)
    real(8), allocatable :: ysec(:), zsec(:), ym(:), zm(:)

contains

    subroutine rw_mesh
        implicit none ! Important ici aussi
        integer :: i, j, k, n, info, n_entities, nelu, iel, ibo
        integer :: i1, i2, i3, i4, j1, j2, j3, ii1, ii2
        real(8) :: x1, y1, z1, x2, y2, z2
        integer :: itmail(100)
        logical :: ficex
    
    
        itmail=0
    
        itmail(15) = 1
    
        itmail(1)  = 2   !bord lin�aire
        itmail(2)  = 3   !triangle lin�aire
    
        itmail(8)  = 3   !bord quadr.
        itmail(9)  = 6   !triangle quadr.
    
        itmail(28) = 6   !bord C0-5
        itmail(25) = 21  !triangle C0-5
    

        inquire( file = 'MAIL.msh' , exist = ficex)
        if(.not.ficex) error stop "fichier MAIL.msh introuvalble"
    
    
        open(unit=10, file='MAIL.msh' ,form='formatted' ,status='old', action='read' , IOSTAT=info)
    
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*)
    
    
    !----------------------------------------------------------
    !   Lecture des contours et de la d�finition de la section
    !----------------------------------------------------------
        read(10,*) noce, nseg, nce
    
        allocate ( ysec(noce), zsec(noce) )
        do i = 1 , noce
            read(10,*)  j, ysec(i), zsec(i)
        enddo
    
        allocate (ind_consec(nseg) , iaseg(nseg), ibseg(nseg) )
        do i = 1 , nseg
            read(10,*)  ind_consec(i) , x1, y1, z1, x2, y2, z2, j1, j2, j3, iaseg(i), ibseg(i)
        
            ibseg(i) = abs(ibseg(i))
        
        enddo
    
    
        allocate ( ics(nce), icseg(100,nce) )
    
        do i = 1 , nce
            read(10,*) j, x1, y1, z1, x2, y2, z2, j1, j2, ics(i), ( icseg(k,i), k=1,ics(i) )
        enddo
    
    
    !----------------------------
    ! Lecture des noeuds :
    !----------------------------
        read(10,*)
        read(10,*)
        read(10,*) n_entities, nno
    

        allocate (  ym(nno), zm(nno) )
        nno = 0
    
        n = 0
        do while ( n .lt. n_entities )
            read(10,*) i1, i2, i3, i4
        
            do i = 1 , i4
                read(10,*)
            enddo
        
            do i = 1 , i4
                read(10,*)  ym(nno+i) , zm(nno+i)
            end do
            nno = nno + i4
        
            n = n + 1
        enddo
    
    
    !---------------------------------------------
    ! Lecture des �l�ments (mailles + bords) :
    !---------------------------------------------
        read(10,*)
        read(10,*)
        read(10,*) n_entities, nelu
    

        if ( nord.eq.2 )  then   
            iel = 6  ;   ibo = 3  
        elseif ( nord.eq.1 )  then
            iel = 3   ;  ibo = 2
        elseif ( nord.eq.5 ) then
            iel = 21  ;  ibo = 6
        endif
        allocate ( ibord(nelu), nodebord(ibo,nelu), ncelem(nelu), nodelem(iel,nelu) )
    
        nelu = 0
        nbord = 0
        nelem = 0
    
        n = 0
        do while ( n .lt. n_entities )
            read(10,*) i1, i2, i3, i4
        
            if ( i3.eq.1 ) then
                do i = 1 , i4
                    nbord = nbord + 1
                    read(10,*) j, nodebord(1,nbord), nodebord(2,nbord)
                enddo
            endif
        
            if ( i3.eq.8 ) then
                do i = 1 , i4
                    nbord = nbord + 1
                    read(10,*) j, nodebord(1,nbord), nodebord(2,nbord), nodebord(3,nbord)
                    ibord(nbord) = i2
                enddo
            endif
        
        
            if ( i3.eq.2 ) then
                do i = 1 , i4
                    nelem = nelem + 1
                    read(10,*) j, (nodelem(k,nelem), k=1,3)
                enddo
            endif
        
            if ( i3.eq.9 ) then
                do i = 1 , i4
                    nelem = nelem + 1
                    read(10,*) j, (nodelem(k,nelem), k=1,6)
                    ncelem(nelem) = i2
                enddo
            endif
        
            nelu = nelu + i4
        
            n = n + 1
        enddo
    
    
        close(10)
    
    
    
    !-------------------------------------------------------------
    !   D�termination des triangles auxquels appartient un bord:
    !-------------------------------------------------------------
        allocate ( nelbord(nbord) ) 
        nelbord = 0
        do j = 1 , nbord
            do i = 1 , nelem
                 ii1=(nodebord(1,j)-nodelem(1,i))*(nodebord(1,j)-nodelem(2,i))*(nodebord(1,j)-nodelem(3,i))
                 ii2=(nodebord(2,j)-nodelem(1,i))*(nodebord(2,j)-nodelem(2,i))*(nodebord(2,j)-nodelem(3,i))
             
                 if ( (ii1==0) .and. (ii2==0) )  nelbord(j) = i

            enddo
        enddo
     
    
    end subroutine
    
    end module