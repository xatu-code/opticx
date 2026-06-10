module parser_optics_xatu_dim
  use constants_math
  use parser_wannier90_tb, &
    only:material_name,R,nRvec !variables
  use parser_input_file, &
    only:xatu_eigval_filepath_in,xatu_states_filepath_in, & !filepaths
      ndim,npointstotal_sq, & !variables
      iflag_xatu,nf,nband_index,norb_ex_cut, & 
      read_line_numbers_int !subroutine
  implicit none

  integer :: nv_ex,nc_ex
  integer :: npointstotal
  integer :: norb_ex
  integer :: norb_ex_band
  integer :: nband_ex
  integer :: naux
  integer :: j

  real(8) G,vcell
  real(8) rkxvector,rkyvector,rkzvector
  real(8) auxr1
  real(8) e_ex
  complex*16 fk_ex
  
  dimension G(3,3)

  allocatable rkxvector(:)
  allocatable rkyvector(:)
  allocatable rkzvector(:)
  allocatable fk_ex(:,:)
  allocatable e_ex(:)
  allocatable auxr1(:)
	  
  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Here we define some BZ variables, either by reading the output
! of Xatu or not
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_optics_xatu_dim()
  implicit none  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,*) '3. Entering parser_optics_xatu_dim'
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Reminder:norb_ex_cut is given as input
  !get band and grid dimensions: from XATU-output or opticx-input
  if (iflag_xatu .eqv. .true.) then
    call get_exciton_dim()
    !nband_ex=nc_ex+nv_ex
  else
    !number of total k-points
    if (ndim==1) npointstotal=npointstotal_sq
    if (ndim==2) npointstotal=npointstotal_sq**2
    if (ndim==3) npointstotal=npointstotal_sq**3
    norb_ex_band=nv_ex*nc_ex
    norb_ex=norb_ex_band*npointstotal   
  end if
  !calculate nv_ex and nc_ex from the array of bands
  nband_ex=size(nband_index,dim=1)
  nv_ex=0
  do j=1,nband_ex
    if (nband_index(j).le.0) nv_ex=nv_ex+1
  end do
  nc_ex=nband_ex-nv_ex

  !change syntax for band counting
  !XATU: ...-1 0 1 2... to explicit band count
  !opticx: ...nf-1,nf,nf+1...
  nband_index(:)=nband_index(:)+nf 

  !allocate grid and exciton arrays
  allocate (rkxvector(npointstotal))
  allocate (rkyvector(npointstotal))
  allocate (rkzvector(npointstotal))
  allocate (e_ex(norb_ex_cut))
  allocate (fk_ex(norb_ex,norb_ex_cut))
  
  !get reciprocal lattice vectors
  call get_reciprocal_vectors()  
  
  !fill exciton and other arrays: from XATU-output or opticx-input
  if (iflag_xatu .eqv. .true.) then
    call get_exciton_data() !get grid and exciton wavefunctions
    write(*,*) '   Exciton data has been read from XATU output'
  else
    call get_grid()
    !get grid and exciton variables set to zero if XATU interface is not requested
    fk_ex=0.0d0
    e_ex=0.0d0
  end if
  write(*,*) "   Grid and band parameters have been set"
  
!   write(*,*) rkxvector(1),rkyvector(1),rkzvector(1)
!   write(*,*) rkxvector(2),rkyvector(2),rkzvector(2)
!   write(*,*) rkxvector(2),rkyvector(2),rkzvector(3)
!   write(*,*) G(1,1),G(1,2),G(1,3)
!   write(*,*) G(2,1),G(2,2),G(2,3)
!   write(*,*) G(3,1),G(3,2),G(3,3)
end subroutine get_optics_xatu_dim   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
! This subroutine prints a part of the exciton wavefunction or 
! the total exciton probability density
subroutine print_exciton_wf(isum,iv,ic,nn)
  implicit none
  integer isum,iv,ic,nn
  integer iv_s,ic_s
  integer iright,i_ex_nn
  integer ibz
  real*8 prob_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  iright=0
  open(10,file='exciton_wf.dat')
  do ibz=1,npointstotal     
    if (isum.eq.1) then
      prob_k=0.0d0
      do iv_s=1,nv_ex
        do ic_s=1,nc_ex
          call get_ex_index_first(nf,nv_ex,nc_ex,iright,ibz,i_ex_nn,ic_s,iv_s)
          prob_k=prob_k+abs(fk_ex(i_ex_nn,nn))**2
        end do
      end do
      write(10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz),prob_k
    else
      call get_ex_index_first(nf,nv_ex,nc_ex,iright,ibz,i_ex_nn,ic,iv)
      write(10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz),abs(fk_ex(i_ex_nn,nn))
    end if
  end do
  close(10)
end subroutine print_exciton_wf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!this subroutine has not been updated in 2025
subroutine get_ex_index_first(nf,nv_ex,nc_ex,iright,ibz,i_ex,ic,iv)
  implicit none
  integer nf,nv_ex,nc_ex,ibz,i_ex,ic,iv
  integer iright 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !get band indeces (respect the fermi level) from A_cv index
    if (iright.eq.1) then
      iv=(nf-nv_ex)+i_ex-int((i_ex-1)/nv_ex)*nv_ex-nf
      ic=(nf+1)+int((i_ex-1)/nv_ex)-nf	
    end if
  !get A_cv index from band indeces
    if (iright.ne.1) then
    i_ex=nc_ex*nv_ex*(ibz-1)+nv_ex*(ic-1)+iv
    end if	  	  
end subroutine get_ex_index_first

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_exciton_dim()
  implicit none
  dimension nband_index_aux1(1000)
  dimension nband_index_aux2(1000)

  integer :: nband_ex_aux
  integer :: nband_index_aux1,nband_index_aux2
  integer :: nband_ex_aux1,nband_ex_aux2
  integer :: iexit
  integer :: i,j,naux,npointstotal_sq
  integer :: hdr1, hdr2, ios
  integer :: nv_ex_local, nc_ex_local, norb_ex_band_local

  real(8) aux1
  character(len=:), allocatable :: file2open
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !This part gets 'nband_ex' and 'nband_index(nband_ex)'
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !save number of valence bands
  file2open=trim(xatu_states_filepath_in)
  open(10,file=file2open)
  read(10,*) 
  do i=1,500
    read(10,*) aux1,aux1,aux1,nband_index_aux1(i) 
    if (i.gt.1) then
      do j=1,i-1
        if (nband_index_aux1(i).eq.nband_index_aux1(j)) then
          nband_ex_aux1=i-1
          goto 128
        end if
      end do
    end if
  end do
  128   continue
  close(10)	 

  !save number of conduction bands
  open(10,file=file2open)
  read(10,*) 
  do i=1,500
    read(10,*) aux1,aux1,aux1,naux,nband_index_aux2(i) 
    if (nband_ex_aux1.gt.1) then
      do j=1,nband_ex_aux1-1
        read(10,*)
      end do
    end if
    if (i.gt.1) then
      do j=1,i-1
        if (nband_index_aux2(i).eq.nband_index_aux2(j)) then
          nband_ex_aux2=(i-1)
          goto 129
        end if
      end do
    end if
  end do
  129   continue
  close(10)	 

  nband_ex_aux=nband_ex_aux1+nband_ex_aux2
  allocate(nband_index(nband_ex_aux))
  do i=1,nband_ex_aux1
    nband_index(i)=nband_index_aux1(i)-nf+1
  end do
  do i=nband_ex_aux1+1,nband_ex_aux
    nband_index(i)=nband_index_aux2(i-nband_ex_aux1)-nf+1
  end do

  nv_ex_local=0
  do i=1,nband_ex_aux
    if (nband_index(i).le.0) nv_ex_local=nv_ex_local+1
  end do
  nc_ex_local=nband_ex_aux-nv_ex_local
  norb_ex_band_local=nv_ex_local*nc_ex_local

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !get nk
  file2open=trim(xatu_eigval_filepath_in)
  open(10,file=file2open)
  read(10,*,iostat=ios) hdr1
  if (ios == 0) then
    read(10,*,iostat=ios) hdr2
    if (ios == 0) then
      npointstotal_sq = hdr1
      naux = hdr2
    else
      ! New xatu Result::writeEigenvalues format: first line = exciton basis size (naux)
      naux = hdr1
      if (norb_ex_band_local > 0) then
        npointstotal = naux / norb_ex_band_local
      else
        npointstotal = 0
      end if
      ! derive npointstotal_sq based on ndim
      if (ndim == 1) then
        npointstotal_sq = npointstotal
      else if (ndim == 2) then
        npointstotal_sq = int(sqrt(dble(npointstotal)) + 0.5d0)
      else if (ndim == 3) then
        npointstotal_sq = int((dble(npointstotal))**(1.0d0/3.0d0) + 0.5d0)
      end if
    end if
  else
    rewind(10)
    read(10,*) npointstotal_sq
    read(10,*) naux
  end if
  close(10)

  !get N_BSE=nv_ex*nc_ex*nk variables
  if (npointstotal == 0) then
    if (ndim==1) npointstotal=npointstotal_sq
    if (ndim==2) npointstotal=npointstotal_sq**2
    if (ndim==3) npointstotal=npointstotal_sq**3
  end if
  if (npointstotal > 0) then
    norb_ex_band = int(naux / npointstotal)
  else
    norb_ex_band = 0
  end if
  norb_ex = norb_ex_band * npointstotal

end subroutine get_exciton_dim
	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_reciprocal_vectors()
  implicit none
  real(8) cx,cy,cz
  real*8 det      
  logical :: active_x, active_y, active_z
  active_x = (NORM2(real(nRvec(:,1))) /= 0.0d0)
  active_y = (NORM2(real(nRvec(:,2))) /= 0.0d0)
  active_z = (NORM2(real(nRvec(:,3))) /= 0.0d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  G=0.0d0

  ! 1D
  if ( ndim == 1 ) then
        
    if (active_z) then
      G(3,3)=2.0d0*pi*(R(3,3))**(-1.0d0)
      vcell=sqrt(R(3,1)**2+R(3,2)**2+R(3,3)**2)
      
    elseif (active_y) then
      G(2,2)=2.0d0*pi*(R(2,2))**(-1.0d0)
      vcell=sqrt(R(2,1)**2+R(2,2)**2+R(2,3)**2)
      
    else
      G(1,1)=2.0d0*pi*(R(1,1))**(-1.0d0)
      vcell=sqrt(R(1,1)**2+R(1,2)**2+R(1,3)**2)
    endif
  
  ! 2D
  elseif ( ndim == 2 ) then
      !write(*,*) "2D, got here"
    if (active_y .and. active_z) then
      !write(*,*) "2D, x zero"
      G(2,2)=2.0d0*pi*(-R(2,3)*R(3,2)+R(2,2)*R(3,3))**(-1.0d0) &
          *(R(3,3))
      G(2,3)=2.0d0*pi*(-R(2,2)*R(3,3)+R(2,3)*R(3,2))**(-1.0d0) &
          *(R(3,2))
      G(3,2)=2.0d0*pi*(-R(2,2)*R(3,3)+R(2,3)*R(3,2))**(-1.0d0) &
          *(R(2,3))
      G(3,3)=2.0d0*pi*(-R(2,3)*R(3,2)+R(2,2)*R(3,3))**(-1.0d0) &
          *(R(2,2))
      call crossproduct(R(3,1),R(3,2),R(3,3),R(2,1),R(2,2),R(2,3),cx,cy,cz)
      
    elseif (active_x .and. active_z) then
      !write(*,*) "2D, y zero"
      G(1,1)=2.0d0*pi*(-R(1,3)*R(3,1)+R(1,1)*R(3,3))**(-1.0d0) &
          *(R(3,3))
      G(1,3)=2.0d0*pi*(-R(1,1)*R(3,3)+R(1,3)*R(3,1))**(-1.0d0) &
          *(R(3,1))
      G(3,1)=2.0d0*pi*(-R(1,1)*R(3,3)+R(1,3)*R(3,1))**(-1.0d0) &
          *(R(1,3))
      G(3,3)=2.0d0*pi*(-R(1,3)*R(3,1)+R(1,1)*R(3,3))**(-1.0d0) &
          *(R(1,1))
      call crossproduct(R(1,1),R(1,2),R(1,3),R(3,1),R(3,2),R(3,3),cx,cy,cz)
      
    else
      !write(*,*) "2D, z zero"
      G(1,1)=2.0d0*pi*(-R(1,1)*R(2,2)+R(1,2)*R(2,1))**(-1.0d0) &
          *(-R(2,2))
      G(1,2)=2.0d0*pi*(-R(1,1)*R(2,2)+R(1,2)*R(2,1))**(-1.0d0) &
          *R(2,1)            
      G(2,1)=2.0d0*pi*(-R(2,1)*R(1,2)+R(2,2)*R(1,1))**(-1.0d0) &
          *(-R(1,2))
      G(2,2)=2.0d0*pi*(-R(2,1)*R(1,2)+R(2,2)*R(1,1))**(-1.0d0) &
          *R(1,1) 
      call crossproduct(R(1,1),R(1,2),R(1,3),R(2,1),R(2,2),R(2,3),cx,cy,cz)
    endif
    
    vcell=sqrt(cx**2+cy**2+cz**2)
    
  ! 3D
  elseif ( ndim == 3 ) then

    ! determinant of 3x3 matrix
    det = R(1,3)*R(2,2)*R(3,1) - R(1,2)*R(2,3)*R(3,1) - R(1,3)*R(2,1)*R(3,2) &
        + R(1,1)*R(2,3)*R(3,2) + R(1,2)*R(2,1)*R(3,3) - R(1,1)*R(2,2)*R(3,3)


    G(1,1)=2.0d0*pi*(det)**(-1.0d0)*( R(2,3)*R(3,2) - R(2,2)*R(3,3) )
    G(1,2)=2.0d0*pi*(det)**(-1.0d0)*(-R(2,3)*R(3,1) + R(2,1)*R(3,3) )
    G(1,3)=2.0d0*pi*(det)**(-1.0d0)*( R(2,2)*R(3,1) - R(2,1)*R(3,2) )
    G(2,1)=2.0d0*pi*(det)**(-1.0d0)*(-R(1,3)*R(3,2) + R(1,2)*R(3,3) )
    G(2,2)=2.0d0*pi*(det)**(-1.0d0)*( R(1,3)*R(3,1) - R(1,1)*R(3,3) )
    G(2,3)=2.0d0*pi*(det)**(-1.0d0)*(-R(1,2)*R(3,1) + R(1,1)*R(3,2) )
    G(3,1)=2.0d0*pi*(det)**(-1.0d0)*( R(1,3)*R(2,2) - R(1,2)*R(2,3) )
    G(3,2)=2.0d0*pi*(det)**(-1.0d0)*(-R(1,3)*R(2,1) + R(1,1)*R(2,3) )
    G(3,3)=2.0d0*pi*(det)**(-1.0d0)*( R(1,2)*R(2,1) - R(1,1)*R(2,2) )
  
    call crossproduct(R(2,1),R(2,2),R(2,3),R(3,1),R(3,2),R(3,3),cx,cy,cz)
    ! norm of triple product
    vcell=abs(R(1,1)*cx+R(1,2)*cy+R(1,3)*cz)
    
  endif
    
end subroutine get_reciprocal_vectors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
subroutine get_exciton_data()
  implicit none
  integer j,nkaka
  integer ib,ibz,ibz_sum,jind  
  real(8) auxr1

  dimension auxr1(2*norb_ex)
  character(len=:), allocatable :: file2open
  integer :: header1, ios
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		  
  !get energies
  file2open=trim(xatu_eigval_filepath_in)
  open(10,file=file2open) 
  read(10,*,iostat=ios) header1
  read(10,*,iostat=ios)
  if (ios == 0) then
    read(10,*,iostat=ios) nkaka
    if (ios == 0) then
      do j=1,norb_ex_cut
        read(10,*,iostat=ios) e_ex(j)
        if (ios /= 0) exit
      end do
    else
      rewind(10)
      read(10,*)
      read(10,*) nkaka,(e_ex(j), j=1,norb_ex_cut)
    end if
  else
    rewind(10)
    read(10,*)
    read(10,*) nkaka,(e_ex(j), j=1,norb_ex_cut)
  end if
  close(10)

  file2open=trim(xatu_states_filepath_in)
    open(10,file=file2open)	  	  
    read(10,*) 

    !reading k-mesh
    do ibz=1,npointstotal
        read(10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz)
        do ib=1,norb_ex_band-1
          read(10,*) 
      end do
    end do
  
  !reading exciton-wf  
  ibz_sum=0
  write(*,*) '   Reading exciton wavefunctions...'
  do ibz=1,norb_ex_cut
    ibz_sum=ibz_sum+1
    if (abs(dble(ibz)/dble(norb_ex_cut))*100.0d0-100.0d0 .lt. 5.0d0) then
      call percentage_index(ibz_sum,norb_ex_cut,nkaka)
    end if
    read(10,*) (auxr1(j),j=1,2*norb_ex)
    do j=1,norb_ex    
      jind=2*j-1
      fk_ex(j,ibz)=complex(auxr1(jind),auxr1(jind+1))
    end do
  end do
  close(10)

  !Please I like to work in atomic units!  	
  e_ex=e_ex/27.211385d0
  rkxvector=rkxvector*0.52917721067121d0 
  rkyvector=rkyvector*0.52917721067121d0 
  rkzvector=rkzvector*0.52917721067121d0 

end subroutine get_exciton_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_grid()
  implicit none
 
  integer :: k,i1,i2,i3
  real(8) :: step,r1,r2,r3
  logical :: active_x, active_y, active_z
  active_x = (NORM2(real(nRvec(:,1))) /= 0.0d0)
  active_y = (NORM2(real(nRvec(:,2))) /= 0.0d0)
  active_z = (NORM2(real(nRvec(:,3))) /= 0.0d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  step=1.0d0/dble(npointstotal_sq-1)

  ! 1D
  if ( ndim == 1 ) then
 
    if ( active_z ) then
      k=1
      do i1=1,npointstotal_sq
        r1=-0.5d0+dble(i1-1)*step
        rkxvector(k)=0.0d0
        rkyvector(k)=0.0d0
        rkzvector(k)=r1*G(3,3)  
        k=k+1       
      end do
      
    elseif ( active_y ) then
      k=1
      do i1=1,npointstotal_sq
        r1=-0.5d0+dble(i1-1)*step
        rkxvector(k)=0.0d0
        rkyvector(k)=r1*G(2,2)  
        rkzvector(k)=0.0d0
        k=k+1       
      end do
      
    else
      k=1
      do i1=1,npointstotal_sq
        r1=-0.5d0+dble(i1-1)*step
        rkxvector(k)=r1*G(1,1) 
        rkyvector(k)=0.0d0 
        rkzvector(k)=0.0d0
        k=k+1       
      end do 
    endif
  
  ! 2D
  ! Loop order convention: the interpolation routine get_fk_ex_k_interp
  ! uses rk1 as the FAST (unit-stride) index and rk2 as the SLOW index:
  !   ibz = nblock + int((rk2+0.5)/slice)*nside
  ! so the grid must be laid out with rk1 varying in the inner loop
  ! and rk2 in the outer loop.
  elseif ( ndim == 2 ) then
            
    if ( active_y .and. active_z ) then
      ! active axes: rk2 (fast) and rk3 (slow)
      k=1
      do i2=1,npointstotal_sq        ! rk3 slow
        r2=-0.5d0+dble(i2-1)*step
        do i1=1,npointstotal_sq      ! rk2 fast
          r1=-0.5d0+dble(i1-1)*step
          rkxvector(k)=0.0d0
          rkyvector(k)=r1*G(2,2)+r2*G(3,2) 
          rkzvector(k)=r1*G(2,3)+r2*G(3,3) 
          k=k+1       
        end do
      end do
      
    elseif ( active_x .and. active_z ) then
      ! active axes: rk1 (fast) and rk3 (slow)
      k=1
      do i2=1,npointstotal_sq        ! rk3 slow
        r2=-0.5d0+dble(i2-1)*step
        do i1=1,npointstotal_sq      ! rk1 fast
          r1=-0.5d0+dble(i1-1)*step
          rkxvector(k)=r1*G(1,1)+r2*G(3,1) 
          rkyvector(k)=0.0d0
          rkzvector(k)=r1*G(1,3)+r2*G(3,3) 
          k=k+1       
        end do
      end do
      
    else
      ! active axes: rk1 (fast) and rk2 (slow)
      k=1
      do i2=1,npointstotal_sq        ! rk2 slow
        r2=-0.5d0+dble(i2-1)*step
        do i1=1,npointstotal_sq      ! rk1 fast
          r1=-0.5d0+dble(i1-1)*step
          rkxvector(k)=r1*G(1,1)+r2*G(2,1) 
          rkyvector(k)=r1*G(1,2)+r2*G(2,2) 
          rkzvector(k)=0.0d0
          k=k+1       
        end do
      end do 
    endif
        
  ! 3D
  ! Interpolation convention:
  !   nblock = int((rk1+0.5)/slice)
  !          + int((rk2+0.5)/slice)*nside
  !          + int((rk3+0.5)/slice)*nside^2 + 1
  ! so rk1 is fastest, rk2 middle, rk3 slowest.
  ! Loop order: i3 outermost, i2 middle, i1 innermost.
  else  
    k=1
    do i3=1,npointstotal_sq          ! rk3 slow
      r3=-0.5d0+dble(i3-1)*step
      do i2=1,npointstotal_sq        ! rk2 middle
        r2=-0.5d0+dble(i2-1)*step
        do i1=1,npointstotal_sq      ! rk1 fast
          r1=-0.5d0+dble(i1-1)*step
          rkxvector(k)=r1*G(1,1)+r2*G(2,1)+r3*G(3,1)
          rkyvector(k)=r1*G(1,2)+r2*G(2,2)+r3*G(3,2)
          rkzvector(k)=r1*G(1,3)+r2*G(2,3)+r3*G(3,3)
          k=k+1      
        end do
      end do
    end do 
    
  endif
  
end subroutine get_grid
 
 
end module parser_optics_xatu_dim
 
