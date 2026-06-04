module parser_optics_xatu_dim
  use constants_math
  use parser_wannier90_tb, &
  only:material_name,R !variables
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
  !TO BE ADAPTED FOR 3D	
  ! Here we define some BZ variables, either by reading the output
  ! of Xatu or not
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine get_optics_xatu_dim()
	    implicit none  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,*) '3. Entering parser_optics_xatu_dim'
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !get reciprocal lattice vectors (to be adapted for 3D!)
      call get_reciprocal_vectors()  
      !write(*,*) (nband_index(j), j=1,size(nband_index,dim=1))
      !pause

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
      !write(*,*) (nband_index(j), j=1,nband_ex)
      !write(*,*) 
      !pause

      !allocate grid and exciton arrays
	    allocate (rkxvector(npointstotal))
	    allocate (rkyvector(npointstotal))
      allocate (rkzvector(npointstotal))
	    allocate (e_ex(norb_ex_cut))
	    allocate (fk_ex(norb_ex,norb_ex_cut))
      
      !fill exciton and other arrays: from XATU-output or opticx-input
      if (iflag_xatu .eqv. .true.) then
        call get_exciton_data() !get grid and exciton wavefunctions
        write(*,*) '   Exciton data has been read from XATU output'
      else
        call get_grid()
        !get frid andexciton variables set to zero if XATU interface is not requested
        fk_ex=0.0d0
        e_ex=0.0d0
      end if
      write(*,*) "   Grid and band parameters have been set"
      
      !write(*,*) rkxvector(1),rkyvector(1),rkzvector(1)
      !write(*,*) G(1,1),G(1,2),G(1,3)
      !write(*,*) G(2,1),G(2,2),G(2,3)
      !pause
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
        !write(*,*) ibz,npointstotal,fk_ex(ibz,1),fk_ex(ibz,1),fk_ex(ibz,3)
        !pause
        !isum=0: print the total probability density
        !isum=1: print the probability density of the exciton with index (iv,ic)
        if (isum.eq.1) then
          prob_k=0.0d0
          do iv_s=1,nv_ex
            do ic_s=1,nc_ex
              call get_ex_index_first(nf,nv_ex,nc_ex,iright,ibz,i_ex_nn,ic_s,iv_s)
              prob_k=prob_k+abs(fk_ex(i_ex_nn,nn))**2
            end do
          end do
          write(10,*) rkxvector(ibz),rkyvector(ibz),prob_k
        else
          call get_ex_index_first(nf,nv_ex,nc_ex,iright,ibz,i_ex_nn,ic,iv)
          write(10,*) rkxvector(ibz),rkyvector(ibz),abs(fk_ex(i_ex_nn,nn)) !,nf,nv_ex,nc_ex,iright,ibz,i_ex_nn,ic,iv
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
      !allocatable :: file2open(:)
      !character :: file2open 
      dimension nband_index_aux1(100)
      dimension nband_index_aux2(100)


      integer :: nband_ex_aux !number of valence bands
      integer :: nband_index_aux1,nband_index_aux2 !auxiliary
      integer :: nband_ex_aux1,nband_ex_aux2 !auxiliary
      integer :: iexit
      integer :: i,j,naux,npointstotal_sq
      integer :: hdr1, hdr2, ios
      
      real(8) aux1
      character(len=:), allocatable :: file2open
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		  

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !UPDATE: 2025-07-09. It is no longer necessary to modify .eigval file to read bandlist
      !This part gets 'nband_ex' and 'nband_index(nband_ex)'
      !Note: nband_ex_aux is used as nband_ex, which will be defined later (historic reasons)
      !      but we should have nband_ex=nband_ex_aux 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !save number of valence bands
      file2open=trim(xatu_states_filepath_in)
      open(10,file=file2open)
      read(10,*) 
      do i=1,50 !we will never use more than 50 bands...
        read(10,*) aux1,aux1,aux1,nband_index_aux1(i) 
        if (i.gt.1) then
          do j=1,i-1
            if (nband_index_aux1(i).eq.nband_index_aux1(j)) then
              nband_ex_aux1=i-1 !save number of valence bands
              goto 128
            end if
          end do
        end if
      end do
128   continue
      close(10)	 
      !write(*,*) '   Number of valence bands:',nband_ex1
      !save number of conduction bands
      open(10,file=file2open)
      read(10,*) 
      do i=1,50 !we will never use more than 50 bands...
        read(10,*) aux1,aux1,aux1,naux,nband_index_aux2(i) 
        !skip similar values
        if (nband_ex_aux1.gt.1) then
          do j=1,nband_ex_aux1-1
            read(10,*)
          end do
        end if
        if (i.gt.1) then
          do j=1,i-1
            if (nband_index_aux2(i).eq.nband_index_aux2(j)) then
              nband_ex_aux2=(i-1) !save number of valence bands
              goto 129
            end if
          end do
        end if
      end do
129   continue
      close(10)	 
      !write(*,*) '   Number of total bands:',nband_ex_aux1,nband_ex_aux2
      nband_ex_aux=nband_ex_aux1+nband_ex_aux2
      allocate(nband_index(nband_ex_aux))
      do i=1,nband_ex_aux1
        nband_index(i)=nband_index_aux1(i)-nf+1
      end do
      do i=nband_ex_aux1+1,nband_ex_aux
        nband_index(i)=nband_index_aux2(i-nband_ex_aux1)-nf+1
      end do
      !old, previous lines replace this call
      !reads a line of numbers and stores them in an allocatable array 'narray'
      !call read_line_numbers_int(nband_index,nband_ex) 
      !write(*,*) nband_ex_aux
      !write(*,*) (nband_index(j), j=1,nband_ex_aux)
      !pause   
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !get nk
      file2open=trim(xatu_eigval_filepath_in)
      open(10,file=file2open)
      ! Try old format: first line npointstotal_sq, second line naux
      read(10,*,iostat=ios) hdr1
      if (ios == 0) then
        read(10,*,iostat=ios) hdr2
        if (ios == 0) then
          npointstotal_sq = hdr1
          naux = hdr2
        else
          ! New xatu Result::writeEigenvalues format: first line = exciton basis size (naux)
          naux = hdr1
          ! compute npointstotal from nv_ex*nc_ex (number of pair combinations per k)
          if (nv_ex*nc_ex > 0) then
            npointstotal = naux / (nv_ex*nc_ex)
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
        ! Fallback: try older format read without iostat
        rewind(10)
        read(10,*) npointstotal_sq
        read(10,*) naux
      end if
      close(10)
      
      !get N_BSE=nv_ex*nc_ex*nk**2 variables
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
    !This routine needs to be adapted for 3D
    subroutine get_reciprocal_vectors()
      implicit none
      real(8) cx,cy,cz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    G=0.0d0
      G(1,1)=2.0d0*pi*(-R(1,1)*R(2,2)+R(1,2)*R(2,1))**(-1.0d0) &
             *(-R(2,2))
      G(1,2)=2.0d0*pi*(-R(1,1)*R(2,2)+R(1,2)*R(2,1))**(-1.0d0) &
             *R(2,1)            
      G(2,1)=2.0d0*pi*(-R(2,1)*R(1,2)+R(2,2)*R(1,1))**(-1.0d0) &
             *(-R(1,2))
      G(2,2)=2.0d0*pi*(-R(2,1)*R(1,2)+R(2,2)*R(1,1))**(-1.0d0) &
             *R(1,1)    
	  		
      call crossproduct(R(1,1),R(1,2),0.0d0,R(2,1), &
      R(2,2),0.0d0,cx,cy,cz)         
      vcell=sqrt(cx**2+cy**2+cz**2)
    
    end subroutine
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
      ! Try to read the format produced by xatu Result::writeEigenvalues:
      ! Line 1: exciton basis size (ignored here)
      ! Line 2: number of eigenvalues (nkaka)
      ! Next lines: one eigenvalue per line
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
          ! fallback to older single-line format
          rewind(10)
          read(10,*)
          read(10,*) nkaka,(e_ex(j), j=1,norb_ex_cut)
        end if
      else
        ! fallback to older single-line format
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
      ibz_sum=0 !counter to display percentage
      write(*,*) '   Reading exciton wavefunctions...'
      do ibz=1,norb_ex_cut
        ibz_sum=ibz_sum+1
        !display percentaje
        if (abs(dble(ibz)/dble(norb_ex_cut))*100.0d0-100.0d0 .lt. 5.0d0) then
          call percentage_index(ibz_sum,norb_ex_cut,nkaka)
		      !write(*,*) '   Reading exciton wf:',int(abs(dble(ibz)/dble(norb_ex_cut))*100.0d0),'%'
		      !write(*,*) '   Reading exciton wf:',abs(dble(ibz)/dble(norb_ex_cut))*100.0d0,'%'
        end if
		    read(10,*) (auxr1(j),j=1,2*norb_ex)
        !write(*,*) (auxr1(j),j=1,2*norb_ex)
        !pause
        do j=1,norb_ex    
          jind=2*j-1
          fk_ex(j,ibz)=complex(auxr1(jind),auxr1(jind+1))
        end do
        !pause
	    end do
      
	    close(10)

      !Please I like to work in atomic units!  	
	    e_ex=e_ex/27.211385d0
      rkxvector=rkxvector*0.52917721067121d0 
      rkyvector=rkyvector*0.52917721067121d0 
	    rkzvector=rkzvector*0.52917721067121d0 

	  end subroutine get_exciton_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !TO BE ADAPTED FOR 3D
    subroutine get_grid()
      implicit none

      integer :: k,i1,i2
      real(8) :: step,r1,r2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      step=1.0d0/dble(npointstotal_sq-1)
      
      !atomic units provided that G are in atomic units
      k=1
      do i1=1,npointstotal_sq
		    r1=-0.5d0+dble(i1-1)*step
        do i2=1,npointstotal_sq
          r2=-0.5d0+dble(i2-1)*step
          rkxvector(k)=r1*G(1,1)+r2*G(2,1) 
          rkyvector(k)=r1*G(1,2)+r2*G(2,2) 
          rkzvector(k)=0.0d0
          k=k+1       
        end do
      end do  
                           
    end subroutine get_grid


end module parser_optics_xatu_dim