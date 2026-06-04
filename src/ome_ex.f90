module ome_ex
  use constants_math
  use parser_input_file, &
  only:nf,e1,e2,eta,nw
  use parser_wannier90_tb, &
  only:material_name,norb
  use parser_optics_xatu_dim, &
  only:npointstotal,vcell, &
  norb_ex,norb_ex_cut,nv_ex,nc_ex,nband_ex,e_ex,fk_ex, &
  get_ex_index_first,print_exciton_wf, & !routines
  rkxvector,rkyvector,rkzvector !k-vectors only used for testing
  use exciton_envelopes, &
  only:fk_ex_der,get_fk_ex_der_k !routines
  implicit none

  !ex-vme
  allocatable :: xme_ex(:,:) !R_{0N}  
  allocatable :: vme_ex(:,:) !V_{0N} 

  !ex-rme
  allocatable :: qme_ex_inter1(:,:,:) !Q_{NN'} (1)
  allocatable :: qme_ex_inter2(:,:,:) !Q_{NN'} (2)
  allocatable :: qme_ex_inter(:,:,:)  !Q_{NN'} 
  allocatable :: yme_ex_inter1(:,:,:) !Y_{NN'} (1)
  allocatable :: yme_ex_inter2(:,:,:) !Y_{NN'} (2)
  allocatable :: yme_ex_inter(:,:,:)  !Y_{NN'} 
  allocatable :: xme_ex_inter(:,:,:)  !R_{NN'}
  allocatable :: vme_ex_inter1(:,:,:) !V_{NN'} (1)
  allocatable :: vme_ex_inter2(:,:,:) !V_{NN'} (2)  
  allocatable :: vme_ex_inter(:,:,:)  !V_{NN'}
  
  complex*16 :: xme_ex,vme_ex
  complex*16 :: qme_ex_inter1,qme_ex_inter2,qme_ex_inter
  complex*16 :: yme_ex_inter1,yme_ex_inter2,yme_ex_inter
  complex*16 :: xme_ex_inter
  complex*16 :: vme_ex_inter1,vme_ex_inter2,vme_ex_inter

  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_ome_ex(iflag_norder)
  implicit none

  integer, intent(in) :: iflag_norder
  integer :: ibz
  integer :: nn, nnp
  integer :: nj

  ! For k-resolved excitonic linear output
  integer :: u_exk
  logical :: do_write_exk
  complex*16, allocatable :: vme_ex_k(:,:)   ! (3, norb_ex_cut)

  ! auxiliary arrays used to evaluate ex-ome
  dimension ek(npointstotal,nband_ex)
  dimension xme_ex_band(npointstotal,3,nband_ex,nband_ex) ! only here! provisional
  dimension vme_ex_band(npointstotal,3,nband_ex,nband_ex)
  dimension berry_eigen_ex_band(npointstotal,3,nband_ex,nband_ex)
  dimension gen_der_ex_band(npointstotal,3,3,nband_ex,nband_ex)
  dimension shift_vector_ex_band(npointstotal,3,3,nband_ex,nband_ex)

  real*8 :: ek
  real*8 :: shift_vector_ex_band
  complex*16 :: xme_ex_band
  complex*16 :: vme_ex_band
  complex*16 :: berry_eigen_ex_band
  complex*16 :: gen_der_ex_band

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,*) '6. Entering ome_ex'

  ! Decide whether to write the k-resolved file:
  ! Only makes sense for linear (iflag_norder == 1)
  do_write_exk = (iflag_norder .eq. 1)

  ! read SP optical matrix elements from file
  write(*,*) '   Reading optical matrix elements (sp)...'
  if (iflag_norder .eq. 1) then
    vme_ex_band = 0.0d0
    ek = 0.0d0
    call read_ome_sp_linear(iflag_norder, npointstotal, nband_ex, vme_ex_band, ek)
  end if

  if (iflag_norder .eq. 2) then
    ek = 0.0d0
    vme_ex_band = 0.0d0
    berry_eigen_ex_band = 0.0d0
    gen_der_ex_band = 0.0d0
    shift_vector_ex_band = 0.0d0
    call read_ome_sp_nonlinear(iflag_norder, npointstotal, nband_ex, berry_eigen_ex_band, &
                               gen_der_ex_band, shift_vector_ex_band, vme_ex_band, ek)

    ! Provisional (19/07/2025): we evaluate here xme_ex_band at a given k-point
    xme_ex_band = 0.0d0
    call get_ome_sp_xme_ex_band(ek, vme_ex_band, xme_ex_band)
  end if

  ! allocate arrays for ex-ome
  ! linear conductivity
  if (iflag_norder .eq. 1 .or. iflag_norder .eq. 2) then
    allocate(vme_ex(3, norb_ex_cut))
    allocate(xme_ex(3, norb_ex_cut))
    vme_ex = 0.0d0
    xme_ex = 0.0d0
  end if

  ! allocate k-resolved buffer + open file (linear only)
  if (do_write_exk) then
    allocate(vme_ex_k(3, norb_ex_cut))
    vme_ex_k = 0.0d0
    u_exk = 77
    call write_ome_ex_linear_kresolved_init(u_exk, material_name, norb_ex_cut)
  end if

  ! second order ones
  if (iflag_norder .eq. 2) then
    allocate(qme_ex_inter1(3, norb_ex_cut, norb_ex_cut))
    allocate(qme_ex_inter2(3, norb_ex_cut, norb_ex_cut))
    allocate(qme_ex_inter(3, norb_ex_cut, norb_ex_cut))
    allocate(yme_ex_inter1(3, norb_ex_cut, norb_ex_cut))
    allocate(yme_ex_inter2(3, norb_ex_cut, norb_ex_cut))
    allocate(yme_ex_inter(3, norb_ex_cut, norb_ex_cut))
    allocate(xme_ex_inter(3, norb_ex_cut, norb_ex_cut))
    allocate(vme_ex_inter1(3, norb_ex_cut, norb_ex_cut))
    allocate(vme_ex_inter2(3, norb_ex_cut, norb_ex_cut))
    allocate(vme_ex_inter(3, norb_ex_cut, norb_ex_cut))

    qme_ex_inter1 = 0.0d0
    qme_ex_inter2 = 0.0d0
    qme_ex_inter  = 0.0d0
    yme_ex_inter1 = 0.0d0
    yme_ex_inter2 = 0.0d0
    yme_ex_inter  = 0.0d0
    xme_ex_inter  = 0.0d0
    vme_ex_inter1 = 0.0d0
    vme_ex_inter2 = 0.0d0
    vme_ex_inter  = 0.0d0

    ! allocate and get derivative of exciton envelope function with respect to k
    allocate(fk_ex_der(3, norb_ex, norb_ex_cut))
    call get_fk_ex_der_k()
  end if

  if (iflag_norder .eq. 1) then
    write(*,*) '   Evaluating excitonic optical matrix elements for linear conductivity...'
  end if
  if (iflag_norder .eq. 2) then
    write(*,*) '   Evaluating excitonic optical matrix elements for nonlinear conductivity...'
  end if

  ! k-space integration of excitonic optical matrix elements
  do ibz = 1, npointstotal
    write(*,*) '   Optical matrix elements (ex): k-point', ibz, '/', npointstotal

    ! Fill V_{0N} and X_{0N} (summed over k) for linear conductivity
    if (iflag_norder .eq. 1 .or. iflag_norder .eq. 2) then
      call get_ome_gs_ex_sum_k(ibz, ek, xme_ex_band, vme_ex_band)
    end if

    ! Also write k-resolved contribution for linear case:
    if (do_write_exk) then
      call get_ome_gs_ex_kresolved(ibz, ek, vme_ex_band, vme_ex_k)
      call write_ome_ex_linear_kresolved_point(u_exk, rkxvector(ibz), rkyvector(ibz), rkzvector(ibz), &
                                               norb_ex_cut, vme_ex_k)
    end if

    if (iflag_norder .eq. 2) then
      call get_ome_inter_ex_sum_k(ibz, ek, xme_ex_band, vme_ex_band, berry_eigen_ex_band)
    end if
  end do

  if (do_write_exk) then
    call write_ome_ex_linear_kresolved_close(u_exk)
    deallocate(vme_ex_k)
    write(*,*) '   k-resolved linear excitonic matrix elements written (omeexk)'
  end if

  if (iflag_norder .eq. 2) then
    ! Sum all terms together for matrix elements (N,N')
    do nn = 1, norb_ex_cut
      do nnp = 1, norb_ex_cut
        do nj = 1, 3
          qme_ex_inter(nj, nn, nnp) = qme_ex_inter1(nj, nn, nnp) + qme_ex_inter2(nj, nn, nnp)
          yme_ex_inter(nj, nn, nnp) = yme_ex_inter1(nj, nn, nnp) + yme_ex_inter2(nj, nn, nnp)
          xme_ex_inter(nj, nn, nnp) = yme_ex_inter(nj, nn, nnp) + qme_ex_inter(nj, nn, nnp)
          vme_ex_inter(nj, nn, nnp) = vme_ex_inter1(nj, nn, nnp) + vme_ex_inter2(nj, nn, nnp)
        end do
      end do
    end do
  end if

  write(*,*) '   Optical matrix elements (ex) have been evaluated'

  ! writing excitonic optical matrix elements (summed over k)
  if (iflag_norder .eq. 1) then
    call write_ome_ex_linear(vme_ex)
  end if
  write(*,*) '   Optical matrix elements (ex, N -> GS) have been written in file'
  write(*,*) '   Optical matrix elements (ex, N -> N^prime) will not be printed in this version'

end subroutine get_ome_ex
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_ome_sp_xme_ex_band(ek,vme_ex_band,xme_ex_band)
    implicit none
    
    !in/out
    dimension ek(npointstotal,nband_ex)
    dimension vme_ex_band(npointstotal,3,nband_ex,nband_ex)
    dimension xme_ex_band(npointstotal,3,nband_ex,nband_ex)
    
    real*8 :: ek
    complex*16 :: vme_ex_band,xme_ex_band

    !here
    integer :: ibz,i,j,nj
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    xme_ex_band=0.0d0
	  do ibz=1,npointstotal
      do nj=1,3
		    do i=1,nband_ex
		      do j=1,nband_ex
			      if (abs(ek(ibz,i)-ek(ibz,j)).lt.1.0d-06) then
			        xme_ex_band(ibz,nj,i,j)=0.0d0
			      else		  
              xme_ex_band(ibz,nj,i,j)=-complex(0.0d0,1.0d0)/(ek(ibz,i)-ek(ibz,j))* &
                                     vme_ex_band(ibz,nj,i,j)
			      end if
			      if (abs(xme_ex_band(ibz,nj,i,j)).gt.20.0d0) then
			        xme_ex_band(ibz,nj,i,j)=0.0d0
			      end if
			    end do
		    end do		 	  
      end do
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end subroutine get_ome_sp_xme_ex_band


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_ome_gs_ex_sum_k(ibz,ek,xme_ex_band,vme_ex_band)
    implicit none
    integer :: ibz

    integer :: nn,ic,iv,iright,nj
    integer :: i_ex_nn
    integer :: i,j
    
    dimension ek(npointstotal,nband_ex)
    dimension vme_ex_band(npointstotal,3,nband_ex,nband_ex)
    dimension xme_ex_band(npointstotal,3,nband_ex,nband_ex)
    
    real*8 :: ek
    complex*16 :: vme_ex_band,xme_ex_band
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    do nn=1,norb_ex_cut		  
      do ic=1,nc_ex
  	    do iv=1,nv_ex
          iright=0
          call get_ex_index_first(nf,nv_ex,nc_ex,iright,ibz,i_ex_nn,ic,iv)			
          do nj=1,3			   
  		      vme_ex(nj,nn)=vme_ex(nj,nn)+fk_ex(i_ex_nn,nn)*vme_ex_band(ibz,nj,iv,nv_ex+ic)	
            xme_ex(nj,nn)=xme_ex(nj,nn)+fk_ex(i_ex_nn,nn)*xme_ex_band(ibz,nj,iv,nv_ex+ic)							  				
          end do			  
        end do	
      end do	
    end do
  end subroutine get_ome_gs_ex_sum_k

  !-----------------------------------------------------------------
  ! Compute k-resolved contribution to V_{0N} for a single k-point
  ! This fills vme_ex_k(:,nn) with the contributions coming from
  ! the given k-point `ibz`, using the same index mapping as
  ! get_ome_gs_ex_sum_k but without summing over k-points.
  !-----------------------------------------------------------------
  subroutine get_ome_gs_ex_kresolved(ibz,ek,vme_ex_band,vme_ex_k)
    implicit none
    integer :: ibz

    integer :: nn,ic,iv,nj
    integer :: i_ex_nn

    dimension ek(npointstotal,nband_ex)
    dimension vme_ex_band(npointstotal,3,nband_ex,nband_ex)
    complex*16 :: vme_ex_band
    complex*16 :: vme_ex_k(3,norb_ex_cut)
    real*8 :: ek

    ! initialize
    vme_ex_k = 0.0d0

    do nn = 1, norb_ex_cut
      do ic = 1, nc_ex
        do iv = 1, nv_ex
          call get_ex_index_first(nf, nv_ex, nc_ex, 0, ibz, i_ex_nn, ic, iv)
          do nj = 1, 3
            vme_ex_k(nj, nn) = vme_ex_k(nj, nn) + fk_ex(i_ex_nn, nn) * vme_ex_band(ibz, nj, iv, nv_ex+ic)
          end do
        end do
      end do
    end do

  end subroutine get_ome_gs_ex_kresolved

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_ome_inter_ex_sum_k(ibz,ek,xme_ex_band,vme_ex_band,berry_eigen_ex_band)
    implicit none

    !in/out
    integer :: ibz

	  dimension ek(npointstotal,nband_ex)
    dimension vme_ex_band(npointstotal,3,nband_ex,nband_ex)
	  dimension berry_eigen_ex_band(npointstotal,3,nband_ex,nband_ex)
	  dimension xme_ex_band(npointstotal,3,nband_ex,nband_ex)

    real*8 :: ek
    complex*16 :: vme_ex_band,berry_eigen_ex_band,xme_ex_band
    complex*16 :: aux2

    !here
    integer :: nn,nnp
    integer :: ic,iv,iright
    integer :: icp,ivp
    integer :: nj
    integer :: i_ex_nn,i_ex_nnp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !N,N' loop. We fill every (N,N') exciton matrix element with the contribution
    !of the given k-point. The N loop is parallelized, while the k-point loop (external to this routine)
    !is in serial. Note that this is different in the SP calculation, where the k-points are parallized.
    !Old tests showed that this is faster than parallelizing the k-point loop.
  
		!!$OMP CRITICAL
		!$OMP PARALLEL DO PRIVATE(nnp,iv,ivp,ic,icp,nj), &  
    !$OMP PRIVATE(i_ex_nn,i_ex_nnp,aux2) 
    do nn=1,norb_ex_cut		
      do nnp=1,norb_ex_cut         
        do ic=1,nc_ex
  	      do iv=1,nv_ex
            iright=0
            call get_ex_index_first(nf,nv_ex,nc_ex,iright,ibz,i_ex_nn,ic,iv)			
            !Dimension loop here is better
            do nj=1,3
              !Here we fill qme_ex_inter with the contribution of every k-point			   
				      qme_ex_inter1(nj,nn,nnp)=qme_ex_inter1(nj,nn,nnp)+ &
                                complex(0.0d0,1.0d0)*conjg(fk_ex(i_ex_nn,nn))*fk_ex_der(nj,i_ex_nn,nnp)             
				      aux2=-complex(0.0d0,1.0d0)*fk_ex(i_ex_nn,nnp) &
	              *(berry_eigen_ex_band(ibz,nj,nv_ex+ic,nv_ex+ic)-berry_eigen_ex_band(ibz,nj,iv,iv))  

	            qme_ex_inter2(nj,nn,nnp)=qme_ex_inter2(nj,nn,nnp)+complex(0.0d0,1.0d0)*conjg(fk_ex(i_ex_nn,nn))*aux2
				      !Extra band loop for vme_ex_inter
              do icp=1,nc_ex
                iright=0
		            call get_ex_index_first(nf,nv_ex,nc_ex,iright,ibz,i_ex_nnp,icp,iv)  
                !Fill vme_ex_inter, term 1			   
		            vme_ex_inter1(nj,nn,nnp)=vme_ex_inter1(nj,nn,nnp) &
				  	                      +conjg(fk_ex(i_ex_nn,nn))*fk_ex(i_ex_nnp,nnp)* &
			  	                        vme_ex_band(ibz,nj,nv_ex+ic,nv_ex+icp)
                if (ic.eq.icp) then
                  continue
                else
                  !Fill yme_ex_inter, term 1. Note that intraband contribution is accounted by qme_ex_inter
                  yme_ex_inter1(nj,nn,nnp)=yme_ex_inter1(nj,nn,nnp) &
					                          +conjg(fk_ex(i_ex_nn,nn))*fk_ex(i_ex_nnp,nnp) &
		                                *xme_ex_band(ibz,nj,nv_ex+ic,nv_ex+icp)
				        end if					          
				      end do

				      do ivp=1,nv_ex
                iright=0
		            call get_ex_index_first(nf,nv_ex,nc_ex,iright,ibz,i_ex_nnp,ic,ivp)   					   
		            !Fill vme_ex_inter, term 2
                vme_ex_inter2(nj,nn,nnp)=vme_ex_inter2(nj,nn,nnp) &
				  	                      -conjg(fk_ex(i_ex_nn,nn))*fk_ex(i_ex_nnp,nnp) &
				  	                      *vme_ex_band(ibz,nj,ivp,iv)
				        if (iv.eq.ivp) then
                  continue
                else
                  !Fill yme_ex_inter, term 2. Note that intraband contribution is accounted by qme_ex_inter
					        yme_ex_inter2(nj,nn,nnp)=yme_ex_inter2(nj,nn,nnp) &
					                          -conjg(fk_ex(i_ex_nn,nn))*fk_ex(i_ex_nnp,nnp) &
		                                *xme_ex_band(ibz,nj,ivp,iv)                 
				        end if
		          end do
				  
            end do
          end do	
        end do	
      
      end do
    end do
    !!$OMP END CRITICAL	
    !$OMP END PARALLEL DO


  end subroutine get_ome_inter_ex_sum_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_ome_ex_linear(vme_ex)
    implicit none
    integer nn,nj
    dimension vme_ex(3,norb_ex_cut)
    complex*16 vme_ex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  open(10,file='ome_linear_ex_'//trim(material_name)//'.omeex') 
    write(10,*) 1    
    do nn=1,norb_ex_cut
      write(10,*) nn,(realpart(vme_ex(nj,nn)),aimag(vme_ex(nj,nn)), nj=1,3)
    end do
    close(10)
  end subroutine write_ome_ex_linear
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine write_ome_ex_linear_kresolved_init(unitno, material_name, norb_ex_cut)
    implicit none
    integer, intent(in) :: unitno, norb_ex_cut
    character(len=*), intent(in) :: material_name

    open(unitno, file='ome_linear_ex_k_'//trim(material_name)//'.omeexk', status='replace')
    ! Simple header:
    ! line 1: tag
    ! line 2: norb_ex_cut
    write(unitno,*) 1
    write(unitno,*) norb_ex_cut
  end subroutine write_ome_ex_linear_kresolved_init


  subroutine write_ome_ex_linear_kresolved_point(unitno, kx, ky, kz, norb_ex_cut, vme_ex_k)
    implicit none
    integer, intent(in) :: unitno, norb_ex_cut
    real*8, intent(in) :: kx, ky, kz
    complex*16, intent(in) :: vme_ex_k(3, norb_ex_cut)

    integer :: nn
    ! For each k-point:
    ! line: kx ky kz
    ! then norb_ex_cut lines:
    !   nn Re(Vx) Im(Vx) Re(Vy) Im(Vy) Re(Vz) Im(Vz)
    write(unitno,*) kx, ky, kz
    do nn = 1, norb_ex_cut
      write(unitno,*) nn, dble(vme_ex_k(1,nn)), dimag(vme_ex_k(1,nn)), &
                          dble(vme_ex_k(2,nn)), dimag(vme_ex_k(2,nn)), &
                          dble(vme_ex_k(3,nn)), dimag(vme_ex_k(3,nn))
    end do
  end subroutine write_ome_ex_linear_kresolved_point


  subroutine write_ome_ex_linear_kresolved_close(unitno)
    implicit none
    integer, intent(in) :: unitno
    close(unitno)
  end subroutine write_ome_ex_linear_kresolved_close


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine read_ome_sp_linear(iflag_norder,npointstotal,nband_ex,vme_ex_band,ek)
    implicit none
    integer iflag_norder
    integer npointstotal,nband_ex
    integer ibz
    integer nj,i,j

    dimension ek(npointstotal,nband_ex)
    dimension vme_ex_band(npointstotal,3,nband_ex,nband_ex)
    
    real*8 :: ek
    real*8 :: a1,a2,a3,b1,b2,b3,b4,b5,b6
    complex*16 vme_ex_band
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  open(10,file='ome_linear_sp_'//trim(material_name)//'.omesp')     
    read(10,*) iflag_norder
    do ibz=1,npointstotal
	    read(10,*) a1,a2,a3,(ek(ibz,j),j=1,nband_ex)	
	    do i=1,nband_ex			
	      do j=1,nband_ex
	        read(10,*) a1,a2,a3,b1,b2,b3,b4,b5,b6
	        vme_ex_band(ibz,1,i,j)=complex(b1,b2)
	        vme_ex_band(ibz,2,i,j)=complex(b3,b4)
	        vme_ex_band(ibz,3,i,j)=complex(b5,b6)
        end do
      end do
    end do

    close(10)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  end subroutine read_ome_sp_linear

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine read_ome_sp_nonlinear(iflag_norder,npointstotal,nband_ex,berry_eigen_ex_band, &
                                   gen_der_ex_band,shift_vector_ex_band,vme_ex_band,ek)
    implicit none
    integer iflag_norder
    integer npointstotal,nband_ex
    integer ibz
    integer nj,i,j

    dimension ek(npointstotal,nband_ex)
    dimension vme_ex_band(npointstotal,3,nband_ex,nband_ex)
    dimension berry_eigen_ex_band(npointstotal,3,nband_ex,nband_ex)
    dimension gen_der_ex_band(npointstotal,3,3,nband_ex,nband_ex)
    dimension shift_vector_ex_band(npointstotal,3,3,nband_ex,nband_ex)

    real*8 :: ek
    real*8 :: a1,a2,a3,b1,b2,b3,b4,b5,b6
    real*8 :: shift_vector_ex_band
    complex*16 vme_ex_band
    complex*16 berry_eigen_ex_band
    complex*16 gen_der_ex_band
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  open(10,file='ome_nonlinear_sp_'//trim(material_name)//'.omesp')     
    read(10,*) iflag_norder
    do ibz=1,npointstotal
      read(10,*) a1,a2,a3,(ek(ibz,j),j=1,nband_ex)	
      do i=1,nband_ex			
        do j=1,nband_ex
	        read(10,*) a1,a2,a3,b1,b2,b3,b4,b5,b6
	        vme_ex_band(ibz,1,i,j)=complex(b1,b2)
	        vme_ex_band(ibz,2,i,j)=complex(b3,b4)
	        vme_ex_band(ibz,3,i,j)=complex(b5,b6)
	        read(10,*) a1,a2,a3,b1,b2,b3,b4,b5,b6
	        berry_eigen_ex_band(ibz,1,i,j)=complex(b1,b2)
	        berry_eigen_ex_band(ibz,2,i,j)=complex(b3,b4)
	        berry_eigen_ex_band(ibz,3,i,j)=complex(b5,b6)
	        do nj=1,3
	          read(10,*) a1,a2,a3,b1,b2,b3
	          shift_vector_ex_band(ibz,nj,1,i,j)=b1
	          shift_vector_ex_band(ibz,nj,2,i,j)=b2
	          shift_vector_ex_band(ibz,nj,3,i,j)=b3
          end do
	        do nj=1,3
	          read(10,*) a1,a2,a3,b1,b2,b3,b4,b5,b6
	          gen_der_ex_band(ibz,nj,1,i,j)=complex(b1,b2)
	          gen_der_ex_band(ibz,nj,2,i,j)=complex(b3,b4)
	          gen_der_ex_band(ibz,nj,3,i,j)=complex(b5,b6)
          end do
	      end do
      end do
    end do	  
    close(10)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  end subroutine read_ome_sp_nonlinear



end module ome_ex