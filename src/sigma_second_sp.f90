module sigma_second_sp
  use constants_math
  use parser_input_file, &
  only:nf,e1,e2,eta,nw,response_text,broadening_type_text
  use parser_wannier90_tb, &
  only:material_name
  use parser_optics_xatu_dim, &
  only:npointstotal,vcell, &
  norb_ex_cut,nv_ex,nc_ex,nband_ex,e_ex,fk_ex, &
  get_ex_index_first,print_exciton_wf, & !routines
  rkxvector,rkyvector,rkzvector !k-vectors only used for testing
  use ome_ex, &
  only:read_ome_sp_nonlinear !routine
  implicit none

  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_sigma_second_sp(nwp,nwq)
    implicit none 
    !in/out
    integer :: nwp,nwq
    
    !here
    integer :: iflag_norder
    integer :: ibz,j
    
    !energies and vme in k-mesh and auxiliary arrays (sp)
    dimension ek(npointstotal,nband_ex)
    dimension vme_ex_band(npointstotal,3,nband_ex,nband_ex)
    dimension berry_eigen_ex_band(npointstotal,3,nband_ex,nband_ex)
    dimension gen_der_ex_band(npointstotal,3,3,nband_ex,nband_ex)
    dimension shift_vector_ex_band(npointstotal,3,3,nband_ex,nband_ex)
    
    !energies and VME (ex)
    dimension wp(nw)
    dimension sigma_w_sp(3,3,nw)

    real*8 ek
    real*8 wp
    real*8 shift_vector_ex_band
    complex*16 vme_ex_band
    complex*16 berry_eigen_ex_band
    complex*16 gen_der_ex_band   
    complex*16 sigma_w_sp 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) '10. Entering sigma_second_sp'

    !read matrix elements from file
    call read_ome_sp_nonlinear(iflag_norder,npointstotal,nband_ex,berry_eigen_ex_band, &
                                   gen_der_ex_band,shift_vector_ex_band,vme_ex_band,ek)
    write(*,*) '    Optical matrix elements (sp) have been read from file'
    
    !compute shift conductivity
    if (nwp.eq.1 .and. nwq.eq.(-1)) then
      call get_sigma_shift_sp(npointstotal,nband_ex,berry_eigen_ex_band, &
                        gen_der_ex_band,shift_vector_ex_band,vme_ex_band,ek)
      !write(*,*) 'The optical response',response_text,'has been evaluated'
    end if
    
    !compute shg susceptibility
    if (nwp.eq.1 .and. nwq.eq.1) then   
      !call get_sigma_shg(npointstotal,nband_ex,berry_eigen_ex_band, &
                        !gen_der_ex_band,shift_vector_ex_band,vme_ex_band,ek)
    end if


  end subroutine get_sigma_second_sp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine get_sigma_shift_sp(npointstotal, nband_ex, berry_eigen_ex_band, &
                               gen_der_ex_band, shift_vector_ex_band, vme_ex_band, ek)
  implicit none

  integer,    intent(in) :: nband_ex, npointstotal
  real*8,     intent(in) :: ek(npointstotal, nband_ex)
  real*8,     intent(in) :: shift_vector_ex_band(npointstotal, 3, 3, nband_ex, nband_ex)
  complex*16, intent(in) :: berry_eigen_ex_band(npointstotal, 3, nband_ex, nband_ex)
  complex*16, intent(in) :: gen_der_ex_band(npointstotal, 3, 3, nband_ex, nband_ex)
  complex*16, intent(in) :: vme_ex_band(npointstotal, 3, nband_ex, nband_ex)

  real*8,     allocatable :: wp(:)
  real*8,     allocatable :: shift_vector_w(:,:,:)
  complex*16, allocatable :: sigma_w_sp(:,:,:,:)
  real*8 :: eta2

  integer :: ibz, i, j, nj, njp

  ! Per-thread private arrays (allocated inside PARALLEL block)
  real*8,     allocatable :: e_nband(:)
  complex*16, allocatable :: vme_nband(:,:,:)
  real*8,     allocatable :: shift_vector_nband(:,:,:,:)
  complex*16, allocatable :: gen_der_nband(:,:,:,:)
  complex*16, allocatable :: vme_der_nband(:,:,:,:)
  real*8,     allocatable :: shift_vector_w_t(:,:,:)
  complex*16, allocatable :: sigma_w_sp_t(:,:,:,:)

  ! ---------------------------------------------------------------------------
  ! Allocate FIRST, then initialize
  ! ---------------------------------------------------------------------------
  allocate(wp(nw))
  allocate(sigma_w_sp(3, 3, 3, nw))
  allocate(shift_vector_w(3, 3, nw))

  call initialize_sigma_second_arrays(nw, wp, eta2, sigma_w_sp)

  shift_vector_w = 0.0d0   ! initialize doesn't touch this one

  write(*,*) '    Evaluating shift conductivity (sp)...'

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP   SHARED(npointstotal, nband_ex, ek, nw, vme_ex_band, &
  !$OMP          shift_vector_ex_band, gen_der_ex_band, &
  !$OMP          wp, eta2, sigma_w_sp, shift_vector_w) &
  !$OMP   PRIVATE(ibz, i, j, nj, njp, &
  !$OMP           e_nband, vme_nband, shift_vector_nband, &
  !$OMP           gen_der_nband, vme_der_nband, &
  !$OMP           shift_vector_w_t, sigma_w_sp_t)

  allocate(e_nband(nband_ex))
  allocate(vme_nband(3, nband_ex, nband_ex))
  allocate(shift_vector_nband(3, 3, nband_ex, nband_ex))
  allocate(gen_der_nband(3, 3, nband_ex, nband_ex))
  allocate(vme_der_nband(3, 3, nband_ex, nband_ex))
  allocate(sigma_w_sp_t(3, 3, 3, nw))
  allocate(shift_vector_w_t(3, 3, nw))

  vme_der_nband    = (0.0d0, 0.0d0)
  sigma_w_sp_t     = (0.0d0, 0.0d0)
  shift_vector_w_t = 0.0d0

  !$OMP DO SCHEDULE(DYNAMIC)
  do ibz = 1, npointstotal

    do i = 1, nband_ex
      e_nband(i) = ek(ibz, i)
      do j = 1, nband_ex
        do nj = 1, 3
          vme_nband(nj, i, j) = vme_ex_band(ibz, nj, i, j)
          do njp = 1, 3
            shift_vector_nband(nj, njp, i, j) = shift_vector_ex_band(ibz, nj, njp, i, j)
            gen_der_nband(nj, njp, i, j)      = gen_der_ex_band(ibz, nj, njp, i, j)
          end do
        end do
      end do
    end do

    call get_shift_intens_sp(nband_ex, nw, e_nband, vme_nband, &
         shift_vector_nband, gen_der_nband, vme_der_nband, &
         shift_vector_w_t, wp, eta2, sigma_w_sp_t)

  end do
  !$OMP END DO

  !$OMP CRITICAL
    sigma_w_sp     = sigma_w_sp     + sigma_w_sp_t
    shift_vector_w = shift_vector_w + shift_vector_w_t
  !$OMP END CRITICAL

  deallocate(e_nband, vme_nband, shift_vector_nband, &
             gen_der_nband, vme_der_nband, &
             sigma_w_sp_t, shift_vector_w_t)

  !$OMP END PARALLEL

  call print_sigma_second_sp(nw, wp, sigma_w_sp, shift_vector_w)
  write(*,*) '    Shift conductivity (sp) has been printed'

  deallocate(wp, sigma_w_sp, shift_vector_w)

end subroutine get_sigma_shift_sp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_shift_intens_sp(nband_ex, nw, e_nband, vme_nband, &
     shift_vector_nband, gen_der_nband, vme_der_nband, &
     shift_vector_w, wp, eta2, sigma_w_sp)
  implicit none
 
  ! ---------------------------------------------------------------------------
  ! Dummy arguments
  ! ---------------------------------------------------------------------------
  integer,    intent(in)    :: nband_ex, nw
  real*8,     intent(in)    :: e_nband(nband_ex)
  complex*16, intent(in)    :: vme_nband(3, nband_ex, nband_ex)
  real*8,     intent(in)    :: shift_vector_nband(3, 3, nband_ex, nband_ex)
  complex*16, intent(in)    :: gen_der_nband(3, 3, nband_ex, nband_ex)
  complex*16, intent(in)    :: vme_der_nband(3, 3, nband_ex, nband_ex)
  real*8,     intent(in)    :: wp(nw)
  real*8,     intent(in)    :: eta2
  real*8,     intent(inout) :: shift_vector_w(3, 3, nw)     ! per-thread accumulator
  complex*16, intent(inout) :: sigma_w_sp(3, 3, 3, nw)      ! per-thread accumulator
 
  ! ---------------------------------------------------------------------------
  ! Local variables  (all stack-allocated, automatically private per call)
  ! ---------------------------------------------------------------------------
  integer    :: iw, nj, njp, njpp, nn, nnp
  real*8     :: factor1, fnn, fnnp, delta_nnp
  complex*16 :: abc(3, nband_ex, nband_ex)
  complex*16 :: shift1, shift2, shift
 
  ! ---------------------------------------------------------------------------
  ! Module-level read-only parameters used below:
  !   nv_ex               -- valence band count
  !   npointstotal        -- total k-points (for normalisation)
  !   vcell               -- unit cell volume
  !   broadening_type_text -- 'gaussian' or 'lorentzian'
  !   response_text        -- 'shift_sumrule', 'shift_shiftvector', etc.
  ! These are read-only globals; reading them from multiple threads is safe.
  ! ---------------------------------------------------------------------------
 
  do iw = 1, nw
    do nn = 1, nband_ex
 
      ! Fermi occupation: valence = 1, conduction = 0
      if (nn .le. nv_ex) then
        fnn = 1.0d0
      else
        fnn = 0.0d0
      end if
 
      do nnp = 1, nband_ex
 
        if (nnp .le. nv_ex) then
          fnnp = 1.0d0
        else
          fnnp = 0.0d0
        end if
 
        factor1 = fnn - fnnp
 
        ! Broadening lineshape
        if (trim(broadening_type_text) == 'gaussian') then
          delta_nnp = 1.0d0/eta2 * 1.0d0/sqrt(2.0d0*pi) * &
            exp(-0.5d0/(eta2**2) * (wp(iw) - e_nband(nn) + e_nband(nnp))**2)
        else if (trim(broadening_type_text) == 'lorentzian') then
          delta_nnp = 1.0d0/pi * aimag(1.0d0 / (wp(iw) - e_nband(nn) + &
            e_nband(nnp) - complex(0.0d0, eta2)))
        else
          ! Default to Gaussian
          delta_nnp = 1.0d0/eta2 * 1.0d0/sqrt(2.0d0*pi) * &
            exp(-0.5d0/(eta2**2) * (wp(iw) - e_nband(nn) + e_nband(nnp))**2)
        end if
 
        do nj = 1, 3
          do njp = 1, 3
 
            ! Shift vector spectral function (2-index, always computed)
            shift_vector_w(nj, njp, iw) = shift_vector_w(nj, njp, iw) + &
              1.0d0/(dble(npointstotal)*vcell) * factor1 * &
              shift_vector_nband(nj, njp, nn, nnp) * delta_nnp
 
            do njpp = 1, 3
 
              if (fnn .eq. fnnp) then
                ! Same occupation: no interband contribution
                shift = (0.0d0, 0.0d0)
              else
 
                if (response_text == 'shift_sumrule') then
                  ! Sum-rule form using generalised derivative
                  shift1 = -complex(0.0d0, 1.0d0) / (e_nband(nn) - e_nband(nnp)) * &
                             vme_nband(njp,  nn, nnp) * gen_der_nband(njpp, nj, nnp, nn)
                  shift2 = -complex(0.0d0, 1.0d0) / (e_nband(nn) - e_nband(nnp)) * &
                             vme_nband(njpp, nn, nnp) * gen_der_nband(njp,  nj, nnp, nn)
                  shift  = -complex(0.0d0, 1.0d0) * (shift1 + shift2)
 
                  sigma_w_sp(nj, njp, njpp, iw) = sigma_w_sp(nj, njp, njpp, iw) + &
                    0.5d0*pi / (dble(npointstotal)*vcell) * factor1 * shift * delta_nnp
 
                end if
 
                if (response_text == 'shift_shiftvector') then
                  ! Nagaosa shift-vector form (TRS)
                  shift = -(shift_vector_nband(nj, njp,  nnp, nn) - &
                             shift_vector_nband(nj, njpp, nn,  nnp)) * &
                           vme_nband(njpp, nn, nnp) * vme_nband(njp, nnp, nn) / &
                           (e_nband(nn) - e_nband(nnp))**2
 
                  sigma_w_sp(nj, njp, njpp, iw) = sigma_w_sp(nj, njp, njpp, iw) + &
                    0.5d0*pi / (dble(npointstotal)*vcell) * factor1 * shift * delta_nnp
 
                end if
 
                if (response_text == 'shift_gender') then
                  ! Placeholder: numerical generalised derivative (Toni's paper)
                  ! TODO: implement
                end if
 
              end if  ! fnn .ne. fnnp
 
            end do  ! njpp
          end do  ! njp
        end do  ! nj
 
      end do  ! nnp
    end do  ! nn
  end do  ! iw
 
end subroutine get_shift_intens_sp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine get_sigma_shg_sp(npointstotal,nband_ex)
    implicit none
    !in/out
    integer :: nband_ex,npointstotal
    !here
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  end subroutine get_sigma_shg_sp


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine initialize_sigma_second_arrays(nw,wp,eta2,sigma_w)
    implicit none
  
    integer :: nw
    integer :: i

    dimension :: wp(nw)
    dimension :: sigma_w(3,3,3,nw)

    real(8) :: wrange,wp,eta2
    complex(8) :: sigma_w
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    wp=0.0d0
    sigma_w=0.0d0
    wrange=e2-e1
    do i=1,nw
      wp(i)=(e1+wrange/dble(nw)*dble(i-1))/27.211385d0
    end do  
    eta2=eta/27.211385d0 !change units to hartree units

  
  end subroutine initialize_sigma_second_arrays

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine print_sigma_second_sp(nw,wp,sigma_w_sp,shift_vector_w)
    implicit none
 
    !in/out
    integer :: iw
    integer :: nw
    dimension :: wp(nw)
    dimension :: sigma_w_sp(3,3,3,nw)
    dimension :: shift_vector_w(3,3,nw)
    
    real*8 :: wp
    complex*16 :: sigma_w_sp
    real*8 :: shift_vector_w
    
    !here
    real*8 :: feps
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    !write frequency dependent conductivity	  
    open(90,file='shift_sp_lengthgauge_'//trim(material_name)//'.dat')
    open(100,file='shift_vector.dat')
    do iw=1,nw
      !feps=-6.623618d-03/(27.21138**2)*1.0d06 !%go from au to (\mu A /V^2)*Angstrongs
      !d=2.6d0 !thickness in angstrongs for MoS2
      !d=3.28d0 !thickness in angstrongs for h-BN
      !feps=feps/(d/0.52917721067121d0) 
      feps=(6.623618d-03)*(1.0d+06)*(27.211386**(-2))*(5.291772d-11)*(1.0d+09) !%go from au to (\mu A /V^2)*nm	
      write(90,*) wp(iw)*27.211385d0,&
                realpart(feps*sigma_w_sp(1,1,1,iw)), &
		        realpart(feps*sigma_w_sp(1,1,2,iw)), &
		        realpart(feps*sigma_w_sp(1,1,3,iw)), &
		        realpart(feps*sigma_w_sp(1,2,1,iw)), &
		  	    realpart(feps*sigma_w_sp(1,2,2,iw)), & 
		  	    realpart(feps*sigma_w_sp(1,2,3,iw)), &
		        realpart(feps*sigma_w_sp(1,3,1,iw)), &
		        realpart(feps*sigma_w_sp(1,3,2,iw)), &
		        realpart(feps*sigma_w_sp(1,3,3,iw)), &
                realpart(feps*sigma_w_sp(2,1,1,iw)), &
		        realpart(feps*sigma_w_sp(2,1,2,iw)), &
		        realpart(feps*sigma_w_sp(2,1,3,iw)), &
		        realpart(feps*sigma_w_sp(2,2,1,iw)), &
		  	    realpart(feps*sigma_w_sp(2,2,2,iw)), & 
		  	    realpart(feps*sigma_w_sp(2,2,3,iw)), &
		        realpart(feps*sigma_w_sp(2,3,1,iw)), &
		        realpart(feps*sigma_w_sp(2,3,2,iw)), &
		        realpart(feps*sigma_w_sp(2,3,3,iw)), &
                realpart(feps*sigma_w_sp(3,1,1,iw)), &
		        realpart(feps*sigma_w_sp(3,1,2,iw)), &
		        realpart(feps*sigma_w_sp(3,1,3,iw)), &
		        realpart(feps*sigma_w_sp(3,2,1,iw)), &
		  	    realpart(feps*sigma_w_sp(3,2,2,iw)), & 
		  	    realpart(feps*sigma_w_sp(3,2,3,iw)), &
		        realpart(feps*sigma_w_sp(3,3,1,iw)), &
		        realpart(feps*sigma_w_sp(3,3,2,iw)), &
		        realpart(feps*sigma_w_sp(3,3,3,iw))
      write(100,*) wp(iw)*27.211385d0,&
            shift_vector_w(1,1,iw),shift_vector_w(1,2,iw), &
            shift_vector_w(2,1,iw),shift_vector_w(2,2,iw)
    end do
    close(90)
    close(95)
    close(100)
  end subroutine print_sigma_second_sp
end module sigma_second_sp

