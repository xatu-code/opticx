module sigma_first_sp
  use constants_math
  use parser_input_file, &
  only:iflag_xatu,nf,e1,e2,eta,nw,broadening_type_text
  use parser_wannier90_tb, &
  only:material_name,norb
  use parser_optics_xatu_dim, &
  only:npointstotal,vcell, &
  norb_ex_cut,nv_ex,nc_ex,nband_ex, &
  rkxvector,rkyvector,rkzvector !k-vectors only used for testing
  use ome_ex, &
  only:read_ome_sp_linear !routine
  use omp_lib  ! <-- OpenMP

  implicit none
  
  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_sigma_first_sp()
    implicit none
    integer :: iflag_norder
    integer :: ibz, j

    ! k-mesh arrays
    real(8)     :: ek(npointstotal, nband_ex)
    complex(8) :: vme_ex_band(npointstotal, 3, nband_ex, nband_ex)

    ! frequency grid and conductivity tensor
    real(8)     :: wp(nw)
    real(8)     :: eta1
    complex(8) :: sigma_w_sp(3, 3, nw)

    ! Per-thread scratch arrays (listed in PRIVATE clause below)
    real(8)     :: e_nband_local(nband_ex)
    complex(8) :: vme_nband_local(3, nband_ex, nband_ex)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) '8. Entering sigma_first_sp'
    !initialize sp arrays
    vme_ex_band=0.0d0
    ek=0.0d0   
    !read optical matrix elements from file
    write(*,*) '   Reading optical matrix elements...'
    call read_ome_sp_linear(iflag_norder,npointstotal,nband_ex,vme_ex_band,ek)
    !allocate conductivity arrays
    call fill_allocate_sigma_arrays(eta1,nw,wp,sigma_w_sp)
    
    write(*,*) '   Evaluating linear conductivity (sp)...'
    !write(*,*) '   Using',omp_get_max_threads(),'OpenMP threads'

    ! ----------------------------------------------------------------
    ! Parallelise over k-points (ibz).
    !
    ! * PRIVATE  : each thread owns its own copy of the scratch arrays
    !              e_nband_local and vme_nband_local so there is no
    !              data race when filling them.
    ! * REDUCTION: sigma_w_sp is accumulated across threads safely.
    ! ----------------------------------------------------------------
    !$OMP PARALLEL DO          &
    !$OMP   DEFAULT(SHARED)    &
    !$OMP   PRIVATE(ibz, e_nband_local, vme_nband_local) &
    !$OMP   REDUCTION(+:sigma_w_sp) &
    !$OMP   SCHEDULE(dynamic)
    do ibz=1,npointstotal
      !fill auxiliary arrays (thread-local copies)
      e_nband_local(:)     = ek(ibz,:)
      vme_nband_local(:,:,:) = vme_ex_band(ibz,:,:,:)
      !fill sigma(w) for a given k point
      call get_kubo_intens_sp(nband_ex,npointstotal,vcell, &
                              e_nband_local,vme_nband_local, &
                              nw,wp,eta1,sigma_w_sp)
    end do
    !$OMP END PARALLEL DO

    !print conductivity tensor
    write(*,*) '   Printing sigma first...'
    call print_sigma_first_sp(nw,wp,sigma_w_sp)

  end subroutine get_sigma_first_sp


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine fill_allocate_sigma_arrays(eta1,nw,wp,sigma_w)
  implicit none

  integer,    intent(in)  :: nw
  integer                 :: i
  real(8),    intent(out) :: wp(nw), eta1
  complex(8), intent(out) :: sigma_w(3,3,nw)
  real(8)                 :: wrange
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  wp=0.0d0
  sigma_w=0.0d0
  wrange=e2-e1
  do i=1,nw
    wp(i)=(e1+wrange/dble(nw)*dble(i-1))/27.211385d0
  end do  
  eta1=eta/27.211385d0 !change units to hartree units
  
  end subroutine fill_allocate_sigma_arrays
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine print_sigma_first_sp(nw,wp,sigma_w_sp)
    implicit none
    integer,    intent(in) :: nw
    integer                :: iw
    real(8),    intent(in) :: wp(nw)
    complex(8),intent(in) :: sigma_w_sp(3,3,nw)
    real(8)                :: feps
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    !write frequency dependent conductivity	  
    open(50,file='sigma_first_sp_real_'//trim(material_name)//'.dat')
    open(55,file='sigma_first_sp_imag_'//trim(material_name)//'.dat')

    do iw=1,nw
      feps=1.0d0 !use atomic units
      write(50,*) wp(iw)*27.211385d0, &
        realpart(feps*sigma_w_sp(1,1,iw)), &
        realpart(feps*sigma_w_sp(1,2,iw)), &
        realpart(feps*sigma_w_sp(1,3,iw)), &
        realpart(feps*sigma_w_sp(2,1,iw)), &
        realpart(feps*sigma_w_sp(2,2,iw)), &
        realpart(feps*sigma_w_sp(2,3,iw)), &
        realpart(feps*sigma_w_sp(3,1,iw)), &
        realpart(feps*sigma_w_sp(3,2,iw)), &
        realpart(feps*sigma_w_sp(3,3,iw))
  
      write(55,*) wp(iw)*27.211385d0,aimag(feps*sigma_w_sp(1,1,iw)), &
          aimag(feps*sigma_w_sp(1,2,iw)), &
          aimag(feps*sigma_w_sp(1,3,iw)), &
          aimag(feps*sigma_w_sp(2,1,iw)), &
          aimag(feps*sigma_w_sp(2,2,iw)), &
          aimag(feps*sigma_w_sp(2,3,iw)), &
          aimag(feps*sigma_w_sp(3,1,iw)), &
          aimag(feps*sigma_w_sp(3,2,iw)), &
          aimag(feps*sigma_w_sp(3,3,iw))	
      
    end do

    close(50)
    close(55)

  end subroutine print_sigma_first_sp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_kubo_intens_sp(nband_ex,npointstotal,vcell,e,vme,nw,wp,eta1,sigma_w_sp)
    implicit none
    integer,    intent(in)    :: nw, nband_ex, npointstotal
    integer                   :: iw, nn, nnp, nj, njp
    real(8),    intent(in)    :: e(nband_ex), wp(nw), eta1, vcell
    complex(8),intent(in)    :: vme(3,nband_ex,nband_ex)
    complex(8),intent(inout) :: sigma_w_sp(3,3,nw)

    real(8)    :: fnn, fnnp, factor1
    complex(8):: delta_nnp, vme_prod
    complex(8):: sigma_local(3,3,nw)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    sigma_local = (0.0d0, 0.0d0)

    ! ----------------------------------------------------------------
    ! Parallelise the band-band loop (nn, nnp).  The iw loop is kept
    ! inside so each (nn,nnp) pair sweeps all frequencies in one shot,
    ! preserving cache locality on wp and sigma_local.
    !
    ! REDUCTION on sigma_local collects partial sums from all threads.
    ! ----------------------------------------------------------------
    !$OMP PARALLEL DO                          &
    !$OMP   DEFAULT(SHARED)                    &
    !$OMP   PRIVATE(nn, nnp, nj, njp,          &
    !$OMP           fnn, fnnp, factor1,        &
    !$OMP           delta_nnp, vme_prod, iw)   &
    !$OMP   REDUCTION(+:sigma_local)           &
    !$OMP   SCHEDULE(static)
    do nn=1,nband_ex
      !fermi distribution
      if (nn.le.nv_ex) then
        fnn=1.0d0
      else
        fnn=0.0d0
      end if

      do nnp=1,nband_ex
        !fermi distribution
        if (nnp.le.nv_ex) then
          fnnp=1.0d0
        else
          fnnp=0.0d0
        end if

        !DECIDE PREFACTOR WITH OCCUPATION
        if (abs(fnn-fnnp).lt.0.1d0) then
          factor1=0.0d0
        else
          factor1=(fnn-fnnp)/(e(nn)-e(nnp))
        end if

        ! Skip entirely if no contribution
        if (factor1 == 0.0d0) cycle

        do iw=1,nw
          ! Broadening: gaussian or lorentzian based on parser input
          if (trim(broadening_type_text) == 'gaussian') then
            delta_nnp = -1.0d0/eta1*1.0d0/sqrt(2.0d0*pi)*&
              exp(-0.5d0/(eta1**2)*(wp(iw)-e(nn)+e(nnp))**2)
          else if (trim(broadening_type_text) == 'lorentzian') then
            delta_nnp = 1.0d0/pi*aimag(1.0d0/(-wp(iw)+e(nn)-e(nnp)+&
              complex(0.0d0,eta1)))
          else
            delta_nnp = -1.0d0/eta1*1.0d0/sqrt(2.0d0*pi)*&
              exp(-0.5d0/(eta1**2)*(wp(iw)-e(nn)+e(nnp))**2)
          end if

          do nj=1,3
            do njp=1,3
              vme_prod = vme(nj,nn,nnp)*vme(njp,nnp,nn)
              sigma_local(nj,njp,iw) = sigma_local(nj,njp,iw) + &
                pi/(dble(npointstotal)*vcell)*factor1*vme_prod*delta_nnp
            end do
          end do

        end do ! iw
      end do   ! nnp
    end do     ! nn
    !$OMP END PARALLEL DO

    ! Accumulate into the caller's array (already protected by the outer
    ! k-point REDUCTION in get_sigma_first_sp)
    sigma_w_sp = sigma_w_sp + sigma_local

  end subroutine get_kubo_intens_sp
end module sigma_first_sp
