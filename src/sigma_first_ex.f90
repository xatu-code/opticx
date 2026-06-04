module sigma_first_ex
  use constants_math
  use parser_input_file, &
  only:iflag_xatu,nf,e1,e2,eta,nw,broadening_type_text
  use parser_wannier90_tb, &
  only:material_name,norb
  use parser_optics_xatu_dim, &
  only:npointstotal,vcell, &
  norb_ex_cut,nv_ex,nc_ex,nband_ex, &
  e_ex,fk_ex, &
  rkxvector,rkyvector,rkzvector !k-vectors only used for testing
  use ome_ex, &
  only:read_ome_sp_linear !routine
  use sigma_first_sp, &
  only:fill_allocate_sigma_arrays
  
  implicit none

  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_sigma_first_ex()
    implicit none 
    integer iflag_norder
    integer :: ibz,j
    
    !energies and vme in k-mesh and auxiliary arrays (sp)
    dimension ek(npointstotal,nband_ex)
    dimension vme_ex_band(npointstotal,3,nband_ex,nband_ex)
    dimension vme_nband(3,nband_ex,nband_ex)
    dimension e_nband(nband_ex)
    
    !energies and VME (ex)
    dimension e_ex(norb_ex_cut)
    dimension vme_ex(3,norb_ex_cut)

    dimension wp(nw)
    dimension sigma_w_sp(3,3,nw),sigma_w_ex(3,3,nw)
    
    real*8 wp,eta1
    real*8 ek,e_nband
    real*8 e_ex
    complex*16 vme_ex_band,vme_nband 
    complex*16 vme_ex
    complex*16 sigma_w_sp,sigma_w_ex
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) '9. Entering sigma_first_ex'
    !initialize ex arrays
    vme_ex=0.0d0
    call read_ome_ex_linear(vme_ex)

    !allocate conductivity arrays
    call fill_allocate_sigma_arrays(eta1,nw,wp,sigma_w_ex)
    
    write(*,*) '   Evaluating linear conductivity (ex)...'
    !get excitonic frequency tensor
    call get_kubo_intens_ex(vme_ex,nw,wp,eta1,sigma_w_ex)

    !print conductivity tensor
    write(*,*) '   Printing sigma first (ex)...'
    call print_sigma_first_ex(nw,wp,sigma_w_ex)

  end subroutine get_sigma_first_ex
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_ome_ex_linear(vme_ex)
    implicit none
    integer nn,nkaka
    dimension vme_ex(3,norb_ex_cut)
    complex*16 vme_ex

    real*8 :: a1,a2,a3,a4,a5,a6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  open(10,file='ome_linear_ex_'//trim(material_name)//'.omeex') 
    read(10,*)     
    do nn=1,norb_ex_cut
      read(10,*) nkaka,a1,a2,a3,a4,a5,a6
      vme_ex(1,nn)=complex(a1,a2)
      vme_ex(2,nn)=complex(a3,a4)
      vme_ex(3,nn)=complex(a5,a6)
    end do
    close(10)
  end subroutine read_ome_ex_linear
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_kubo_intens_ex(vme_ex,nw,wp,eta1,sigma_w_ex)
  implicit none
  dimension skubo_ex_int(3,norb_ex_cut,norb_ex_cut)
  dimension wp(nw),sigma_w_ex(3,3,nw)
  dimension vme_ex(3,norb_ex_cut)
  
  integer nw
  integer iw,nn,nj,njp
  real*8 delta_n_ex
  real*8 wp,eta1
  
  complex*16 :: vme_ex
  complex*16 :: skubo_ex_int, sigma_w_ex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
    skubo_ex_int=0.0d0
    sigma_w_ex=0.0d0

    do iw=1,nw  
      do nn=1,norb_ex_cut
	      do nj=1,3
	        do njp=1,3
            
            !N integrand
            skubo_ex_int(nj,njp,nn)=pi/(dble(npointstotal)*vcell) &
            *conjg(vme_ex(nj,nn))*vme_ex(njp,nn)/e_ex(nn)   !pick the correct order of operators
            
            !at a given frequency
            !delta function
            if (trim(broadening_type_text) == 'gaussian') then
              delta_n_ex = pi*1.0d0/eta1*1.0d0/sqrt(2.0d0*pi)*&
                exp(-0.5d0/(eta1**2)*(wp(iw)-e_ex(nn))**2)
            else if (trim(broadening_type_text) == 'lorentzian') then
              delta_n_ex = 1.0d0/pi*aimag(1.0d0/(wp(iw)-e_ex(nn)-&
                complex(0.0d0,eta1)))
            else
              delta_n_ex = pi*1.0d0/eta1*1.0d0/sqrt(2.0d0*pi)*&
                exp(-0.5d0/(eta1**2)*(wp(iw)-e_ex(nn))**2)
            end if
            !sigma_w
			      sigma_w_ex(nj,njp,iw)=sigma_w_ex(nj,njp,iw)+skubo_ex_int(nj,njp,nn)*delta_n_ex
          
	        end do
	      end do
      end do
    end do  

  end subroutine get_kubo_intens_ex

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine print_sigma_first_ex(nw,wp,sigma_w_ex)
    implicit none
    integer :: iw
    integer :: nw
    dimension :: wp(nw)
    dimension :: sigma_w_ex(3,3,nw)
    
    real*8 :: wp,feps
    complex*16 :: sigma_w_ex
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    !write frequency dependent conductivity	  
    open(50,file='sigma_first_ex_real_'//trim(material_name)//'.dat')
    open(55,file='sigma_first_ex_imag_'//trim(material_name)//'.dat')

    do iw=1,nw
      feps=1.0d0 !use atomic units
      write(50,*) wp(iw)*27.211385d0, &
        realpart(feps*sigma_w_ex(1,1,iw)), &
        realpart(feps*sigma_w_ex(1,2,iw)), &
        realpart(feps*sigma_w_ex(1,3,iw)), &
        realpart(feps*sigma_w_ex(2,1,iw)), &
        realpart(feps*sigma_w_ex(2,2,iw)), &
        realpart(feps*sigma_w_ex(2,3,iw)), &
        realpart(feps*sigma_w_ex(3,1,iw)), &
        realpart(feps*sigma_w_ex(3,2,iw)), &
        realpart(feps*sigma_w_ex(3,3,iw))
  
      write(55,*) wp(iw)*27.211385d0, &
          aimag(feps*sigma_w_ex(1,1,iw)), &
          aimag(feps*sigma_w_ex(1,2,iw)), &
          aimag(feps*sigma_w_ex(1,3,iw)), &
          aimag(feps*sigma_w_ex(2,1,iw)), &
          aimag(feps*sigma_w_ex(2,2,iw)), &
          aimag(feps*sigma_w_ex(2,3,iw)), &
          aimag(feps*sigma_w_ex(3,1,iw)), &
          aimag(feps*sigma_w_ex(3,2,iw)), &
          aimag(feps*sigma_w_ex(3,3,iw))	
      
    end do

    close(50)
    close(55)

  end subroutine print_sigma_first_ex

end module sigma_first_ex

