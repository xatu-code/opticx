module bands
  use constants_math
  use parser_wannier90_tb, &
  only:material_name,nR,nRvec,norb,R,shop,hhop
  use parser_optics_xatu_dim, &
  only:G
  implicit none
  
  allocatable rkxvector_path(:),rkyvector_path(:),rkzvector_path(:)
  allocatable rklengthvector_path(:)
  allocatable b1(:),b2(:) !only a path in a 2D BZ is allowed
  
  real*8 rkxvector_path,rkyvector_path,rkzvector_path,rklengthvector_path
  real*8 b1,b2
  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_energy_bands()
	  implicit none
    integer npointstotal_path,npaths
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  	  
	  npointstotal_path=200
	  !npointstotal_path=200
    npaths=2 !hBn, hexagonal
    !npaths=2         !GeS, square
    allocate (rkxvector_path(npointstotal_path))
    allocate (rkyvector_path(npointstotal_path))
    allocate (rkzvector_path(npointstotal_path))
    allocate (rklengthvector_path(npointstotal_path))
    allocate (b1(npaths+1))
    allocate (b2(npaths+1))
  	b1(npaths+1)=0.0d0
    b2(npaths+1)=0.0d0

    !Default path: from 0.5*\bold{G}_1 to 0.5*\bold{G}_2
    b1(1)=0.5d0
    b2(1)=0.0d0
    b1(2)=0.0d0
    b2(2)=0.0d0
    b1(3)=0.0d0
    b2(3)=0.5d0

    !introduce path manually for now
		!hBN_Pedersen, hexagonal BZ
    !b1(1)=0.0d0
		!b2(1)=0.0d0
		!b1(2)=0.5d0
		!b2(2)=0.5d0
		!b1(3)=2.0d0/3.0d0
		!b2(3)=1.0d0/3.0d0		
		!b1(4)=0.0d0
		!b2(4)=0.0d0
    !MoS2 wannier maria, hexagonal bz
  	!b1(1)=0.0d0
    !b2(1)=0.0d0
    !b1(2)=0.5d0
    !b2(2)=0.0d0
    !b1(3)=1.0d0/3.0d0
    !b2(3)=1.0d0/3.0d0
    !b1(4)=0.0d0
    !b2(4)=0.0d0 	

    !get reciprocal space paths  
    call get_path(npointstotal_path,npaths)
    call get_eigenenergies(npointstotal_path)
    !pause
	end subroutine get_energy_bands
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_path(npointstotal_path,npaths)
    implicit none
    integer npointstotal_path,npaths
    integer npp
    integer index  
    integer i,ii,j
    real*8 r1x,r1y,dx,dy
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
    !reciprocal lattice vectors    
    npp=npointstotal_path/npaths !this has to be divisible by npaths
    do i=1,npaths             
      r1x=(b1(i+1)-b1(i))*G(1,1)+(b2(i+1)-b2(i))*G(2,1)
      r1y=(b1(i+1)-b1(i))*G(1,2)+(b2(i+1)-b2(i))*G(2,2)
      dx=r1x/dble(npp)
      dy=r1y/dble(npp)
      do j=1,npp
        index=(i-1)*npp+j
        rkxvector_path(index)=b1(i)*G(1,1)+b2(i)*G(2,1)+dble(j)*dx
        rkyvector_path(index)=b1(i)*G(1,2)+b2(i)*G(2,2)+dble(j)*dy
      end do     
    end do  

    !accumulated k-path
    do j=1,npointstotal_path
      rklengthvector_path(j)=0.0d0
      do ii=1,j
        if (ii.gt.1) then
          rklengthvector_path(j)=rklengthvector_path(j)+ &
                      sqrt((rkxvector_path(ii)-rkxvector_path(ii-1))**2+ &
                      (rkyvector_path(ii)-rkyvector_path(ii-1))**2)
        !else  
          !rklengthvector_path(j)=sqrt(rkxvector_path(1)**2+rkyvector_path(1)**2)
        end if
      end do
    end do
  end subroutine get_path
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_eigenenergies(npointstotal_path)
    implicit none
    
    integer npointstotal_path
    integer ialpha
    integer ialphap
    integer iRp
    integer ibz
    integer j

    dimension skernel(norb,norb)
    dimension hkernel(norb,norb)
    
    dimension e(norb)
    
    real*8 e
    real*8 Rx,Ry,Rz
    real*8 :: rkx,rky,rkz

    complex*16 skernel,hkernel
    complex*16 phase,factor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  
    open(10,file='bands_'//trim(material_name)//'.dat')	   	  
	  do ibz=1,npointstotal_path
        !write(*,*) 'point:',ibz,npointstotal_path        
        rkx=rkxvector_path(ibz)
        rky=rkyvector_path(ibz)
        rkz=rkzvector_path(ibz)
        hkernel=0.0d0
        skernel=0.0d0
        e=0.0d0
        do ialpha=1,norb
          do ialphap=1,ialpha		  
            do iRp=1,nR
              Rx=dble(nRvec(iRp,1))*R(1,1)+dble(nRvec(iRp,2))*R(2,1)+dble(nRvec(iRp,3))*R(3,1)
              Ry=dble(nRvec(iRp,1))*R(1,2)+dble(nRvec(iRp,2))*R(2,2)+dble(nRvec(iRp,3))*R(3,2)
              Rz=dble(nRvec(iRp,1))*R(1,3)+dble(nRvec(iRp,2))*R(2,3)+dble(nRvec(iRp,3))*R(3,3)
              phase=complex(0.0d0,rkx*Rx+rky*Ry+rkz*Rz)
              factor=exp(phase)                       
              hkernel(ialpha,ialphap)=hkernel(ialpha,ialphap)+ &
              factor*hhop(iRp,ialpha,ialphap)                
              skernel(ialpha,ialphap)=skernel(ialpha,ialphap)+ &
              factor*shop(iRp,ialpha,ialphap)          
            end do
            hkernel(ialphap,ialpha)=conjg(hkernel(ialpha,ialphap))
            skernel(ialphap,ialpha)=conjg(skernel(ialpha,ialphap))    
          end do
        end do 
        call diagoz(norb,e,hkernel)  

        write(10,*) rkx,rky,rkz,rklengthvector_path(ibz),(e(j)*27.211385d0,j=1,norb)
        !if (rkx.eq.0.0d0 .and. rky.eq.0.0d0) then
          !write(*,*) 'eigenenergies:',e(60)*27.211385d0,e(61)*27.211385d0,(e(61)-e(60))*27.211385d0
        !end if
	  end do
    close(10)


  end subroutine get_eigenenergies
      
end module bands