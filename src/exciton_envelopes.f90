module exciton_envelopes
  use constants_math
  use parser_optics_xatu_dim, &
	only:G,npointstotal,rkxvector,rkyvector,rkzvector, &
       norb_ex,norb_ex_cut,norb_ex_band,fk_ex
  use parser_input_file, &
	only:ndim
  use parser_wannier90_tb, &
    only:nRvec
  implicit none
 
  allocatable fk_ex_der(:,:,:)
  complex*16 fk_ex_der
 
  contains
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_fk_ex_der_k()
	implicit none
 
	allocatable :: fk_ex_aux(:,:),fk_ex_der_aux(:,:,:)
 
	!here
	dimension rk1vector(npointstotal),rk2vector(npointstotal),rk3vector(npointstotal)
 
	integer :: ibz
	integer :: nn,j,j_aux
	integer :: nj
 
	real*8 :: rk1vector,rk2vector,rk3vector
	real*8 :: xc,yc,zc,xp_bz,yp_bz,zp_bz,xc_bz,yc_bz,zc_bz
	real*8 :: rk1,rk2,rk3,rkxp,rkyp,rkzp
	real*8 :: rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz
 
	complex*16 :: fk_ex_aux,fk_ex_der_aux
	complex*16 :: fk_ex_k_interp,fk_ex_k_interp_back,fk_ex_k_interp_for
	complex*16 :: f_grid
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	logical :: active_x, active_y, active_z
	active_x = (NORM2(real(nRvec(:,1))) /= 0.0d0)
	active_y = (NORM2(real(nRvec(:,2))) /= 0.0d0)
	active_z = (NORM2(real(nRvec(:,3))) /= 0.0d0)
 
	!auxiliary arrays within this subroutine
	allocate(fk_ex_aux(npointstotal,norb_ex_cut))
	allocate(fk_ex_der_aux(3,npointstotal,norb_ex_cut))
	fk_ex_aux=0.0d0
	fk_ex_der_aux=0.0d0
 
	!get_k_kc gives the crystal coordinates of the k-point.
	!we first fill a vector with all kc points. Note that it is trivial (-.5 to 0.5 ...)
	!but we can test the routine with this simple task
	do ibz=1,npointstotal
		call get_k_kc(G,rkxvector(ibz),rkyvector(ibz),rkzvector(ibz),rk1vector(ibz), &
		rk2vector(ibz),rk3vector(ibz),xp_bz,yp_bz,zp_bz,xc_bz,yc_bz,zc_bz)
	end do
 
	!We evaluate the k-derivative of the exciton envelope function by finite differences with interpolation
	write(*,*) '   Evaluating exciton envelope function derivative with respect to k...'
 
	!Sum over e-h pairs
	do j=1,norb_ex_band
 
		!fill a matrix with vector with A_vc(npointstotal,norb_ex)
		do ibz=1,npointstotal
			j_aux=(ibz-1)*norb_ex_band+j
			do nn=1,norb_ex_cut
				fk_ex_aux(ibz,nn)=fk_ex(j_aux,nn)
			end do
		end do
 
		!Sum over exciton states
		do nn=1,norb_ex_cut
			do ibz=1,npointstotal
 
				rkxp=rkxvector(ibz)
				rkyp=rkyvector(ibz)
				rkzp=rkzvector(ibz)
				call get_k_kc(G,rkxp,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
 
				! Default: all derivatives zero (border or inactive directions)
				fk_ex_der_aux(1,ibz,nn)=0.0d0
				fk_ex_der_aux(2,ibz,nn)=0.0d0
				fk_ex_der_aux(3,ibz,nn)=0.0d0
 
				! 1D: only one direction is active; finite difference only along that axis
				if (ndim == 1) then
 
					if ( active_z ) then
						! active axis: k3
						if (abs(abs(rk3)-0.5d0) .lt. 0.00001d0) cycle
 
						rkxp=rkxvector(ibz)
						rkyp=rkyvector(ibz)
						rkzp=rkzvector(ibz)+dk
						call get_k_kc(G,rkxp,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
						call get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector,rk3vector,&
							norb_ex,norb_ex_cut,fk_ex_aux,rk1,rk2,rk3,nn,fk_ex_k_interp_for)
 
						rkzp=rkzvector(ibz)-dk
						call get_k_kc(G,rkxp,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
						call get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector,rk3vector,&
							norb_ex,norb_ex_cut,fk_ex_aux,rk1,rk2,rk3,nn,fk_ex_k_interp_back)
 
						fk_ex_der_aux(3,ibz,nn)=(fk_ex_k_interp_for-fk_ex_k_interp_back)/(2.0d0*dk)
 
					elseif ( active_y ) then
						! active axis: k2
						if (abs(abs(rk2)-0.5d0) .lt. 0.00001d0) cycle
 
						rkxp=rkxvector(ibz)
						rkyp=rkyvector(ibz)+dk
						rkzp=rkzvector(ibz)
						call get_k_kc(G,rkxp,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
						call get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector,rk3vector,&
							norb_ex,norb_ex_cut,fk_ex_aux,rk1,rk2,rk3,nn,fk_ex_k_interp_for)
 
						rkyp=rkyvector(ibz)-dk
						call get_k_kc(G,rkxp,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
						call get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector,rk3vector,&
							norb_ex,norb_ex_cut,fk_ex_aux,rk1,rk2,rk3,nn,fk_ex_k_interp_back)
 
						fk_ex_der_aux(2,ibz,nn)=(fk_ex_k_interp_for-fk_ex_k_interp_back)/(2.0d0*dk)
 
					else
						! active axis: k1
						if (abs(abs(rk1)-0.5d0) .lt. 0.00001d0) cycle
 
						rkxp=rkxvector(ibz)+dk
						rkyp=rkyvector(ibz)
						rkzp=rkzvector(ibz)
						call get_k_kc(G,rkxp,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
						call get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector,rk3vector,&
							norb_ex,norb_ex_cut,fk_ex_aux,rk1,rk2,rk3,nn,fk_ex_k_interp_for)
 
						rkxp=rkxvector(ibz)-dk
						call get_k_kc(G,rkxp,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
						call get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector,rk3vector,&
							norb_ex,norb_ex_cut,fk_ex_aux,rk1,rk2,rk3,nn,fk_ex_k_interp_back)
 
						fk_ex_der_aux(1,ibz,nn)=(fk_ex_k_interp_for-fk_ex_k_interp_back)/(2.0d0*dk)
 
					endif
 
				! 2D: two directions are active; finite difference along both active axes
				elseif (ndim == 2) then
 
					if ( active_x .and. active_y ) then
						! active axes: k1, k2
						if (abs(abs(rk1)-0.5d0) .lt. 0.00001d0 .or. &
							abs(abs(rk2)-0.5d0) .lt. 0.00001d0) cycle
 
						! d/dk1
						rkxp=rkxvector(ibz)+dk
						rkyp=rkyvector(ibz)
						rkzp=rkzvector(ibz)
						call get_k_kc(G,rkxp,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
						call get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector,rk3vector,&
							norb_ex,norb_ex_cut,fk_ex_aux,rk1,rk2,rk3,nn,fk_ex_k_interp_for)
 
						rkxp=rkxvector(ibz)-dk
						call get_k_kc(G,rkxp,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
						call get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector,rk3vector,&
							norb_ex,norb_ex_cut,fk_ex_aux,rk1,rk2,rk3,nn,fk_ex_k_interp_back)
 
						fk_ex_der_aux(1,ibz,nn)=(fk_ex_k_interp_for-fk_ex_k_interp_back)/(2.0d0*dk)
 
						! d/dk2
						rkxp=rkxvector(ibz)
						rkyp=rkyvector(ibz)+dk
						rkzp=rkzvector(ibz)
						call get_k_kc(G,rkxp,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
						call get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector,rk3vector,&
							norb_ex,norb_ex_cut,fk_ex_aux,rk1,rk2,rk3,nn,fk_ex_k_interp_for)
 
						rkyp=rkyvector(ibz)-dk
						call get_k_kc(G,rkxp,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
						call get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector,rk3vector,&
							norb_ex,norb_ex_cut,fk_ex_aux,rk1,rk2,rk3,nn,fk_ex_k_interp_back)
 
						fk_ex_der_aux(2,ibz,nn)=(fk_ex_k_interp_for-fk_ex_k_interp_back)/(2.0d0*dk)
 
					elseif ( active_x .and. active_z ) then
						! active axes: k1, k3
						if (abs(abs(rk1)-0.5d0) .lt. 0.00001d0 .or. &
							abs(abs(rk3)-0.5d0) .lt. 0.00001d0) cycle
 
						! d/dk1
						rkxp=rkxvector(ibz)+dk
						rkyp=rkyvector(ibz)
						rkzp=rkzvector(ibz)
						call get_k_kc(G,rkxp,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
						call get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector,rk3vector,&
							norb_ex,norb_ex_cut,fk_ex_aux,rk1,rk2,rk3,nn,fk_ex_k_interp_for)
 
						rkxp=rkxvector(ibz)-dk
						call get_k_kc(G,rkxp,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
						call get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector,rk3vector,&
							norb_ex,norb_ex_cut,fk_ex_aux,rk1,rk2,rk3,nn,fk_ex_k_interp_back)
 
						fk_ex_der_aux(1,ibz,nn)=(fk_ex_k_interp_for-fk_ex_k_interp_back)/(2.0d0*dk)
 
						! d/dk3
						rkxp=rkxvector(ibz)
						rkyp=rkyvector(ibz)
						rkzp=rkzvector(ibz)+dk
						call get_k_kc(G,rkxp,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
						call get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector,rk3vector,&
							norb_ex,norb_ex_cut,fk_ex_aux,rk1,rk2,rk3,nn,fk_ex_k_interp_for)
 
						rkzp=rkzvector(ibz)-dk
						call get_k_kc(G,rkxp,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
						call get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector,rk3vector,&
							norb_ex,norb_ex_cut,fk_ex_aux,rk1,rk2,rk3,nn,fk_ex_k_interp_back)
 
						fk_ex_der_aux(3,ibz,nn)=(fk_ex_k_interp_for-fk_ex_k_interp_back)/(2.0d0*dk)
 
					else
						! active axes: k2, k3
						if (abs(abs(rk2)-0.5d0) .lt. 0.00001d0 .or. &
							abs(abs(rk3)-0.5d0) .lt. 0.00001d0) cycle
 
						! d/dk2
						rkxp=rkxvector(ibz)
						rkyp=rkyvector(ibz)+dk
						rkzp=rkzvector(ibz)
						call get_k_kc(G,rkxp,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
						call get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector,rk3vector,&
							norb_ex,norb_ex_cut,fk_ex_aux,rk1,rk2,rk3,nn,fk_ex_k_interp_for)
 
						rkyp=rkyvector(ibz)-dk
						call get_k_kc(G,rkxp,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
						call get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector,rk3vector,&
							norb_ex,norb_ex_cut,fk_ex_aux,rk1,rk2,rk3,nn,fk_ex_k_interp_back)
 
						fk_ex_der_aux(2,ibz,nn)=(fk_ex_k_interp_for-fk_ex_k_interp_back)/(2.0d0*dk)
 
						! d/dk3
						rkxp=rkxvector(ibz)
						rkyp=rkyvector(ibz)
						rkzp=rkzvector(ibz)+dk
						call get_k_kc(G,rkxp,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
						call get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector,rk3vector,&
							norb_ex,norb_ex_cut,fk_ex_aux,rk1,rk2,rk3,nn,fk_ex_k_interp_for)
 
						rkzp=rkzvector(ibz)-dk
						call get_k_kc(G,rkxp,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
						call get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector,rk3vector,&
							norb_ex,norb_ex_cut,fk_ex_aux,rk1,rk2,rk3,nn,fk_ex_k_interp_back)
 
						fk_ex_der_aux(3,ibz,nn)=(fk_ex_k_interp_for-fk_ex_k_interp_back)/(2.0d0*dk)
 
					endif
 
				! 3D: all three directions active
				else
 
					if (abs(abs(rk1)-0.5d0) .lt. 0.00001d0 .or. &
						abs(abs(rk2)-0.5d0) .lt. 0.00001d0 .or. &
						abs(abs(rk3)-0.5d0) .lt. 0.00001d0) cycle
 
					! d/dk1
					rkxp=rkxvector(ibz)+dk
					rkyp=rkyvector(ibz)
					rkzp=rkzvector(ibz)
					call get_k_kc(G,rkxp,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
					call get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector,rk3vector,&
						norb_ex,norb_ex_cut,fk_ex_aux,rk1,rk2,rk3,nn,fk_ex_k_interp_for)
 
					rkxp=rkxvector(ibz)-dk
					call get_k_kc(G,rkxp,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
					call get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector,rk3vector,&
						norb_ex,norb_ex_cut,fk_ex_aux,rk1,rk2,rk3,nn,fk_ex_k_interp_back)
 
					fk_ex_der_aux(1,ibz,nn)=(fk_ex_k_interp_for-fk_ex_k_interp_back)/(2.0d0*dk)
 
					! d/dk2
					rkxp=rkxvector(ibz)
					rkyp=rkyvector(ibz)+dk
					rkzp=rkzvector(ibz)
					call get_k_kc(G,rkxp,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
					call get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector,rk3vector,&
						norb_ex,norb_ex_cut,fk_ex_aux,rk1,rk2,rk3,nn,fk_ex_k_interp_for)
 
					rkyp=rkyvector(ibz)-dk
					call get_k_kc(G,rkxp,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
					call get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector,rk3vector,&
						norb_ex,norb_ex_cut,fk_ex_aux,rk1,rk2,rk3,nn,fk_ex_k_interp_back)
 
					fk_ex_der_aux(2,ibz,nn)=(fk_ex_k_interp_for-fk_ex_k_interp_back)/(2.0d0*dk)
 
					! d/dk3
					rkxp=rkxvector(ibz)
					rkyp=rkyvector(ibz)
					rkzp=rkzvector(ibz)+dk
					call get_k_kc(G,rkxp,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
					call get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector,rk3vector,&
						norb_ex,norb_ex_cut,fk_ex_aux,rk1,rk2,rk3,nn,fk_ex_k_interp_for)
 
					rkzp=rkzvector(ibz)-dk
					call get_k_kc(G,rkxp,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
					call get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector,rk3vector,&
						norb_ex,norb_ex_cut,fk_ex_aux,rk1,rk2,rk3,nn,fk_ex_k_interp_back)
 
					fk_ex_der_aux(3,ibz,nn)=(fk_ex_k_interp_for-fk_ex_k_interp_back)/(2.0d0*dk)
 
				endif
 
			end do
		end do
 
		!Save e-h pair wf
		do ibz=1,npointstotal
			j_aux=(ibz-1)*norb_ex_band+j
			do nn=1,norb_ex_cut
				do nj=1,3
					fk_ex_der(nj,j_aux,nn)=fk_ex_der_aux(nj,ibz,nn)
				end do
			end do
		end do
 
	end do
 
end subroutine get_fk_ex_der_k
 

 


 
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
subroutine get_fk_ex_k_interp(npointstotal,rk1vector,rk2vector,rk3vector,&
			norb_ex,norb_ex_cut,fk_ex,rk1,rk2,rk3,nn,fk_ex_k_interp)
	implicit none 
 
	!in/out
	integer :: npointstotal,norb_ex,norb_ex_cut,nn
 
	dimension rk1vector(npointstotal),rk2vector(npointstotal),rk3vector(npointstotal)
	dimension fk_ex(npointstotal,norb_ex_cut)
 
	real*8 :: rk1vector,rk2vector,rk3vector
	real*8 :: rk1,rk2,rk3
	complex*16 :: fk_ex,fk_ex_k_interp
 
	!here
	integer :: nside,nblock
	integer :: ibz_q11,ibz_q21,ibz_q31,ibz_q12,ibz_q22,ibz_q32,ibz_q13,ibz_q23,ibz_q33
	integer :: ibz_q111,ibz_q211,ibz_q121,ibz_q221,ibz_q112,ibz_q212,ibz_q122,ibz_q222
	real*8 :: slice
	real*8 :: x1,x2,x3,y1,y2,y3,z1,z2,z3,x,y,z
	complex*16 f_q11,f_q12,f_q13,f_q21,f_q22,f_q23,f_q31,f_q32,f_q33
	complex*16 f_q111,f_q211,f_q121,f_q221,f_q112,f_q212,f_q122,f_q222
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	logical :: active_x, active_y, active_z
	active_x = (NORM2(real(nRvec(:,1))) /= 0.0d0)
	active_y = (NORM2(real(nRvec(:,2))) /= 0.0d0)
	active_z = (NORM2(real(nRvec(:,3))) /= 0.0d0)
	
	! 1D
	if (ndim == 1) then
		!identify the block in which the point is
		nside=nint(dble(npointstotal))
		slice=1.0d0/dble(nside-1)
 
		if ( active_z ) then
			! only k3 is active
			nblock=int((rk3+0.5d0)/slice)+1
			ibz_q13=nblock
			ibz_q33=nblock+1
 
			if (ibz_q13.gt.npointstotal .or. &
				ibz_q33.gt.npointstotal) then
				write(*,*) 'interpolation went wrong'
				fk_ex_k_interp=0.0d0
				pause
			end if
			if (ibz_q13.lt.1 .or. &
				ibz_q33.lt.1) then
				write(*,*) 'interpolation went wrong'
				fk_ex_k_interp=0.0d0
				pause
			end if
 
			!interpolate
			z1=rk3vector(ibz_q13)
			z3=rk3vector(ibz_q33)
			z=rk3
 
			f_q13=fk_ex(ibz_q13,nn)
			f_q33=fk_ex(ibz_q33,nn)
 
			fk_ex_k_interp=f_q13*(z3-z)/(z3-z1) + f_q33*(z-z1)/(z3-z1)
 
		elseif ( active_y ) then
			! only k2 is active
			nblock=int((rk2+0.5d0)/slice)+1
			ibz_q12=nblock
			ibz_q22=nblock+1
 
			if (ibz_q12.gt.npointstotal .or. &
				ibz_q22.gt.npointstotal) then
				write(*,*) 'interpolation went wrong'
				fk_ex_k_interp=0.0d0
				pause
			end if
			if (ibz_q12.lt.1 .or. &
				ibz_q22.lt.1) then
				write(*,*) 'interpolation went wrong'
				fk_ex_k_interp=0.0d0
				pause
			end if
 
			!interpolate
			y1=rk2vector(ibz_q12)
			y2=rk2vector(ibz_q22)
			y=rk2
 
			f_q12=fk_ex(ibz_q12,nn)
			f_q22=fk_ex(ibz_q22,nn)
 
			fk_ex_k_interp=f_q12*(y2-y)/(y2-y1) + f_q22*(y-y1)/(y2-y1)
 
		else
			! only k1 is active
			nblock=int((rk1+0.5d0)/slice)+1
			ibz_q11=nblock
			ibz_q21=nblock+1
 
			if (ibz_q11.gt.npointstotal .or. &
				ibz_q21.gt.npointstotal) then
				write(*,*) 'interpolation went wrong'
				fk_ex_k_interp=0.0d0
				pause
			end if
			if (ibz_q11.lt.1 .or. &
				ibz_q21.lt.1) then
				write(*,*) 'interpolation went wrong'
				fk_ex_k_interp=0.0d0
				pause
			end if
 
			!interpolate
			x1=rk1vector(ibz_q11)
			x2=rk1vector(ibz_q21)
			x=rk1
 
			f_q11=fk_ex(ibz_q11,nn)
			f_q21=fk_ex(ibz_q21,nn)
 
			fk_ex_k_interp=f_q11*(x2-x)/(x2-x1) + f_q21*(x-x1)/(x2-x1)
 
		endif
 
	! 2D
	elseif ( ndim == 2) then
		!identify the block in which the point is
		nside=nint(sqrt(dble(npointstotal)))
		slice=1.0d0/dble(nside-1)
		
		if ( active_x .and. active_y ) then
			!get the corners
			nblock=int((rk1+0.5d0)/slice)+1
			ibz_q11=nblock+int((rk2+0.5d0)/slice)*nside
			ibz_q21=ibz_q11+1
			
			ibz_q12=ibz_q11+nside
			ibz_q22=ibz_q21+nside  
			
			if (ibz_q11.gt.npointstotal .or. &
				ibz_q21.gt.npointstotal .or. &
				ibz_q12.gt.npointstotal .or. &
				ibz_q22.gt.npointstotal) then
				write(*,*) 'interpolation went wrong'
				fk_ex_k_interp=0.0d0
				pause
			end if
			if (ibz_q11.lt.1 .or. &
				ibz_q21.lt.1 .or. &
				ibz_q12.lt.1 .or. &
				ibz_q22.lt.1) then
				write(*,*) 'interpolation went wrong'
				fk_ex_k_interp=0.0d0
				pause
			end if
 
			!interpolate
			x1=rk1vector(ibz_q11)
			x2=rk1vector(ibz_q21)
			
			y1=rk2vector(ibz_q11)
			y2=rk2vector(ibz_q12)
			
			x=rk1
			y=rk2
 
			f_q11=fk_ex(ibz_q11,nn)
			f_q21=fk_ex(ibz_q21,nn)
			
			f_q12=fk_ex(ibz_q12,nn)
			f_q22=fk_ex(ibz_q22,nn)
			
			fk_ex_k_interp=1.0d0/((x2-x1)*(y2-y1))* &
				(f_q11*(x2-x)*(y2-y)+f_q21*(x-x1)*(y2-y) &
				+f_q12*(x2-x)*(y-y1)+f_q22*(x-x1)*(y-y1))
				
		elseif ( active_x .and. active_z ) then
			! y->z
			!get the corners
			nblock=int((rk1+0.5d0)/slice)+1
			ibz_q11=nblock+int((rk3+0.5d0)/slice)*nside
			ibz_q31=ibz_q11+1
			
			ibz_q13=ibz_q11+nside
			ibz_q33=ibz_q31+nside  
			
			if (ibz_q11.gt.npointstotal .or. &
				ibz_q31.gt.npointstotal .or. &
				ibz_q13.gt.npointstotal .or. &
				ibz_q33.gt.npointstotal) then
				write(*,*) 'interpolation went wrong'
				fk_ex_k_interp=0.0d0
				pause
			end if
			if (ibz_q11.lt.1 .or. &
				ibz_q31.lt.1 .or. &
				ibz_q13.lt.1 .or. &
				ibz_q33.lt.1) then
				write(*,*) 'interpolation went wrong'
				fk_ex_k_interp=0.0d0
				pause
			end if
 
			!interpolate
			x1=rk1vector(ibz_q11)
			x3=rk1vector(ibz_q31)
			
			z1=rk3vector(ibz_q11)
			z3=rk3vector(ibz_q13)
			
			x=rk1
			z=rk3
 
			f_q11=fk_ex(ibz_q11,nn)
			f_q31=fk_ex(ibz_q31,nn)
			
			f_q13=fk_ex(ibz_q13,nn)
			f_q33=fk_ex(ibz_q33,nn)
			
			fk_ex_k_interp=1.0d0/((x3-x1)*(z3-z1))* &
				(f_q11*(x3-x)*(z3-z)+f_q31*(x-x1)*(z3-z) &
				+f_q13*(x3-x)*(z-z1)+f_q33*(x-x1)*(z-z1))
		else
			! x->y
			!get the corners
			nblock=int((rk2+0.5d0)/slice)+1
			ibz_q22=nblock+int((rk3+0.5d0)/slice)*nside
			ibz_q32=ibz_q22+1
			
			ibz_q23=ibz_q22+nside
			ibz_q33=ibz_q32+nside  
			
			if (ibz_q22.gt.npointstotal .or. &
				ibz_q32.gt.npointstotal .or. &
				ibz_q23.gt.npointstotal .or. &
				ibz_q33.gt.npointstotal) then
				write(*,*) 'interpolation went wrong'
				fk_ex_k_interp=0.0d0
				pause
			end if
			if (ibz_q22.lt.1 .or. &
				ibz_q32.lt.1 .or. &
				ibz_q23.lt.1 .or. &
				ibz_q33.lt.1) then
				write(*,*) 'interpolation went wrong'
				fk_ex_k_interp=0.0d0
				pause
			end if
 
			!interpolate
			y2=rk2vector(ibz_q22)
			y3=rk2vector(ibz_q32)
			
			z2=rk3vector(ibz_q22)
			z3=rk3vector(ibz_q23)
			
			y=rk2
			z=rk3
 
			f_q22=fk_ex(ibz_q22,nn)
			f_q32=fk_ex(ibz_q32,nn)
			
			f_q23=fk_ex(ibz_q23,nn)
			f_q33=fk_ex(ibz_q33,nn)
			
			fk_ex_k_interp=1.0d0/((y3-y2)*(z3-z2))* &
				(f_q22*(y3-y)*(z3-z)+f_q32*(y-y2)*(z3-z) &
				+f_q23*(y3-y)*(z-z2)+f_q33*(y-y2)*(z-z2))
		endif
 
	! 3D
	else
		!identify the block in which the point is
		nside=nint(dble(npointstotal)**(1.0d0/3.0d0))
		slice=1.0d0/dble(nside-1)
 
		! indices follow: ibz_qXYZ where X=k1 corner (1=lo,2=hi),
		!                                   Y=k2 corner (1=lo,2=hi),
		!                                   Z=k3 corner (1=lo,2=hi)
		! flat layout: index = ix + iy*nside + iz*nside^2  (1-based)
		nblock=    int((rk1+0.5d0)/slice) &
		        +  int((rk2+0.5d0)/slice)*nside &
		        +  int((rk3+0.5d0)/slice)*nside*nside + 1
 
		ibz_q111=nblock
		ibz_q211=nblock+1
		ibz_q121=nblock+nside
		ibz_q221=nblock+nside+1
		ibz_q112=nblock+nside*nside
		ibz_q212=nblock+nside*nside+1
		ibz_q122=nblock+nside*nside+nside
		ibz_q222=nblock+nside*nside+nside+1
 
		if (ibz_q111.gt.npointstotal .or. &
			ibz_q211.gt.npointstotal .or. &
			ibz_q121.gt.npointstotal .or. &
			ibz_q221.gt.npointstotal .or. &
			ibz_q112.gt.npointstotal .or. &
			ibz_q212.gt.npointstotal .or. &
			ibz_q122.gt.npointstotal .or. &
			ibz_q222.gt.npointstotal) then
			write(*,*) 'interpolation went wrong'
			fk_ex_k_interp=0.0d0
			pause
		end if
		if (ibz_q111.lt.1 .or. &
			ibz_q211.lt.1 .or. &
			ibz_q121.lt.1 .or. &
			ibz_q221.lt.1 .or. &
			ibz_q112.lt.1 .or. &
			ibz_q212.lt.1 .or. &
			ibz_q122.lt.1 .or. &
			ibz_q222.lt.1) then
			write(*,*) 'interpolation went wrong'
			fk_ex_k_interp=0.0d0
			pause
		end if
 
		!interpolate
		x1=rk1vector(ibz_q111)
		x2=rk1vector(ibz_q211)
 
		y1=rk2vector(ibz_q111)
		y2=rk2vector(ibz_q121)
 
		z1=rk3vector(ibz_q111)
		z2=rk3vector(ibz_q112)
 
		x=rk1
		y=rk2
		z=rk3
 
		f_q111=fk_ex(ibz_q111,nn)
		f_q211=fk_ex(ibz_q211,nn)
		f_q121=fk_ex(ibz_q121,nn)
		f_q221=fk_ex(ibz_q221,nn)
		f_q112=fk_ex(ibz_q112,nn)
		f_q212=fk_ex(ibz_q212,nn)
		f_q122=fk_ex(ibz_q122,nn)
		f_q222=fk_ex(ibz_q222,nn)
 
		fk_ex_k_interp=1.0d0/((x2-x1)*(y2-y1)*(z2-z1))* &
			(f_q111*(x2-x)*(y2-y)*(z2-z) &
			+f_q211*(x-x1)*(y2-y)*(z2-z) &
			+f_q121*(x2-x)*(y-y1)*(z2-z) &
			+f_q221*(x-x1)*(y-y1)*(z2-z) &
			+f_q112*(x2-x)*(y2-y)*(z-z1) &
			+f_q212*(x-x1)*(y2-y)*(z-z1) &
			+f_q122*(x2-x)*(y-y1)*(z-z1) &
			+f_q222*(x-x1)*(y-y1)*(z-z1))
 
	endif
end subroutine get_fk_ex_k_interp
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
subroutine get_k_kc(G,xp,yp,zp,xc,yc,zc,xp_bz,yp_bz,zp_bz,xc_bz,yc_bz,zc_bz)
	implicit none
 
	dimension G(3,3)
 
	real*8 :: G
	real*8 :: xp,yp,zp,xc,yc,zc,xp_bz,yp_bz,zp_bz,xc_bz,yc_bz,zc_bz
	real*8 :: det
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
	logical :: active_x, active_y, active_z
	active_x = (NORM2(real(nRvec(:,1))) /= 0.0d0)
	active_y = (NORM2(real(nRvec(:,2))) /= 0.0d0)
	active_z = (NORM2(real(nRvec(:,3))) /= 0.0d0)
	! 1D
	if (ndim == 1) then
		if ( active_z ) then
			xc=0.0d0
			yc=0.0d0
			zc=zp/(G(3,3))
			
			xc_bz=0.0d0
			yc_bz=0.0d0
			zc_bz=zc-dble(int(zc/0.5d0))
			
			xp_bz=0.0d0
			yp_bz=0.0d0
			zp_bz=zc_bz*G(3,3)
			
		elseif ( active_y ) then
			xc=0.0d0
			yc=yp/(G(2,2))
			zc=0.0d0
			
			xc_bz=0.0d0
			yc_bz=yc-dble(int(yc/0.5d0))
			zc_bz=0.0d0
			
			xp_bz=0.0d0
			yp_bz=yc_bz*G(2,2)
			zp_bz=0.0d0
		else
			xc=xp/(G(1,1))
			yc=0.0d0
			zc=0.0d0
			
			xc_bz=xc-dble(int(xc/0.5d0))
			yc_bz=0.0d0
			zc_bz=0.0d0
			
			xp_bz=xc_bz*G(1,1)
			yp_bz=0.0d0
			zp_bz=0.0d0
		endif
	
	! 2D
	elseif ( ndim == 2) then
	
		if ( active_x .and. active_y ) then
			xc=(xp*G(2,2)-G(2,1)*yp)/(G(1,1)*G(2,2)-G(2,1)*G(1,2))
			yc=(G(1,1)*yp-xp*G(1,2))/(G(1,1)*G(2,2)-G(2,1)*G(1,2))
			zc=0.0d0
			
			xc_bz=xc-dble(int(xc/0.5d0))
			yc_bz=yc-dble(int(yc/0.5d0))
			zc_bz=0.0d0
			
			xp_bz=xc_bz*G(1,1)+yc_bz*G(2,1)
			yp_bz=xc_bz*G(1,2)+yc_bz*G(2,2)
			zp_bz=0.0d0
		elseif ( active_x .and. active_z ) then
			xc=(xp*G(3,3)-G(3,1)*yp)/(G(1,1)*G(3,3)-G(3,1)*G(1,3))
			yc=0.0d0
			zc=(G(1,1)*yp-xp*G(1,3))/(G(1,1)*G(3,3)-G(3,1)*G(1,3))
			
			xc_bz=xc-dble(int(xc/0.5d0))
			yc_bz=0.0d0
			zc_bz=zc-dble(int(zc/0.5d0))
			
			xp_bz=xc_bz*G(1,1)+yc_bz*G(3,1)
			yp_bz=0.0d0
			zp_bz=xc_bz*G(1,3)+yc_bz*G(3,3)
		else
			xc=0.0d0
			yc=(yp*G(3,3)-G(3,2)*zp)/(G(2,2)*G(3,3)-G(3,2)*G(2,3))
			zc=(G(2,2)*zp-yp*G(2,3))/(G(2,2)*G(3,3)-G(3,2)*G(2,3))
			
			xc_bz=0.0d0
			yc_bz=yc-dble(int(yc/0.5d0))
			zc_bz=zc-dble(int(zc/0.5d0))
			
			xp_bz=0.0d0
			yp_bz=yc_bz*G(2,2)+zc_bz*G(3,2)
			zp_bz=yc_bz*G(2,3)+zc_bz*G(3,3)
		endif
		
	! 3D
	else
		det=(- G(1,3)*G(2,2)*G(3,1) + G(1,2)*G(2,3)*G(3,1) + G(1,3)*G(2,1)*G(3,2) &
        - G(1,1)*G(2,3)*G(3,2) - G(1,2)*G(2,1)*G(3,3) + G(1,1)*G(2,2)*G(3,3) )
		
		xc=( xp*(G(2,2)*G(3,3)-G(2,3)*G(3,2)) + yp*(G(1,3)*G(3,2)-G(1,2)*G(3,3)) + zp*(G(1,2)*G(2,3)-G(1,3)*G(2,2)) )/(det)
		yc=( xp*(G(2,3)*G(3,1)-G(2,1)*G(3,3)) + yp*(G(1,1)*G(3,3)-G(1,3)*G(3,1)) + zp*(G(1,3)*G(2,1)-G(1,1)*G(2,3)) )/(det)
		zc=( xp*(G(2,1)*G(3,2)-G(2,2)*G(3,1)) + yp*(G(1,2)*G(3,1)-G(1,1)*G(3,2)) + zp*(G(1,1)*G(2,2)-G(1,2)*G(2,1)) )/(det)
 
		xc_bz=xc-dble(int(xc/0.5d0))
		yc_bz=yc-dble(int(yc/0.5d0))
		zc_bz=zc-dble(int(zc/0.5d0))
 
		xp_bz=xc_bz*G(1,1)+yc_bz*G(2,1)+zc_bz*G(3,1)
		yp_bz=xc_bz*G(1,2)+yc_bz*G(2,2)+zc_bz*G(3,2)
		zp_bz=xc_bz*G(1,3)+yc_bz*G(2,3)+zc_bz*G(3,3)
	endif
end subroutine get_k_kc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
end module exciton_envelopes
