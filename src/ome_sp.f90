module ome_sp
   use constants_math
   use parser_wannier90_tb, &
      only:material_name,nR,nRvec,norb,R,shop,hhop,rhop_c
   use parser_optics_xatu_dim, &
      only:npointstotal,rkxvector,rkyvector,rkzvector, &
      nband_ex,nband_index,nv_ex,nc_ex

   implicit none

   allocatable vme_ex_band(:,:,:,:)
   allocatable ek(:,:)
   allocatable gen_der_ex_band(:,:,:,:,:)
   allocatable shift_vector_ex_band(:,:,:,:,:)
   allocatable berry_eigen_ex_band(:,:,:,:)

   real*8 ek
   real*8 shift_vector_ex_band
   complex*16 vme_ex_band
   complex*16 gen_der_ex_band
   complex*16 berry_eigen_ex_band
contains
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_ome_sp(iflag_norder)
      implicit none

      integer iflag_norder
      integer ibz,ibz_sum
      integer i,j,ii,jj,nj
      integer naux1

      allocatable :: skernel(:,:)
      allocatable :: hkernel(:,:)
      allocatable :: sderkernel(:,:,:)
      allocatable :: hderkernel(:,:,:)
      allocatable :: akernel(:,:,:)

      allocatable :: e(:)
      allocatable :: ecomplex(:)
      allocatable :: hk_ev(:,:)
      allocatable :: vme(:,:,:)

      allocatable :: abc(:,:,:)
      allocatable :: gd1(:,:,:,:)
      allocatable :: gd2(:,:,:,:)
      allocatable :: gd3(:,:,:,:)
      allocatable :: gen_der(:,:,:,:)
      allocatable :: vme_der(:,:,:,:)
      allocatable :: shift_vector(:,:,:,:)
      allocatable :: berry_eigen1(:,:,:)
      allocatable :: berry_eigen2(:,:,:)
      allocatable :: berry_eigen(:,:,:)

      real(8) rkx,rky,rkz
      real(8) :: e
      real(8) :: shift_vector

      complex*16 :: skernel,hkernel,sderkernel,hderkernel,akernel
      complex*16 :: hk_ev
      complex*16 :: ecomplex
      complex*16 :: vjseudoa,vjseudob
      complex*16 :: vme
      complex*16 :: vme_der
      complex*16 :: berry_eigen1,berry_eigen2,berry_eigen

      complex*16 abc,gd1,gd2,gd3,gen_der
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,*) '5. Entering ome_sp'
      !allocate to save k-dependent vme and energies between involved band
      allocate(vme_ex_band(npointstotal,3,nband_ex,nband_ex))
      allocate(ek(npointstotal,nband_ex))
      allocate(berry_eigen_ex_band(npointstotal,3,nband_ex,nband_ex))
      allocate(gen_der_ex_band(npointstotal,3,3,nband_ex,nband_ex))
      allocate(shift_vector_ex_band(npointstotal,3,3,nband_ex,nband_ex))
      gen_der_ex_band=0.0d0
      shift_vector_ex_band=0.0d0
      berry_eigen_ex_band=0.0d0
      vme_ex_band=0.0d0
      ek=0.0d0


      allocate(skernel(norb,norb))
      allocate(hkernel(norb,norb))
      allocate(sderkernel(3,norb,norb))
      allocate(hderkernel(3,norb,norb))
      allocate(akernel(3,norb,norb))

      allocate(e(norb))
      allocate(ecomplex(norb))
      allocate(hk_ev(norb,norb))
      allocate(vme(3,norb,norb))

      allocate(abc(3,norb,norb))
      allocate(gd1(3,3,norb,norb))
      allocate(gd2(3,3,norb,norb))
      allocate(gd3(3,3,norb,norb))
      allocate(gen_der(3,3,norb,norb))
      allocate(vme_der(3,3,norb,norb))
      allocate(shift_vector(3,3,norb,norb))
      allocate(berry_eigen1(3,norb,norb))
      allocate(berry_eigen2(3,norb,norb))
      allocate(berry_eigen(3,norb,norb))

      !if (iflag_norder.eq.2) then
      !allocate (gen_der_ex_band(npointstotal,3,3,nband_ex,nband_ex))
      !allocate (shift_vector_ex_band(npointstotal,3,3,nband_ex,nband_ex))
      !gen_der_ex_band=0.0d0
      !shift_vector_ex_band=0.0d0
      !end if
      !open(50,file='coefs_new.dat')
      !Brillouin zone sampling - parallelization

      !write
      !only:material_name,nR,nRvec,norb,R,shop,hhop,rhop_c
      !use parser_optics_xatu_dim, &
      !only:npointstotal,rkxvector,rkyvector,rkzvector, &
      !nband_ex,nband_index,nv_ex,nc_ex
      !write(*,*) material_name
      !write(*,*) nR
      !do i=1,nR
      !write(*,*) nRvec(i,1),nRvec(i,2),nRvec(i,3)
      !end do
      !do i=1,npointstotal
      !write(*,*) rkxvector(i),rkyvector(i),rkzvector(i)
      !end do
      !pause

      write(*,*) '   Calculating optical matrix elements (sp): sampling BZ...'
      ibz_sum=0 !counter for the number of k points in the BZ

      !initializing variables
      !$OMP PARALLEL DO PRIVATE(rkx,rky,rkz), &
      !$OMP PRIVATE(hkernel,skernel,sderkernel,hderkernel,akernel), &
      !$OMP PRIVATE(hk_ev,e,vme), &
      !$OMP PRIVATE(abc,gen_der,gd1,gd2,gd3), &
      !$OMP PRIVATE(vme_der,shift_vector,berry_eigen1,berry_eigen2,berry_eigen)
      do ibz=1,npointstotal
         !!$OMP CRITICAL
         !ibz_sum=ibz_sum+1
         !call percentage_index(ibz_sum,npointstotal,naux1)
         !!$OMP END CRITICAL
         write(*,*) '   Optical matrix elements (sp): k-point',ibz,'/',npointstotal
         !pause
         rkx=rkxvector(ibz)
         rky=rkyvector(ibz)
         rkz=rkzvector(ibz)

         !pause

         !get matrices in the \alpha, \alpha' basis (orbitals,k)
         call get_vme_kernels_ome(rkx,rky,rkz,norb,skernel,sderkernel, &
            hkernel,hderkernel,akernel)
         call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
            hk_ev,e,vme)

         !matrix elements for second-order calculation
         if (iflag_norder.eq.2) then
            call get_gen_der_sumrule(norb,vme,e,abc,gen_der,gd1,gd2,gd3)
            call get_berry_eigen_fourpoint(rkx,rky,rkz,norb,vme_der, &
               shift_vector,berry_eigen1,berry_eigen2,berry_eigen)
         end if

         !write(*,*) ibz,rkx,rky,shift_vector(1,1,2,3),berry_eigen1(1,2,3),berry_eigen2(1,2,3)

         !saving eigenvalues and optical matrix elements at this k point for the bandlist
         do i=1,nband_ex
            ii=nband_index(i)
            ek(ibz,i)=e(ii)
            do nj=1,3
               do j=1,nband_ex
                  jj=nband_index(j)
                  vme_ex_band(ibz,nj,i,j)=vme(nj,ii,jj)

                  if (iflag_norder.eq.2) then
                     shift_vector_ex_band(ibz,nj,1,i,j)=shift_vector(nj,1,ii,jj)
                     shift_vector_ex_band(ibz,nj,2,i,j)=shift_vector(nj,2,ii,jj)
                     shift_vector_ex_band(ibz,nj,3,i,j)=shift_vector(nj,3,ii,jj)

                     gen_der_ex_band(ibz,nj,1,i,j)=gen_der(nj,1,ii,jj)
                     gen_der_ex_band(ibz,nj,2,i,j)=gen_der(nj,2,ii,jj)
                     gen_der_ex_band(ibz,nj,3,i,j)=gen_der(nj,3,ii,jj)

                     berry_eigen_ex_band(ibz,nj,i,j)=berry_eigen(nj,ii,jj)
                     berry_eigen_ex_band(ibz,nj,i,j)=berry_eigen(nj,ii,jj)
                     berry_eigen_ex_band(ibz,nj,i,j)=berry_eigen(nj,ii,jj)

                  end if

               end do
            end do
         end do

         !pause
      end do
      !$OMP END PARALLEL DO
      !close(50)

      !write matrix elements into file
      write(*,*) '   Writing optical matrix elements (sp) into file'
      if (iflag_norder.eq.1) then
         call write_ome_sp_linear(iflag_norder,npointstotal,nband_ex,vme_ex_band,ek)
      end if
      if (iflag_norder.eq.2) then
         call write_ome_sp_nonlinear(iflag_norder,npointstotal,nband_ex,vme_ex_band,ek, &
            gen_der_ex_band,shift_vector_ex_band,berry_eigen_ex_band)
      end if
      write(*,*) '   Optical matrix elements (sp) have been written in file'
      !pause
      !deallocate(vme_ex_band)
      !deallocate(ek)
   end subroutine get_ome_sp
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !you have to addapt this to 3D
   subroutine get_berry_eigen_fourpoint(rkx,rky,rkz,norb,vme_der, &
      shift_vector,berry_eigen1,berry_eigen2,berry_eigen)
      implicit none

      integer nR,norb
      integer nn,nnp
      integer ialpha,ialphap
      integer nj,njp

      dimension hkernel(norb,norb),skernel(norb,norb)
      dimension sderkernel(3,norb,norb),hderkernel(3,norb,norb)

      dimension hk_alpha(norb,norb),hk_ev(norb,norb),e(norb)
      dimension akernel(3,norb,norb)

      dimension vjseudoa(3,norb,norb),vjseudob(3,norb,norb),vme(3,norb,norb)
      dimension berry_eigen1(3,norb,norb),berry_eigen2(3,norb,norb)
      dimension berry_eigen(3,norb,norb)
      dimension vme_der(3,3,norb,norb)
      dimension vme_der_phase(3,3,norb,norb)
      dimension hk_ev_neigh(5,norb,norb)
      dimension vme_neigh(5,3,norb,norb)

      dimension shift_vector(3,3,norb,norb)

      real*8 rkx,rky,rkz,rkx_neigh,rky_neigh,rkz_neigh
      real*8 e
      real*8 vme_der_phase
      real*8 shift_vector
      real*8 ph1,ph2,ph3,ph4

      complex*16 hkernel,akernel,skernel,sderkernel,hderkernel
      complex*16 hk_alpha,hk_ev,vjseudoa,vjseudob,vme
      complex*16 amu,amup,aux1,aux2,aux3,aux4,factor
      complex*16 vme_der
      complex*16 berry_eigen1,berry_eigen2,berry_eigen
      complex*16 hk_ev_neigh,vme_neigh
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      vme=0.0d0

      hk_ev_neigh=0.0d0
      vme_neigh=0.0d0

      berry_eigen1=0.0d0
      berry_eigen2=0.0d0
      berry_eigen=0.0d0
      vme_der=0.0d0
      shift_vector=0.0d0

      rkx_neigh=rkx+dk
      rky_neigh=rky
      rkz_neigh=rkz
      call get_vme_kernels_ome(rkx_neigh,rky_neigh,rkz_neigh,norb,skernel, &
         sderkernel,hkernel,hderkernel,akernel)
      !velocity matrix elements
      call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
         hk_ev_neigh(1,:,:),e,vme_neigh(1,:,:,:))

      !lets try to do only a n,n' loop
      rkx_neigh=rkx
      rky_neigh=rky-dk
      rkz_neigh=rkz
      call get_vme_kernels_ome(rkx_neigh,rky_neigh,rkz_neigh,norb,skernel, &
         sderkernel,hkernel,hderkernel,akernel)
      !velocity matrix elements
      call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
         hk_ev_neigh(2,:,:),e,vme_neigh(2,:,:,:))

      !lets try to do only a n,n' loop
      rkx_neigh=rkx-dk
      rky_neigh=rky
      rkz_neigh=rkz
      call get_vme_kernels_ome(rkx_neigh,rky_neigh,rkz_neigh,norb,skernel, &
         sderkernel,hkernel,hderkernel,akernel)
      !velocity matrix elements
      call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
         hk_ev_neigh(3,:,:),e,vme_neigh(3,:,:,:))

      !lets try to do only a n,n' loop
      rkx_neigh=rkx
      rky_neigh=rky+dk
      rkz_neigh=rkz
      call get_vme_kernels_ome(rkx_neigh,rky_neigh,rkz_neigh,norb,skernel, &
         sderkernel,hkernel,hderkernel,akernel)
      !velocity matrix elements
      call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
         hk_ev_neigh(4,:,:),e,vme_neigh(4,:,:,:))

      !lets try to do only a n,n' loop
      rkx_neigh=rkx
      rky_neigh=rky
      rkz_neigh=rkz
      call get_vme_kernels_ome(rkx_neigh,rky_neigh,rkz_neigh,norb,skernel, &
         sderkernel,hkernel,hderkernel,akernel)
      !velocity matrix elements
      call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
         hk_ev_neigh(5,:,:),e,vme_neigh(5,:,:,:))



      do nn=1,norb
         do nnp=1,norb
            !Direct evaluation of Berry connection. See Esteve-Paredes and Palacios, Scipost 2023
            do ialpha=1,norb
               do ialphap=1,norb
                  !x-dir
                  aux1=(hk_ev_neigh(1,ialphap,nnp)-hk_ev_neigh(3,ialphap,nnp))/(2.0d0*dk)
                  berry_eigen1(1,nn,nnp)=berry_eigen1(1,nn,nnp)+ &
                     complex(0.0d0,1.0d0)*conjg(hk_ev_neigh(5,ialpha,nn))*skernel(ialpha,ialphap)*aux1
                  berry_eigen2(1,nn,nnp)=berry_eigen2(1,nn,nnp)+ &
                     conjg(hk_ev_neigh(5,ialpha,nn))*hk_ev_neigh(5,ialphap,nnp)*akernel(1,ialpha,ialphap)

                  !if (nn.eq.1 .and. nnp.eq.1 .and. ialpha.eq.1 .and. ialphap.eq.1) then
                  !write(*,*) ialpha,ialphap,berry_eigen1(1,nn,nnp),berry_eigen2(1,nn,nnp),akernel(1,ialpha,ialphap)
                  !pause
                  !end if

                  !y-dir
                  aux1=(hk_ev_neigh(4,ialphap,nnp)-hk_ev_neigh(2,ialphap,nnp))/(2.0d0*dk)
                  berry_eigen1(2,nn,nnp)=berry_eigen1(2,nn,nnp)+ &
                     complex(0.0d0,1.0d0)*conjg(hk_ev_neigh(5,ialpha,nn))*skernel(ialpha,ialphap)*aux1
                  berry_eigen2(2,nn,nnp)=berry_eigen2(2,nn,nnp)+ &
                     conjg(hk_ev_neigh(5,ialpha,nn))*hk_ev_neigh(5,ialphap,nnp)*akernel(2,ialpha,ialphap)

                  !z-dir
                  berry_eigen1(3,nn,nnp)=0.0d0
                  berry_eigen2(3,nn,nnp)=berry_eigen2(3,nn,nnp)+ &
                     conjg(hk_ev_neigh(5,ialpha,nn))*hk_ev_neigh(5,ialphap,nnp)*akernel(3,ialpha,ialphap)
               end do
            end do
            do nj=1,3
               !if BC's is diverging at degenrate points, set it to zero
               berry_eigen(nj,nn,nnp)=berry_eigen1(nj,nn,nnp)+berry_eigen2(nj,nn,nnp)
               if (abs(berry_eigen(nj,nn,nnp)).gt.50.0d0) then
                  berry_eigen(nj,nn,nnp)=0.0d0
               end if

               aux1=vme_neigh(1,nj,nn,nnp)
               aux3=vme_neigh(3,nj,nn,nnp)
               aux4=vme_neigh(4,nj,nn,nnp)
               aux2=vme_neigh(2,nj,nn,nnp)
               vme_der(1,nj,nn,nnp)=(aux1-aux3)/(2.0d0*dk)
               vme_der(2,nj,nn,nnp)=(aux4-aux2)/(2.0d0*dk)

               call get_phase(vme_neigh(1,nj,nn,nnp),ph1)
               call get_phase(vme_neigh(3,nj,nn,nnp),ph3)
               call get_phase(vme_neigh(4,nj,nn,nnp),ph4)
               call get_phase(vme_neigh(2,nj,nn,nnp),ph2)
               vme_der_phase(1,nj,nn,nnp)=(ph1-ph3)/(2.0d0*dk)
               vme_der_phase(2,nj,nn,nnp)=(ph4-ph2)/(2.0d0*dk)
            end do
         end do
      end do

      !shift vector
      do nn=1,norb
         do nnp=1,norb
            do nj=1,3
               do njp=1,3
                  shift_vector(nj,njp,nn,nnp)=-vme_der_phase(nj,njp,nn,nnp) &
                     +(realpart(berry_eigen(nj,nn,nn))-realpart(berry_eigen(nj,nnp,nnp)))
                  if (abs(shift_vector(nj,njp,nn,nnp)).gt.50.0d0) then
                     shift_vector(nj,njp,nn,nnp)=0.0d0
                  end if
                  vme_der(nj,njp,nn,nnp)=vme_der(nj,njp,nn,nnp) &
                     -complex(0.0d0,1.0d0)*vme_neigh(5,njp,nn,nnp) &
                     *(realpart(berry_eigen(nj,nn,nn))-realpart(berry_eigen(nj,nnp,nnp)))
               end do
            end do
         end do
      end do

   end subroutine get_berry_eigen_fourpoint
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_phase(aux1,ph)
      real*8 ph
      complex*16 aux1
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (abs(aux1).lt.10d-6) then
         ph=0.0d0
      else
         ph=aimag(log(aux1))
      end if
   end subroutine get_phase
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_gen_der_sumrule(norb,vme,e,abc,gen_der,gd1,gd2,gd3)
      implicit real*8 (a-h,o-z)

      integer norb,norb_inter_cut
      integer nn,nnp,nnpp
      integer nj,njp

      dimension e(norb)
      dimension vme(3,norb,norb)
      dimension gen_der(3,3,norb,norb)
      dimension gd1(3,3,norb,norb)
      dimension gd2(3,3,norb,norb)
      dimension gd3(3,3,norb,norb)
      dimension abc(3,norb,norb)

      real*8 e
      complex*16 vme,gen_der,abc,gd1,gd2,gd3,aux1,aux2,aux3
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      abc=0.0d0
      gd1=0.0d0
      gd2=0.0d0
      gd3=0.0d0
      gen_der=0.0d0
      norb_inter_cut=norb
      !write(*,*) 'computing berry connection, generalized derivative'
      do nn=1,norb
         do nnp=1,norb
            do nj=1,3
               !if (nn.eq.nnp) then
               if (abs(e(nn)-e(nnp)).lt.1.0d-05) then
                  abc(nj,nn,nnp)=0.0d0
               else
                  abc(nj,nn,nnp)=-complex(0.0d0,1.0d0)*vme(nj,nn,nnp)/(e(nn)-e(nnp))
                  !abc(nj,nnp,nn)=conjg(abc(nj,nn,nnp))
               end if
               do njp=1,3
                  !if (nn.eq.nnp) then
                  if (abs(e(nn)-e(nnp)).lt.1.0d-05) then
                     gd1(nj,njp,nn,nnp)=0.0d0
                     gd2(nj,njp,nn,nnp)=0.0d0
                     gd3(nj,njp,nn,nnp)=0.0d0
                  else
                     gd1(nj,njp,nn,nnp)=complex(0.0d0,1.0d0)/((e(nn)-e(nnp))**2)* &
                        (vme(nj,nn,nnp)*(vme(njp,nn,nn)-vme(njp,nnp,nnp))+ &
                        vme(njp,nn,nnp)*(vme(nj,nn,nn)-vme(nj,nnp,nnp)))

                     gd2(nj,njp,nn,nnp)=0.0d0
                     gd3(nj,njp,nn,nnp)=0.0d0
                     do nnpp=1,norb_inter_cut
                        !if (nnpp.eq.nn .or. nnpp.eq.nnp) then
                        if (abs(e(nnpp)-e(nn)).lt.1.0d-05 .or. &
                           abs(e(nnpp)-e(nnp)).lt.1.0d-05) then
                           gd2(nj,njp,nn,nnp)=gd2(nj,njp,nn,nnp)+0.0d0
                           gd3(nj,njp,nn,nnp)=gd3(nj,njp,nn,nnp)+0.0d0
                        else
                           gd2(nj,njp,nn,nnp)=gd2(nj,njp,nn,nnp)+ &
                              complex(0.0d0,1.0d0)/(e(nn)-e(nnp))* &
                              (vme(nj,nn,nnpp)*vme(njp,nnpp,nnp)/(e(nnpp)-e(nnp)))

                           gd3(nj,njp,nn,nnp)=gd3(nj,njp,nn,nnp)+ &
                              complex(0.0d0,1.0d0)/(e(nn)-e(nnp))* &
                              (-vme(njp,nn,nnpp)*vme(nj,nnpp,nnp)/(e(nn)-e(nnpp)))
                        end if
                     end do
                     !momentums and A and B term
                  end if
                  gen_der(nj,njp,nn,nnp)=gd1(nj,njp,nn,nnp)+gd2(nj,njp,nn,nnp) &
                     +gd3(nj,njp,nn,nnp)
                  !gen_der(nj,njp,nnp,nn)=conjg(gen_der(nj,njp,nn,nnp))
               end do
            end do
         end do
      end do

   end subroutine get_gen_der_sumrule

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_ome_sp_linear(iflag_norder,npointstotal,nband_ex,vme_ex_band,ek)
      implicit none
      integer iflag_norder
      integer npointstotal,nband_ex
      integer ibz
      integer nj,i,j

      dimension ek(npointstotal,nband_ex)
      dimension vme_ex_band(npointstotal,3,nband_ex,nband_ex)

      real*8 ek
      complex*16 vme_ex_band
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      open(10,file='ome_linear_sp_'//trim(material_name)//'.omesp')
      write(10,*) iflag_norder
      do ibz=1,npointstotal
         write(10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz),(ek(ibz,j),j=1,nband_ex)
         do i=1,nband_ex
            do j=1,nband_ex
               write(10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz), &
                  (realpart(vme_ex_band(ibz,nj,i,j)),aimag(vme_ex_band(ibz,nj,i,j)), nj=1,3)
            end do
         end do
      end do
      close(10)
   end subroutine write_ome_sp_linear
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_ome_sp_nonlinear(iflag_norder,npointstotal,nband_ex,vme_ex_band,ek,gen_der_ex_band, &
      shift_vector_ex_band,berry_eigen_ex_band)
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

      real*8 ek
      real*8 shift_vector_ex_band
      complex*16 vme_ex_band
      complex*16 berry_eigen_ex_band
      complex*16 gen_der_ex_band
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      open(10,file='ome_nonlinear_sp_'//trim(material_name)//'.omesp')
      write(10,*) iflag_norder
      do ibz=1,npointstotal
         write(10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz),(ek(ibz,j),j=1,nband_ex)
         do i=1,nband_ex
            do j=1,nband_ex
               write(10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz), &
                  (realpart(vme_ex_band(ibz,nj,i,j)),aimag(vme_ex_band(ibz,nj,i,j)), nj=1,3)
               write(10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz),(realpart(berry_eigen_ex_band(ibz,nj,i,j)), &
                  aimag(berry_eigen_ex_band(ibz,nj,i,j)), nj=1,3)
               write(10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz),(shift_vector_ex_band(ibz,1,nj,i,j), nj=1,3)
               write(10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz),(shift_vector_ex_band(ibz,2,nj,i,j), nj=1,3)
               write(10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz),(shift_vector_ex_band(ibz,3,nj,i,j), nj=1,3)
               write(10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz),(realpart(gen_der_ex_band(ibz,1,nj,i,j)), &
                  aimag(gen_der_ex_band(ibz,1,nj,i,j)), nj=1,3)
               write(10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz),(realpart(gen_der_ex_band(ibz,2,nj,i,j)), &
                  aimag(gen_der_ex_band(ibz,2,nj,i,j)), nj=1,3)
               write(10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz),(realpart(gen_der_ex_band(ibz,3,nj,i,j)), &
                  aimag(gen_der_ex_band(ibz,3,nj,i,j)), nj=1,3)
            end do
         end do
      end do
      close(10)
   end subroutine write_ome_sp_nonlinear
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! This routine evaluates the alpha,alpha' matrices
   subroutine get_vme_kernels_ome(rkx,rky,rkz,norb,skernel,sderkernel, &
      hkernel,hderkernel,akernel)
      implicit none

      integer norb
      integer ialpha
      integer ialphap
      integer iRp
      integer nj

      dimension skernel(norb,norb)
      dimension hkernel(norb,norb)
      dimension sderkernel(3,norb,norb)
      dimension hderkernel(3,norb,norb)
      dimension akernel(3,norb,norb)

      dimension hderhop(3,nR,norb,norb)
      dimension sderhop(3,nR,norb,norb)

      real(8) Rx,Ry,Rz
      real(8) rkx,rky,rkz


      complex*16 skernel,sderkernel,hkernel,hderkernel,akernel
      complex*16 phase,factor
      complex*16 hderhop,sderhop,rhop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      hkernel=0.0d0
      hderkernel=0.0d0
      skernel=0.0d0
      sderkernel=0.0d0
      akernel=0.0d0

      hderhop=0.0d0
      sderhop=0.0d0

      do ialpha=1,norb
         do ialphap=1,ialpha
            !write(*,*) akernel(1,1,1)
            !pause
            do iRp=1,nR
               Rx=dble(nRvec(iRp,1))*R(1,1)+dble(nRvec(iRp,2))*R(2,1)
               Ry=dble(nRvec(iRp,1))*R(1,2)+dble(nRvec(iRp,2))*R(2,2)
               Rz=0.0d0
               phase=complex(0.0d0,rkx*Rx+rky*Ry+rkz*Rz)
               factor=exp(phase)

               hkernel(ialpha,ialphap)=hkernel(ialpha,ialphap)+ &
                  factor*hhop(iRp,ialpha,ialphap)
               skernel(ialpha,ialphap)=skernel(ialpha,ialphap)+ &
                  factor*shop(iRp,ialpha,ialphap)

               hderhop(1,iRp,ialpha,ialphap)=complex(0.0d0,Rx)*hhop(iRp,ialpha,ialphap)
               hderhop(2,iRp,ialpha,ialphap)=complex(0.0d0,Ry)*hhop(iRp,ialpha,ialphap)
               hderhop(3,iRp,ialpha,ialphap)=complex(0.0d0,Rz)*hhop(iRp,ialpha,ialphap)
               sderhop(1,iRp,ialpha,ialphap)=complex(0.0d0,Rx)*shop(iRp,ialpha,ialphap)
               sderhop(2,iRp,ialpha,ialphap)=complex(0.0d0,Ry)*shop(iRp,ialpha,ialphap)
               sderhop(3,iRp,ialpha,ialphap)=complex(0.0d0,Rz)*shop(iRp,ialpha,ialphap)
               do nj=1,3
                  sderkernel(nj,ialpha,ialphap)=sderkernel(nj,ialpha,ialphap)+ &
                     factor*sderhop(nj,iRp,ialpha,ialphap)
                  hderkernel(nj,ialpha,ialphap)=hderkernel(nj,ialpha,ialphap)+ &
                     factor*hderhop(nj,iRp,ialpha,ialphap)
                  akernel(nj,ialpha,ialphap)=akernel(nj,ialpha,ialphap)+ &
                     factor*(rhop_c(nj,iRp,ialpha,ialphap)+ &
                     complex(0.0d0,1.0d0)*sderhop(nj,iRp,ialpha,ialphap))
               end do
               !if (ialpha.eq.1 .and. ialphap.eq.1) then
               !write(*,*) iRp,akernel(1,1,1),factor,rhop_c(nj,iRp,ialpha,ialphap), &
               !sderhop(nj,iRp,ialpha,ialphap)
               !end if
            end do
            !write(*,*) akernel(1,1,1)
            !pause
            !complex conjugate
            do nj=1,3
               hkernel(ialphap,ialpha)=conjg(hkernel(ialpha,ialphap))
               skernel(ialphap,ialpha)=conjg(skernel(ialpha,ialphap))
               sderkernel(nj,ialphap,ialpha)=conjg(sderkernel(nj,ialpha,ialphap))
               hderkernel(nj,ialphap,ialpha)=conjg(hderkernel(nj,ialpha,ialphap))
               akernel(nj,ialphap,ialpha)=conjg(akernel(nj,ialpha,ialphap))+ &
                  complex(0.0d0,1.0d0)*conjg(sderkernel(nj,ialpha,ialphap))
            end do

         end do
      end do
   end subroutine get_vme_kernels_ome


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! This routine evaluates the alpha,alpha' matrices
   subroutine get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
      hk_ev,e,vme)
      implicit none

      integer :: norb
      integer :: ialpha
      integer :: ialphap
      integer :: iRp
      integer :: nj
      integer :: i,j,ii,jj,nn,nnp

      dimension skernel(norb,norb)
      dimension hkernel(norb,norb)
      dimension sderkernel(3,norb,norb)
      dimension hderkernel(3,norb,norb)
      dimension akernel(3,norb,norb)

      dimension vjseudoa(3,norb,norb)
      dimension vjseudob(3,norb,norb)
      dimension e(norb)
      dimension hk_ev(norb,norb)
      dimension vme(3,norb,norb)

      real*8 e
      complex*16 skernel,sderkernel,hkernel,hderkernel,akernel
      complex*16 hk_ev,vjseudoa,vjseudob,vme
      complex*16 amu,amup,aux1,factor

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !diagonalization
      e=0.0d0
      call diagoz(norb,e,hkernel)
      hk_ev(:,:)=hkernel(:,:)
      !Multiply the eigenvectors (C_{\alpha=1 n}(k_0), C_{\alpha=1 n}(k_0), ...)
      !by a phase phi_n(k_0), this is new eigenvectors are
      !exp(-i*phi_n(k_0))*(C_{\alpha=1 n}(k_0), C_{\alpha=1 n}(k_0), ...)
      !
      !Note: right now this phase is used to give a locally smooth phase in the BZ
      call phase_eigvec_nk(norb,hk_ev)
      vme=0.0d0
      vjseudoa=0.0d0
      vjseudob=0.0d0

      do nn=1,norb
         do nnp=1,nn
            !momentums and A and B term
            do ialpha=1,norb
               do ialphap=1,norb
                  amu=hk_ev(ialpha,nn)
                  amup=hk_ev(ialphap,nnp)
                  do nj=1,3
                     vjseudoa(nj,nn,nnp)=vjseudoa(nj,nn,nnp)+ &
                        conjg(amu)*amup*hderkernel(nj,ialpha,ialphap)
                     vjseudob(nj,nn,nnp)=vjseudob(nj,nn,nnp)+conjg(amu)*amup* &
                        (e(nn)*akernel(nj,ialpha,ialphap)-e(nnp)*conjg(akernel(nj,ialphap,ialpha)))* &
                        complex(0.0d0,1.0d0)
                  end do
               end do
            end do
            !call cpu_time(time2)
            !write(*,*) 'k-sampling time',norb,time2-time1,'s'
            !pause
            do nj=1,3
               vme(nj,nn,nnp)=vjseudoa(nj,nn,nnp)+vjseudob(nj,nn,nnp)
               vme(nj,nnp,nn)=conjg(vme(nj,nn,nnp))
               vjseudoa(nj,nnp,nn)=conjg(vjseudoa(nj,nn,nnp))
               vjseudob(nj,nnp,nn)=conjg(vjseudob(nj,nn,nnp))
            end do
         end do
      end do

   end subroutine get_vme_eigen_ome

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine phase_eigvec_nk(norb,hk_ev)
      implicit none

      integer norb
      integer i,j,ii

      dimension hk_ev(norb,norb)

      real*8 :: arg
      complex*16 :: aux1,hk_ev,factor
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !phase election: this is done to smooth the gauge
      do j=1,norb
         aux1=0.0d0
         do i=1,norb
            aux1=aux1+hk_ev(i,j)
         end do
         !argument of the sym
         arg=atan2(aimag(aux1),realpart(aux1))
         factor=exp(complex(0.0d0,-arg))
         !write(*,*) 'sum is now:',aux1*factor
         do ii=1,norb
            hk_ev(ii,j)=hk_ev(ii,j)*factor
            !hk_ev(ii,j)=hk_ev(ii,j)*1.0d0
         end do
      end do

   end subroutine
end module ome_sp
