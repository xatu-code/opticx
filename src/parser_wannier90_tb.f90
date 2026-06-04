module parser_wannier90_tb
   use parser_input_file, only:read_line_numbers_int
   implicit none
   private
   public :: norb,nR,n1,n2,n3,nRvec
   public :: R,shop,hhop,rhop_c
   public :: material_name
   public :: wannier90_get

   integer norb
   integer nR
   integer n1,n2,n3
   integer nRvec

   real(8) R,shop
   complex*16 hhop,rhop_c
   character(1000) material_name

   dimension   R(3,3)
   integer, allocatable :: Degen(:)
   allocatable nRvec(:,:)
   allocatable hhop(:,:,:)
   allocatable rhop_c(:,:,:,:)
   allocatable shop(:,:,:)

contains

!--------------------------------------------------------------------
   pure function to_lower(str) result(out)
!> Convert a string to lower‑case (portable; no compiler extension)
      character(len=*), intent(in) :: str
      character(len=len(str))      :: out
      integer :: i
      out = str
      do i = 1, len(str)
         if (out(i:i) >= 'A' .and. out(i:i) <= 'Z') &
            out(i:i) = char(iachar(out(i:i)) + 32)
         end do
      end function to_lower
      !--------------------------------------------------------------------
      subroutine wannier90_get(material_name_in)
         implicit none
         ! -----------------------------------------------------------------
         character(len=*), intent(in) :: material_name_in    ! full path read from input.txt
         ! -----------------------------------------------------------------
         integer          :: fp, iR, i, ialpha, ialphap
         integer          :: nkk1, nkk2, nRzero
         real(8)          :: a1,a2,a3,a4,a5,a6
         character(len=:), allocatable :: file2open, basename
         integer          :: p, ext_pos
         integer :: num_chunks
         integer :: rem
         ! -----------------------------------------------------------------
         write(*,*) '2. Entering parser_wannier90_tb'

      ! === 1.  Use the path exactly as supplied ========================
      !write(*,*) "MATERIAL NAME PARSED:", material_name_in
      file2open = trim(material_name_in)

      ! === 2.  Derive clean material name (no dir, no _tb.dat) =========
      p = max( index(file2open,'/',back=.true.),  &
         index(file2open,'\',back=.true.) )   ! works on Win/Linux
      if (p == 0) then
         basename = file2open
      else
         basename = file2open(p+1:)
      end if

      ext_pos = len_trim(basename) - len('_tb.dat') + 1
      if ( ext_pos > 0 .and. to_lower(basename(ext_pos:)) == '_tb.dat' ) then
         basename = basename(:ext_pos-1)
      end if
      material_name = adjustl(basename)

      ! === 3.  Open Wannier90 TB file ==================================
      open(unit=fp, file=file2open, action='read', status='old')
      read(fp,*)
      read(fp,*) R(1,1),R(1,2),R(1,3)
      read(fp,*) R(2,1),R(2,2),R(2,3)
      read(fp,*) R(3,1),R(3,2),R(3,3)
      read(fp,*) norb
      read(fp,*) nR

      !allocate nR, h-,r-, and s-hoppings
      allocate (nRvec(nR,3))
      allocate (hhop(nR,norb,norb))
      allocate (shop(nR,norb,norb))
      allocate (rhop_c(3,nR,norb,norb))
      allocate (Degen(nR))

      num_chunks = nR / 15
      rem = MOD(nR, 15)
      do i = 1, num_chunks
         read(fp, *) Degen((i - 1) * 15 + 1:(i - 1) * 15 + 15)
      end do
      if (rem > 0) then
         read(fp, *) Degen(num_chunks * 15 + 1:num_chunks * 15 + rem)
      end if

      !get the hopping matrices
      do iR=1,nR
         read(fp,*) nRvec(iR,:)
         do ialphap=1,norb
            do ialpha=1,norb
               read(fp,*) nkk1,nkk2,a1,a2
               hhop(iR,nkk1,nkk2)=complex(a1,a2)
            end do
         end do
         read(fp,*)
      end do

      !get rhoppings
      do iR=1,nR
         read(fp,*) !nRvec is already strored
         do ialphap=1,norb
            do ialpha=1,norb
               read(fp,*) nkk1,nkk2,a1,a2,a3,a4,a5,a6
               rhop_c(1,iR,nkk1,nkk2)=complex(a1,a2)
               rhop_c(2,iR,nkk1,nkk2)=complex(a3,a4)
               rhop_c(3,iR,nkk1,nkk2)=complex(a5,a6)
            end do
         end do
         if (iR /= nR) read(fp,*) !blank line
      end do
      close(fp)
      !get orthogonal overlap: this variable is a reminiscent
      !of the interface with the original crystal interface.
      !I maintain the overlap matrix though

      !locate the (0,0,0) element of nRvec
      do iR=1,nR
         if (nRvec(iR,1)==0 .and. nRvec(iR,2)==0 .and. nRvec(iR,3)==0) then
            nRzero=iR
            exit
         end if
      end do
      !wannier functions are orthonormal
      shop=0.0d0
      do ialpha=1,norb
         shop(nRzero,ialpha,ialpha)=1.0d0
      end do


      !APPLY BIAS BY HAND
      !do iR=1,nR
      !do ialpha=1,norb
      !do ialphap=1,norb
      !hhop(iR,ialpha,ialphap)=hhop(iR,ialpha,ialphap)-0.02d0*rhop_c(3,iR,ialpha,ialphap)
      !end do
      !end do
      !end do
      !do ialpha=1,norb
      !hhop(nRzero,ialpha,ialpha)=hhop(nRzero,ialpha,ialpha)-0.1d0*rhop_c(3,nRzero,ialpha,ialpha)
      !end do

      !convert units: to Hartree and bohrs
      hhop=hhop/27.211385d0
      rhop_c=rhop_c/0.52917721067121d0
      R=R/0.52917721067121d0
      write(*,*) '   Wannier hamiltonian has been read'
   end subroutine wannier90_get


end module parser_wannier90_tb
