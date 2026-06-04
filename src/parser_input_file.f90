module parser_input_file
  implicit none
  private
  public :: material_name_in
  public :: filename_input
  public :: xatu_eigval_filepath_in
  public :: xatu_states_filepath_in
  public :: iflag_xatu_text
  public :: iflag_ome_sp_text
  public :: iflag_ome_ex_text
  public :: response_text
  public :: iflag_xatu
  public :: iflag_ome_sp
  public :: iflag_ome_ex
  public :: ndim,nf,npointstotal_sq
  public :: e1,e2,eta,nw
  public :: get_input_file
  public :: nband_index
  public :: norb_ex_cut
  public :: broadening_type_text
  public :: read_line_numbers_int !subroutine

  character(len=1000) :: material_name_in
  character(len=100) :: filename_input
  character(len=100) :: iflag_xatu_text
  character(len=100) :: iflag_ome_sp_text
  character(len=100) :: iflag_ome_ex_text
  character(len=100) :: broadening_type_text
  character(len=1000) :: xatu_eigval_filepath_in
  character(len=1000) :: xatu_states_filepath_in
  character(len=100) :: response_text

  logical :: iflag_xatu
  logical :: iflag_ome_sp
  logical :: iflag_ome_ex

  integer :: ndim
  integer :: nf
  integer :: npointstotal_sq
  integer :: norb_ex_cut
  integer :: nband_index
  integer :: nw
  real(8) :: e1,e2,eta

  allocatable :: nband_index(:)

  contains
    function to_lower(str) result(lower_str)
      implicit none
      character(len=*), intent(in) :: str
      character(len=len(str)) :: lower_str
      integer :: i, ic
      
      do i = 1, len(str)
        ic = iachar(str(i:i))
        if (ic >= iachar('A') .and. ic <= iachar('Z')) then
          lower_str(i:i) = achar(ic + 32)
        else
          lower_str(i:i) = str(i:i)
        end if
      end do
    end function to_lower
    
    subroutine get_input_file()
      implicit none
      integer, allocatable :: narray(:) 
      integer :: num_values, ios
      character(len=1000) :: line
      character(len=100) :: param_name
      logical :: ndim_found, material_found, xatu_found, bandlist_found
      logical :: ncells_found, nfermi_found, ome_sp_found, ome_ex_found
      logical :: response_found, energy_found, exciton_found
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,*) '1. Entering parser_input_file'
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      ! Initialize flags
      ndim_found = .false.
      material_found = .false.
      xatu_found = .false.
      bandlist_found = .false.
      ncells_found = .false.
      nfermi_found = .false.
      ome_sp_found = .false.
      ome_ex_found = .false.
      response_found = .false.
      energy_found = .false.
      exciton_found = .false.
      ! default broadening
      broadening_type_text = 'gaussian'
      
      call get_command_argument(1,filename_input)
      open(10,file=adjustl(filename_input))
      
      ! Read file sequentially and process parameters based on their labels
      do
        read(10,'(A)',iostat=ios) line
        if (ios /= 0) exit  ! End of file
        
        line = adjustl(line)
        
        ! Check if this is a comment line (parameter label)
        if (line(1:1) == '#') then
          ! Extract parameter name and read corresponding value
          param_name = adjustl(line(3:))  ! Remove "# " prefix
          
          if (index(param_name, 'Periodic dimensions') > 0) then
            read(10,*) ndim
            ndim_found = .true.
            
          else if (index(param_name, 'Wannier90_filename') > 0) then
            read(10,'(A)') material_name_in
            material_found = .true.
            
          else if (index(param_name, 'Xatu_interface') > 0) then
            read(10,*) iflag_xatu_text
            xatu_found = .true.
            
            if (iflag_xatu_text == 'true') then
              iflag_xatu = .true.
              ! Read the eigval and states file paths that follow
              read(10,'(A)') xatu_eigval_filepath_in
              read(10,'(A)') xatu_states_filepath_in
            else if (iflag_xatu_text == 'false') then
              iflag_xatu = .false.
            else
              write(*,*) 'Error: Invalid value in Xatu_interface. Expected "true" or "false".'
              stop
            end if
            
          else if (index(param_name, 'Exciton_cutoff') > 0) then
            read(10,*) norb_ex_cut
            exciton_found = .true.
            
          else if (index(param_name, 'Bandlist') > 0) then
            call read_line_numbers_int(narray, num_values)
            bandlist_found = .true.
            
          else if (index(param_name, 'Ncells') > 0) then
            read(10,*) npointstotal_sq
            ncells_found = .true.
            
          else if (index(param_name, 'Nfermi') > 0) then
            read(10,*) nf
            nfermi_found = .true.
            
          else if (index(param_name, 'OME_sp') > 0 .or. index(param_name, 'OME_SP') > 0) then
            read(10,*) iflag_ome_sp_text
            ome_sp_found = .true.
            
          else if (index(param_name, 'OME_ex') > 0 .or. index(param_name, 'OME_EX') > 0) then
            read(10,*) iflag_ome_ex_text
            ome_ex_found = .true.
            
          else if (index(param_name, 'Response') > 0) then
            read(10,*) response_text
            response_found = .true.
            
          else if (index(param_name, 'Energy_variables') > 0) then
            read(10,*) e1, e2, eta, nw
            energy_found = .true.
          else if (index(param_name, 'Broadening_type') > 0 .or. index(param_name,'Broadening')>0) then
            read(10,'(A)') broadening_type_text
            broadening_type_text = adjustl(broadening_type_text)
            broadening_type_text = to_lower(broadening_type_text)
          
          end if
        end if
      end do
      
      close(10)
      
      ! Handle bandlist case: allocate nband_index if bandlist was found
      if (bandlist_found) then
        allocate(nband_index(num_values))
        nband_index(:) = narray(:)
      end if
      
      ! Set npointstotal_sq to 0 if using xatu interface
      if (iflag_xatu) then
        npointstotal_sq = 0
      end if
      
      ! Declare flags from text strings
      if (iflag_ome_sp_text == 'true') then
        iflag_ome_sp = .true.
      else
        iflag_ome_sp = .false.
      end if
      
      if (iflag_ome_ex_text == 'true') then
        iflag_ome_ex = .true.
      else
        iflag_ome_ex = .false.
      end if
      
      write(*,*) '   Input file has been read'
    end subroutine get_input_file
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !This routine reads a line of numbers into an array
    subroutine read_line_numbers_int(narray,num_values)
    implicit none
    allocatable :: narray(:)

    integer :: ios,ncount,i
    integer :: narray
    integer :: num_values
    integer :: temp_num 
    integer :: istart,iposition
    character(len=1000) :: line
 
    !Read the line of text
    read(10,'(A)', iostat=ios) line
    if (ios /= 0) then
      print *, 'Error reading file'
    stop
    end if

    ncount=0
    do i=1,len_trim(line)
      if (line(i:i) == ' ') ncount=ncount+1
    end do
    ncount=ncount+1  ! One more than the number of spaces

    !Allocate the array based on the number of values
    allocate(narray(ncount))

    !Reset the number of values counter
    num_values=0
    istart=1
    !Now, sequentially extract numbers from the line
    do i=1,ncount
      !Find the next number in the line
      read(line(istart:),*,iostat=ios) temp_num
      if (ios == 0) then
        num_values=num_values+1
        narray(num_values)=temp_num
        !Move the starting  position to the next number
        iposition=scan(line(istart:),' ')
        if (iposition>0) then
          istart=istart+iposition
        end if
      end if
    end do

    end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module parser_input_file