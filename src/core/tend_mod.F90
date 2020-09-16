module tend_mod

  use flogger
  use const_mod
  use mesh_mod
  use namelist_mod
  use allocator_mod

  implicit none

  private

  public tend_type
  public tends
  public tend_init_root

  type tend_type
    type(mesh_type), pointer :: mesh => null()
    real(r8), allocatable, dimension(:,:) :: du
    real(r8), allocatable, dimension(:,:) :: dv
    real(r8), allocatable, dimension(:,:) :: dgd
    ! Individual tendencies
    real(r8), allocatable, dimension(:,:) :: u_adv_lon
    real(r8), allocatable, dimension(:,:) :: u_adv_lat
    real(r8), allocatable, dimension(:,:) :: fv
    real(r8), allocatable, dimension(:,:) :: u_pgf
    real(r8), allocatable, dimension(:,:) :: v_adv_lon
    real(r8), allocatable, dimension(:,:) :: v_adv_lat
    real(r8), allocatable, dimension(:,:) :: v_sinuv
    real(r8), allocatable, dimension(:,:) :: fu
    real(r8), allocatable, dimension(:,:) :: v_pgf
    real(r8), allocatable, dimension(:,:) :: dmfdlon
    real(r8), allocatable, dimension(:,:) :: dmfdlat
  contains
    procedure :: init => tend_init
    procedure :: clear => tend_clear
    final :: tend_final
  end type tend_type
  
  type(tend_type), allocatable, target :: tends(:)

contains 

  subroutine tend_init_root()
    
    integer i

    if (.not. allocated(tends)) then
      select case (trim(time_scheme))
      case ('pc2')
        allocate(tends(3))
      case ('rk3')
        allocate(tends(4))
      case ('rk4')
        allocate(tends(5))
      case default
        call log_error('Wrong time scheme!')
      end select
      do i = lbound(tends, 1), ubound(tends, 1)
        call tends(i)%init(mesh)
      end do
      call log_notice('tend module is initialized.')
    end if

  end subroutine tend_init_root

  subroutine tend_init(this, mesh)
    
    class(tend_type), intent(inout)         :: this
    type(mesh_type ), intent(in   ), target :: mesh

    call this%clear()

    this%mesh => mesh

    call allocate_array(mesh, this%du        , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%dv        , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%dgd       , full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%u_adv_lon , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%u_adv_lat , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%fv        , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%u_pgf     , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%v_adv_lon , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%v_adv_lat , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%v_sinuv   , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%fu        , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%v_pgf     , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%dmfdlon   , full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%dmfdlat   , full_lon=.true., full_lat=.true.)
  end subroutine tend_init

  subroutine tend_clear(this)
    
    class(tend_type), intent(inout) :: this

    if (allocated(this%du        )) deallocate(this%du        )
    if (allocated(this%dv        )) deallocate(this%dv        )
    if (allocated(this%dgd       )) deallocate(this%dgd       )
    if (allocated(this%u_adv_lon )) deallocate(this%u_adv_lon )
    if (allocated(this%u_adv_lat )) deallocate(this%u_adv_lat )
    if (allocated(this%fv        )) deallocate(this%fv        )
    if (allocated(this%u_pgf     )) deallocate(this%u_pgf     )
    if (allocated(this%v_adv_lon )) deallocate(this%v_adv_lon )
    if (allocated(this%v_adv_lat )) deallocate(this%v_adv_lat )
    if (allocated(this%v_sinuv   )) deallocate(this%v_sinuv   )
    if (allocated(this%fu        )) deallocate(this%fu        )
    if (allocated(this%v_pgf     )) deallocate(this%v_pgf     )
    if (allocated(this%dmfdlon   )) deallocate(this%dmfdlon   )
    if (allocated(this%dmfdlat   )) deallocate(this%dmfdlat   )

  end subroutine tend_clear

  subroutine tend_final(this)
    
    type(tend_type), intent(inout) :: this

    call this%clear()

  end subroutine tend_final
end module tend_mod
