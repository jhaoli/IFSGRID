module operators_mod

  use flogger
  use const_mod
  use namelist_mod
  use mesh_mod
  use static_mod
  use state_mod
  use tend_mod
  use parallel_mod

  implicit none

  private

  public operators_prepare
  public calc_uv_adv
  public calc_v_sinuv
  public calc_fv_fu
  public calc_uv_pgf
  public calc_dmfdlon_dmfdlat

contains

  subroutine operators_prepare(state)

    type(state_type), intent(inout) :: state

    call calc_m_lon_m_lat(state)
    call calc_m_vtx(state)
    call calc_vor(state)
    call calc_mf_lon_n_mf_lat_n(state)
    call calc_u_lat_v_lon(state)

  end subroutine operators_prepare

  subroutine calc_m_lon_m_lat(state)
    
    type(state_type), intent(inout) :: state
    type(mesh_type), pointer :: mesh
    integer i, j

    mesh => state%mesh

    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        state%m_lon(i,j) = (state%gd(i,j) + state%gd(i+1,j)) * 0.5_r8
      end do
    end do 

    do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
#ifdef V_POLE
        state%m_lat(i,j) = (state%gd(i,j) + state%gd(i,j-1)) * 0.5_r8
#else
        state%m_lat(i,j) = (state%gd(i,j) + state%gd(i,j+1)) * 0.5_r8
#endif
      end do
    end do

  end subroutine calc_m_lon_m_lat
 
  subroutine calc_m_vtx(state)
    
    type(state_type), intent(inout) :: state
    type(mesh_type), pointer :: mesh

    integer i, j
    real(r8) pole

    mesh => state%mesh

    do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
#ifdef V_POLE
        state%m_vtx(i,j) = (state%gd(i,j  ) + state%gd(i+1,j  ) + &
                            state%gd(i,j-1) + state%gd(i+1,j-1)) * 0.25_r8
#else
        state%m_vtx(i,j) = (state%gd(i,j  ) + state%gd(i+1,j  ) + &
                            state%gd(i,j+1) + state%gd(i+1,j+1)) * 0.25_r8 
#endif
      end do
    end do
#ifdef V_POLE
    if (mesh%has_south_pole()) then
      j = mesh%half_lat_start_idx
      pole = 0.0_r8
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        pole = pole + state%gd(i,j)
      end do
      pole = pole / mesh%num_half_lon
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        state%m_vtx(i,j) = pole
      end do
    endif
    if (mesh%has_north_pole()) then
      j = mesh%half_lat_end_idx
      pole = 0.0_r8
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        pole = pole + state%gd(i,j-1)
      end do
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        state%m_vtx(i,j) = pole
      end do
    end if
#endif

  end subroutine calc_m_vtx
  
  subroutine calc_vor(state)
    
    type(state_type), intent(inout) :: state
    type(mesh_type), pointer :: mesh
    integer i, j
    real(r8) pole

    mesh => state%mesh

    do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
#ifdef V_POLE
        state%vor(i,j) = ((state%v(i+1,j) - state%v(i,j  )) / mesh%dlon / mesh%half_cos_lat(j) - &
                          (state%u(i,  j) - state%u(i,j-1)) / mesh%dlat) / (radius * mesh%half_cos_lat(j))
#else
        state%vor(i,j) = ((state%v(i+1,j) - state%v(i,j)) / mesh%dlon / mesh%half_cos_lat(j) -&
                          (state%u(i,j+1) - state%u(i,j)) /mesh%dlat) / (radius * mesh%half_cos_lat(j))
#endif
      end do
    end do
#ifdef V_POLE
  if (mesh%has_south_pole()) then
    j = mesh%half_lat_start_idx
    pole = 0.0_r8
    do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
      pole = pole - state%u(i,j)
    end do
    pole = pole / mesh%num_half_lon / (radius * mesh%dlat * 0.5_r8)
    do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
      state%vor(i,j) = pole
    end do
  end if
  if (mesh%has_north_pole()) then
    j = mesh%half_lat_end_idx
    pole = 0.0_r8
    do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
      pole = pole + state%u(i,j-1)
    end do
    pole = pole / mesh%num_half_lon / (radius * mesh%dlat * 0.5_r8)
    do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
      state%vor(i,j) = pole
    end do
  end if
#endif

  end subroutine calc_vor

  subroutine calc_mf_lon_n_mf_lat_n(state)
  
    type(state_type), intent(inout) :: state
    type(mesh_type), pointer :: mesh
    integer i, j
    
    mesh => state%mesh

    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        state%mf_lon_n(i,j) = state%m_lon(i,j) * state%u(i,j) / mesh%full_cos_lat(j)
      end do
    end do
    call parallel_fill_halo(state%mesh, state%mf_lon_n)

    do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        state%mf_lat_n(i,j) = state%m_lat(i,j) * state%v(i,j) / mesh%half_cos_lat(j)
      end do
    end do
    call parallel_fill_halo(state%mesh, state%mf_lat_n)

  end subroutine calc_mf_lon_n_mf_lat_n
  
  subroutine calc_u_lat_v_lon(state)

    type(state_type), intent(inout) :: state
    type(mesh_type), pointer :: mesh
    integer i, j

    mesh => state%mesh

    do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
#ifdef V_POLE
      state%u_lat(i,j) = (state%u(i-1,j) + state%u(i,j) + state%u(i-1,j-1) + state%u(i,j-1)) * 0.25_r8
#else
      state%u_lat(i,j) = (state%u(i-1,j) + state%u(i,j) + state%u(i-1,j+1) + state%u(i,j+1)) * 0.25_r8 
#endif
      end do
    end do
    
    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
#ifdef V_POLE
        state%v_lon(i,j) = ((state%v(i,j  ) + state%v(i+1,j  )) * mesh%half_cos_lat(j  ) + &
                            (state%v(i,j+1) + state%v(i+1,j+1)) * mesh%half_cos_lat(j+1)) * 0.25_r8 / mesh%full_cos_lat(j)

#else
        state%v_lon(i,j) = ((state%v(i,j  ) + state%v(i+1,j  )) * mesh%half_cos_lat(j  ) + &
                            (state%v(i,j-1) + state%v(i+1,j-1)) * mesh%half_cos_lat(j-1)) * 0.25_r8 / mesh%full_cos_lat(j)
#endif
      end do
    end do

  end subroutine calc_u_lat_v_lon
 
  subroutine calc_uv_adv(state, tend, dt)

    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend
    real(r8), intent(in) :: dt

    type(mesh_type), pointer :: mesh
    integer i, j, ip

    mesh => state%mesh

    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        tend%u_adv_lon(i,j) = state%u(i,j) * (state%u(i+1,j) - state%u(i-1,j)) / &
                         (2.0_r8 * mesh%dlon) / radius / mesh%full_cos_lat(j)**2
      end do
    end do

    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
#ifdef V_POLE
        ! if (j == mesh%full_lat_start_idx_no_pole) then
        !   tend%u_adv_lat(i,j) = state%v_lon(i,j) * mesh%full_cos_lat(j) * &
        !                         (state%u(i,j+1) - state%u(i,j)) / mesh%dlat / &
        !                         radius / mesh%full_cos_lat(j)**2
        ! else if (j == mesh%full_lat_end_idx_no_pole) then
        !   tend%u_adv_lat(i,j) = state%v_lon(i,j) * mesh%full_cos_lat(j) * &
        !                         (state%u(i,j) - state%u(i,j-1)) / mesh%dlat / &
        !                         radius / mesh%full_cos_lat(j)**2
        ! else 
        ip = i + mesh%num_half_lon / 2
        if (ip > mesh%num_half_lon / 2) then
          ip = ip - mesh%num_half_lon /2
        end if
        if (j == mesh%full_lat_start_idx_no_pole) then
          state%u(i,j-1) = state%u(ip,j)
        else if (j == mesh%full_lat_end_idx_no_pole) then 
          state%u(i,j+1) = state%u(ip,j)
        end if
          tend%u_adv_lat(i,j) = state%v_lon(i,j) * mesh%full_cos_lat(j) * &
                                (state%u(i,j+1) - state%u(i,j-1)) / (2.0_r8 * mesh%dlat) / &
                                radius / mesh%full_cos_lat(j)**2
        ! end if
#else
        tend%u_adv_lat(i,j) = state%v_lon(i,j) * mesh%full_cos_lat(j) * &
                              (state%u(i,j+1) - state%u(i,j-1)) / (2.0_r8 * mesh%dlat) / &
                              radius / mesh%full_cos_lat(j)**2
#endif
      end do
    end do

    do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        tend%v_adv_lon(i,j) = state%u_lat(i,j) * &
                              (state%v(i+1,j) - state%v(i-1,j)) / (2.0_r8 * mesh%dlon) / &
                              radius / mesh%half_cos_lat(j)**2
      end do
    end do

    do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
#ifdef V_POLE
        tend%v_adv_lat(i,j) = state%v(i,j) * mesh%half_cos_lat(j) * &
                              (state%v(i,j+1) - state%v(i,j-1)) / (2.0_r8 * mesh%dlat) / &
                              radius / mesh%half_cos_lat(j)**2
#else
        ! if (j == mesh%half_lat_start_idx_no_pole) then
        !   tend%v_adv_lat(i,j) = state%v(i,j) * mesh%half_cos_lat(j) * (state%v(i,j+1) - state%v(i,j)) / &
        !                        mesh%dlat / radius / mesh%half_cos_lat(j)**2
        ! else if (j == mesh%half_lat_end_idx_no_pole) then
        !   tend%v_adv_lat(i,j) = state%v(i,j) * mesh%half_cos_lat(j) * (state%v(i,j) - state%v(i,j-1)) / &
        !                        mesh%dlat / radius / mesh%half_cos_lat(j)**2
        ! else
!        ip = i + mesh%num_full_lon
!        if (ip > mesh%num_full_lon / 2) then
!          ip = ip - mesh%num_full_lon / 2
!        end if
!        if (j == mesh%half_lat_start_idx_no_pole) then 
!          state%v(i,j-1) = state%v(ip,j)
!        else if (j == mesh%half_lat_end_idx_no_pole) then 
!          state%v(i,j+1) = state%v(ip,j)
!        end if
          tend%v_adv_lat(i,j) = state%v(i,j) * mesh%half_cos_lat(j) * &
                                (state%v(i,j+1) - state%v(i,j-1)) / (2.0_r8 * mesh%dlat) / &
                                radius / mesh%half_cos_lat(j)**2
        ! end if
#endif
      end do
    end do

  end subroutine calc_uv_adv
  
  subroutine calc_v_sinuv(state, tend, dt)
  
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    real(r8), intent(in) :: dt

    type(mesh_type), pointer :: mesh
    integer i, j

    mesh => state%mesh

    do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        tend%v_sinuv(i,j) = mesh%half_sin_lat(j) * (state%u_lat(i,j)**2 + state%v(i,j)**2) / &
                                                    radius / mesh%half_cos_lat(j)**2
      end do
    end do
  
  end subroutine calc_v_sinuv

  subroutine calc_fv_fu(state, tend, dt)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    real(r8), intent(in) :: dt

    type(mesh_type), pointer :: mesh
    integer i, j

    mesh => state%mesh

    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        tend%fv(i,j) = mesh%full_f(j) * state%v_lon(i,j)    
      end do
    end do

    do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
#ifdef V_POLE
       tend%fu(i,j) = (mesh%full_f(j  ) * (state%u(i-1,j  ) + state%u(i,j   )) + &
                       mesh%full_f(j-1) * (state%u(i-1,j-1) + state%u(i,j-1 ))) * 0.25_r8
        ! tend%fu(i,j) = mesh%half_f(j) * state%u_lat(i,j)
#else
       tend%fu(i,j) = (mesh%full_f(j  ) * (state%u(i-1,j  ) + state%u(i,j   )) + &
                       mesh%full_f(j+1) * (state%u(i-1,j+1) + state%u(i,j+1 ))) * 0.25_r8
#endif
        ! tend%fu(i,j) = mesh%half_f(j) * state%u_lat(i,j)
      end do
    end do

  end subroutine calc_fv_fu

  subroutine calc_uv_pgf(state, static, tend, dt)

    type(state_type), intent(in) :: state
    type(static_type), intent(in) :: static
    type(tend_type), intent(inout) :: tend
    real(r8), intent(in) :: dt

    type(mesh_type), pointer :: mesh
    integer i, j

    mesh => state%mesh

    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        tend%u_pgf(i,j) = (state%gd(i+1,j) + static%ghs(i+1,j) - &
                           state%gd(i  ,j) - static%ghs(i  ,j)) / mesh%dlon / radius
      end do
    end do

    do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
#ifdef V_POLE
        tend%v_pgf(i,j) = mesh%half_cos_lat(j) * &
                         (state%gd(i,j  ) + static%ghs(i,j  ) - &
                          state%gd(i,j-1) - static%ghs(i,j-1)) / mesh%dlat / radius
#else
        tend%v_pgf(i,j) = mesh%half_cos_lat(j) * &
                         (state%gd(i,j+1) + static%ghs(i,j+1) - &
                          state%gd(i,j  ) - static%ghs(i,j  )) / mesh%dlat / radius
#endif
      end do
    end do

  end subroutine calc_uv_pgf

  subroutine calc_dmfdlon_dmfdlat(state, tend, dt)
  
    type(state_type), intent(inout) :: state
    type(tend_type ), intent(inout) :: tend
    real(r8)        , intent(in   ) :: dt 

    type(mesh_type), pointer :: mesh
    integer i, j
    real(r8) pole

    mesh => state%mesh

    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        tend%dmfdlon(i,j) = (state%mf_lon_n(i,j) - state%mf_lon_n(i-1,j)) / &
                            radius / mesh%full_cos_lat(j) / mesh%dlon
      end do
    end do

    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
#ifdef V_POLE
        tend%dmfdlat(i,j) = (state%mf_lat_n(i,j+1) * mesh%half_cos_lat(j+1) - &
                             state%mf_lat_n(i,j  ) * mesh%half_cos_lat(j  )) /& 
                             radius / mesh%full_cos_lat(j) / mesh%dlat
#else
        tend%dmfdlat(i,j) = (state%mf_lat_n(i,j  ) * mesh%half_cos_lat(j  ) - &
                             state%mf_lat_n(i,j-1) * mesh%half_cos_lat(j-1)) /&
                             radius / mesh%full_cos_lat(j) / mesh%dlat
#endif
      end do
    end do
#ifndef V_POLE
    if (mesh%has_south_pole()) then
      j = mesh%full_lat_start_idx
      pole = 0.0_r8
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        pole = pole + state%mf_lat_n(i,j)
      end do
      pole = pole * 4.0_r8 / mesh%num_full_lon / (radius * mesh%dlat)
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        tend%dmfdlat(i,j) = pole
      end do
    end if

    if (mesh%has_north_pole()) then
      j = mesh%full_lat_end_idx
      pole = 0.0_r8
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        pole = pole + state%mf_lat_n(i,j-1)
      end do
      pole = -pole * 4.0_r8 / mesh%num_full_lon / (radius * mesh%dlat)
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        tend%dmfdlat(i,j) = pole
      end do
    end if
#endif
  end subroutine calc_dmfdlon_dmfdlat

end module operators_mod
