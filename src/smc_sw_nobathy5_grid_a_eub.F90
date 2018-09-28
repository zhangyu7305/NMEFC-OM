program smc_sw_nopolar_nobathy
  use netcdf
  implicit none
  real(kind=8), parameter :: dt = 10.0d0
  real(kind=8), parameter :: pie=3.141592654d0, grav = 9.8d0, almu = 100000.0d0, avm = 1.0d-5, rho0 = 1035.0d0
  real(kind=8), parameter :: eps = 1.0d-8, rearth=6.37122d6, omega = 7.2921d-5, gama = 0.3d0
  integer(kind=4), parameter :: ta = 10000
  integer(kind=4), parameter :: neig = 2
  integer(kind=4), parameter :: rn_shlat = 0
  integer(kind=4) :: status, dimid, varid
  integer(kind=4) :: ncid_geo, ncid_sbc, ncid_ini, ncid_out
  integer(kind=4) :: nglo_b_dimid, t_dimid, dimids2d(2)
  integer(kind=4) :: lon_id, lat_id, time_id, h_id, u_id, v_id, vor_id, div_id
  integer(kind=4) :: i, j, j1, j2, t, tout, idx, idx_b, i_poln, i_pols, i_st, i_en
  integer(kind=4) :: nglo_b, nbdy_b
  integer(kind=4) :: tn_1, tn_2, ts_1, ts_2, tw_1, tw_2, te_1, te_2
  integer(kind=4) :: clock_count1, clock_count2, clock_rate, clock_max
  integer(kind=4), allocatable :: tn(:,:), ts(:,:), tw(:,:), te(:,:), idx_b_m(:), idx_m_b(:)
  character(len=nf90_max_name) :: fn_out, fn_grid, fn_ini, fn_sbc
  real(kind=8) :: time(1)
  real(kind=8) :: rad, nl_temp, nl_temp1, nl_temp2
  real(kind=8) :: u_n,u_s,u_w,u_e,u_c
  real(kind=8) :: v_n,v_s,v_w,v_e,v_c
  real(kind=8) :: h_n,h_s,h_w,h_e,h_c
  real(kind=8) :: b_n,b_s,b_w,b_e,b_c
  real(kind=8) :: u_ne,u_se,u_we,u_ee
  real(kind=8) :: v_ne,v_se,v_we,v_ee
  real(kind=8) :: h_ne,h_se,h_we,h_ee
  real(kind=8) :: t_n,t_s,t_w,t_e
  real(kind=8) :: ke_n,ke_s,ke_w,ke_e
  real(kind=8) :: u_f_pol, v_f_pol, u_pol, v_pol
  real(kind=8) :: vor_n,vor_s,vor_w,vor_e,vor_c
  real(kind=8) :: tmask_n, tmask_s, tmask_e, tmask_w
  real(kind=8) :: umask_n, umask_s, umask_e, umask_w
  real(kind=8) :: vmask_n, vmask_s, vmask_e, vmask_w
  real(kind=8) :: fmask_n, fmask_s, fmask_e, fmask_w
  real(kind=8) :: lat_polex1, lat_polex2, lat_polex3, lat_polex4
  real(kind=8) :: u_arc_tra_n, u_arc_trb_n, u_arc_tra_e, u_arc_trb_e
  real(kind=8) :: v_arc_tra_n, v_arc_trb_n, v_arc_tra_e, v_arc_trb_e
  real(kind=8) :: u_arc_itra_n, u_arc_itrb_n, u_arc_itra_e, u_arc_itrb_e
  real(kind=8) :: v_arc_itra_n, v_arc_itrb_n, v_arc_itra_e, v_arc_itrb_e
  real(kind=8) :: temps1, temps2, temps3
  
  real(kind=8) :: u0
  real(kind=8) :: phi0, phi1, lat_rad, lon_rad, ddd0, lam0, lam1, lam2, lam3, eta_dis, eta_dis_p2
  
  real(kind=8), allocatable :: lat_t(:), lon_t(:), ff(:), bathy(:)
  real(kind=8), allocatable :: hdfc1(:), hdfc2(:), alh(:), udfc1(:), udfc2(:), vdfc1(:), vdfc2(:)
  real(kind=8), allocatable :: lat_u(:), lon_u(:), lat_v(:), lon_v(:)
  real(kind=8), allocatable :: e1t(:), e1u(:), e1v(:), e1f(:)
  real(kind=8), allocatable :: e2t(:), e2u(:), e2v(:), e2f(:)
  real(kind=8), allocatable :: ubar_bef(:), vbar_bef(:), eta_bef(:), vor_bef(:), div_bef(:), divbar_unror(:)
  real(kind=8), allocatable :: ubar_unror(:), vbar_unror(:)
  real(kind=8), allocatable :: ubar_now(:), vbar_now(:), eta_now(:), vor_now(:), div_now(:)
  real(kind=8), allocatable :: uw(:), vw(:), hw(:), vorw(:), divw(:)
  real(kind=8), allocatable :: t_u(:), t_v(:), u2_i(:), v2_j(:)
  real(kind=8), allocatable :: vbar_f(:), ubar_f(:), temp(:)
  real(kind=8), allocatable :: temp1(:), temp2(:), temp3(:), temp4(:), temp5(:), temp6(:)
  real(kind=8), allocatable :: tmask(:), umask(:), vmask(:), fmask(:)
  real(kind=8), allocatable :: u_sbc(:), v_sbc(:)
  real(kind=8), allocatable :: u_arc_tra(:), u_arc_trb(:), v_arc_tra(:), v_arc_trb(:)
  real(kind=8), allocatable :: u_arc_itra(:), u_arc_itrb(:), v_arc_itra(:), v_arc_itrb(:)
  real(kind=8), allocatable :: u_arc_tra2(:), u_arc_trb2(:), v_arc_tra2(:), v_arc_trb2(:)
  real(kind=8), allocatable :: u_arc_itra2(:), u_arc_itrb2(:), v_arc_itra2(:), v_arc_itrb2(:)
  real(kind=8), allocatable :: vor_part_u(:), keg_part_u(:), hpg_part_u(:), diff_part_u(:)
  real(kind=8), allocatable :: vor_part_v(:), keg_part_v(:), hpg_part_v(:), diff_part_v(:)
  real(kind=8), allocatable :: sbc_part_u(:), bot_fric_u(:), adv_part_t(:)
  real(kind=8), allocatable :: sbc_part_v(:), bot_fric_v(:)
  logical :: euler = .true.
  character (len = 10):: cdate, ctime
  
!=====start set basic parameters==================
  fn_grid = "smc_zhangyu.nc"
  fn_ini = "smc_ini.nc"
  fn_sbc = "smc_sbc.nc"
  fn_out = "smc_sw_nobathy_grid_a_out_test.nc"
  rad = pie/180.0d0
!=====end set basic parameters==================

!=====start read basic grid infomation============
  status = nf90_open(trim(fn_grid), nf90_nowrite, ncid_geo)
  status = nf90_inq_dimid(ncid_geo, 'nglo_b', dimid)
  status = nf90_inquire_dimension(ncid_geo, dimid, len = nglo_b)
  allocate(lat_t(nglo_b), lon_t(nglo_b))
  allocate(lat_u(nglo_b), lon_u(nglo_b), lat_v(nglo_b), lon_v(nglo_b))
  allocate(e1t(nglo_b), e1u(nglo_b), e1v(nglo_b), e1f(nglo_b))
  allocate(e2t(nglo_b), e2u(nglo_b), e2v(nglo_b), e2f(nglo_b))
  allocate(tn(nglo_b,neig), ts(nglo_b,neig), tw(nglo_b,neig), te(nglo_b,neig))
  status = nf90_get_att(ncid_geo, nf90_global, "lat_polex1", lat_polex1)
  status = nf90_get_att(ncid_geo, nf90_global, "lat_polex2", lat_polex2)
  status = nf90_get_att(ncid_geo, nf90_global, "lat_polex3", lat_polex3)
  status = nf90_get_att(ncid_geo, nf90_global, "lat_polex4", lat_polex4)
  status = nf90_inq_varid(ncid_geo, 'lon_t_b', varid)
  status = nf90_get_var(ncid_geo, varid, lon_t)
  status = nf90_inq_varid(ncid_geo, 'lat_t_b', varid)
  status = nf90_get_var(ncid_geo, varid, lat_t)
  status = nf90_inq_varid(ncid_geo, 'lon_u_b', varid)
  status = nf90_get_var(ncid_geo, varid, lon_u)
  status = nf90_inq_varid(ncid_geo, 'lat_u_b', varid)
  status = nf90_get_var(ncid_geo, varid, lat_u)
  status = nf90_inq_varid(ncid_geo, 'lon_v_b', varid)
  status = nf90_get_var(ncid_geo, varid, lon_v)
  status = nf90_inq_varid(ncid_geo, 'lat_v_b', varid)
  status = nf90_get_var(ncid_geo, varid, lat_v)
  status = nf90_inq_varid(ncid_geo, 'e1_t_b', varid)
  status = nf90_get_var(ncid_geo, varid, e1t)
  status = nf90_inq_varid(ncid_geo, 'e1_u_b', varid)
  status = nf90_get_var(ncid_geo, varid, e1u)
  status = nf90_inq_varid(ncid_geo, 'e1_v_b', varid)
  status = nf90_get_var(ncid_geo, varid, e1v)
  status = nf90_inq_varid(ncid_geo, 'e1_f_b', varid)
  status = nf90_get_var(ncid_geo, varid, e1f)
  status = nf90_inq_varid(ncid_geo, 'e2_t_b', varid)
  status = nf90_get_var(ncid_geo, varid, e2t)
  status = nf90_inq_varid(ncid_geo, 'e2_u_b', varid)
  status = nf90_get_var(ncid_geo, varid, e2u)
  status = nf90_inq_varid(ncid_geo, 'e2_v_b', varid)
  status = nf90_get_var(ncid_geo, varid, e2v)
  status = nf90_inq_varid(ncid_geo, 'e2_f_b', varid)
  status = nf90_get_var(ncid_geo, varid, e2f)
  status = nf90_inq_varid(ncid_geo, 'tgrid_b_n', varid)
  status = nf90_get_var(ncid_geo, varid, tn)
  status = nf90_inq_varid(ncid_geo, 'tgrid_b_s', varid)
  status = nf90_get_var(ncid_geo, varid, ts)
  status = nf90_inq_varid(ncid_geo, 'tgrid_b_w', varid)
  status = nf90_get_var(ncid_geo, varid, tw)
  status = nf90_inq_varid(ncid_geo, 'tgrid_b_e', varid)
  status = nf90_get_var(ncid_geo, varid, te)
  status = nf90_close(ncid_geo)
!=====end read basic grid infomation============

!======start define initial field======================
  
  if (lat_t(1) == -90.0d0) then
    i_poln = 1
    i_st = 2
  else
    i_st = 1
    i_poln = 0
  end if
  if (lat_t(nglo_b) == 90.0d0) then
    i_en = nglo_b - 1
    i_pols = nglo_b
  else
    i_en = nglo_b
    i_pols = 0
  end if

  nbdy_b = 0
  idx_b = 1
  do i=1,nglo_b
    if (lat_t(i) == lat_polex1) then
      where(ts == i) ts = nglo_b+idx_b
      nbdy_b = nbdy_b + 1
      idx_b = idx_b + 1
      !write(*,"(1x,i5,a,f8.3)") i, " lat = ", lat_t(i)
    end if
    if (lat_t(i) == lat_polex2) then
      where(tn == i) tn = nglo_b+idx_b
      nbdy_b = nbdy_b + 1
      idx_b = idx_b + 1
      !write(*,"(1x,i5,a,f8.3)") i, " lat = ", lat_t(i)
    end if
    if (lat_t(i) == lat_polex3) then
      where(ts == i) ts = nglo_b+idx_b
      nbdy_b = nbdy_b + 1
      idx_b = idx_b + 1
      !write(*,"(1x,i5,a,f8.3)") i, " lat = ", lat_t(i)
    end if
    if (lat_t(i) == lat_polex4) then
      where(tn == i) tn = nglo_b+idx_b
      nbdy_b = nbdy_b + 1
      idx_b = idx_b + 1
      !write(*,"(1x,i5,a,f8.3)") i, " lat = ", lat_t(i)
    end if
  end do
  write(*,"(1x,a,f8.3)") "lat_polex1 = ", lat_polex1
  write(*,"(1x,a,f8.3)") "lat_polex2 = ", lat_polex2
  write(*,"(1x,a,f8.3)") "lat_polex3 = ", lat_polex3
  write(*,"(1x,a,f8.3)") "lat_polex4 = ", lat_polex4
  write(*,"(1x,a,i8)") "nbdy_b = ", nbdy_b
  
  allocate(idx_b_m(nbdy_b), idx_m_b(nglo_b))
  allocate(ubar_bef(nglo_b+nbdy_b+1), vbar_bef(nglo_b+nbdy_b+1), eta_bef(nglo_b+nbdy_b+1))
  allocate(vor_bef(nglo_b+nbdy_b+1), div_bef(nglo_b+nbdy_b+1), divbar_unror(nglo_b+nbdy_b+1))
  allocate(ubar_now(nglo_b+nbdy_b+1), vbar_now(nglo_b+nbdy_b+1), eta_now(nglo_b+nbdy_b+1))
  allocate(vor_now(nglo_b+nbdy_b+1), div_now(nglo_b+nbdy_b+1))
  allocate(vbar_f(nglo_b+nbdy_b+1), ubar_f(nglo_b+nbdy_b+1), ubar_unror(nglo_b+nbdy_b+1), vbar_unror(nglo_b+nbdy_b+1))
  allocate(t_u(nglo_b+nbdy_b+1), t_v(nglo_b+nbdy_b+1), u2_i(nglo_b+nbdy_b+1), v2_j(nglo_b+nbdy_b+1))
  allocate(u_arc_tra(nglo_b+nbdy_b+1), u_arc_trb(nglo_b+nbdy_b+1), v_arc_tra(nglo_b+nbdy_b+1), v_arc_trb(nglo_b+nbdy_b+1))
  allocate(u_arc_itra(nglo_b+nbdy_b+1), u_arc_itrb(nglo_b+nbdy_b+1), v_arc_itra(nglo_b+nbdy_b+1), v_arc_itrb(nglo_b+nbdy_b+1))
  allocate(u_arc_tra2(nglo_b+1), u_arc_trb2(nglo_b+1), v_arc_tra2(nglo_b+1), v_arc_trb2(nglo_b+1))
  allocate(u_arc_itra2(nglo_b+1), u_arc_itrb2(nglo_b+1), v_arc_itra2(nglo_b+1), v_arc_itrb2(nglo_b+1))
  allocate(vor_part_u(nglo_b+1), keg_part_u(nglo_b+1), hpg_part_u(nglo_b+1), diff_part_u(nglo_b+1))
  allocate(vor_part_v(nglo_b+1), keg_part_v(nglo_b+1), hpg_part_v(nglo_b+1), diff_part_v(nglo_b+1))
  allocate(sbc_part_u(nglo_b+1), bot_fric_u(nglo_b+1), adv_part_t(nglo_b+1))
  allocate(sbc_part_v(nglo_b+1), bot_fric_v(nglo_b+1))
  allocate(tmask(nglo_b+nbdy_b+1),umask(nglo_b+nbdy_b+1),vmask(nglo_b+nbdy_b+1),fmask(nglo_b+nbdy_b+1),ff(nglo_b+nbdy_b+1))
  allocate(bathy(nglo_b+nbdy_b+1))
  allocate(hdfc1(nglo_b+nbdy_b+1), hdfc2(nglo_b+nbdy_b+1), alh(nglo_b+nbdy_b+1), udfc1(nglo_b+nbdy_b+1), udfc2(nglo_b+nbdy_b+1))
  allocate(vdfc1(nglo_b+nbdy_b+1), vdfc2(nglo_b+nbdy_b+1))
  allocate(temp1(nglo_b+nbdy_b+1),temp2(nglo_b+nbdy_b+1),temp3(nglo_b+nbdy_b+1),temp4(nglo_b+nbdy_b+1),temp5(nglo_b+nbdy_b+1))
  allocate(temp6(nglo_b+nbdy_b+1))
  allocate(u_sbc(nglo_b+1), v_sbc(nglo_b+1))
  allocate(temp(nglo_b))
  
  status = nf90_open(trim(fn_grid), nf90_nowrite, ncid_geo)
  status = nf90_inq_varid(ncid_geo, 'smc_tgrid_bathy_b', varid)
  status = nf90_get_var(ncid_geo, varid, temp)
  bathy(1:nglo_b) = temp(1:nglo_b)
  bathy(nglo_b+nbdy_b+1) = 0.0d0
  status = nf90_inq_varid(ncid_geo, 'ff_f_b', varid)
  status = nf90_get_var(ncid_geo, varid, temp)
  ff(1:nglo_b) = temp(1:nglo_b)
  ff(nglo_b+nbdy_b+1) = 0.0d0
  status = nf90_inq_varid(ncid_geo, 'smc_tgrid_b_mask', varid)
  status = nf90_get_var(ncid_geo, varid, temp)
  tmask(1:nglo_b) = temp(1:nglo_b)
  tmask(nglo_b+nbdy_b+1) = 0.0d0
  status = nf90_inq_varid(ncid_geo, 'smc_ugrid_b_mask', varid)
  status = nf90_get_var(ncid_geo, varid, temp)
  umask(1:nglo_b) = temp(1:nglo_b)
  umask(nglo_b+nbdy_b+1) = 0.0d0
  status = nf90_inq_varid(ncid_geo, 'smc_vgrid_b_mask', varid)
  status = nf90_get_var(ncid_geo, varid, temp)
  vmask(1:nglo_b) = temp(1:nglo_b)
  vmask(nglo_b+nbdy_b+1) = 0.0d0
  status = nf90_inq_varid(ncid_geo, 'smc_fgrid_b_mask', varid)
  status = nf90_get_var(ncid_geo, varid, temp)
  fmask(1:nglo_b) = temp(1:nglo_b)
  fmask(nglo_b+nbdy_b+1) = 3.0d0
  status = nf90_close(ncid_geo)
  where(tn<0) tn = nglo_b+nbdy_b+1
  where(ts<0) ts = nglo_b+nbdy_b+1
  where(tw<0) tw = nglo_b+nbdy_b+1
  where(te<0) te = nglo_b+nbdy_b+1
  
  idx_b_m = 0
  idx_m_b = 0
  ubar_bef = 0.0d0
  vbar_bef = 0.0d0
  eta_bef = 0.0d0
  vor_bef = 0.0d0
  div_bef = 0.0d0
  ubar_now = 0.0d0
  vbar_now = 0.0d0
  eta_now = 0.0d0
  vor_now = 0.0d0
  div_now = 0.0d0
  vbar_f = 0.0d0
  ubar_f = 0.0d0
  ubar_unror = 0.0d0
  vbar_unror = 0.0d0
  divbar_unror = 0.0d0
  t_u = 0.0d0
  t_v = 0.0d0
  u2_i = 0.0d0
  v2_j = 0.0d0
  u_arc_tra = 0.0d0
  u_arc_trb = 0.0d0
  v_arc_tra = 0.0d0
  v_arc_trb = 0.0d0
  u_arc_itra = 0.0d0
  u_arc_itrb = 0.0d0
  v_arc_itra = 0.0d0
  v_arc_itrb = 0.0d0
  u_arc_tra2 = 0.0d0
  u_arc_trb2 = 0.0d0
  v_arc_tra2 = 0.0d0
  v_arc_trb2 = 0.0d0
  u_arc_itra2 = 0.0d0
  u_arc_itrb2 = 0.0d0
  v_arc_itra2 = 0.0d0
  v_arc_itrb2 = 0.0d0
  vor_part_u = 0.0d0
  vor_part_v = 0.0d0
  keg_part_u = 0.0d0
  keg_part_v = 0.0d0
  hpg_part_u = 0.0d0
  hpg_part_v = 0.0d0
  diff_part_u = 0.0d0
  diff_part_v = 0.0d0
  sbc_part_u = 0.0d0
  sbc_part_v = 0.0d0
  bot_fric_u = 0.0d0
  bot_fric_v = 0.0d0
  adv_part_t = 0.0d0
  u_sbc = 0.0d0
  v_sbc = 0.0d0
  temp1 = 0.0d0
  temp2 = 0.0d0
  temp3 = 0.0d0
  temp4 = 0.0d0
  temp5 = 0.0d0
  temp6 = 0.0d0
  alh = almu
  hdfc1 = 0.0d0
  hdfc2 = 0.0d0
  udfc1 = 0.0d0
  udfc2 = 0.0d0
  vdfc1 = 0.0d0
  vdfc2 = 0.0d0
      
  idx_b = 1
  do i=1,nglo_b
    if (lat_t(i) == lat_polex1) then
      idx_b_m(idx_b) = i
      idx_m_b(i) = idx_b + nglo_b
      idx_b = idx_b + 1
    end if
    if (lat_t(i) == lat_polex2) then
      idx_b_m(idx_b) = i
      idx_m_b(i) = idx_b + nglo_b
      idx_b = idx_b + 1
    end if
    if (lat_t(i) == lat_polex3) then
      idx_b_m(idx_b) = i
      idx_m_b(i) = idx_b + nglo_b
      idx_b = idx_b + 1
    end if
    if (lat_t(i) == lat_polex4) then
      idx_b_m(idx_b) = i
      idx_m_b(i) = idx_b + nglo_b
      idx_b = idx_b + 1
    end if
    lat_rad = lat_t(i)*rad
    alh(i) = almu*(1.0d0 - gama + gama*sin(lat_rad)**2)
  end do
  do i=nglo_b+1,nglo_b+nbdy_b
    j = idx_b_m(i-nglo_b)
    tmask(i) = tmask(j)
  	umask(i) = umask(j)
  	vmask(i) = vmask(j)
  	fmask(i) = fmask(j)
  	ff(i) = ff(j)
  	bathy(i) = bathy(j)
  	alh(i) = alh(j)
  end do
  
  do i=1,nglo_b
    if ((lat_t(i) <= lat_polex1) .or. (lat_t(i) >= lat_polex4)) then
  	  u_arc_tra(i) = cos(rad*lon_t(i))*sin(rad*lat_t(i))/sqrt(1-(cos(rad*lon_t(i))*cos(rad*lat_t(i)))**2)
  	  u_arc_trb(i) = -sin(rad*lon_t(i))/sqrt(1-(cos(rad*lon_t(i))*cos(rad*lat_t(i)))**2)
  	  v_arc_tra(i) = sin(rad*lon_t(i))/sqrt(1-(cos(rad*lon_t(i))*cos(rad*lat_t(i)))**2)
  	  v_arc_trb(i) = cos(rad*lon_t(i))*sin(rad*lat_t(i))/sqrt(1-(cos(rad*lon_t(i))*cos(rad*lat_t(i)))**2)
  	
  	  u_arc_itra(i) = cos(rad*lon_t(i))*sin(rad*lat_t(i))/sqrt(1-(cos(rad*lon_t(i))*cos(rad*lat_t(i)))**2)
  	  u_arc_itrb(i) = sin(rad*lon_t(i))/sqrt(1-(cos(rad*lon_t(i))*cos(rad*lat_t(i)))**2)
  	  v_arc_itra(i) = -sin(rad*lon_t(i))/sqrt(1-(cos(rad*lon_t(i))*cos(rad*lat_t(i)))**2)
  	  v_arc_itrb(i) = cos(rad*lon_t(i))*sin(rad*lat_t(i))/sqrt(1-(cos(rad*lon_t(i))*cos(rad*lat_t(i)))**2)
  	else
  	  u_arc_tra(i) = 1.0d0
  	  u_arc_trb(i) = 0.0d0
  	  v_arc_tra(i) = 0.0d0
  	  v_arc_trb(i) = 1.0d0
  	
  	  u_arc_itra(i) = 1.0d0
  	  u_arc_itrb(i) = 0.0d0
  	  v_arc_itra(i) = 0.0d0
  	  v_arc_itrb(i) = 1.0d0
    end if
  	u_arc_tra2(i) = cos(rad*lon_t(i))*sin(rad*lat_t(i))/sqrt(1-(cos(rad*lon_t(i))*cos(rad*lat_t(i)))**2)
    u_arc_trb2(i) = -(rad*lon_t(i))/sqrt(1-(cos(rad*lon_t(i))*cos(rad*lat_t(i)))**2)
  	v_arc_tra2(i) = sin(rad*lon_t(i))/sqrt(1-(cos(rad*lon_t(i))*cos(rad*lat_t(i)))**2)
  	v_arc_trb2(i) = cos(rad*lon_t(i))*sin(rad*lat_t(i))/sqrt(1-(cos(rad*lon_t(i))*cos(rad*lat_t(i)))**2)
  	
  	u_arc_itra2(i) = cos(rad*lon_t(i))*sin(rad*lat_t(i))/sqrt(1-(cos(rad*lon_t(i))*cos(rad*lat_t(i)))**2)
  	u_arc_itrb2(i) = sin(rad*lon_t(i))/sqrt(1-(cos(rad*lon_t(i))*cos(rad*lat_t(i)))**2)
  	v_arc_itra2(i) = -sin(rad*lon_t(i))/sqrt(1-(cos(rad*lon_t(i))*cos(rad*lat_t(i)))**2)
  	v_arc_itrb2(i) = cos(rad*lon_t(i))*sin(rad*lat_t(i))/sqrt(1-(cos(rad*lon_t(i))*cos(rad*lat_t(i)))**2)
  end do
  do i=nglo_b+1,nglo_b+nbdy_b
    j = idx_b_m(i-nglo_b)
    if (lat_t(j) == lat_polex1) then
      u_arc_tra(i) = 1.0d0
  	  u_arc_trb(i) = 0.0d0
  	  v_arc_tra(i) = 0.0d0
  	  v_arc_trb(i) = 1.0d0
  	  u_arc_itra(i) = cos(rad*lon_t(j))*sin(rad*lat_t(j))/sqrt(1-(cos(rad*lon_t(j))*cos(rad*lat_t(j)))**2)
  	  u_arc_itrb(i) = sin(rad*lon_t(j))/sqrt(1-(cos(rad*lon_t(j))*cos(rad*lat_t(j)))**2)
  	  v_arc_itra(i) = -sin(rad*lon_t(j))/sqrt(1-(cos(rad*lon_t(j))*cos(rad*lat_t(j)))**2)
  	  v_arc_itrb(i) = cos(rad*lon_t(j))*sin(rad*lat_t(j))/sqrt(1-(cos(rad*lon_t(j))*cos(rad*lat_t(j)))**2)
    end if
    if (lat_t(j) == lat_polex2) then
      u_arc_tra(i) = cos(rad*lon_t(j))*sin(rad*lat_t(j))/sqrt(1-(cos(rad*lon_t(j))*cos(rad*lat_t(j)))**2)
  	  u_arc_trb(i) = -sin(rad*lon_t(j))/sqrt(1-(cos(rad*lon_t(j))*cos(rad*lat_t(j)))**2)
  	  v_arc_tra(i) = sin(rad*lon_t(j))/sqrt(1-(cos(rad*lon_t(j))*cos(rad*lat_t(j)))**2)
  	  v_arc_trb(i) = cos(rad*lon_t(j))*sin(rad*lat_t(j))/sqrt(1-(cos(rad*lon_t(j))*cos(rad*lat_t(j)))**2)
  	  u_arc_itra(i) = u_arc_tra(i)
  	  u_arc_itrb(i) = u_arc_trb(i)
  	  v_arc_itra(i) = v_arc_tra(i)
  	  v_arc_itrb(i) = v_arc_trb(i)
    end if
    if (lat_t(j) == lat_polex3) then
      u_arc_tra(i) = cos(rad*lon_t(j))*sin(rad*lat_t(j))/sqrt(1-(cos(rad*lon_t(j))*cos(rad*lat_t(j)))**2)
  	  u_arc_trb(i) = -sin(rad*lon_t(j))/sqrt(1-(cos(rad*lon_t(j))*cos(rad*lat_t(j)))**2)
  	  v_arc_tra(i) = sin(rad*lon_t(j))/sqrt(1-(cos(rad*lon_t(j))*cos(rad*lat_t(j)))**2)
  	  v_arc_trb(i) = cos(rad*lon_t(j))*sin(rad*lat_t(j))/sqrt(1-(cos(rad*lon_t(j))*cos(rad*lat_t(j)))**2)
  	  u_arc_itra(i) = u_arc_tra(i)
  	  u_arc_itrb(i) = u_arc_trb(i)
  	  v_arc_itra(i) = v_arc_tra(i)
  	  v_arc_itrb(i) = v_arc_trb(i)
    end if
    if (lat_t(j) == lat_polex4) then
      u_arc_tra(i) = 1.0d0
  	  u_arc_trb(i) = 0.0d0
  	  v_arc_tra(i) = 0.0d0
  	  v_arc_trb(i) = 1.0d0
  	  u_arc_itra(i) = cos(rad*lon_t(j))*sin(rad*lat_t(j))/sqrt(1-(cos(rad*lon_t(j))*cos(rad*lat_t(j)))**2)
  	  u_arc_itrb(i) = sin(rad*lon_t(j))/sqrt(1-(cos(rad*lon_t(j))*cos(rad*lat_t(j)))**2)
  	  v_arc_itra(i) = -sin(rad*lon_t(j))/sqrt(1-(cos(rad*lon_t(j))*cos(rad*lat_t(j)))**2)
  	  v_arc_itrb(i) = cos(rad*lon_t(j))*sin(rad*lat_t(j))/sqrt(1-(cos(rad*lon_t(j))*cos(rad*lat_t(j)))**2)
    end if
  end do
  
!!---just for test-----
  u_arc_tra = 1.0d0
  u_arc_trb = 0.0d0
  v_arc_tra = 0.0d0
  v_arc_trb = 1.0d0
  u_arc_itra = 1.0d0
  u_arc_itrb = 0.0d0
  v_arc_itra = 0.0d0
  v_arc_itrb = 1.0d0
  u_arc_tra2 = 1.0d0
  u_arc_trb2 = 0.0d0
  v_arc_tra2 = 0.0d0
  v_arc_trb2 = 1.0d0
  u_arc_itra2 = 1.0d0
  u_arc_itrb2 = 0.0d0
  v_arc_itra2 = 0.0d0
  v_arc_itrb2 = 1.0d0
!!!!

  status = nf90_open(trim(fn_ini), nf90_nowrite, ncid_ini)
  status = nf90_inq_varid(ncid_ini, 'h_ini', varid)
  status = nf90_get_var(ncid_ini, varid, temp)
  eta_bef(1:nglo_b) = temp(1:nglo_b)
  status = nf90_inq_varid(ncid_ini, 'u_ini', varid)
  status = nf90_get_var(ncid_ini, varid, temp)
  ubar_bef(1:nglo_b) = temp(1:nglo_b)
  status = nf90_close(ncid_ini)

!!! simple test
!  eta_bef = 0.0d0
!  ubar_bef = 0.0d0
!  eta_bef(70698) = 20.0d0
  
!!! test2
!  u0 = omega*rearth/12.0d0
!  do i=1,nglo_b
!    lat_rad = lat_t(i)*rad
!    eta_bef(i) = (3000.0d0*grav-u0*(omega*rearth+0.5*u0)*sin(lat_rad)**2)/grav
!    ubar_bef(i) = u0*cos(lat_rad)
!    vbar_bef(i) = 0.0d0
!  end do

!!! test 6  
!  phi0 = pie/7.0d0
!  phi1 = pie/2.0d0 - pie/7.0d0
!  lam0 = 1.5d0*pie
!  ddd0 = 120.0d0
!  do i=1,nglo_b
!    lat_rad = lat_t(i)*rad
!    lon_rad = lon_t(i)*rad
!    if ((lat_rad > phi0) .and. (lat_rad < phi1)) then
!      ubar_bef(i) = 80.0d0*exp(4.0d0/(phi1 - phi0)**2-1.0d0/((lat_rad - phi0)*(phi1 - lat_rad)))
!    end if
!    lam1 = (lon_rad - 2.0d0*pie - lam0)**2
!    lam2 = (lon_rad - lam0)**2
!    lam3 = (lon_rad + 2.0d0*pie - lam0)**2
!    eta_dis_p2 = min(lam1,lam2)
!    eta_dis_p2 = min(eta_dis_p2,lam3)
!    eta_dis = ddd0*cos(lat_rad)*exp(-9.0d0*eta_dis_p2-225.0d0*(lat_rad-pie/4.0d0)**2)
!    eta_bef(i) = eta_dis
!  end do
  
  do i=nglo_b+1,nglo_b+nbdy_b
    j = idx_b_m(i-nglo_b)
    ubar_bef(i) = ubar_bef(j)
    vbar_bef(i) = vbar_bef(j)
    eta_bef(i) = eta_bef(j)
  end do
  do i=1,nglo_b+nbdy_b
    ubar_unror(i) = ubar_bef(i)
    vbar_unror(i) = vbar_bef(i)
    u_pol = u_arc_tra(i)*ubar_bef(i) + u_arc_trb(i)*vbar_bef(i)
    v_pol = v_arc_tra(i)*ubar_bef(i) + v_arc_trb(i)*vbar_bef(i)
    ubar_bef(i) = u_pol
    vbar_bef(i) = v_pol
  end do
  
!  status = nf90_open(trim(fn_sbc), nf90_nowrite, ncid_sbc)
!  status = nf90_inq_varid(ncid_sbc, 'u_sbc', varid)
!  status = nf90_get_var(ncid_sbc, varid, temp)
!  do i=1,nglo_b
!    u_sbc(i) = temp(i)
!  end do
!  status = nf90_inq_varid(ncid_sbc, 'v_sbc', varid)
!  status = nf90_get_var(ncid_sbc, varid, temp)
!  do i=1,nglo_b
!    v_sbc(i) = temp(i)
!  end do
!  do i=1,nglo_b
!    u_pol = u_arc_tra(i)*u_sbc(i) + u_arc_trb(i)*v_sbc(i)
!    v_pol = v_arc_tra(i)*u_sbc(i) + v_arc_trb(i)*v_sbc(i)
!    u_sbc(i) = u_pol
!    v_sbc(i) = v_pol
!  end do
!  status = nf90_close(ncid_sbc)
!======end define initial field======================

  tout = 0
!  ff = 0
  
!======start define output nc file====================== 
  allocate(uw(nglo_b), vw(nglo_b), hw(nglo_b), vorw(nglo_b), divw(nglo_b))
  status = nf90_create(trim(fn_out), nf90_64bit_offset, ncid_out)
  status = nf90_def_dim(ncid_out, "nglo_b", nglo_b, nglo_b_dimid)
  status = nf90_def_dim(ncid_out, "time", nf90_unlimited, t_dimid)
  dimids2d = (/ nglo_b_dimid, t_dimid/)
  
  status = nf90_def_var(ncid_out, "lon", nf90_double, nglo_b_dimid, lon_id)
  status = nf90_def_var(ncid_out, "lat", nf90_double, nglo_b_dimid, lat_id)
  status = nf90_def_var(ncid_out, "time",nf90_double, t_dimid, time_id)
  status = nf90_def_var(ncid_out, "h", nf90_double, dimids2d, h_id)
  status = nf90_def_var(ncid_out, "ubar", nf90_double, dimids2d, u_id)
  status = nf90_def_var(ncid_out, "vbar", nf90_double, dimids2d, v_id)
  status = nf90_def_var(ncid_out, "vor", nf90_double, dimids2d, vor_id)
  status = nf90_def_var(ncid_out, "div", nf90_double, dimids2d, div_id)
  status = nf90_put_att(ncid_out, nf90_global, "Author", 'Yu Zhang')
  status = nf90_put_att(ncid_out, nf90_global, "Organization", 'National Marine Environmental Forecasting Center')
  call date_and_time(cdate, ctime)
  status = nf90_put_att(ncid_out, nf90_global, "date", trim(cdate)//', '//trim(ctime(1:2))//':'//trim(ctime(3:4)) &
            & //':'//trim(ctime(5:6)))
  status = nf90_put_att(ncid_out, nf90_global, "about", 'SMC experiment output file')
  status = nf90_enddef(ncid_out)
  status = nf90_put_var(ncid_out, lon_id, lon_t )
  status = nf90_put_var(ncid_out, lat_id, lat_t )
  status = nf90_close(ncid_out)
  
  tout = tout + 1
  time(1) = tout
  hw(:) = eta_bef(1:nglo_b)
  uw(:) = ubar_bef(1:nglo_b)
  vw(:) = vbar_bef(1:nglo_b)
  vorw(:) = vor_bef(1:nglo_b) - ff(1:nglo_b)
  divw(:) = divbar_unror(1:nglo_b)
  status = nf90_open(trim(fn_out), nf90_write, ncid_out)
  status = nf90_inq_varid(ncid_out, 'time', time_id)
  status = nf90_put_var(ncid_out, time_id, time, (/tout/), (/1/))
  status = nf90_inq_varid(ncid_out, 'h', h_id)
  status = nf90_put_var(ncid_out, h_id, hw, (/1,tout/), (/nglo_b,1/))
  status = nf90_inq_varid(ncid_out, 'ubar', u_id)
  status = nf90_put_var(ncid_out, u_id, uw, (/1,tout/), (/nglo_b,1/))
  status = nf90_inq_varid(ncid_out, 'vbar', v_id)
  status = nf90_put_var(ncid_out, v_id, vw, (/1,tout/), (/nglo_b,1/))
  status = nf90_inq_varid(ncid_out, 'vor', vor_id)
  status = nf90_put_var(ncid_out, vor_id, vorw, (/1,tout/), (/nglo_b,1/))
  status = nf90_inq_varid(ncid_out, 'div', div_id)
  status = nf90_put_var(ncid_out, div_id, divw, (/1,tout/), (/nglo_b,1/))
  status = nf90_close(ncid_out)
  
!======end define output nc file======================
  
  call system_clock ( clock_count1, clock_rate, clock_max )
  
  !$acc data copyin(tmask) copyin(umask) copyin(vmask) copyin(fmask) &
  !$acc      copyin(e1t) copyin(e1u) copyin(e1v) copyin(e1f) copyin(lat_t) &
  !$acc      copyin(e2t) copyin(e2u) copyin(e2v) copyin(e2f) copyin(ff) copyin(bathy) &
  !$acc      copyin(tn) copyin(ts) copyin(tw) copyin(te) copyin(idx_b_m) copyin(idx_m_b) &
  !$acc      copyin(u_sbc) copyin(v_sbc) copyin(ubar_f) copyin(vbar_f) copyin(t_u) copyin(t_v) &
  !$acc      copyin(u2_i) copyin(v2_j) copyin(ubar_unror) copyin(vbar_unror) copyin(divbar_unror) &
  !$acc      copyin(eta_bef) copyin(vor_bef) copyin(div_bef) copyin(ubar_bef) copyin(vbar_bef) &
  !$acc      copyin(eta_now) copyin(vor_now) copyin(div_now) copyin(ubar_now) copyin(vbar_now) &
  !$acc      copyin(u_arc_tra) copyin(u_arc_trb) copyin(v_arc_tra) copyin(v_arc_trb) &
  !$acc      copyin(u_arc_tra2) copyin(u_arc_trb2) copyin(v_arc_tra2) copyin(v_arc_trb2) &
  !$acc      copyin(u_arc_itra) copyin(u_arc_itrb) copyin(v_arc_itra) copyin(v_arc_itrb) &
  !$acc      copyin(u_arc_itra2) copyin(u_arc_itrb2) copyin(v_arc_itra2) copyin(v_arc_itrb2) &
  !$acc      copyin(vor_part_u) copyin(keg_part_u) copyin(hpg_part_u) copyin(diff_part_u) copyin(sbc_part_u) &
  !$acc      copyin(vor_part_v) copyin(keg_part_v) copyin(hpg_part_v) copyin(diff_part_v) copyin(sbc_part_v) &
  !$acc      copyin(bot_fric_u) copyin(bot_fric_v) copyin(adv_part_t) &
  !$acc      copyin(hdfc1) copyin(hdfc2) copyin(udfc1) copyin(vdfc1) copyin(udfc2) copyin(vdfc2) copyin(alh) &
  !$acc      copyin(temp1) copyin(temp2) copyin(temp3) copyin(temp4) copyin(temp5) copyin(temp6)
  
  do t=1,ta
!=================================================
!=============start integral======================
!=================================================
    !$acc kernels present(tmask, umask, vmask, fmask, &
    !$acc                 e1t, e1u, e1v, e1f, lat_t, &
    !$acc                 e2t, e2u, e2v, e2f, ff, bathy, &
    !$acc                 tn, ts, tw, te, idx_b_m, idx_m_b, &
    !$acc                 u_sbc, v_sbc, ubar_f, vbar_f, t_u, t_v, &
    !$acc                 u2_i, v2_j, ubar_unror, vbar_unror, divbar_unror, &
    !$acc                 eta_bef, vor_bef, div_bef, ubar_bef, vbar_bef, &
    !$acc                 eta_now, vor_now, div_now, ubar_now, vbar_now, &
    !$acc                 u_arc_tra, u_arc_trb, v_arc_tra, v_arc_trb, &
    !$acc                 u_arc_tra2, u_arc_trb2, v_arc_tra2, v_arc_trb2, &
    !$acc                 u_arc_itra, u_arc_itrb, v_arc_itra, v_arc_itrb, &
    !$acc                 u_arc_itra2, u_arc_itrb2, v_arc_itra2, v_arc_itrb2, &
    !$acc                 vor_part_u, keg_part_u, hpg_part_u, diff_part_u, sbc_part_u, &
    !$acc                 vor_part_v, keg_part_v, hpg_part_v, diff_part_v, sbc_part_v, &
    !$acc                 bot_fric_u, bot_fric_v, adv_part_t, &
    !$acc                 hdfc1, hdfc2, udfc1, vdfc1, udfc2, vdfc2, alh, &
    !$acc                 temp1, temp2, temp3, temp4, temp5, temp6)

    
      
!==============First Euler====================     

!=============================================
!--------compute udfc, vdfc, hdfc-------------
!============================================= 
    do i=i_st,i_en
      tn_1 = tn(i,1)
      tn_2 = tn(i,2)
      ts_1 = ts(i,1)
      ts_2 = ts(i,2)
      tw_1 = tw(i,1)
      tw_2 = tw(i,2)
      te_1 = te(i,1)
      te_2 = te(i,2)
      tmask_w = tmask(tw_1)*tmask(tw_2)
      tmask_s = tmask(ts_1)*tmask(ts_2)
      tmask_n = tmask(tn_1)*tmask(tn_2)
      tmask_e = tmask(te_1)*tmask(te_2)
      u_w = (ubar_bef(tw_1)+ubar_bef(tw_2))*0.5d0*tmask_w
      u_e = (ubar_bef(te_1)+ubar_bef(te_2))*0.5d0*tmask_e
      u_s = (ubar_bef(ts_1)+ubar_bef(ts_2))*0.5d0*tmask_s
      u_n = (ubar_bef(tn_1)+ubar_bef(tn_2))*0.5d0*tmask_n
      v_w = (vbar_bef(tw_1)+vbar_bef(tw_2))*0.5d0*tmask_w
      v_e = (vbar_bef(te_1)+vbar_bef(te_2))*0.5d0*tmask_e
      v_s = (vbar_bef(ts_1)+vbar_bef(ts_2))*0.5d0*tmask_s
      v_n = (vbar_bef(tn_1)+vbar_bef(tn_2))*0.5d0*tmask_n
      h_w = (eta_bef(tw_1)+eta_bef(tw_2))*0.5d0*tmask_w
      h_e = (eta_bef(te_1)+eta_bef(te_2))*0.5d0*tmask_e
      h_s = (eta_bef(ts_1)+eta_bef(ts_2))*0.5d0*tmask_s
      h_n = (eta_bef(tn_1)+eta_bef(tn_2))*0.5d0*tmask_n
      b_w = (bathy(tw_1)+bathy(tw_2))*0.5d0*tmask_w
      b_e = (bathy(te_1)+bathy(te_2))*0.5d0*tmask_e
      b_s = (bathy(ts_1)+bathy(ts_2))*0.5d0*tmask_s
      b_n = (bathy(tn_1)+bathy(tn_2))*0.5d0*tmask_n
      u_c = ubar_bef(i)*tmask(i)
      v_c = vbar_bef(i)*tmask(i)
      h_c = eta_bef(i)*tmask(i)
      b_c = bathy(i)*tmask(i)
      if ((tmask_e) .eq. 0.0d0) then
        h_e = h_c
        b_e = b_c
      end if
      if ((tmask_w) .eq. 0.0d0) then
        h_w = h_c
        b_w = b_c
      end if
      if ((tmask_n) .eq. 0.0d0) then
        h_n = h_c
        b_n = b_c
      end if
      if ((tmask_s) .eq. 0.0d0) then
        h_s = h_c
        b_s = b_c
      end if
      if (rn_shlat == 0) then
        if ((tmask_e) .eq. 0.0d0) then
          v_e = v_c
        end if
        if ((tmask_w) .eq. 0.0d0) then
          v_w = v_c
        end if
        if ((tmask_n) .eq. 0.0d0) then
          u_n = u_c
        end if
        if ((tmask_s) .eq. 0.0d0) then
          u_s = u_c
        end if
      else if (rn_shlat == 1) then
        if ((tmask_e) .eq. 0.0d0) then
          v_e = 0.0d0
        end if
        if ((tmask_w) .eq. 0.0d0) then
          v_w = 0.0d0
        end if
        if ((tmask_n) .eq. 0.0d0) then
          u_n = 0.0d0
        end if
        if ((tmask_s) .eq. 0.0d0) then
          u_s = 0.0d0
        end if
      else if (rn_shlat == 2) then
        if ((tmask_e) .eq. 0.0d0) then
          v_e = -1.0d0*v_c
        end if
        if ((tmask_w) .eq. 0.0d0) then
          v_w = -1.0d0*v_c
        end if
        if ((tmask_n) .eq. 0.0d0) then
          u_n = -1.0d0*u_c
        end if
        if ((tmask_s) .eq. 0.0d0) then
          u_s = -1.0d0*u_c
        end if
      end if
      u_we = tmask_w*(u_c*(h_c+b_c)+u_w*(h_w+b_w))/(h_c+b_c+h_w+b_w+eps)
      u_ee = tmask_e*(u_c*(h_c+b_c)+u_e*(h_e+b_e))/(h_c+b_c+h_e+b_e+eps)
      u_ne = tmask_n*(u_c*(h_c+b_c)+u_n*(h_n+b_n))/(h_c+b_c+h_n+b_n+eps)
      u_se = tmask_s*(u_c*(h_c+b_c)+u_s*(h_s+b_s))/(h_c+b_c+h_s+b_s+eps)
      v_we = tmask_w*(v_c*(h_c+b_c)+v_w*(h_w+b_w))/(h_c+b_c+h_w+b_w+eps)
      v_ee = tmask_e*(v_c*(h_c+b_c)+v_e*(h_e+b_e))/(h_c+b_c+h_e+b_e+eps)
      v_ne = tmask_n*(v_c*(h_c+b_c)+v_n*(h_n+b_n))/(h_c+b_c+h_n+b_n+eps)
      v_se = tmask_s*(v_c*(h_c+b_c)+v_s*(h_s+b_s))/(h_c+b_c+h_s+b_s+eps)
      h_we = (h_c+h_w)/2.0d0
      h_ee = (h_c+h_e)/2.0d0
      h_ne = (h_c+h_n)/2.0d0
      h_se = (h_c+h_s)/2.0d0
      udfc1(i) = alh(i) * ((u_ee-u_we)*e2t(i) + (u_ne-u_se)*e1t(i)) / (e1t(i)*e2t(i))
      vdfc1(i) = alh(i) * ((v_ee-v_we)*e2t(i) + (v_ne-v_se)*e1t(i)) / (e1t(i)*e2t(i))
      hdfc1(i) = alh(i) * ((h_ee-h_we)*e2t(i) + (h_ne-h_se)*e1t(i)) / (e1t(i)*e2t(i))
    end do
    !$acc loop independent
    do i=nglo_b+1,nglo_b+nbdy_b
      j = idx_b_m(i-nglo_b)
  	  hdfc1(i) = hdfc1(j)
  	  udfc1(i) = udfc1(j)
  	  vdfc1(i) = vdfc1(j)
    end do
    if (i_st > 1) then
      hdfc1(1) = (hdfc1(2) + hdfc1(3) + hdfc1(4) + hdfc1(5) + hdfc1(6) + hdfc1(7) &
                 & + hdfc1(8) + hdfc1(9))/8.0d0
      udfc1(1) = (udfc1(2) + udfc1(3) + udfc1(4) + udfc1(5) + udfc1(6) + udfc1(7) &
                 & + udfc1(8) + udfc1(9))/8.0d0
      vdfc1(1) = (vdfc1(2) + vdfc1(3) + vdfc1(4) + vdfc1(5) + vdfc1(6) + vdfc1(7) &
                 & + vdfc1(8) + vdfc1(9))/8.0d0
    end if
    if (i_en < nglo_b) then
      hdfc1(nglo_b) = (hdfc1(nglo_b-1) + hdfc1(nglo_b-2) + hdfc1(nglo_b-3) + hdfc1(nglo_b-4) &
                      & + hdfc1(nglo_b-5) + hdfc1(nglo_b-6) + hdfc1(nglo_b-7) + hdfc1(nglo_b-8))/8.0d0
      udfc1(nglo_b) = (udfc1(nglo_b-1) + udfc1(nglo_b-2) + udfc1(nglo_b-3) + udfc1(nglo_b-4) &
                      & + udfc1(nglo_b-5) + udfc1(nglo_b-6) + udfc1(nglo_b-7) + udfc1(nglo_b-8))/8.0d0
      vdfc1(nglo_b) = (vdfc1(nglo_b-1) + vdfc1(nglo_b-2) + vdfc1(nglo_b-3) + vdfc1(nglo_b-4) &
                      & + vdfc1(nglo_b-5) + vdfc1(nglo_b-6) + vdfc1(nglo_b-7) + vdfc1(nglo_b-8))/8.0d0
    end if
    do i=i_st,i_en
      tn_1 = tn(i,1)
      tn_2 = tn(i,2)
      ts_1 = ts(i,1)
      ts_2 = ts(i,2)
      tw_1 = tw(i,1)
      tw_2 = tw(i,2)
      te_1 = te(i,1)
      te_2 = te(i,2)
      tmask_w = tmask(tw_1)*tmask(tw_2)
      tmask_s = tmask(ts_1)*tmask(ts_2)
      tmask_n = tmask(tn_1)*tmask(tn_2)
      tmask_e = tmask(te_1)*tmask(te_2)
      h_w = (hdfc1(tw_1)+hdfc1(tw_2))*0.5d0*tmask_w
      h_e = (hdfc1(te_1)+hdfc1(te_2))*0.5d0*tmask_e
      h_s = (hdfc1(ts_1)+hdfc1(ts_2))*0.5d0*tmask_s
      h_n = (hdfc1(tn_1)+hdfc1(tn_2))*0.5d0*tmask_n
      h_c = hdfc1(i)*tmask(i)
      u_w = (udfc1(tw_1)+udfc1(tw_2))*0.5d0*tmask_w
      u_e = (udfc1(te_1)+udfc1(te_2))*0.5d0*tmask_e
      u_s = (udfc1(ts_1)+udfc1(ts_2))*0.5d0*tmask_s
      u_n = (udfc1(tn_1)+udfc1(tn_2))*0.5d0*tmask_n
      u_c = udfc1(i)*tmask(i)
      v_w = (vdfc1(tw_1)+vdfc1(tw_2))*0.5d0*tmask_w
      v_e = (vdfc1(te_1)+vdfc1(te_2))*0.5d0*tmask_e
      v_s = (vdfc1(ts_1)+vdfc1(ts_2))*0.5d0*tmask_s
      v_n = (vdfc1(tn_1)+vdfc1(tn_2))*0.5d0*tmask_n
      v_c = vdfc1(i)*tmask(i)
      if ((tmask_e) .eq. 0.0d0) then
        h_e = h_c
        u_e = u_c
        v_e = v_c
      end if
      if ((tmask_w) .eq. 0.0d0) then
        h_w = h_c
        u_w = u_c
        v_w = v_c
      end if
      if ((tmask_n) .eq. 0.0d0) then
        h_n = h_c
        u_n = u_c
        v_n = v_c
      end if
      if ((tmask_s) .eq. 0.0d0) then
        h_s = h_c
        u_s = u_c
        v_s = v_c
      end if
      h_we = (h_c+h_w)/2.0d0
      h_ee = (h_c+h_e)/2.0d0
      h_ne = (h_c+h_n)/2.0d0
      h_se = (h_c+h_s)/2.0d0
      u_we = (u_c+u_w)/2.0d0
      u_ee = (u_c+u_e)/2.0d0
      u_ne = (u_c+u_n)/2.0d0
      u_se = (u_c+u_s)/2.0d0
      v_we = (v_c+v_w)/2.0d0
      v_ee = (v_c+v_e)/2.0d0
      v_ne = (v_c+v_n)/2.0d0
      v_se = (v_c+v_s)/2.0d0
      hdfc2(i) = ((h_ee-h_we)*e2t(i) + (h_ne-h_se)*e1t(i)) / (e1t(i)*e2t(i))
      udfc2(i) = ((u_ee-u_we)*e2t(i) + (u_ne-u_se)*e1t(i)) / (e1t(i)*e2t(i))
      vdfc2(i) = ((v_ee-v_we)*e2t(i) + (v_ne-v_se)*e1t(i)) / (e1t(i)*e2t(i))
    end do


!=============================================
!--------compute eta_now, div, adv_t, vor-----
!============================================= 
    do i=i_st,i_en
      tn_1 = tn(i,1)
      tn_2 = tn(i,2)
      ts_1 = ts(i,1)
      ts_2 = ts(i,2)
      tw_1 = tw(i,1)
      tw_2 = tw(i,2)
      te_1 = te(i,1)
      te_2 = te(i,2)
      tmask_w = tmask(tw_1)*tmask(tw_2)
      tmask_s = tmask(ts_1)*tmask(ts_2)
      tmask_n = tmask(tn_1)*tmask(tn_2)
      tmask_e = tmask(te_1)*tmask(te_2)
      u_w = (ubar_bef(tw_1)+ubar_bef(tw_2))*0.5d0*tmask_w
      u_e = (ubar_bef(te_1)+ubar_bef(te_2))*0.5d0*tmask_e
      u_s = (ubar_bef(ts_1)+ubar_bef(ts_2))*0.5d0*tmask_s
      u_n = (ubar_bef(tn_1)+ubar_bef(tn_2))*0.5d0*tmask_n
      v_w = (vbar_bef(tw_1)+vbar_bef(tw_2))*0.5d0*tmask_w
      v_e = (vbar_bef(te_1)+vbar_bef(te_2))*0.5d0*tmask_e
      v_s = (vbar_bef(ts_1)+vbar_bef(ts_2))*0.5d0*tmask_s
      v_n = (vbar_bef(tn_1)+vbar_bef(tn_2))*0.5d0*tmask_n
      h_w = (eta_bef(tw_1)+eta_bef(tw_2))*0.5d0*tmask_w
      h_e = (eta_bef(te_1)+eta_bef(te_2))*0.5d0*tmask_e
      h_s = (eta_bef(ts_1)+eta_bef(ts_2))*0.5d0*tmask_s
      h_n = (eta_bef(tn_1)+eta_bef(tn_2))*0.5d0*tmask_n
      b_w = (bathy(tw_1)+bathy(tw_2))*0.5d0*tmask_w
      b_e = (bathy(te_1)+bathy(te_2))*0.5d0*tmask_e
      b_s = (bathy(ts_1)+bathy(ts_2))*0.5d0*tmask_s
      b_n = (bathy(tn_1)+bathy(tn_2))*0.5d0*tmask_n
      u_c = ubar_bef(i)*tmask(i)
      v_c = vbar_bef(i)*tmask(i)
      h_c = eta_bef(i)*tmask(i)
      b_c = bathy(i)*tmask(i)
      if ((tmask_e) .eq. 0.0d0) then
        h_e = h_c
        b_e = b_c
      end if
      if ((tmask_w) .eq. 0.0d0) then
        h_w = h_c
        b_w = b_c
      end if
      if ((tmask_n) .eq. 0.0d0) then
        h_n = h_c
        b_n = b_c
      end if
      if ((tmask_s) .eq. 0.0d0) then
        h_s = h_c
        b_s = b_c
      end if
      if (rn_shlat == 0) then
        if ((tmask_e) .eq. 0.0d0) then
          v_e = v_c
        end if
        if ((tmask_w) .eq. 0.0d0) then
          v_w = v_c
        end if
        if ((tmask_n) .eq. 0.0d0) then
          u_n = u_c
        end if
        if ((tmask_s) .eq. 0.0d0) then
          u_s = u_c
        end if
      else if (rn_shlat == 1) then
        if ((tmask_e) .eq. 0.0d0) then
          v_e = 0.0d0
        end if
        if ((tmask_w) .eq. 0.0d0) then
          v_w = 0.0d0
        end if
        if ((tmask_n) .eq. 0.0d0) then
          u_n = 0.0d0
        end if
        if ((tmask_s) .eq. 0.0d0) then
          u_s = 0.0d0
        end if
      else if (rn_shlat == 2) then
        if ((tmask_e) .eq. 0.0d0) then
          v_e = -1.0d0*v_c
        end if
        if ((tmask_w) .eq. 0.0d0) then
          v_w = -1.0d0*v_c
        end if
        if ((tmask_n) .eq. 0.0d0) then
          u_n = -1.0d0*u_c
        end if
        if ((tmask_s) .eq. 0.0d0) then
          u_s = -1.0d0*u_c
        end if
      end if
      u_we = tmask_w*(u_c*(h_c+b_c)+u_w*(h_w+b_w))/(h_c+b_c+h_w+b_w+eps)
      u_ee = tmask_e*(u_c*(h_c+b_c)+u_e*(h_e+b_e))/(h_c+b_c+h_e+b_e+eps)
      u_ne = tmask_n*(u_c*(h_c+b_c)+u_n*(h_n+b_n))/(h_c+b_c+h_n+b_n+eps)
      u_se = tmask_s*(u_c*(h_c+b_c)+u_s*(h_s+b_s))/(h_c+b_c+h_s+b_s+eps)
      v_we = tmask_w*(v_c*(h_c+b_c)+v_w*(h_w+b_w))/(h_c+b_c+h_w+b_w+eps)
      v_ee = tmask_e*(v_c*(h_c+b_c)+v_e*(h_e+b_e))/(h_c+b_c+h_e+b_e+eps)
      v_ne = tmask_n*(v_c*(h_c+b_c)+v_n*(h_n+b_n))/(h_c+b_c+h_n+b_n+eps)
      v_se = tmask_s*(v_c*(h_c+b_c)+v_s*(h_s+b_s))/(h_c+b_c+h_s+b_s+eps)
      h_we = (h_c+h_w)/2.0d0
      h_ee = (h_c+h_e)/2.0d0
      h_ne = (h_c+h_n)/2.0d0
      h_se = (h_c+h_s)/2.0d0
      ke_w = u_we**2+v_we**2
      ke_e = u_ee**2+v_ee**2
      ke_n = u_ne**2+v_ne**2
      ke_s = u_se**2+v_se**2
      
      div_bef(i) = ((u_ee-u_we)*e2t(i) + (v_ne-v_se)*e1t(i)) / (e1t(i)*e2t(i))
      vor_bef(i) = ((v_ee-v_we)*e2t(i) - (u_ne-u_se)*e1t(i)) / (e1t(i)*e2t(i)) + ff(i)
      vor_part_u(i) = v_c*vor_bef(i)*tmask(i)
      vor_part_v(i) =-u_c*vor_bef(i)*tmask(i)
      keg_part_u(i) = -0.5d0*tmask(i)*(ke_e - ke_w) / e1t(i)
      keg_part_v(i) = -0.5d0*tmask(i)*(ke_n - ke_s) / e2t(i)
      if ((lat_t(i) >= lat_polex2) .and. (lat_t(i) <= lat_polex3)) then
        divbar_unror(i) = div_bef(i)
        adv_part_t(i) = ((h_ee-h_we)*u_c*e2t(i) + (h_ne-h_se)*v_c*e1t(i)) / (e1t(i)*e2t(i))
        eta_now(i) = (hdfc2(i) - adv_part_t(i) - divbar_unror(i)*h_c)*dt*tmask(i) + h_c
      else
        u_w = (ubar_unror(tw_1)+ubar_unror(tw_2))*0.5d0*tmask_w
        u_e = (ubar_unror(te_1)+ubar_unror(te_2))*0.5d0*tmask_e
        v_s = (vbar_unror(ts_1)+vbar_unror(ts_2))*0.5d0*tmask_s
        v_n = (vbar_unror(tn_1)+vbar_unror(tn_2))*0.5d0*tmask_n
        u_c = ubar_unror(i)*tmask(i)
        v_c = vbar_unror(i)*tmask(i)
        u_we = tmask_w*(u_c*(h_c+b_c)+u_w*(h_w+b_w))/(h_c+b_c+h_w+b_w+eps)
        u_ee = tmask_e*(u_c*(h_c+b_c)+u_e*(h_e+b_e))/(h_c+b_c+h_e+b_e+eps)
        v_ne = tmask_n*(v_c*(h_c+b_c)+v_n*(h_n+b_n))/(h_c+b_c+h_n+b_n+eps)
        v_se = tmask_s*(v_c*(h_c+b_c)+v_s*(h_s+b_s))/(h_c+b_c+h_s+b_s+eps)
        divbar_unror(i) = ((u_ee-u_we)*e2t(i) + (v_ne-v_se)*e1t(i)) / (e1t(i)*e2t(i))
        adv_part_t(i) = ((h_ee-h_we)*u_c*e2t(i) + (h_ne-h_se)*v_c*e1t(i)) / (e1t(i)*e2t(i))
        eta_now(i) = (hdfc2(i) - adv_part_t(i) - divbar_unror(i)*h_c)*dt*tmask(i) + h_c
      end if 
    end do
    !$acc loop independent
    do i=nglo_b+1,nglo_b+nbdy_b
      j = idx_b_m(i-nglo_b)
  	  eta_now(i) = eta_now(j)
    end do
    if (i_st > 1) then
      eta_now(1) = (eta_now(2) + eta_now(3) + eta_now(4) + eta_now(5) + eta_now(6) + eta_now(7) &
                 & + eta_now(8) + eta_now(9))/8.0d0
    end if
    if (i_en < nglo_b) then
      eta_now(nglo_b) = (eta_now(nglo_b-1) + eta_now(nglo_b-2) + eta_now(nglo_b-3) + eta_now(nglo_b-4) &
                      & + eta_now(nglo_b-5) + eta_now(nglo_b-6) + eta_now(nglo_b-7) + eta_now(nglo_b-8))/8.0d0
    end if
    
!------------------------------------------------
!--------compute hpg_part ubar vbar--------------
!------------------------------------------------
    do i=i_st,i_en
      tn_1 = tn(i,1)
      tn_2 = tn(i,2)
      ts_1 = ts(i,1)
      ts_2 = ts(i,2)
      tw_1 = tw(i,1)
      tw_2 = tw(i,2)
      te_1 = te(i,1)
      te_2 = te(i,2)
      tmask_w = tmask(tw_1)*tmask(tw_2)
      tmask_s = tmask(ts_1)*tmask(ts_2)
      tmask_n = tmask(tn_1)*tmask(tn_2)
      tmask_e = tmask(te_1)*tmask(te_2)
      h_w = (eta_bef(tw_1)+eta_bef(tw_2))*0.5d0*tmask_w
      h_e = (eta_bef(te_1)+eta_bef(te_2))*0.5d0*tmask_e
      h_s = (eta_bef(ts_1)+eta_bef(ts_2))*0.5d0*tmask_s
      h_n = (eta_bef(tn_1)+eta_bef(tn_2))*0.5d0*tmask_n
      b_w = (bathy(tw_1)+bathy(tw_2))*0.5d0*tmask_w
      b_e = (bathy(te_1)+bathy(te_2))*0.5d0*tmask_e
      b_s = (bathy(ts_1)+bathy(ts_2))*0.5d0*tmask_s
      b_n = (bathy(tn_1)+bathy(tn_2))*0.5d0*tmask_n
      h_c = eta_bef(i)*tmask(i)
      b_c = bathy(i)*tmask(i)
      if ((tmask_e) .eq. 0.0d0) then
        h_e = h_c
        b_e = b_c
      end if
      if ((tmask_w) .eq. 0.0d0) then
        h_w = h_c
        b_w = b_c
      end if
      if ((tmask_n) .eq. 0.0d0) then
        h_n = h_c
        b_n = b_c
      end if
      if ((tmask_s) .eq. 0.0d0) then
        h_s = h_c
        b_s = b_c
      end if
      h_we = (h_c+b_c+h_w+b_w)/2.0d0
      h_ee = (h_c+b_c+h_e+b_e)/2.0d0
      h_ne = (h_c+b_c+h_n+b_n)/2.0d0
      h_se = (h_c+b_c+h_s+b_s)/2.0d0
      hpg_part_u(i) = -grav*tmask(i)*(h_ee-h_we)/e1t(i)
      hpg_part_v(i) = -grav*tmask(i)*(h_ne-h_se)/e2t(i)
      u_pol = u_arc_tra(i)*hpg_part_u(i) + u_arc_trb(i)*hpg_part_v(i)
    	v_pol = v_arc_tra(i)*hpg_part_u(i) + v_arc_trb(i)*hpg_part_v(i)
    	hpg_part_u(i) = u_pol
    	hpg_part_v(i) = v_pol
    	ubar_now(i) = (vor_part_u(i) + keg_part_u(i) + hpg_part_u(i) + udfc2(i))*dt*tmask(i) + ubar_bef(i)
    	vbar_now(i) = (vor_part_v(i) + keg_part_v(i) + hpg_part_v(i) + vdfc2(i))*dt*tmask(i) + vbar_bef(i)
    	!ubar_now(i) = (keg_part_u(i) + hpg_part_u(i))*dt*tmask(i) + ubar_bef(i)
    	!vbar_now(i) = (keg_part_v(i) + hpg_part_v(i))*dt*tmask(i) + vbar_bef(i)
    end do
    do i=i_st,i_en
      u_pol = u_arc_itra(i)*ubar_now(i) + u_arc_itrb(i)*vbar_now(i)
    	v_pol = v_arc_itra(i)*ubar_now(i) + v_arc_itrb(i)*vbar_now(i)
    	ubar_unror(i) = u_pol
    	vbar_unror(i) = v_pol
    end do
    do i=nglo_b+1,nglo_b+nbdy_b
      j = idx_b_m(i-nglo_b)
      temp1(i) = ubar_now(j)
  	  temp2(i) = vbar_now(j)
  	  u_pol = u_arc_itra(i)*temp1(i) + u_arc_itrb(i)*temp2(i)
    	v_pol = v_arc_itra(i)*temp1(i) + v_arc_itrb(i)*temp2(i)
    	temp1(i) = u_pol
    	temp2(i) = v_pol
    	temp3(i) = ubar_unror(j)
  	  temp4(i) = vbar_unror(j)
    end do
    do i=nglo_b+1,nglo_b+nbdy_b
      ubar_now(i) = temp1(i)
      vbar_now(i) = temp2(i)
      ubar_unror(i) = temp3(i)
  	  vbar_unror(i) = temp4(i)
    end do
    if (i_st > 1) then
      ubar_now(1) = (ubar_now(2) + ubar_now(3) + ubar_now(4) + ubar_now(5) + ubar_now(6) + ubar_now(7) &
                 & + ubar_now(8) + ubar_now(9))/8.0d0
      vbar_now(1) = (vbar_now(2) + vbar_now(3) + vbar_now(4) + vbar_now(5) + vbar_now(6) + vbar_now(7) &
                 & + vbar_now(8) + vbar_now(9))/8.0d0
      ubar_unror(1) = (ubar_unror(2) + ubar_unror(3) + ubar_unror(4) + ubar_unror(5) + ubar_unror(6) + ubar_unror(7) &
                 & + ubar_unror(8) + ubar_unror(9))/8.0d0
      vbar_unror(1) = (vbar_unror(2) + vbar_unror(3) + vbar_unror(4) + vbar_unror(5) + vbar_unror(6) + vbar_unror(7) &
                 & + vbar_unror(8) + vbar_unror(9))/8.0d0
    end if
    if (i_en < nglo_b) then
      ubar_now(nglo_b) = (ubar_now(nglo_b-1) + ubar_now(nglo_b-2) + ubar_now(nglo_b-3) + ubar_now(nglo_b-4) &
                      & + ubar_now(nglo_b-5) + ubar_now(nglo_b-6) + ubar_now(nglo_b-7) + ubar_now(nglo_b-8))/8.0d0
      vbar_now(nglo_b) = (vbar_now(nglo_b-1) + vbar_now(nglo_b-2) + vbar_now(nglo_b-3) + vbar_now(nglo_b-4) &
                      & + vbar_now(nglo_b-5) + vbar_now(nglo_b-6) + vbar_now(nglo_b-7) + vbar_now(nglo_b-8))/8.0d0
      ubar_unror(nglo_b) = (ubar_unror(nglo_b-1) + ubar_unror(nglo_b-2) + ubar_unror(nglo_b-3) + ubar_unror(nglo_b-4) &
                      & + ubar_unror(nglo_b-5) + ubar_unror(nglo_b-6) + ubar_unror(nglo_b-7) + ubar_unror(nglo_b-8))/8.0d0
      vbar_unror(nglo_b) = (vbar_unror(nglo_b-1) + vbar_unror(nglo_b-2) + vbar_unror(nglo_b-3) + vbar_unror(nglo_b-4) &
                      & + vbar_unror(nglo_b-5) + vbar_unror(nglo_b-6) + vbar_unror(nglo_b-7) + vbar_unror(nglo_b-8))/8.0d0
    end if

!==============Second Back====================     
!=============================================
!--------compute udfc, vdfc, hdfc-------------
!============================================= 
    do i=i_st,i_en
      tn_1 = tn(i,1)
      tn_2 = tn(i,2)
      ts_1 = ts(i,1)
      ts_2 = ts(i,2)
      tw_1 = tw(i,1)
      tw_2 = tw(i,2)
      te_1 = te(i,1)
      te_2 = te(i,2)
      tmask_w = tmask(tw_1)*tmask(tw_2)
      tmask_s = tmask(ts_1)*tmask(ts_2)
      tmask_n = tmask(tn_1)*tmask(tn_2)
      tmask_e = tmask(te_1)*tmask(te_2)
      u_w = (ubar_now(tw_1)+ubar_now(tw_2))*0.5d0*tmask_w
      u_e = (ubar_now(te_1)+ubar_now(te_2))*0.5d0*tmask_e
      u_s = (ubar_now(ts_1)+ubar_now(ts_2))*0.5d0*tmask_s
      u_n = (ubar_now(tn_1)+ubar_now(tn_2))*0.5d0*tmask_n
      v_w = (vbar_now(tw_1)+vbar_now(tw_2))*0.5d0*tmask_w
      v_e = (vbar_now(te_1)+vbar_now(te_2))*0.5d0*tmask_e
      v_s = (vbar_now(ts_1)+vbar_now(ts_2))*0.5d0*tmask_s
      v_n = (vbar_now(tn_1)+vbar_now(tn_2))*0.5d0*tmask_n
      h_w = (eta_now(tw_1)+eta_now(tw_2))*0.5d0*tmask_w
      h_e = (eta_now(te_1)+eta_now(te_2))*0.5d0*tmask_e
      h_s = (eta_now(ts_1)+eta_now(ts_2))*0.5d0*tmask_s
      h_n = (eta_now(tn_1)+eta_now(tn_2))*0.5d0*tmask_n
      b_w = (bathy(tw_1)+bathy(tw_2))*0.5d0*tmask_w
      b_e = (bathy(te_1)+bathy(te_2))*0.5d0*tmask_e
      b_s = (bathy(ts_1)+bathy(ts_2))*0.5d0*tmask_s
      b_n = (bathy(tn_1)+bathy(tn_2))*0.5d0*tmask_n
      u_c = ubar_now(i)*tmask(i)
      v_c = vbar_now(i)*tmask(i)
      h_c = eta_now(i)*tmask(i)
      b_c = bathy(i)*tmask(i)
      if ((tmask_e) .eq. 0.0d0) then
        h_e = h_c
        b_e = b_c
      end if
      if ((tmask_w) .eq. 0.0d0) then
        h_w = h_c
        b_w = b_c
      end if
      if ((tmask_n) .eq. 0.0d0) then
        h_n = h_c
        b_n = b_c
      end if
      if ((tmask_s) .eq. 0.0d0) then
        h_s = h_c
        b_s = b_c
      end if
      if (rn_shlat == 0) then
        if ((tmask_e) .eq. 0.0d0) then
          v_e = v_c
        end if
        if ((tmask_w) .eq. 0.0d0) then
          v_w = v_c
        end if
        if ((tmask_n) .eq. 0.0d0) then
          u_n = u_c
        end if
        if ((tmask_s) .eq. 0.0d0) then
          u_s = u_c
        end if
      else if (rn_shlat == 1) then
        if ((tmask_e) .eq. 0.0d0) then
          v_e = 0.0d0
        end if
        if ((tmask_w) .eq. 0.0d0) then
          v_w = 0.0d0
        end if
        if ((tmask_n) .eq. 0.0d0) then
          u_n = 0.0d0
        end if
        if ((tmask_s) .eq. 0.0d0) then
          u_s = 0.0d0
        end if
      else if (rn_shlat == 2) then
        if ((tmask_e) .eq. 0.0d0) then
          v_e = -1.0d0*v_c
        end if
        if ((tmask_w) .eq. 0.0d0) then
          v_w = -1.0d0*v_c
        end if
        if ((tmask_n) .eq. 0.0d0) then
          u_n = -1.0d0*u_c
        end if
        if ((tmask_s) .eq. 0.0d0) then
          u_s = -1.0d0*u_c
        end if
      end if
      u_we = tmask_w*(u_c*(h_c+b_c)+u_w*(h_w+b_w))/(h_c+b_c+h_w+b_w+eps)
      u_ee = tmask_e*(u_c*(h_c+b_c)+u_e*(h_e+b_e))/(h_c+b_c+h_e+b_e+eps)
      u_ne = tmask_n*(u_c*(h_c+b_c)+u_n*(h_n+b_n))/(h_c+b_c+h_n+b_n+eps)
      u_se = tmask_s*(u_c*(h_c+b_c)+u_s*(h_s+b_s))/(h_c+b_c+h_s+b_s+eps)
      v_we = tmask_w*(v_c*(h_c+b_c)+v_w*(h_w+b_w))/(h_c+b_c+h_w+b_w+eps)
      v_ee = tmask_e*(v_c*(h_c+b_c)+v_e*(h_e+b_e))/(h_c+b_c+h_e+b_e+eps)
      v_ne = tmask_n*(v_c*(h_c+b_c)+v_n*(h_n+b_n))/(h_c+b_c+h_n+b_n+eps)
      v_se = tmask_s*(v_c*(h_c+b_c)+v_s*(h_s+b_s))/(h_c+b_c+h_s+b_s+eps)
      h_we = (h_c+h_w)/2.0d0
      h_ee = (h_c+h_e)/2.0d0
      h_ne = (h_c+h_n)/2.0d0
      h_se = (h_c+h_s)/2.0d0
      udfc1(i) = alh(i) * ((u_ee-u_we)*e2t(i) + (u_ne-u_se)*e1t(i)) / (e1t(i)*e2t(i))
      vdfc1(i) = alh(i) * ((v_ee-v_we)*e2t(i) + (v_ne-v_se)*e1t(i)) / (e1t(i)*e2t(i))
      hdfc1(i) = alh(i) * ((h_ee-h_we)*e2t(i) + (h_ne-h_se)*e1t(i)) / (e1t(i)*e2t(i))
    end do
    !$acc loop independent
    do i=nglo_b+1,nglo_b+nbdy_b
      j = idx_b_m(i-nglo_b)
  	  hdfc1(i) = hdfc1(j)
  	  udfc1(i) = udfc1(j)
  	  vdfc1(i) = vdfc1(j)
    end do
    if (i_st > 1) then
      hdfc1(1) = (hdfc1(2) + hdfc1(3) + hdfc1(4) + hdfc1(5) + hdfc1(6) + hdfc1(7) &
                 & + hdfc1(8) + hdfc1(9))/8.0d0
      udfc1(1) = (udfc1(2) + udfc1(3) + udfc1(4) + udfc1(5) + udfc1(6) + udfc1(7) &
                 & + udfc1(8) + udfc1(9))/8.0d0
      vdfc1(1) = (vdfc1(2) + vdfc1(3) + vdfc1(4) + vdfc1(5) + vdfc1(6) + vdfc1(7) &
                 & + vdfc1(8) + vdfc1(9))/8.0d0
    end if
    if (i_en < nglo_b) then
      hdfc1(nglo_b) = (hdfc1(nglo_b-1) + hdfc1(nglo_b-2) + hdfc1(nglo_b-3) + hdfc1(nglo_b-4) &
                      & + hdfc1(nglo_b-5) + hdfc1(nglo_b-6) + hdfc1(nglo_b-7) + hdfc1(nglo_b-8))/8.0d0
      udfc1(nglo_b) = (udfc1(nglo_b-1) + udfc1(nglo_b-2) + udfc1(nglo_b-3) + udfc1(nglo_b-4) &
                      & + udfc1(nglo_b-5) + udfc1(nglo_b-6) + udfc1(nglo_b-7) + udfc1(nglo_b-8))/8.0d0
      vdfc1(nglo_b) = (vdfc1(nglo_b-1) + vdfc1(nglo_b-2) + vdfc1(nglo_b-3) + vdfc1(nglo_b-4) &
                      & + vdfc1(nglo_b-5) + vdfc1(nglo_b-6) + vdfc1(nglo_b-7) + vdfc1(nglo_b-8))/8.0d0
    end if
    do i=i_st,i_en
      tn_1 = tn(i,1)
      tn_2 = tn(i,2)
      ts_1 = ts(i,1)
      ts_2 = ts(i,2)
      tw_1 = tw(i,1)
      tw_2 = tw(i,2)
      te_1 = te(i,1)
      te_2 = te(i,2)
      tmask_w = tmask(tw_1)*tmask(tw_2)
      tmask_s = tmask(ts_1)*tmask(ts_2)
      tmask_n = tmask(tn_1)*tmask(tn_2)
      tmask_e = tmask(te_1)*tmask(te_2)
      h_w = (hdfc1(tw_1)+hdfc1(tw_2))*0.5d0*tmask_w
      h_e = (hdfc1(te_1)+hdfc1(te_2))*0.5d0*tmask_e
      h_s = (hdfc1(ts_1)+hdfc1(ts_2))*0.5d0*tmask_s
      h_n = (hdfc1(tn_1)+hdfc1(tn_2))*0.5d0*tmask_n
      h_c = hdfc1(i)*tmask(i)
      u_w = (udfc1(tw_1)+udfc1(tw_2))*0.5d0*tmask_w
      u_e = (udfc1(te_1)+udfc1(te_2))*0.5d0*tmask_e
      u_s = (udfc1(ts_1)+udfc1(ts_2))*0.5d0*tmask_s
      u_n = (udfc1(tn_1)+udfc1(tn_2))*0.5d0*tmask_n
      u_c = udfc1(i)*tmask(i)
      v_w = (vdfc1(tw_1)+vdfc1(tw_2))*0.5d0*tmask_w
      v_e = (vdfc1(te_1)+vdfc1(te_2))*0.5d0*tmask_e
      v_s = (vdfc1(ts_1)+vdfc1(ts_2))*0.5d0*tmask_s
      v_n = (vdfc1(tn_1)+vdfc1(tn_2))*0.5d0*tmask_n
      v_c = vdfc1(i)*tmask(i)
      if ((tmask_e) .eq. 0.0d0) then
        h_e = h_c
        u_e = u_c
        v_e = v_c
      end if
      if ((tmask_w) .eq. 0.0d0) then
        h_w = h_c
        u_w = u_c
        v_w = v_c
      end if
      if ((tmask_n) .eq. 0.0d0) then
        h_n = h_c
        u_n = u_c
        v_n = v_c
      end if
      if ((tmask_s) .eq. 0.0d0) then
        h_s = h_c
        u_s = u_c
        v_s = v_c
      end if
      h_we = (h_c+h_w)/2.0d0
      h_ee = (h_c+h_e)/2.0d0
      h_ne = (h_c+h_n)/2.0d0
      h_se = (h_c+h_s)/2.0d0
      u_we = (u_c+u_w)/2.0d0
      u_ee = (u_c+u_e)/2.0d0
      u_ne = (u_c+u_n)/2.0d0
      u_se = (u_c+u_s)/2.0d0
      v_we = (v_c+v_w)/2.0d0
      v_ee = (v_c+v_e)/2.0d0
      v_ne = (v_c+v_n)/2.0d0
      v_se = (v_c+v_s)/2.0d0
      hdfc2(i) = ((h_ee-h_we)*e2t(i) + (h_ne-h_se)*e1t(i)) / (e1t(i)*e2t(i))
      udfc2(i) = ((u_ee-u_we)*e2t(i) + (u_ne-u_se)*e1t(i)) / (e1t(i)*e2t(i))
      vdfc2(i) = ((v_ee-v_we)*e2t(i) + (v_ne-v_se)*e1t(i)) / (e1t(i)*e2t(i))
    end do
!=============================================
!--------compute eta_bef, div, adv_t, vor-----
!============================================= 
    do i=i_st,i_en
      tn_1 = tn(i,1)
      tn_2 = tn(i,2)
      ts_1 = ts(i,1)
      ts_2 = ts(i,2)
      tw_1 = tw(i,1)
      tw_2 = tw(i,2)
      te_1 = te(i,1)
      te_2 = te(i,2)
      tmask_w = tmask(tw_1)*tmask(tw_2)
      tmask_s = tmask(ts_1)*tmask(ts_2)
      tmask_n = tmask(tn_1)*tmask(tn_2)
      tmask_e = tmask(te_1)*tmask(te_2)
      u_w = (ubar_now(tw_1)+ubar_now(tw_2))*0.5d0*tmask_w
      u_e = (ubar_now(te_1)+ubar_now(te_2))*0.5d0*tmask_e
      u_s = (ubar_now(ts_1)+ubar_now(ts_2))*0.5d0*tmask_s
      u_n = (ubar_now(tn_1)+ubar_now(tn_2))*0.5d0*tmask_n
      v_w = (vbar_now(tw_1)+vbar_now(tw_2))*0.5d0*tmask_w
      v_e = (vbar_now(te_1)+vbar_now(te_2))*0.5d0*tmask_e
      v_s = (vbar_now(ts_1)+vbar_now(ts_2))*0.5d0*tmask_s
      v_n = (vbar_now(tn_1)+vbar_now(tn_2))*0.5d0*tmask_n
      h_w = (eta_now(tw_1)+eta_now(tw_2))*0.5d0*tmask_w
      h_e = (eta_now(te_1)+eta_now(te_2))*0.5d0*tmask_e
      h_s = (eta_now(ts_1)+eta_now(ts_2))*0.5d0*tmask_s
      h_n = (eta_now(tn_1)+eta_now(tn_2))*0.5d0*tmask_n
      b_w = (bathy(tw_1)+bathy(tw_2))*0.5d0*tmask_w
      b_e = (bathy(te_1)+bathy(te_2))*0.5d0*tmask_e
      b_s = (bathy(ts_1)+bathy(ts_2))*0.5d0*tmask_s
      b_n = (bathy(tn_1)+bathy(tn_2))*0.5d0*tmask_n
      u_c = ubar_now(i)*tmask(i)
      v_c = vbar_now(i)*tmask(i)
      h_c = eta_now(i)*tmask(i)
      b_c = bathy(i)*tmask(i)
      if ((tmask_e) .eq. 0.0d0) then
        h_e = h_c
        b_e = b_c
      end if
      if ((tmask_w) .eq. 0.0d0) then
        h_w = h_c
        b_w = b_c
      end if
      if ((tmask_n) .eq. 0.0d0) then
        h_n = h_c
        b_n = b_c
      end if
      if ((tmask_s) .eq. 0.0d0) then
        h_s = h_c
        b_s = b_c
      end if
      if (rn_shlat == 0) then
        if ((tmask_e) .eq. 0.0d0) then
          v_e = v_c
        end if
        if ((tmask_w) .eq. 0.0d0) then
          v_w = v_c
        end if
        if ((tmask_n) .eq. 0.0d0) then
          u_n = u_c
        end if
        if ((tmask_s) .eq. 0.0d0) then
          u_s = u_c
        end if
      else if (rn_shlat == 1) then
        if ((tmask_e) .eq. 0.0d0) then
          v_e = 0.0d0
        end if
        if ((tmask_w) .eq. 0.0d0) then
          v_w = 0.0d0
        end if
        if ((tmask_n) .eq. 0.0d0) then
          u_n = 0.0d0
        end if
        if ((tmask_s) .eq. 0.0d0) then
          u_s = 0.0d0
        end if
      else if (rn_shlat == 2) then
        if ((tmask_e) .eq. 0.0d0) then
          v_e = -1.0d0*v_c
        end if
        if ((tmask_w) .eq. 0.0d0) then
          v_w = -1.0d0*v_c
        end if
        if ((tmask_n) .eq. 0.0d0) then
          u_n = -1.0d0*u_c
        end if
        if ((tmask_s) .eq. 0.0d0) then
          u_s = -1.0d0*u_c
        end if
      end if
      u_we = tmask_w*(u_c*(h_c+b_c)+u_w*(h_w+b_w))/(h_c+b_c+h_w+b_w+eps)
      u_ee = tmask_e*(u_c*(h_c+b_c)+u_e*(h_e+b_e))/(h_c+b_c+h_e+b_e+eps)
      u_ne = tmask_n*(u_c*(h_c+b_c)+u_n*(h_n+b_n))/(h_c+b_c+h_n+b_n+eps)
      u_se = tmask_s*(u_c*(h_c+b_c)+u_s*(h_s+b_s))/(h_c+b_c+h_s+b_s+eps)
      v_we = tmask_w*(v_c*(h_c+b_c)+v_w*(h_w+b_w))/(h_c+b_c+h_w+b_w+eps)
      v_ee = tmask_e*(v_c*(h_c+b_c)+v_e*(h_e+b_e))/(h_c+b_c+h_e+b_e+eps)
      v_ne = tmask_n*(v_c*(h_c+b_c)+v_n*(h_n+b_n))/(h_c+b_c+h_n+b_n+eps)
      v_se = tmask_s*(v_c*(h_c+b_c)+v_s*(h_s+b_s))/(h_c+b_c+h_s+b_s+eps)
      h_we = (h_c+h_w)/2.0d0
      h_ee = (h_c+h_e)/2.0d0
      h_ne = (h_c+h_n)/2.0d0
      h_se = (h_c+h_s)/2.0d0
      ke_w = u_we**2+v_we**2
      ke_e = u_ee**2+v_ee**2
      ke_n = u_ne**2+v_ne**2
      ke_s = u_se**2+v_se**2
      
      div_bef(i) = ((u_ee-u_we)*e2t(i) + (v_ne-v_se)*e1t(i)) / (e1t(i)*e2t(i))
      vor_bef(i) = ((v_ee-v_we)*e2t(i) - (u_ne-u_se)*e1t(i)) / (e1t(i)*e2t(i)) + ff(i)
      vor_part_u(i) = v_c*vor_bef(i)*tmask(i)
      vor_part_v(i) =-u_c*vor_bef(i)*tmask(i)
      keg_part_u(i) = -0.5d0*tmask(i)*(ke_e - ke_w) / e1t(i)
      keg_part_v(i) = -0.5d0*tmask(i)*(ke_n - ke_s) / e2t(i)
      if ((lat_t(i) >= lat_polex2) .and. (lat_t(i) <= lat_polex3)) then
        divbar_unror(i) = div_bef(i)
        adv_part_t(i) = ((h_ee-h_we)*u_c*e2t(i) + (h_ne-h_se)*v_c*e1t(i)) / (e1t(i)*e2t(i))
        eta_bef(i) = (hdfc2(i) - adv_part_t(i) - divbar_unror(i)*h_c)*dt*tmask(i) + eta_bef(i)
      else
        u_w = (ubar_unror(tw_1)+ubar_unror(tw_2))*0.5d0*tmask_w
        u_e = (ubar_unror(te_1)+ubar_unror(te_2))*0.5d0*tmask_e
        v_s = (vbar_unror(ts_1)+vbar_unror(ts_2))*0.5d0*tmask_s
        v_n = (vbar_unror(tn_1)+vbar_unror(tn_2))*0.5d0*tmask_n
        u_c = ubar_unror(i)*tmask(i)
        v_c = vbar_unror(i)*tmask(i)
        u_we = tmask_w*(u_c*(h_c+b_c)+u_w*(h_w+b_w))/(h_c+b_c+h_w+b_w+eps)
        u_ee = tmask_e*(u_c*(h_c+b_c)+u_e*(h_e+b_e))/(h_c+b_c+h_e+b_e+eps)
        v_ne = tmask_n*(v_c*(h_c+b_c)+v_n*(h_n+b_n))/(h_c+b_c+h_n+b_n+eps)
        v_se = tmask_s*(v_c*(h_c+b_c)+v_s*(h_s+b_s))/(h_c+b_c+h_s+b_s+eps)
        divbar_unror(i) = ((u_ee-u_we)*e2t(i) + (v_ne-v_se)*e1t(i)) / (e1t(i)*e2t(i))
        adv_part_t(i) = ((h_ee-h_we)*u_c*e2t(i) + (h_ne-h_se)*v_c*e1t(i)) / (e1t(i)*e2t(i))
        eta_bef(i) = (hdfc2(i) - adv_part_t(i) - divbar_unror(i)*h_c)*dt*tmask(i) + eta_bef(i)
      end if 
    end do
    !$acc loop independent
    do i=nglo_b+1,nglo_b+nbdy_b
      j = idx_b_m(i-nglo_b)
  	  eta_bef(i) = eta_bef(j)
    end do
    if (i_st > 1) then
      eta_bef(1) = (eta_bef(2) + eta_bef(3) + eta_bef(4) + eta_bef(5) + eta_bef(6) + eta_bef(7) &
                 & + eta_bef(8) + eta_bef(9))/8.0d0
    end if
    if (i_en < nglo_b) then
      eta_bef(nglo_b) = (eta_bef(nglo_b-1) + eta_bef(nglo_b-2) + eta_bef(nglo_b-3) + eta_bef(nglo_b-4) &
                      & + eta_bef(nglo_b-5) + eta_bef(nglo_b-6) + eta_bef(nglo_b-7) + eta_bef(nglo_b-8))/8.0d0
    end if
    
!------------------------------------------------
!--------compute hpg_part ubar vbar--------------
!------------------------------------------------
    do i=i_st,i_en
      tn_1 = tn(i,1)
      tn_2 = tn(i,2)
      ts_1 = ts(i,1)
      ts_2 = ts(i,2)
      tw_1 = tw(i,1)
      tw_2 = tw(i,2)
      te_1 = te(i,1)
      te_2 = te(i,2)
      tmask_w = tmask(tw_1)*tmask(tw_2)
      tmask_s = tmask(ts_1)*tmask(ts_2)
      tmask_n = tmask(tn_1)*tmask(tn_2)
      tmask_e = tmask(te_1)*tmask(te_2)
      h_w = (eta_now(tw_1)+eta_now(tw_2))*0.5d0*tmask_w
      h_e = (eta_now(te_1)+eta_now(te_2))*0.5d0*tmask_e
      h_s = (eta_now(ts_1)+eta_now(ts_2))*0.5d0*tmask_s
      h_n = (eta_now(tn_1)+eta_now(tn_2))*0.5d0*tmask_n
      b_w = (bathy(tw_1)+bathy(tw_2))*0.5d0*tmask_w
      b_e = (bathy(te_1)+bathy(te_2))*0.5d0*tmask_e
      b_s = (bathy(ts_1)+bathy(ts_2))*0.5d0*tmask_s
      b_n = (bathy(tn_1)+bathy(tn_2))*0.5d0*tmask_n
      h_c = eta_now(i)*tmask(i)
      b_c = bathy(i)*tmask(i)
      if ((tmask_e) .eq. 0.0d0) then
        h_e = h_c
        b_e = b_c
      end if
      if ((tmask_w) .eq. 0.0d0) then
        h_w = h_c
        b_w = b_c
      end if
      if ((tmask_n) .eq. 0.0d0) then
        h_n = h_c
        b_n = b_c
      end if
      if ((tmask_s) .eq. 0.0d0) then
        h_s = h_c
        b_s = b_c
      end if
      h_we = (h_c+b_c+h_w+b_w)/2.0d0
      h_ee = (h_c+b_c+h_e+b_e)/2.0d0
      h_ne = (h_c+b_c+h_n+b_n)/2.0d0
      h_se = (h_c+b_c+h_s+b_s)/2.0d0
      hpg_part_u(i) = -grav*tmask(i)*(h_ee-h_we)/e1t(i)
      hpg_part_v(i) = -grav*tmask(i)*(h_ne-h_se)/e2t(i)
      u_pol = u_arc_tra(i)*hpg_part_u(i) + u_arc_trb(i)*hpg_part_v(i)
    	v_pol = v_arc_tra(i)*hpg_part_u(i) + v_arc_trb(i)*hpg_part_v(i)
    	hpg_part_u(i) = u_pol
    	hpg_part_v(i) = v_pol
    	ubar_bef(i) = (vor_part_u(i) + keg_part_u(i) + hpg_part_u(i) + udfc2(i))*dt*tmask(i) + ubar_bef(i)
    	vbar_bef(i) = (vor_part_v(i) + keg_part_v(i) + hpg_part_v(i) + vdfc2(i))*dt*tmask(i) + vbar_bef(i)
    end do
    do i=i_st,i_en
      u_pol = u_arc_itra(i)*ubar_bef(i) + u_arc_itrb(i)*vbar_bef(i)
    	v_pol = v_arc_itra(i)*ubar_bef(i) + v_arc_itrb(i)*vbar_bef(i)
    	ubar_unror(i) = u_pol
    	vbar_unror(i) = v_pol
    end do
    do i=nglo_b+1,nglo_b+nbdy_b
      j = idx_b_m(i-nglo_b)
      temp1(i) = ubar_bef(j)
  	  temp2(i) = vbar_bef(j)
  	  u_pol = u_arc_itra(i)*temp1(i) + u_arc_itrb(i)*temp2(i)
    	v_pol = v_arc_itra(i)*temp1(i) + v_arc_itrb(i)*temp2(i)
    	temp1(i) = u_pol
    	temp2(i) = v_pol
    	temp3(i) = ubar_unror(j)
  	  temp4(i) = vbar_unror(j)
    end do
    do i=nglo_b+1,nglo_b+nbdy_b
      ubar_bef(i) = temp1(i)
      vbar_bef(i) = temp2(i)
      ubar_unror(i) = temp3(i)
  	  vbar_unror(i) = temp4(i)
    end do
    if (i_st > 1) then
      ubar_bef(1) = (ubar_bef(2) + ubar_bef(3) + ubar_bef(4) + ubar_bef(5) + ubar_bef(6) + ubar_bef(7) &
                 & + ubar_bef(8) + ubar_bef(9))/8.0d0
      vbar_bef(1) = (vbar_bef(2) + vbar_bef(3) + vbar_bef(4) + vbar_bef(5) + vbar_bef(6) + vbar_bef(7) &
                 & + vbar_bef(8) + vbar_bef(9))/8.0d0
      ubar_unror(1) = (ubar_unror(2) + ubar_unror(3) + ubar_unror(4) + ubar_unror(5) + ubar_unror(6) + ubar_unror(7) &
                 & + ubar_unror(8) + ubar_unror(9))/8.0d0
      vbar_unror(1) = (vbar_unror(2) + vbar_unror(3) + vbar_unror(4) + vbar_unror(5) + vbar_unror(6) + vbar_unror(7) &
                 & + vbar_unror(8) + vbar_unror(9))/8.0d0
    end if
    if (i_en < nglo_b) then
      ubar_bef(nglo_b) = (ubar_bef(nglo_b-1) + ubar_bef(nglo_b-2) + ubar_bef(nglo_b-3) + ubar_bef(nglo_b-4) &
                      & + ubar_bef(nglo_b-5) + ubar_bef(nglo_b-6) + ubar_bef(nglo_b-7) + ubar_bef(nglo_b-8))/8.0d0
      vbar_bef(nglo_b) = (vbar_bef(nglo_b-1) + vbar_bef(nglo_b-2) + vbar_bef(nglo_b-3) + vbar_bef(nglo_b-4) &
                      & + vbar_bef(nglo_b-5) + vbar_bef(nglo_b-6) + vbar_bef(nglo_b-7) + vbar_bef(nglo_b-8))/8.0d0
      ubar_unror(nglo_b) = (ubar_unror(nglo_b-1) + ubar_unror(nglo_b-2) + ubar_unror(nglo_b-3) + ubar_unror(nglo_b-4) &
                      & + ubar_unror(nglo_b-5) + ubar_unror(nglo_b-6) + ubar_unror(nglo_b-7) + ubar_unror(nglo_b-8))/8.0d0
      vbar_unror(nglo_b) = (vbar_unror(nglo_b-1) + vbar_unror(nglo_b-2) + vbar_unror(nglo_b-3) + vbar_unror(nglo_b-4) &
                      & + vbar_unror(nglo_b-5) + vbar_unror(nglo_b-6) + vbar_unror(nglo_b-7) + vbar_unror(nglo_b-8))/8.0d0
    end if
    
    !$acc end kernels
        
    if (euler) then
      euler = .false.
    end if
       
!-------------start write result------------------------
    idx = mod(t,1000)
    if (idx == 0) then
      write(*,"(1x,a,i8,a,f9.0)") "write model compute reach time step = ", t," time = ", t*dt
    end if
    idx = mod(t,1000)
    if (t == 1) then
      idx = 0
    end if
    !idx = 0
    !$acc update self(eta_bef,vor_bef,divbar_unror,ubar_unror,vbar_unror) if (idx == 0)
    if (idx == 0) then
      write(*,"(1x,a,i8,a,f9.0)") "write model output at time step = ", t," time = ", t*dt
      tout = tout + 1
      time(1) = tout
      hw(:) = eta_bef(1:nglo_b)
      uw(:) = ubar_unror(1:nglo_b)
      vw(:) = vbar_unror(1:nglo_b)
      vorw(:) = vor_bef(1:nglo_b) - ff(1:nglo_b)
      divw(:) = divbar_unror(1:nglo_b)
      status = nf90_open(trim(fn_out), nf90_write, ncid_out)
      status = nf90_inq_varid(ncid_out, 'time', time_id)
      status = nf90_put_var(ncid_out, time_id, time, (/tout/), (/1/))
      status = nf90_inq_varid(ncid_out, 'h', h_id)
      status = nf90_put_var(ncid_out, h_id, hw, (/1,tout/), (/nglo_b,1/))
      status = nf90_inq_varid(ncid_out, 'ubar', u_id)
      status = nf90_put_var(ncid_out, u_id, uw, (/1,tout/), (/nglo_b,1/))
      status = nf90_inq_varid(ncid_out, 'vbar', v_id)
      status = nf90_put_var(ncid_out, v_id, vw, (/1,tout/), (/nglo_b,1/))
      status = nf90_inq_varid(ncid_out, 'vor', vor_id)
      status = nf90_put_var(ncid_out, vor_id, vorw, (/1,tout/), (/nglo_b,1/))
      status = nf90_inq_varid(ncid_out, 'div', div_id)
      status = nf90_put_var(ncid_out, div_id, divw, (/1,tout/), (/nglo_b,1/))
      status = nf90_close(ncid_out)
    end if
!-------------end write result------------------------
!=============end integral======================
  end do
  !$acc end data
  
  deallocate(idx_b_m, idx_m_b)
  deallocate(uw, vw, hw, vorw, divw)
  deallocate(vbar_f, ubar_f, ubar_unror, vbar_unror)
  deallocate(t_u, t_v)
  deallocate(u2_i, v2_j)
  deallocate(lat_t, lon_t, ff, bathy)
  deallocate(hdfc1, hdfc2, alh, udfc1, udfc2, vdfc1, vdfc2)
  deallocate(lat_u, lon_u, lat_v, lon_v)
  deallocate(u_arc_tra, u_arc_trb, v_arc_tra, v_arc_trb)
  deallocate(u_arc_itra, u_arc_itrb, v_arc_itra, v_arc_itrb)
  deallocate(u_arc_tra2, u_arc_trb2, v_arc_tra2, v_arc_trb2)
  deallocate(u_arc_itra2, u_arc_itrb2, v_arc_itra2, v_arc_itrb2)
  deallocate(e1t, e1u, e1v, e1f)
  deallocate(e2t, e2u, e2v, e2f)
  deallocate(tn, ts, tw, te)
  deallocate(ubar_bef, vbar_bef, eta_bef, vor_bef, div_bef, divbar_unror)
  deallocate(ubar_now, vbar_now, eta_now, vor_now, div_now)
  deallocate(temp, tmask, umask, vmask, fmask)
  deallocate(u_sbc, v_sbc)
  deallocate(vor_part_u, keg_part_u, hpg_part_u, diff_part_u)
  deallocate(vor_part_v, keg_part_v, hpg_part_v, diff_part_v)
  deallocate(sbc_part_u, bot_fric_u, adv_part_t)
  deallocate(sbc_part_v, bot_fric_v)
  deallocate(temp1, temp2, temp3, temp4, temp5, temp6)
  
  call system_clock ( clock_count2, clock_rate, clock_max )
  write (*,"(1x,a,f10.3,a)") 'Elapsed real time = ', real(clock_count2-clock_count1)/real(clock_rate),'s'
  
end program smc_sw_nopolar_nobathy
