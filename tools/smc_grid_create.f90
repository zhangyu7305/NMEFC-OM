program smc_grid_create
  use netcdf
  implicit none
  integer(kind=4), parameter :: fdiv=8, level2=1,level3=1
  real(kind=8), parameter :: stlon=0.0d0, enlon=360.0d0, stlat=-90.0d0, enlat=90.0d0, dlat=0.5d0
  real(kind=8), parameter :: rearth=6.37122d6, pie=3.141592654d0, omega = 7.2921d-5
  real(kind=8), parameter :: ratio_min=0.5d0, ratio_max=0.9d0, dratio_lim=0.01d0
  real(kind=8), parameter :: pol_s=-70.0d0, pol_n=65.0d0
    
  character(len=nf90_max_name) :: temp_str1, temp_str2
  integer(kind=4) :: i,j,k,k_start,nx,ny,ng_l1,ng_l1_bathy,ng_l2,ng_l3,nw_l1,nw_l2,nw_l3,f
  integer(kind=4) :: status, ncid, varid
  integer(kind=4) :: umask_n, vmask_e, umask_n_ls, vmask_e_ls
  real(kind=8) :: lat,dx,dy,dis
  real(kind=8) :: ratio_lim_temp, maxratio, minratio, avgratio, varratio
  real(kind=8) :: ratio_lim_select, varratio_select, avgratio_select
  real(kind=8) :: rupdown, rnm, rem, rwm, rsm
  real(kind=8) :: rad
  real(kind=8) :: dlat_temp, dlon_temp
  real(kind=8) :: lat_polex1, lat_polex2, lat_polex3, lat_polex4, min_polex2, min_polex3, min_polex
  logical :: first_run, alive
  
  integer(kind=4), dimension(:), allocatable :: smc_tgrid_back_l1_i, smc_tgrid_back_l1_j
  integer(kind=4), dimension(:), allocatable :: par_id, idx_bathy
  integer(kind=4), dimension(:), allocatable :: smc_tgrid_l1_mask, smc_ugrid_l1_mask, smc_vgrid_l1_mask, smc_fgrid_l1_mask
  integer(kind=4), dimension(:), allocatable :: smc_ugrid_l1_mask_ls, smc_vgrid_l1_mask_ls
  integer(kind=4), dimension(:,:), allocatable :: smc_tgrid_back_l1, smc_tgrid_back_l2, smc_tgrid_back_l3
  integer(kind=4), dimension(:,:), allocatable :: smc_tgrid_l1_n, smc_tgrid_l1_s, smc_tgrid_l1_e, smc_tgrid_l1_w
  integer(kind=4), dimension(:,:), allocatable :: smc_tgrid_l1_n_temp, smc_tgrid_l1_s_temp, smc_tgrid_l1_e_temp, smc_tgrid_l1_w_temp
  integer(kind=4), dimension(:,:), allocatable :: smc_tgrid_l2_n, smc_tgrid_l2_s, smc_tgrid_l2_e, smc_tgrid_l2_w
  integer(kind=4), dimension(:,:), allocatable :: smc_tgrid_l3_n, smc_tgrid_l3_s, smc_tgrid_l3_e, smc_tgrid_l3_w
  integer(kind=4), dimension(:,:), allocatable :: smc_tgrid_b_l1_n, smc_tgrid_b_l1_s, smc_tgrid_b_l1_e, smc_tgrid_b_l1_w
  real(kind=8), dimension(:), allocatable :: ratio,dlon,lat_t_l2,lon_t_l2,lat_t_l3,lon_t_l3
  real(kind=8), dimension(:), allocatable :: lat_t_l1,lat_u_l1,lat_v_l1,lat_f_l1
  real(kind=8), dimension(:), allocatable :: lon_t_l1,lon_u_l1,lon_v_l1,lon_f_l1
  real(kind=8), dimension(:), allocatable :: dlon_t_l1, dlat_t_l1, dlon_t_l2, dlat_t_l2, dlon_t_l3, dlat_t_l3
  real(kind=8), dimension(:), allocatable :: lon_t_b_l1,lon_u_b_l1,lon_v_b_l1,lon_f_b_l1
  real(kind=8), dimension(:), allocatable :: lat_t_b_l1,lat_u_b_l1,lat_v_b_l1,lat_f_b_l1,ff_f_b_l1
  real(kind=8), dimension(:), allocatable :: e1_t_b_l1,e1_u_b_l1,e1_v_b_l1,e1_f_b_l1
  real(kind=8), dimension(:), allocatable :: e2_t_b_l1,e2_u_b_l1,e2_v_b_l1,e2_f_b_l1
  real(kind=8), dimension(:), allocatable :: e1_t_l1,e1_u_l1,e1_v_l1,e1_f_l1
  real(kind=8), dimension(:), allocatable :: e2_t_l1,e2_u_l1,e2_v_l1,e2_f_l1
  real(kind=8), dimension(:), allocatable :: dlat_t_b_l1,dlon_t_b_l1
  real(kind=8), dimension(:), allocatable :: bathy_t_l1, bathy_t_b_l1
  real(kind=8), dimension(:), allocatable :: smc_tgrid_b_l1_mask, smc_ugrid_b_l1_mask, smc_vgrid_b_l1_mask, smc_fgrid_b_l1_mask
  real(kind=8), dimension(:,:), allocatable :: bathy_t_l1_2d
 
!=====define global mother grid meridional ny and dy info, these unchange 
  ny = int((enlat-stlat)/dlat)
  allocate(ratio(ny/2),dlon(ny+1))
  dy = dlat*2*pie*rearth/360.0d0
  rad = pie/180.0d0
  
!=====find the best smc divide parameters, the mean grid aspect ratio close to 1, 
!=====and decide on which latitude we divide one in two grid
  varratio_select = 999.9d0
  ratio_lim_temp = ratio_min
  do while(ratio_lim_temp<ratio_max)
    dlon(1) = 360.0d0/fdiv
    first_run = .true.
    i = 1
    do j=2,ny/2+1
      lat = stlat + (j-1.0d0)*dlat
      dx = dlon(1)*2*pie*rearth*cos(lat*pie/180.0d0)/360.0d0
      ratio(i) = dy/dx
      if (first_run) then
        if (ratio(i)<ratio_lim_temp) then
          write(*,*) "ratio limit too large! reduce ratio_lim or increase fdiv"
          exit
        end if
        first_run = .false.
      end if
      do while(ratio(i)<ratio_lim_temp)
        dlon(1) = dlon(1)/2
        dx = dlon(1)*2*pie*rearth*cos(lat*pie/180.0d0)/360.0d0
        ratio(i) = dy/dx
      end do
      i = i + 1
      !write(*,*) lat,dx,dy,dy/dx,dlon(1)
    end do
    
    call cal_var(varratio,ratio,ny/2)
    maxratio = maxval(ratio)
    minratio = minval(ratio)
    avgratio = sum(ratio)/(max(1,size(ratio)))
    
    if(abs(avgratio-1)<0.1) then
      if (varratio<varratio_select) then
        varratio_select = varratio
        ratio_lim_select = ratio_lim_temp
        avgratio_select = avgratio
      end if
    end if
    
    ratio_lim_temp = ratio_lim_temp + dratio_lim
  end do
  !write(*,*) ratio_lim_select,varratio_select,avgratio_select
  
!=====divide global zonal again use the best parameters
  dlon(1) = 360.0d0
  dlon(ny+1) = dlon(1)
  dlon(2) = 360.0d0/fdiv
  i = 1
  do j=2,ny/2+1
    lat = stlat + (j-1.0d0)*dlat
    dx = dlon(j)*2*pie*rearth*cos(lat*pie/180.0d0)/360.0d0
    ratio(i) = dy/dx
    do while(ratio(i)<ratio_lim_select)
      dlon(j) = dlon(j)/2
      dx = dlon(j)*2*pie*rearth*cos(lat*pie/180.0d0)/360.0d0
      ratio(i) = dy/dx
    end do
    dlon(ny+1-j+1) = dlon(j)
    !write(*,*) lat,dx,dy,ratio(i),dlon(j)
    dlon(j+1) = dlon(j)
    i = i + 1
  end do
  
!=====calculate and distribute number global smc grid level 1 into a 1D array
  ng_l1 = 0
  do j=1,ny+1
    ng_l1 = ng_l1 + int(360.0d0/dlon(j))
  end do
  allocate(dlat_t_l1(ng_l1),dlon_t_l1(ng_l1))
  allocate(lat_t_l1(ng_l1),lat_u_l1(ng_l1),lat_v_l1(ng_l1),lat_f_l1(ng_l1))
  allocate(lon_t_l1(ng_l1),lon_u_l1(ng_l1),lon_v_l1(ng_l1),lon_f_l1(ng_l1))
  lat_t_l1 = 99999.9d0
  lon_t_l1 = 99999.9d0
  
  nx = int(360.0d0/minval(dlon))
  allocate(smc_tgrid_back_l1(nx,ny+1))
  allocate(smc_tgrid_back_l1_i(ng_l1),smc_tgrid_back_l1_j(ng_l1))
  smc_tgrid_back_l1 = -9999
  
  lat_polex2 = -70.0d0
  lat_polex3 = 65.0d0
  min_polex2 = 100000.0d0
  min_polex3 = 100000.0d0
  k = 1
  do j=1,ny+1
    nx = int(360.0d0/dlon(j))
    lat = stlat + (j-1.0d0)*dlat
    do i=1,nx
      lat_t_l1(k) = lat
      lon_t_l1(k) = stlon + dlon(j)/2 + (i-1)*dlon(j)
      dlat_t_l1(k) = dlat
      dlon_t_l1(k) = dlon(j)
      smc_tgrid_back_l1(i,j) = k
      smc_tgrid_back_l1_i(k) = i
      smc_tgrid_back_l1_j(k) = j
      k = k + 1
    end do
    write(*,"(1x,i8,f8.3,1x,f8.3)") nx,lat_t_l1(k-1),lon_t_l1(k-1)
    min_polex = abs(lat_t_l1(k-1)-pol_s)
    if (min_polex<=min_polex2) then
      min_polex2 = min_polex
      lat_polex2 = lat_t_l1(k-1)
      lat_polex1 = lat_polex2 - dlat
    end if
    min_polex = abs(lat_t_l1(k-1)-pol_n)
    if (min_polex<=min_polex3) then
      min_polex3 = min_polex
      lat_polex3 = lat_t_l1(k-1)
      lat_polex4 = lat_polex3 + dlat
    end if
  end do
  
  write(*,"(1x,a,4f8.3)") "polar exchange zone position:", lat_polex1, lat_polex2, lat_polex3, lat_polex4
  
  i = 1
  lat_u_l1(i) = lat_t_l1(i)
  lon_u_l1(i) = lon_t_l1(i)
  lat_v_l1(i) = lat_t_l1(i)
  lon_v_l1(i) = lon_t_l1(i)
  lat_f_l1(i) = lat_t_l1(i)
  lon_f_l1(i) = lon_t_l1(i)
  do i=2,ng_l1-1
    lat_u_l1(i) = lat_t_l1(i)
  	lon_u_l1(i) = lon_t_l1(i) + dlon_t_l1(i)/2
  	if (lon_u_l1(i) >= 360.0) then
  	  lon_u_l1(i) = lon_u_l1(i) - 360.0
  	end if

  	lat_v_l1(i) = lat_t_l1(i) + dlat_t_l1(i)/2
  	lon_v_l1(i) = lon_t_l1(i)

  	lat_f_l1(i) = lat_v_l1(i)
  	lon_f_l1(i) = lon_u_l1(i)
  end do
  i = ng_l1
  lat_u_l1(i) = lat_t_l1(i)
  lon_u_l1(i) = lon_t_l1(i)
  lat_v_l1(i) = lat_t_l1(i)
  lon_v_l1(i) = lon_t_l1(i)
  lat_f_l1(i) = lat_t_l1(i)
  lon_f_l1(i) = lon_t_l1(i)
  
!=====find neighbors of every smc grid level 1
  allocate(smc_tgrid_l1_n(ng_l1,2), smc_tgrid_l1_s(ng_l1,2), smc_tgrid_l1_e(ng_l1,2), smc_tgrid_l1_w(ng_l1,2))
  allocate(smc_tgrid_l1_n_temp(ng_l1,2), smc_tgrid_l1_s_temp(ng_l1,2), smc_tgrid_l1_e_temp(ng_l1,2), smc_tgrid_l1_w_temp(ng_l1,2))
  smc_tgrid_l1_n = -9999
  smc_tgrid_l1_s = -9999
  smc_tgrid_l1_e = -9999
  smc_tgrid_l1_w = -9999
  smc_tgrid_l1_n(1,1) = smc_tgrid_back_l1(1,2)
  smc_tgrid_l1_n(1,2) = smc_tgrid_back_l1(fdiv,2)
  smc_tgrid_l1_s(1,1) = smc_tgrid_back_l1(1,2)
  smc_tgrid_l1_s(1,2) = smc_tgrid_back_l1(fdiv,2)
  smc_tgrid_l1_e(1,1) = smc_tgrid_back_l1(1,1)
  smc_tgrid_l1_e(1,2) = smc_tgrid_back_l1(1,1)
  smc_tgrid_l1_w(1,1) = smc_tgrid_back_l1(1,1)
  smc_tgrid_l1_w(1,2) = smc_tgrid_back_l1(1,1)
  smc_tgrid_l1_n(ng_l1,1) = smc_tgrid_back_l1(1,ny)
  smc_tgrid_l1_n(ng_l1,2) = smc_tgrid_back_l1(fdiv,ny)
  smc_tgrid_l1_s(ng_l1,1) = smc_tgrid_back_l1(1,ny)
  smc_tgrid_l1_s(ng_l1,2) = smc_tgrid_back_l1(fdiv,ny)
  smc_tgrid_l1_e(ng_l1,1) = smc_tgrid_back_l1(1,ny+1)
  smc_tgrid_l1_e(ng_l1,2) = smc_tgrid_back_l1(1,ny+1)
  smc_tgrid_l1_w(ng_l1,1) = smc_tgrid_back_l1(1,ny+1)
  smc_tgrid_l1_w(ng_l1,2) = smc_tgrid_back_l1(1,ny+1)
  
  !===find neighbor north====
  do j=2,ny-1
    nx = int(360.0d0/dlon(j))
    rupdown = dlon(j+1)/dlon(j)
    do i=1,nx
      if (rupdown<0.6d0) then
        smc_tgrid_l1_n(smc_tgrid_back_l1(i,j),1) = smc_tgrid_back_l1(i*2-1,j+1)
        smc_tgrid_l1_n(smc_tgrid_back_l1(i,j),2) = smc_tgrid_back_l1(i*2  ,j+1)
      else if (rupdown>1.9d0) then
        smc_tgrid_l1_n(smc_tgrid_back_l1(i,j),1) = smc_tgrid_back_l1((i+1)/2,j+1)
        smc_tgrid_l1_n(smc_tgrid_back_l1(i,j),2) = smc_tgrid_back_l1((i+1)/2,j+1)
      else
        smc_tgrid_l1_n(smc_tgrid_back_l1(i,j),1) = smc_tgrid_back_l1(i,j+1)
        smc_tgrid_l1_n(smc_tgrid_back_l1(i,j),2) = smc_tgrid_back_l1(i,j+1)
      end if
    end do
  end do
  j=ny
  nx = int(360.0d0/dlon(j))
  do i=1,nx
    smc_tgrid_l1_n(smc_tgrid_back_l1(i,j),1) = smc_tgrid_back_l1(1,j+1)
    smc_tgrid_l1_n(smc_tgrid_back_l1(i,j),2) = smc_tgrid_back_l1(1,j+1)
  end do
  
  write(*,*) maxval(smc_tgrid_l1_n),minval(smc_tgrid_l1_n)
  
  !===find neighbor south====
  do j=3,ny
    nx = int(360.0d0/dlon(j))
    rupdown = dlon(j)/dlon(j-1)
    do i=1,nx
      if (rupdown<0.6d0) then
        smc_tgrid_l1_s(smc_tgrid_back_l1(i,j),1) = smc_tgrid_back_l1((i+1)/2,j-1)
        smc_tgrid_l1_s(smc_tgrid_back_l1(i,j),2) = smc_tgrid_back_l1((i+1)/2,j-1)
      else if (rupdown>1.9d0) then
        smc_tgrid_l1_s(smc_tgrid_back_l1(i,j),1) = smc_tgrid_back_l1(i*2-1,j-1)
        smc_tgrid_l1_s(smc_tgrid_back_l1(i,j),2) = smc_tgrid_back_l1(i*2  ,j-1)
      else
        smc_tgrid_l1_s(smc_tgrid_back_l1(i,j),1) = smc_tgrid_back_l1(i,j-1)
        smc_tgrid_l1_s(smc_tgrid_back_l1(i,j),2) = smc_tgrid_back_l1(i,j-1)
      end if
    end do
  end do
  j=2
  nx = int(360.0d0/dlon(j))
  do i=1,nx
    smc_tgrid_l1_s(smc_tgrid_back_l1(i,j),1) = smc_tgrid_back_l1(1,j-1)
    smc_tgrid_l1_s(smc_tgrid_back_l1(i,j),2) = smc_tgrid_back_l1(1,j-1)
  end do
  
  write(*,*) maxval(smc_tgrid_l1_s),minval(smc_tgrid_l1_s)
  
  !===find neighbor east====
  do j=2,ny
    nx = int(360.0d0/dlon(j))
    do i=1,nx-1
      smc_tgrid_l1_e(smc_tgrid_back_l1(i,j),1) = smc_tgrid_back_l1(i+1,j)
      smc_tgrid_l1_e(smc_tgrid_back_l1(i,j),2) = smc_tgrid_back_l1(i+1,j)
    end do
    i=nx
    smc_tgrid_l1_e(smc_tgrid_back_l1(i,j),1) = smc_tgrid_back_l1(1,j)
    smc_tgrid_l1_e(smc_tgrid_back_l1(i,j),2) = smc_tgrid_back_l1(1,j)
  end do
  
  write(*,*) maxval(smc_tgrid_l1_e),minval(smc_tgrid_l1_e)
  
  !===find neighbor west====
  do j=2,ny
    nx = int(360.0d0/dlon(j))
    do i=2,nx
      smc_tgrid_l1_w(smc_tgrid_back_l1(i,j),1) = smc_tgrid_back_l1(i-1,j)
      smc_tgrid_l1_w(smc_tgrid_back_l1(i,j),2) = smc_tgrid_back_l1(i-1,j)
    end do
    i=1
    smc_tgrid_l1_w(smc_tgrid_back_l1(i,j),1) = smc_tgrid_back_l1(nx,j)
    smc_tgrid_l1_w(smc_tgrid_back_l1(i,j),2) = smc_tgrid_back_l1(nx,j)
  end do
  
  write(*,*) maxval(smc_tgrid_l1_w),minval(smc_tgrid_l1_w)

!=====calculate and distribute number global smc grid level 2 into a 1D array
  if (level2 == 1) then
    nx = int(360.0d0/minval(dlon))
    allocate(smc_tgrid_back_l2(2*nx,2*(ny+1-2)+2))
    smc_tgrid_back_l2 = -9999
    do k=2, ng_l1
      dis = sqrt((lon_t_l1(k)-120.0d0)**2+(lat_t_l1(k)-30.0d0)**2)
      if (dis<20) then
        rnm = dlon_t_l1(smc_tgrid_l1_n(k,1))/dlon_t_l1(k)
        rsm = dlon_t_l1(smc_tgrid_l1_s(k,1))/dlon_t_l1(k)
        rem = dlat_t_l1(smc_tgrid_l1_e(k,1))/dlat_t_l1(k)
        rwm = dlat_t_l1(smc_tgrid_l1_w(k,1))/dlat_t_l1(k)
        if ((rnm<1.5d0) .and. (rsm<1.5d0) .and. (rem<1.5d0) .and. (rwm<1.5d0)) then
          smc_tgrid_back_l1(smc_tgrid_back_l1_i(k),smc_tgrid_back_l1_j(k)) = -4
          smc_tgrid_back_l2(smc_tgrid_back_l1_i(k)*2,smc_tgrid_back_l1_j(k)*2-1) = 1
          smc_tgrid_back_l2(smc_tgrid_back_l1_i(k)*2,smc_tgrid_back_l1_j(k)*2-2) = 1
          smc_tgrid_back_l2(smc_tgrid_back_l1_i(k)*2-1,smc_tgrid_back_l1_j(k)*2-1) = 1
          smc_tgrid_back_l2(smc_tgrid_back_l1_i(k)*2-1,smc_tgrid_back_l1_j(k)*2-2) = 1
        end if
      end if
    end do
  end if
  
  ng_l2 = 0
  do j=1, 2*(ny+1-2)+2
    do i=1, 2*nx
      if (smc_tgrid_back_l2(i,j) == 1) then
        ng_l2 = ng_l2 + 1
        smc_tgrid_back_l2(i,j) = ng_l2
      end if  
    end do
  end do

!=====calculate e1 and e2 for global
  allocate(e1_t_l1(ng_l1),e1_u_l1(ng_l1),e1_v_l1(ng_l1),e1_f_l1(ng_l1))
  allocate(e2_t_l1(ng_l1),e2_u_l1(ng_l1),e2_v_l1(ng_l1),e2_f_l1(ng_l1))
  do i=1,ng_l1
    dlat_temp = abs(lat_t_l1(i)-lat_t_l1(smc_tgrid_l1_s(i,1)))
    dlon_temp = abs(lon_t_l1(i)-lon_t_l1(smc_tgrid_l1_w(i,1)))
    if (dlon_temp > 180) dlon_temp = abs(dlon_temp - 360)
    e1_t_l1(i) = rearth*dlon_temp*rad*cos(rad*lat_t_l1(i))
    e2_t_l1(i) = rearth*dlat_temp*rad
    
    dlat_temp = abs(lat_f_l1(i)-lat_f_l1(smc_tgrid_l1_s(i,1)))
    dlon_temp = abs(lon_t_l1(i)-lon_t_l1(smc_tgrid_l1_e(i,1)))
    if (dlon_temp > 180) dlon_temp = abs(dlon_temp - 360)
    e1_u_l1(i) = rearth*dlon_temp*rad*cos(rad*lat_t_l1(i))
    e2_u_l1(i) = rearth*dlat_temp*rad
    
    dlat_temp = abs(lat_t_l1(i)-lat_t_l1(smc_tgrid_l1_n(i,1)))
    dlon_temp = abs(lon_f_l1(i)-lon_f_l1(smc_tgrid_l1_w(i,1)))
    if (dlon_temp > 180) dlon_temp = abs(dlon_temp - 360)
    e1_v_l1(i) = rearth*dlon_temp*rad*cos(rad*lat_t_l1(i))
    e2_v_l1(i) = rearth*dlat_temp*rad
    
    dlat_temp = abs(lat_u_l1(i)-lat_u_l1(smc_tgrid_l1_n(i,1)))
    dlon_temp = abs(lon_v_l1(i)-lon_v_l1(smc_tgrid_l1_e(i,1)))
    if (dlon_temp > 180) dlon_temp = abs(dlon_temp - 360)
    e1_f_l1(i) = rearth*dlon_temp*rad*cos(rad*lat_t_l1(i))
    e2_f_l1(i) = rearth*dlat_temp*rad
  end do

!=====calculate and distribute number global smc grid level 3 into a 1D array
  if (level3 == 1) then
    nx = int(360.0d0/minval(dlon))
    allocate(smc_tgrid_back_l3(4*nx,4*(ny+1-2)+2))
    smc_tgrid_back_l3 = -9999
    
  end if

!=====modification with bathy
  allocate(bathy_t_l1_2d(ng_l1,2),bathy_t_l1(ng_l1))
  allocate(smc_tgrid_l1_mask(ng_l1),smc_ugrid_l1_mask(ng_l1),smc_vgrid_l1_mask(ng_l1),smc_fgrid_l1_mask(ng_l1))
  allocate(smc_ugrid_l1_mask_ls(ng_l1),smc_vgrid_l1_mask_ls(ng_l1))
  inquire(file="smc_bathy.nc", exist=alive)
  if (alive) then
    status = nf90_open("smc_bathy.nc", nf90_nowrite, ncid)
    status = nf90_inq_varid(ncid, "smc_bathy", varid)
    status = nf90_get_var(ncid, varid, bathy_t_l1_2d, (/1,1,1/), (/ng_l1,2,1/))
    bathy_t_l1(:) = -1*bathy_t_l1_2d(:,1)
  else
  	bathy_t_l1 = 5000.0d0
  	do i=1,ng_l1
  	  !if (abs(lat_t_l1(i)) > 60.0) then
  	  !	bathy_t_l1(i) = 0.0d0
  	  !end if
  	  !if (abs(lon_t_l1(i)) > 300.0) then
  	  !  bathy_t_l1(i) = 0.0d0
  	  !end if
  	  !if (abs(lon_t_l1(i)) < 90.0) then
  	  !  bathy_t_l1(i) = 0.0d0
  	  !end if
  	  !if (lat_t_l1(i) > 70.0) then
  	  !	bathy_t_l1(i) = 0.0d0
  	  !end if
  	  !if (lat_t_l1(i) < 30.0) then
  	  !	bathy_t_l1(i) = 0.0d0
  	  !end if
  	end do
  end if
  smc_tgrid_l1_mask(:) = 0
  smc_ugrid_l1_mask(:) = 0
  smc_vgrid_l1_mask(:) = 0
  smc_fgrid_l1_mask(:) = 0
!-----set the first water grid  
  smc_tgrid_l1_mask(3*ng_l1/4) = 1
!-----------------------------
  do i=1,ng_l1
    do j=1,ng_l1
      if (smc_tgrid_l1_mask(j)>0) then
        do k=1,2
          if (smc_tgrid_l1_n(j,k)>0) then
            if (bathy_t_l1(smc_tgrid_l1_n(j,k))>0) then
              smc_tgrid_l1_mask(smc_tgrid_l1_n(j,k)) = 1
            !else
            !  smc_tgrid_l1_n(j,k) = -99999
            end if
          end if
          if (smc_tgrid_l1_s(j,k)>0) then
            if (bathy_t_l1(smc_tgrid_l1_s(j,k))>0) then
              smc_tgrid_l1_mask(smc_tgrid_l1_s(j,k)) = 1
            !else
            !  smc_tgrid_l1_s(j,k) = -99999
            end if
          end if
          if (smc_tgrid_l1_w(j,k)>0) then
            if (bathy_t_l1(smc_tgrid_l1_w(j,k))>0) then
              smc_tgrid_l1_mask(smc_tgrid_l1_w(j,k)) = 1
            !else
            !  smc_tgrid_l1_w(j,k) = -99999
            end if
          end if
          if (smc_tgrid_l1_e(j,k)>0) then
            if (bathy_t_l1(smc_tgrid_l1_e(j,k))>0) then
              smc_tgrid_l1_mask(smc_tgrid_l1_e(j,k)) = 1
            !else
            !  smc_tgrid_l1_e(j,k) = -99999
            end if
          end if
        end do
      end if
    end do
  end do
  
  do i=1,ng_l1
    if ((smc_tgrid_l1_mask(i)*smc_tgrid_l1_mask(smc_tgrid_l1_e(i,1))*smc_tgrid_l1_mask(smc_tgrid_l1_e(i,2))) .eq. 0) then
      smc_ugrid_l1_mask(i) = 0
    else
      smc_ugrid_l1_mask(i) = 1
    end if
    if ((smc_tgrid_l1_mask(i)*smc_tgrid_l1_mask(smc_tgrid_l1_n(i,1))*smc_tgrid_l1_mask(smc_tgrid_l1_n(i,2))) .eq. 0) then
      smc_vgrid_l1_mask(i) = 0
    else
      smc_vgrid_l1_mask(i) = 1
    end if
  end do
  
  do i=1,ng_l1
    if ((smc_tgrid_l1_mask(i)+smc_tgrid_l1_mask(smc_tgrid_l1_e(i,1))*smc_tgrid_l1_mask(smc_tgrid_l1_e(i,2))) .gt. 0) then
      smc_ugrid_l1_mask_ls(i) = 1
    else
      smc_ugrid_l1_mask_ls(i) = 0
    end if
    if ((smc_tgrid_l1_mask(i)+smc_tgrid_l1_mask(smc_tgrid_l1_n(i,1))*smc_tgrid_l1_mask(smc_tgrid_l1_n(i,2))) .gt. 0) then
      smc_vgrid_l1_mask_ls(i) = 1
    else
      smc_vgrid_l1_mask_ls(i) = 0
    end if
  end do
  
  do i=1,ng_l1
    umask_n = smc_ugrid_l1_mask(smc_tgrid_l1_n(i,1))*smc_ugrid_l1_mask(smc_tgrid_l1_n(i,2))
    vmask_e = smc_vgrid_l1_mask(smc_tgrid_l1_e(i,1))*smc_vgrid_l1_mask(smc_tgrid_l1_e(i,2))
    umask_n_ls = smc_ugrid_l1_mask_ls(smc_tgrid_l1_n(i,1))*smc_ugrid_l1_mask_ls(smc_tgrid_l1_n(i,2))
    vmask_e_ls = smc_vgrid_l1_mask_ls(smc_tgrid_l1_e(i,1))*smc_vgrid_l1_mask_ls(smc_tgrid_l1_e(i,2))
    if ((smc_ugrid_l1_mask(i)*umask_n*vmask_e) .eq. 0) then
      if ((umask_n_ls + vmask_e_ls + smc_ugrid_l1_mask_ls(i) + smc_vgrid_l1_mask_ls(i)) .eq. 0) then
        smc_fgrid_l1_mask(i) = 0
      else
        smc_fgrid_l1_mask(i) = 2
      end if
    else
      smc_fgrid_l1_mask(i) = 1
    end if
  end do
  
  ng_l1_bathy = 0
  do i=1,ng_l1
    if (smc_fgrid_l1_mask(i)>0) then
      ng_l1_bathy = ng_l1_bathy + 1
    end if
  end do
  
  allocate(idx_bathy(ng_l1))
  idx_bathy = -99999
  j = 1
  do i=1,ng_l1
    if (smc_fgrid_l1_mask(i)>0) then
      idx_bathy(i) = j
      j = j + 1
    end if
  end do

  smc_tgrid_l1_n_temp(:,:) = smc_tgrid_l1_n(:,:)
  smc_tgrid_l1_s_temp(:,:) = smc_tgrid_l1_s(:,:)
  smc_tgrid_l1_w_temp(:,:) = smc_tgrid_l1_w(:,:)
  smc_tgrid_l1_e_temp(:,:) = smc_tgrid_l1_e(:,:)

  do i=1,ng_l1
    do k=1,2
      !if (smc_tgrid_l1_n(i,k) > 0) then
        smc_tgrid_l1_n_temp(i,k) = idx_bathy(smc_tgrid_l1_n(i,k))
      !end if
      !if (smc_tgrid_l1_s(i,k) > 0) then
        smc_tgrid_l1_s_temp(i,k) = idx_bathy(smc_tgrid_l1_s(i,k))
      !end if
      !if (smc_tgrid_l1_w(i,k) > 0) then
        smc_tgrid_l1_w_temp(i,k) = idx_bathy(smc_tgrid_l1_w(i,k))
      !end if
      !if (smc_tgrid_l1_e(i,k) > 0) then
        smc_tgrid_l1_e_temp(i,k) = idx_bathy(smc_tgrid_l1_e(i,k))
      !end if
    end do
  end do

!=====arrange grid info (lat lon and e1 e2)===========

  allocate(smc_tgrid_b_l1_n(ng_l1_bathy,2), smc_tgrid_b_l1_s(ng_l1_bathy,2), & 
         & smc_tgrid_b_l1_e(ng_l1_bathy,2), smc_tgrid_b_l1_w(ng_l1_bathy,2))
  allocate(dlat_t_b_l1(ng_l1_bathy),dlon_t_b_l1(ng_l1_bathy))
  allocate(lat_t_b_l1(ng_l1_bathy),lat_u_b_l1(ng_l1_bathy),lat_v_b_l1(ng_l1_bathy),lat_f_b_l1(ng_l1_bathy))
  allocate(lon_t_b_l1(ng_l1_bathy),lon_u_b_l1(ng_l1_bathy),lon_v_b_l1(ng_l1_bathy),lon_f_b_l1(ng_l1_bathy))
  allocate(smc_tgrid_b_l1_mask(ng_l1_bathy),smc_ugrid_b_l1_mask(ng_l1_bathy), &
         & smc_vgrid_b_l1_mask(ng_l1_bathy),smc_fgrid_b_l1_mask(ng_l1_bathy))
  allocate(bathy_t_b_l1(ng_l1_bathy),ff_f_b_l1(ng_l1_bathy))
  allocate(e1_t_b_l1(ng_l1_bathy),e2_t_b_l1(ng_l1_bathy))
  allocate(e1_u_b_l1(ng_l1_bathy),e2_u_b_l1(ng_l1_bathy))
  allocate(e1_v_b_l1(ng_l1_bathy),e2_v_b_l1(ng_l1_bathy))
  allocate(e1_f_b_l1(ng_l1_bathy),e2_f_b_l1(ng_l1_bathy))
  j = 1
  do i=1,ng_l1
    if (smc_fgrid_l1_mask(i)>0) then
      lat_t_b_l1(j) = lat_t_l1(i)
      lon_t_b_l1(j) = lon_t_l1(i)
      lat_u_b_l1(j) = lat_u_l1(i)
      lon_u_b_l1(j) = lon_u_l1(i)
      lat_v_b_l1(j) = lat_v_l1(i)
      lon_v_b_l1(j) = lon_v_l1(i)
      lat_f_b_l1(j) = lat_f_l1(i)
      lon_f_b_l1(j) = lon_f_l1(i)
      
      e1_t_b_l1(j) = e1_t_l1(i)
      e2_t_b_l1(j) = e2_t_l1(i)
      e1_u_b_l1(j) = e1_u_l1(i)
      e2_u_b_l1(j) = e2_u_l1(i)
      e1_v_b_l1(j) = e1_v_l1(i)
      e2_v_b_l1(j) = e2_v_l1(i)
      e1_f_b_l1(j) = e1_f_l1(i)
      e2_f_b_l1(j) = e2_f_l1(i)
      
      dlat_t_b_l1(j) = dlat_t_l1(i)
      dlon_t_b_l1(j) = dlon_t_l1(i)
      bathy_t_b_l1(j) = bathy_t_l1(i)
      do k=1,2
        smc_tgrid_b_l1_n(j,k) = smc_tgrid_l1_n_temp(i,k)
        smc_tgrid_b_l1_s(j,k) = smc_tgrid_l1_s_temp(i,k)
        smc_tgrid_b_l1_w(j,k) = smc_tgrid_l1_w_temp(i,k)
        smc_tgrid_b_l1_e(j,k) = smc_tgrid_l1_e_temp(i,k)
      end do
      smc_tgrid_b_l1_mask(j) = dble(smc_tgrid_l1_mask(i))
      smc_ugrid_b_l1_mask(j) = dble(smc_ugrid_l1_mask(i))
      smc_vgrid_b_l1_mask(j) = dble(smc_vgrid_l1_mask(i))
      smc_fgrid_b_l1_mask(j) = dble(smc_fgrid_l1_mask(i))   
      j = j + 1
    end if
  end do
  
  do i=1,ng_l1_bathy
  	ff_f_b_l1(i) = 2.0d0*omega*sin(lat_t_b_l1(i)*rad)
  end do
  
!=====creat graph format file for parallel computing
  f = 11
  nw_l1 = 0
  if (smc_fgrid_l1_mask(1) .eq. 0) then
    k_start = 1
  else
    do k=2,fdiv+1
      if (smc_fgrid_l1_mask(k) > 0 ) then
        nw_l1 = nw_l1 + 1
      end if
    end do
    k_start = 2
  end if
  
  do k=k_start,ng_l1_bathy
    i=1
    if ((k<smc_tgrid_b_l1_n(k,i)) .and. (smc_tgrid_b_l1_n(k,i)>0)) then
      nw_l1 = nw_l1 + 1
    end if
    if ((k<smc_tgrid_b_l1_s(k,i)) .and. (smc_tgrid_b_l1_s(k,i)>0)) then
      nw_l1 = nw_l1 + 1
    end if
    if ((k<smc_tgrid_b_l1_e(k,i)) .and. (smc_tgrid_b_l1_e(k,i)>0)) then
      nw_l1 = nw_l1 + 1
    end if
    if ((k<smc_tgrid_b_l1_w(k,i)) .and. (smc_tgrid_b_l1_w(k,i)>0)) then
      nw_l1 = nw_l1 + 1
    end if
    i=2
    if ((k<smc_tgrid_b_l1_n(k,i)) .and. (smc_tgrid_b_l1_n(k,1) .ne. smc_tgrid_b_l1_n(k,2)) .and. (smc_tgrid_b_l1_n(k,i)>0)) then
      nw_l1 = nw_l1 + 1
    end if
    if ((k<smc_tgrid_b_l1_s(k,i)) .and. (smc_tgrid_b_l1_s(k,1) .ne. smc_tgrid_b_l1_s(k,2)) .and. (smc_tgrid_b_l1_s(k,i)>0)) then
      nw_l1 = nw_l1 + 1
    end if
    if ((k<smc_tgrid_b_l1_e(k,i)) .and. (smc_tgrid_b_l1_e(k,1) .ne. smc_tgrid_b_l1_e(k,2)) .and. (smc_tgrid_b_l1_e(k,i)>0)) then
      nw_l1 = nw_l1 + 1
    end if
    if ((k<smc_tgrid_b_l1_w(k,i)) .and. (smc_tgrid_b_l1_w(k,1) .ne. smc_tgrid_b_l1_w(k,2)) .and. (smc_tgrid_b_l1_w(k,i)>0)) then
      nw_l1 = nw_l1 + 1
    end if
  end do
  write(*,*) "total edge is",nw_l1
  open(16,file="smc_grid.graph")
  write(temp_str1,'(i9)') ng_l1_bathy
  write(temp_str2,'(i9)') nw_l1
  temp_str1 = trim(adjustl(temp_str1))//" "//trim(adjustl(temp_str2))
  write(temp_str2,'(i9)') f
  temp_str1 = trim(adjustl(temp_str1))//" "//trim(adjustl(temp_str2))
  write(16,*) trim(adjustl(temp_str1))
  
  if (smc_fgrid_l1_mask(1) .eq. 0) then
    k_start = 1
  else
  	temp_str1 = "1"
    do k=2,fdiv+1
      if (smc_fgrid_l1_mask(k) > 0 ) then
        write(temp_str2,'(i9)') idx_bathy(k)
        temp_str1 = trim(adjustl(temp_str1))//" "//trim(adjustl(temp_str2))//" 1"
      end if
    end do
    write(16,*) trim(adjustl(temp_str1))
    k_start = 2
  end if
  do k=k_start,ng_l1_bathy-1
    temp_str1 = "1"
    if (smc_tgrid_b_l1_n(k,1) .eq. smc_tgrid_b_l1_n(k,2)) then
    	if (smc_tgrid_b_l1_n(k,1) > 0) then
        write(temp_str2,'(i9)') smc_tgrid_b_l1_n(k,1)
        temp_str1 = trim(adjustl(temp_str1))//" "//trim(adjustl(temp_str2))//" 1"
      end if
    else
    	if (smc_tgrid_b_l1_n(k,1) > 0) then
        write(temp_str2,'(i9)') smc_tgrid_b_l1_n(k,1)
        temp_str1 = trim(adjustl(temp_str1))//" "//trim(adjustl(temp_str2))//" 1"
      end if
      if (smc_tgrid_b_l1_n(k,2) > 0) then
        write(temp_str2,'(i9)') smc_tgrid_b_l1_n(k,2)
        temp_str1 = trim(adjustl(temp_str1))//" "//trim(adjustl(temp_str2))//" 1"
      end if
    end if
    
    if (smc_tgrid_b_l1_w(k,1) .eq. smc_tgrid_b_l1_w(k,2)) then
    	if (smc_tgrid_b_l1_w(k,1) > 0) then
        write(temp_str2,'(i9)') smc_tgrid_b_l1_w(k,1)
        temp_str1 = trim(adjustl(temp_str1))//" "//trim(adjustl(temp_str2))//" 1"
      end if
    else
    	if (smc_tgrid_b_l1_w(k,1) > 0) then
        write(temp_str2,'(i9)') smc_tgrid_b_l1_w(k,1)
        temp_str1 = trim(adjustl(temp_str1))//" "//trim(adjustl(temp_str2))//" 1"
      end if
      if (smc_tgrid_b_l1_w(k,2) > 0) then
        write(temp_str2,'(i9)') smc_tgrid_b_l1_w(k,2)
        temp_str1 = trim(adjustl(temp_str1))//" "//trim(adjustl(temp_str2))//" 1"
      end if
    end if
    
    if (smc_tgrid_b_l1_s(k,1) .eq. smc_tgrid_b_l1_s(k,2)) then
    	if (smc_tgrid_b_l1_s(k,1) > 0) then
        write(temp_str2,'(i9)') smc_tgrid_b_l1_s(k,1)
        temp_str1 = trim(adjustl(temp_str1))//" "//trim(adjustl(temp_str2))//" 1"
      end if
    else
    	if (smc_tgrid_b_l1_s(k,1) > 0) then
        write(temp_str2,'(i9)') smc_tgrid_b_l1_s(k,1)
        temp_str1 = trim(adjustl(temp_str1))//" "//trim(adjustl(temp_str2))//" 1"
      end if
      if (smc_tgrid_b_l1_s(k,2) > 0) then
        write(temp_str2,'(i9)') smc_tgrid_b_l1_s(k,2)
        temp_str1 = trim(adjustl(temp_str1))//" "//trim(adjustl(temp_str2))//" 1"
      end if
    end if
    
    if (smc_tgrid_b_l1_e(k,1) .eq. smc_tgrid_b_l1_e(k,2)) then
    	if (smc_tgrid_b_l1_e(k,1) > 0) then
        write(temp_str2,'(i9)') smc_tgrid_b_l1_e(k,1)
        temp_str1 = trim(adjustl(temp_str1))//" "//trim(adjustl(temp_str2))//" 1"
      end if
    else
    	if (smc_tgrid_b_l1_e(k,1) > 0) then
        write(temp_str2,'(i9)') smc_tgrid_b_l1_e(k,1)
        temp_str1 = trim(adjustl(temp_str1))//" "//trim(adjustl(temp_str2))//" 1"
      end if
      if (smc_tgrid_b_l1_e(k,2) > 0) then
        write(temp_str2,'(i9)') smc_tgrid_b_l1_e(k,2)
        temp_str1 = trim(adjustl(temp_str1))//" "//trim(adjustl(temp_str2))//" 1"
      end if
    end if
    write(16,*) trim(adjustl(temp_str1))
  end do
  temp_str1 = "1"
  do k=ng_l1-1,ng_l1-fdiv,-1
  	if (smc_fgrid_l1_mask(k) > 0 ) then
  		write(temp_str2,'(i9)') idx_bathy(k)
  		temp_str1 = trim(adjustl(temp_str1))//" "//trim(adjustl(temp_str2))//" 1"
  	end if
  end do
  write(16,*) trim(adjustl(temp_str1))
  close(16)
  
  allocate(par_id(ng_l1_bathy))
  inquire(file="tmppartition", exist=alive)
  if (alive) then
    open(17,file="tmppartition",status='old')
    do k=1,ng_l1_bathy
      read(17,*) par_id(k)
    end do
    close(17)
  else
    par_id = 0
  end if

  call netcdf_write_smc('smc_zhangyu.nc',ng_l1_bathy,ng_l1,lat_t_l1,lon_t_l1, &
    & lat_t_b_l1,lon_t_b_l1,lat_u_b_l1,lon_u_b_l1,lat_v_b_l1,lon_v_b_l1,lat_f_b_l1,lon_f_b_l1, &
    & ff_f_b_l1, &
    & smc_tgrid_l1_n,smc_tgrid_l1_s,smc_tgrid_l1_e,smc_tgrid_l1_w, &
    & smc_tgrid_b_l1_n,smc_tgrid_b_l1_s,smc_tgrid_b_l1_e,smc_tgrid_b_l1_w,par_id, &
    & idx_bathy, &
    & smc_tgrid_l1_mask,smc_ugrid_l1_mask,smc_vgrid_l1_mask,smc_fgrid_l1_mask, &
    & smc_tgrid_b_l1_mask,smc_ugrid_b_l1_mask,smc_vgrid_b_l1_mask,smc_fgrid_b_l1_mask, &
    & bathy_t_l1,bathy_t_b_l1, &
    & e1_t_b_l1,e1_u_b_l1,e1_v_b_l1,e1_f_b_l1,e2_t_b_l1,e2_u_b_l1,e2_v_b_l1,e2_f_b_l1, &
    & lat_polex1, lat_polex2, lat_polex3, lat_polex4)
  
  deallocate(ratio,dlon)
  deallocate(dlat_t_l1,dlon_t_l1)
  deallocate(smc_tgrid_back_l1)
  deallocate(smc_tgrid_back_l1_i,smc_tgrid_back_l1_j)
  deallocate(smc_tgrid_l1_n,smc_tgrid_l1_s,smc_tgrid_l1_e,smc_tgrid_l1_w)
  deallocate(smc_tgrid_l1_n_temp,smc_tgrid_l1_s_temp,smc_tgrid_l1_e_temp,smc_tgrid_l1_w_temp)
  deallocate(smc_tgrid_back_l2)
  deallocate(smc_tgrid_back_l3)
  deallocate(par_id)
  deallocate(bathy_t_l1_2d,bathy_t_l1)
  deallocate(smc_tgrid_l1_mask,smc_ugrid_l1_mask,smc_vgrid_l1_mask,smc_fgrid_l1_mask)
  deallocate(smc_tgrid_b_l1_mask,smc_ugrid_b_l1_mask,smc_vgrid_b_l1_mask,smc_fgrid_b_l1_mask)
  deallocate(smc_ugrid_l1_mask_ls,smc_vgrid_l1_mask_ls)
  deallocate(idx_bathy)
  deallocate(smc_tgrid_b_l1_n,smc_tgrid_b_l1_s,smc_tgrid_b_l1_e,smc_tgrid_b_l1_w)
  deallocate(dlat_t_b_l1,dlon_t_b_l1)
  deallocate(ff_f_b_l1)
  deallocate(lat_t_b_l1,lat_u_b_l1,lat_v_b_l1,lat_f_b_l1)
  deallocate(lon_t_b_l1,lon_u_b_l1,lon_v_b_l1,lon_f_b_l1)
  deallocate(lat_t_l1,lat_u_l1,lat_v_l1,lat_f_l1)
  deallocate(lon_t_l1,lon_u_l1,lon_v_l1,lon_f_l1)
  deallocate(e1_t_b_l1,e1_u_b_l1,e1_v_b_l1,e1_f_b_l1)
  deallocate(e2_t_b_l1,e2_u_b_l1,e2_v_b_l1,e2_f_b_l1)
  deallocate(e1_t_l1,e1_u_l1,e1_v_l1,e1_f_l1)
  deallocate(e2_t_l1,e2_u_l1,e2_v_l1,e2_f_l1)
end program smc_grid_create

!--------------------------------------------------------
subroutine cal_var(varratio,ratio,ny)
!--------------------------------------------------------
  implicit none
  integer(kind=4) :: i,ny
  real(kind=8) :: ratio(ny),diff(ny),varratio
  
  do i=1,ny
    diff(i) = (ratio(i)-1)*(ratio(i)-1)
  end do
  varratio = sum(diff)/(max(1,size(diff)))
  
  return
end subroutine

!--------------------------------------------------------
subroutine check(status)
!--------------------------------------------------------
  use netcdf
  integer(kind=4), intent ( in) :: status

  if(status /= nf90_noerr) then 
     print *, trim(nf90_strerror(status))
     stop "Stopped"
  end if
end subroutine check

!--------------------------------------------------------
subroutine netcdf_write_smc(ncname,nglo_b,nglo,lat,lon, &
    & lat_t_b,lon_t_b,lat_u_b,lon_u_b,lat_v_b,lon_v_b,lat_f_b,lon_f_b, &
    & ff_f_b, &
    & smc_tgrid_n,smc_tgrid_s,smc_tgrid_e,smc_tgrid_w, &
    & smc_tgrid_b_n,smc_tgrid_b_s,smc_tgrid_b_e,smc_tgrid_b_w,smc_tgrid_par, &
    & idx_b, &
    & smc_tgrid_mask,smc_ugrid_mask,smc_vgrid_mask,smc_fgrid_mask, &
    & smc_tgrid_b_mask,smc_ugrid_b_mask,smc_vgrid_b_mask,smc_fgrid_b_mask, &
    & bathy,bathy_b, &
    & e1_t_b,e1_u_b,e1_v_b,e1_f_b,e2_t_b,e2_u_b,e2_v_b,e2_f_b, &
    & lat_polex1, lat_polex2, lat_polex3, lat_polex4)
!--------------------------------------------------------
use netcdf
implicit none
character (len = *) :: ncname
character (len = *), parameter :: nglo_name = "nglo"
character (len = *), parameter :: nglo_b_name = "nglo_b"
character (len = *), parameter :: neig_name = "neig"
character (len = 10):: cdate, ctime
integer(kind=4) :: ncid
integer(kind=4) :: nglo, nglo_b
integer(kind=4) :: nglo_dimid,nglo_b_dimid,neig_dimid,dimid_2d(2),dimid_b_2d(2)
integer(kind=4) :: lat_id,lon_id,lat_2d_id,lon_2d_id
integer(kind=4) :: lat_t_b_id,lat_u_b_id,lat_v_b_id,lat_f_b_id,ff_f_b_id
integer(kind=4) :: lon_t_b_id,lon_u_b_id,lon_v_b_id,lon_f_b_id
integer(kind=4) :: e1_t_b_id,e1_u_b_id,e1_v_b_id,e1_f_b_id
integer(kind=4) :: e2_t_b_id,e2_u_b_id,e2_v_b_id,e2_f_b_id
integer(kind=4) :: tgrid_n_id,tgrid_s_id,tgrid_e_id,tgrid_w_id,tgrid_par_id
integer(kind=4) :: tgrid_b_n_id,tgrid_b_s_id,tgrid_b_e_id,tgrid_b_w_id
integer(kind=4) :: idx_b_id,bathy_id,bathy_b_id
integer(kind=4) :: smc_tgrid_mask_id,smc_ugrid_mask_id,smc_vgrid_mask_id,smc_fgrid_mask_id
integer(kind=4) :: smc_tgrid_b_mask_id,smc_ugrid_b_mask_id,smc_vgrid_b_mask_id,smc_fgrid_b_mask_id
integer(kind=4) :: smc_tgrid_n(nglo,2),smc_tgrid_s(nglo,2),smc_tgrid_e(nglo,2),smc_tgrid_w(nglo,2)
integer(kind=4) :: smc_tgrid_b_n(nglo_b,2),smc_tgrid_b_s(nglo_b,2),smc_tgrid_b_e(nglo_b,2),smc_tgrid_b_w(nglo_b,2)
integer(kind=4) :: smc_tgrid_par(nglo_b),idx_b(nglo)
integer(kind=4) :: smc_tgrid_mask(nglo),smc_ugrid_mask(nglo),smc_vgrid_mask(nglo),smc_fgrid_mask(nglo)
real(kind=8) :: smc_tgrid_b_mask(nglo_b),smc_ugrid_b_mask(nglo_b),smc_vgrid_b_mask(nglo_b),smc_fgrid_b_mask(nglo_b)
real(kind=8) :: lat(nglo),lon(nglo),lat_2d(nglo,2),lon_2d(nglo,2),bathy(nglo),bathy_b(nglo_b)
real(kind=8) :: lat_t_b(nglo_b),lat_u_b(nglo_b),lat_v_b(nglo_b),lat_f_b(nglo_b),ff_f_b(nglo_b)
real(kind=8) :: lon_t_b(nglo_b),lon_u_b(nglo_b),lon_v_b(nglo_b),lon_f_b(nglo_b)
real(kind=8) :: e1_t_b(nglo_b),e1_u_b(nglo_b),e1_v_b(nglo_b),e1_f_b(nglo_b)
real(kind=8) :: e2_t_b(nglo_b),e2_u_b(nglo_b),e2_v_b(nglo_b),e2_f_b(nglo_b)
real(kind=8) :: lat_polex1, lat_polex2, lat_polex3, lat_polex4

lat_2d(:,1) = lat(:)
lat_2d(:,2) = lat(:)
lon_2d(:,1) = lon(:)
lon_2d(:,2) = lon(:)

call date_and_time(cdate, ctime)
call check( nf90_create(ncname, nf90_64bit_offset, ncid))

call check( nf90_def_dim(ncid, nglo_name, nglo, nglo_dimid))
call check( nf90_def_dim(ncid, nglo_b_name, nglo_b, nglo_b_dimid))
call check( nf90_def_dim(ncid, neig_name,    2, neig_dimid))
dimid_2d = (/nglo_dimid, neig_dimid/)
dimid_b_2d = (/nglo_b_dimid, neig_dimid/)

call check( nf90_def_var(ncid, "lon", nf90_double, nglo_dimid, lon_id))
call check( nf90_put_att(ncid, lon_id, "max", maxval(lon)))
call check( nf90_put_att(ncid, lon_id, "min", minval(lon)))

call check( nf90_def_var(ncid, "lat", nf90_double, nglo_dimid, lat_id))
call check( nf90_put_att(ncid, lat_id, "max", maxval(lat)))
call check( nf90_put_att(ncid, lat_id, "min", minval(lat)))

call check( nf90_def_var(ncid, "lon_t_b", nf90_double, nglo_b_dimid, lon_t_b_id))
call check( nf90_put_att(ncid, lon_t_b_id, "max", maxval(lon_t_b)))
call check( nf90_put_att(ncid, lon_t_b_id, "min", minval(lon_t_b)))

call check( nf90_def_var(ncid, "lat_t_b", nf90_double, nglo_b_dimid, lat_t_b_id))
call check( nf90_put_att(ncid, lat_t_b_id, "max", maxval(lat_t_b)))
call check( nf90_put_att(ncid, lat_t_b_id, "min", minval(lat_t_b)))

call check( nf90_def_var(ncid, "lon_u_b", nf90_double, nglo_b_dimid, lon_u_b_id))
call check( nf90_put_att(ncid, lon_u_b_id, "max", maxval(lon_u_b)))
call check( nf90_put_att(ncid, lon_u_b_id, "min", minval(lon_u_b)))

call check( nf90_def_var(ncid, "lat_u_b", nf90_double, nglo_b_dimid, lat_u_b_id))
call check( nf90_put_att(ncid, lat_u_b_id, "max", maxval(lat_u_b)))
call check( nf90_put_att(ncid, lat_u_b_id, "min", minval(lat_u_b)))

call check( nf90_def_var(ncid, "lon_v_b", nf90_double, nglo_b_dimid, lon_v_b_id))
call check( nf90_put_att(ncid, lon_v_b_id, "max", maxval(lon_v_b)))
call check( nf90_put_att(ncid, lon_v_b_id, "min", minval(lon_v_b)))

call check( nf90_def_var(ncid, "lat_v_b", nf90_double, nglo_b_dimid, lat_v_b_id))
call check( nf90_put_att(ncid, lat_v_b_id, "max", maxval(lat_v_b)))
call check( nf90_put_att(ncid, lat_v_b_id, "min", minval(lat_v_b)))

call check( nf90_def_var(ncid, "lon_f_b", nf90_double, nglo_b_dimid, lon_f_b_id))
call check( nf90_put_att(ncid, lon_f_b_id, "max", maxval(lon_f_b)))
call check( nf90_put_att(ncid, lon_f_b_id, "min", minval(lon_f_b)))

call check( nf90_def_var(ncid, "lat_f_b", nf90_double, nglo_b_dimid, lat_f_b_id))
call check( nf90_put_att(ncid, lat_f_b_id, "max", maxval(lat_f_b)))
call check( nf90_put_att(ncid, lat_f_b_id, "min", minval(lat_f_b)))

call check( nf90_def_var(ncid, "ff_f_b", nf90_double, nglo_b_dimid, ff_f_b_id))
call check( nf90_put_att(ncid, ff_f_b_id, "max", maxval(ff_f_b)))
call check( nf90_put_att(ncid, ff_f_b_id, "min", minval(ff_f_b)))

call check( nf90_def_var(ncid, "e1_t_b", nf90_double, nglo_b_dimid, e1_t_b_id))
call check( nf90_put_att(ncid, e1_t_b_id, "max", maxval(e1_t_b)))
call check( nf90_put_att(ncid, e1_t_b_id, "min", minval(e1_t_b)))

call check( nf90_def_var(ncid, "e2_t_b", nf90_double, nglo_b_dimid, e2_t_b_id))
call check( nf90_put_att(ncid, e2_t_b_id, "max", maxval(e2_t_b)))
call check( nf90_put_att(ncid, e2_t_b_id, "min", minval(e2_t_b)))

call check( nf90_def_var(ncid, "e1_u_b", nf90_double, nglo_b_dimid, e1_u_b_id))
call check( nf90_put_att(ncid, e1_u_b_id, "max", maxval(e1_u_b)))
call check( nf90_put_att(ncid, e1_u_b_id, "min", minval(e1_u_b)))

call check( nf90_def_var(ncid, "e2_u_b", nf90_double, nglo_b_dimid, e2_u_b_id))
call check( nf90_put_att(ncid, e2_u_b_id, "max", maxval(e2_u_b)))
call check( nf90_put_att(ncid, e2_u_b_id, "min", minval(e2_u_b)))

call check( nf90_def_var(ncid, "e1_v_b", nf90_double, nglo_b_dimid, e1_v_b_id))
call check( nf90_put_att(ncid, e1_v_b_id, "max", maxval(e1_v_b)))
call check( nf90_put_att(ncid, e1_v_b_id, "min", minval(e1_v_b)))

call check( nf90_def_var(ncid, "e2_v_b", nf90_double, nglo_b_dimid, e2_v_b_id))
call check( nf90_put_att(ncid, e2_v_b_id, "max", maxval(e2_v_b)))
call check( nf90_put_att(ncid, e2_v_b_id, "min", minval(e2_v_b)))

call check( nf90_def_var(ncid, "e1_f_b", nf90_double, nglo_b_dimid, e1_f_b_id))
call check( nf90_put_att(ncid, e1_f_b_id, "max", maxval(e1_f_b)))
call check( nf90_put_att(ncid, e1_f_b_id, "min", minval(e1_f_b)))

call check( nf90_def_var(ncid, "e2_f_b", nf90_double, nglo_b_dimid, e2_f_b_id))
call check( nf90_put_att(ncid, e2_f_b_id, "max", maxval(e2_f_b)))
call check( nf90_put_att(ncid, e2_f_b_id, "min", minval(e2_f_b)))

call check( nf90_def_var(ncid, "lon_2d", nf90_double, dimid_2d, lon_2d_id))
call check( nf90_put_att(ncid, lon_2d_id, "max", maxval(lon_2d)))
call check( nf90_put_att(ncid, lon_2d_id, "min", minval(lon_2d)))

call check( nf90_def_var(ncid, "lat_2d", nf90_double, dimid_2d, lat_2d_id))
call check( nf90_put_att(ncid, lat_2d_id, "max", maxval(lat_2d)))
call check( nf90_put_att(ncid, lat_2d_id, "min", minval(lat_2d)))

call check( nf90_def_var(ncid, "tgrid_n", nf90_int, dimid_2d, tgrid_n_id))
call check( nf90_put_att(ncid, tgrid_n_id, "max", maxval(smc_tgrid_n)))
call check( nf90_put_att(ncid, tgrid_n_id, "min", minval(smc_tgrid_n)))

call check( nf90_def_var(ncid, "tgrid_s", nf90_int, dimid_2d, tgrid_s_id))
call check( nf90_put_att(ncid, tgrid_s_id, "max", maxval(smc_tgrid_s)))
call check( nf90_put_att(ncid, tgrid_s_id, "min", minval(smc_tgrid_s)))

call check( nf90_def_var(ncid, "tgrid_e", nf90_int, dimid_2d, tgrid_e_id))
call check( nf90_put_att(ncid, tgrid_e_id, "max", maxval(smc_tgrid_e)))
call check( nf90_put_att(ncid, tgrid_e_id, "min", minval(smc_tgrid_e)))

call check( nf90_def_var(ncid, "tgrid_w", nf90_int, dimid_2d, tgrid_w_id))
call check( nf90_put_att(ncid, tgrid_w_id, "max", maxval(smc_tgrid_w)))
call check( nf90_put_att(ncid, tgrid_w_id, "min", minval(smc_tgrid_w)))

call check( nf90_def_var(ncid, "tgrid_b_n", nf90_int, dimid_b_2d, tgrid_b_n_id))
call check( nf90_put_att(ncid, tgrid_b_n_id, "max", maxval(smc_tgrid_b_n)))
call check( nf90_put_att(ncid, tgrid_b_n_id, "min", minval(smc_tgrid_b_n)))

call check( nf90_def_var(ncid, "tgrid_b_s", nf90_int, dimid_b_2d, tgrid_b_s_id))
call check( nf90_put_att(ncid, tgrid_b_s_id, "max", maxval(smc_tgrid_b_s)))
call check( nf90_put_att(ncid, tgrid_b_s_id, "min", minval(smc_tgrid_b_s)))

call check( nf90_def_var(ncid, "tgrid_b_e", nf90_int, dimid_b_2d, tgrid_b_e_id))
call check( nf90_put_att(ncid, tgrid_b_e_id, "max", maxval(smc_tgrid_b_e)))
call check( nf90_put_att(ncid, tgrid_b_e_id, "min", minval(smc_tgrid_b_e)))

call check( nf90_def_var(ncid, "tgrid_b_w", nf90_int, dimid_b_2d, tgrid_b_w_id))
call check( nf90_put_att(ncid, tgrid_b_w_id, "max", maxval(smc_tgrid_b_w)))
call check( nf90_put_att(ncid, tgrid_b_w_id, "min", minval(smc_tgrid_b_w)))

call check( nf90_def_var(ncid, "smc_tgrid_par", nf90_int, nglo_b_dimid, tgrid_par_id))
call check( nf90_put_att(ncid, tgrid_par_id, "max", maxval(smc_tgrid_par)))
call check( nf90_put_att(ncid, tgrid_par_id, "min", minval(smc_tgrid_par)))

call check( nf90_def_var(ncid, "idx_bathy", nf90_int, nglo_dimid, idx_b_id))
call check( nf90_put_att(ncid, idx_b_id, "max", maxval(idx_b)))
call check( nf90_put_att(ncid, idx_b_id, "min", minval(idx_b)))

call check( nf90_def_var(ncid, "smc_tgrid_mask", nf90_int, nglo_dimid, smc_tgrid_mask_id))
call check( nf90_put_att(ncid, smc_tgrid_mask_id, "max", maxval(smc_tgrid_mask)))
call check( nf90_put_att(ncid, smc_tgrid_mask_id, "min", minval(smc_tgrid_mask)))

call check( nf90_def_var(ncid, "smc_ugrid_mask", nf90_int, nglo_dimid, smc_ugrid_mask_id))
call check( nf90_put_att(ncid, smc_ugrid_mask_id, "max", maxval(smc_ugrid_mask)))
call check( nf90_put_att(ncid, smc_ugrid_mask_id, "min", minval(smc_ugrid_mask)))

call check( nf90_def_var(ncid, "smc_vgrid_mask", nf90_int, nglo_dimid, smc_vgrid_mask_id))
call check( nf90_put_att(ncid, smc_vgrid_mask_id, "max", maxval(smc_vgrid_mask)))
call check( nf90_put_att(ncid, smc_vgrid_mask_id, "min", minval(smc_vgrid_mask)))

call check( nf90_def_var(ncid, "smc_fgrid_mask", nf90_int, nglo_dimid, smc_fgrid_mask_id))
call check( nf90_put_att(ncid, smc_fgrid_mask_id, "max", maxval(smc_fgrid_mask)))
call check( nf90_put_att(ncid, smc_fgrid_mask_id, "min", minval(smc_fgrid_mask)))

call check( nf90_def_var(ncid, "smc_tgrid_b_mask", nf90_double, nglo_b_dimid, smc_tgrid_b_mask_id))
call check( nf90_put_att(ncid, smc_tgrid_b_mask_id, "max", maxval(smc_tgrid_b_mask)))
call check( nf90_put_att(ncid, smc_tgrid_b_mask_id, "min", minval(smc_tgrid_b_mask)))

call check( nf90_def_var(ncid, "smc_ugrid_b_mask", nf90_double, nglo_b_dimid, smc_ugrid_b_mask_id))
call check( nf90_put_att(ncid, smc_ugrid_b_mask_id, "max", maxval(smc_ugrid_b_mask)))
call check( nf90_put_att(ncid, smc_ugrid_b_mask_id, "min", minval(smc_ugrid_b_mask)))

call check( nf90_def_var(ncid, "smc_vgrid_b_mask", nf90_double, nglo_b_dimid, smc_vgrid_b_mask_id))
call check( nf90_put_att(ncid, smc_vgrid_b_mask_id, "max", maxval(smc_vgrid_b_mask)))
call check( nf90_put_att(ncid, smc_vgrid_b_mask_id, "min", minval(smc_vgrid_b_mask)))

call check( nf90_def_var(ncid, "smc_fgrid_b_mask", nf90_double, nglo_b_dimid, smc_fgrid_b_mask_id))
call check( nf90_put_att(ncid, smc_fgrid_b_mask_id, "max", maxval(smc_fgrid_b_mask)))
call check( nf90_put_att(ncid, smc_fgrid_b_mask_id, "min", minval(smc_fgrid_b_mask)))

call check( nf90_def_var(ncid, "smc_tgrid_bathy", nf90_double, nglo_dimid, bathy_id))
call check( nf90_put_att(ncid, bathy_id, "max", maxval(bathy)))
call check( nf90_put_att(ncid, bathy_id, "min", minval(bathy)))

call check( nf90_def_var(ncid, "smc_tgrid_bathy_b", nf90_double, nglo_b_dimid, bathy_b_id))
call check( nf90_put_att(ncid, bathy_b_id, "max", maxval(bathy_b)))
call check( nf90_put_att(ncid, bathy_b_id, "min", minval(bathy_b)))

call check( nf90_put_att(ncid, nf90_global, "Author", 'Yu Zhang'))
call check( nf90_put_att(ncid, nf90_global, "Organization", 'National Marine Environmental Forecasting Center'))
call check( nf90_put_att(ncid, nf90_global, "date", trim(cdate)//', '//trim(ctime(1:2))//':'//trim(ctime(3:4)) &
            & //':'//trim(ctime(5:6))))
call check( nf90_put_att(ncid, nf90_global, "about", 'SMC grid file'))
call check( nf90_put_att(ncid, nf90_global, "lat_polex1", lat_polex1))
call check( nf90_put_att(ncid, nf90_global, "lat_polex2", lat_polex2))
call check( nf90_put_att(ncid, nf90_global, "lat_polex3", lat_polex3))
call check( nf90_put_att(ncid, nf90_global, "lat_polex4", lat_polex4))
call check( nf90_enddef(ncid) )

call check( nf90_put_var(ncid, lon_id, lon ) )
call check( nf90_put_var(ncid, lat_id, lat ) )
call check( nf90_put_var(ncid, lon_t_b_id, lon_t_b ) )
call check( nf90_put_var(ncid, lat_t_b_id, lat_t_b ) )
call check( nf90_put_var(ncid, lon_u_b_id, lon_u_b ) )
call check( nf90_put_var(ncid, lat_u_b_id, lat_u_b ) )
call check( nf90_put_var(ncid, lon_v_b_id, lon_v_b ) )
call check( nf90_put_var(ncid, lat_v_b_id, lat_v_b ) )
call check( nf90_put_var(ncid, lon_f_b_id, lon_f_b ) )
call check( nf90_put_var(ncid, lat_f_b_id, lat_f_b ) )
call check( nf90_put_var(ncid, ff_f_b_id, ff_f_b ) )
call check( nf90_put_var(ncid, e1_t_b_id, e1_t_b ) )
call check( nf90_put_var(ncid, e2_t_b_id, e2_t_b ) )
call check( nf90_put_var(ncid, e1_u_b_id, e1_u_b ) )
call check( nf90_put_var(ncid, e2_u_b_id, e2_u_b ) )
call check( nf90_put_var(ncid, e1_v_b_id, e1_v_b ) )
call check( nf90_put_var(ncid, e2_v_b_id, e2_v_b ) )
call check( nf90_put_var(ncid, e1_f_b_id, e1_f_b ) )
call check( nf90_put_var(ncid, e2_f_b_id, e2_f_b ) )
call check( nf90_put_var(ncid, lon_2d_id, lon_2d ) )
call check( nf90_put_var(ncid, lat_2d_id, lat_2d ) )
call check( nf90_put_var(ncid, tgrid_n_id, smc_tgrid_n ) )
call check( nf90_put_var(ncid, tgrid_s_id, smc_tgrid_s ) )
call check( nf90_put_var(ncid, tgrid_e_id, smc_tgrid_e ) )
call check( nf90_put_var(ncid, tgrid_w_id, smc_tgrid_w ) )
call check( nf90_put_var(ncid, tgrid_b_n_id, smc_tgrid_b_n ) )
call check( nf90_put_var(ncid, tgrid_b_s_id, smc_tgrid_b_s ) )
call check( nf90_put_var(ncid, tgrid_b_e_id, smc_tgrid_b_e ) )
call check( nf90_put_var(ncid, tgrid_b_w_id, smc_tgrid_b_w ) )
call check( nf90_put_var(ncid, tgrid_par_id, smc_tgrid_par ) )
call check( nf90_put_var(ncid, idx_b_id, idx_b ) )
call check( nf90_put_var(ncid, smc_tgrid_mask_id, smc_tgrid_mask ) )
call check( nf90_put_var(ncid, smc_ugrid_mask_id, smc_ugrid_mask ) )
call check( nf90_put_var(ncid, smc_vgrid_mask_id, smc_vgrid_mask ) )
call check( nf90_put_var(ncid, smc_fgrid_mask_id, smc_fgrid_mask ) )
call check( nf90_put_var(ncid, smc_tgrid_b_mask_id, smc_tgrid_b_mask ) )
call check( nf90_put_var(ncid, smc_ugrid_b_mask_id, smc_ugrid_b_mask ) )
call check( nf90_put_var(ncid, smc_vgrid_b_mask_id, smc_vgrid_b_mask ) )
call check( nf90_put_var(ncid, smc_fgrid_b_mask_id, smc_fgrid_b_mask ) )
call check( nf90_put_var(ncid, bathy_id, bathy ) )
call check( nf90_put_var(ncid, bathy_b_id, bathy_b ) )
call check( nf90_close(ncid) )

print *,"*** SUCCESS write smc grid file ",ncname," !"

return
end subroutine netcdf_write_smc