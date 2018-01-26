c
c ---------------------------------------------------------------------
      module cominput 
c     last edited by Bill Sacks 01.13.10
c ---------------------------------------------------------------------
c
c This file defines parameters that describe a climate input data set.
c For more details on the parameters contained here, see the file:
c input_descriptions/README.
c
      implicit none
      save

c The following are parameters describing the characteristics of the input data set:
c
      logical
     >     read_daily_directly, ! true -> read daily climate data directly, rather than applying daily anomalies to monthly values
     >     one_file_per_year,   ! true -> separate input file for each year
     >     read_tmin_tmax,      ! true -> read tmin and tmax directly, rather than tmean
     >     read_radiation,      ! true -> read incoming radiation directly, rather than cloud cover
     >     correct_sphumid,     ! true -> do corrections of specific humidity based on covariances between RH and temperature
     >     conserve_tdew        ! true -> diurnal cycle of humidity conserves daily mean dewpoint temp. rather than daily mean specific humidity

      integer
     >     istyrm,              ! first year monthly values exist in monthly timeseries files
     >     nanom,               ! number of years in monthly timeseries files
     >     istyrd,              ! first year monthly values exist in daily timeseries files
     >     nanomd               ! number of years in daily timeseries files

      character*8 level_var     ! name of level variable (if any) in input files

c
c The following are the names of the files used with this input data set.
c Note that these file names do NOT contain the trailing '.nc';
c for one_file_per_year = TRUE, they also do not contain the _YYYY suffix
c (so, e.g., a file name of zedx/daily/tmax could stand for zedx/daily/tmax_1948.nc, zedx/daily/tmax_1949.nc, etc.)
c     *_mon_clim: monthly climatologies; 
c     *_mon: monthly timeseries, if provided;
c     *_daily: daily timeseries [possibly specified as anomalies]
c
      character*256
     >     fn_cloud_mon_clim,
     >     fn_cloud_mon,
     >     fn_cloud_daily,
     >     fn_rads_mon_clim,
     >     fn_rads_mon,
     >     fn_rads_daily,
     >     fn_temp_mon_clim,
     >     fn_temp_mon,
     >     fn_temp_daily,
     >     fn_dtr_mon_clim,
     >     fn_dtr_mon,
     >     fn_dtr_daily,
     >     fn_tmin_mon_clim,
     >     fn_tmin_mon,
     >     fn_tmin_daily,
     >     fn_tmax_mon_clim,
     >     fn_tmax_mon,
     >     fn_tmax_daily,
     >     fn_prec_mon_clim,
     >     fn_prec_mon,
     >     fn_prec_daily,
     >     fn_rh_mon_clim,
     >     fn_rh_mon,
     >     fn_rh_daily,
     >     fn_wspd_mon_clim,
     >     fn_wspd_mon,
     >     fn_wspd_daily,
     >     fn_wetd_mon_clim,
     >     fn_wetd_mon,            ! note that there is no fn_wetd_daily, since wetd is only used in the weather generator
     >     fn_tminavgann


c
c The following are the names of the netCDF variables in the files in this data set.
c There should be one variable name for each file name given above (the fn_* variables).
c
      character*32
     >     var_cloud_mon_clim,
     >     var_cloud_mon,
     >     var_cloud_daily,
     >     var_rads_mon_clim,
     >     var_rads_mon,
     >     var_rads_daily,
     >     var_temp_mon_clim,
     >     var_temp_mon,
     >     var_temp_daily,
     >     var_dtr_mon_clim,
     >     var_dtr_mon,
     >     var_dtr_daily,
     >     var_tmin_mon_clim,
     >     var_tmin_mon,
     >     var_tmin_daily,
     >     var_tmax_mon_clim,
     >     var_tmax_mon,
     >     var_tmax_daily,
     >     var_prec_mon_clim,
     >     var_prec_mon,
     >     var_prec_daily,
     >     var_rh_mon_clim,
     >     var_rh_mon,
     >     var_rh_daily,
     >     var_wspd_mon_clim,
     >     var_wspd_mon,
     >     var_wspd_daily,
     >     var_wetd_mon_clim,
     >     var_wetd_mon,            ! note that there is no var_wetd_daily, since wetd is only used in the weather generator      
     >     var_tminavgann

c
c The following are variables set in the code (i.e., not read in from the parameter file)
c
      integer 
     >     overlap_start,     ! start year of overlap between monthly & daily files
     >     overlap_end        ! end year of overlap between monthly & daily files

c Available years of climate data, divided into leap years and non-leap years
c (Used for spin-up when weather_generator is FALSE)
c As of 02.16.10, these are allocated and set by clim_file_utils :
c climate_years_avail, called from read_input_description
      integer, dimension(:), allocatable ::
     >        avail_leaps,
     >        avail_nonleaps



      end module cominput
