c Last modified by Bill Sacks May, 2010
c 
c program crops.f
c ----------------------------------------------------------------------------
      subroutine initialcrop
c ----------------------------------------------------------------------------
c
c called from main.f or physiology.f
c
c subroutine to initialize crop related variables that are dependent  
c on a crop growing season basis.  Variables are initialized at the end
c of the growing season after harvest date and after data is written out (io.f)  
c
c additionally, this routine needs to be called at (1) first time
c crops are grown in a model run (a restart), or (2) when
c crops are replacing natural vegetation (which could also be a restart)
c
c uses:
c
      use comgrid
      use compar
      use comsum
      use comveg
      use comcrop
      use comnitr
c
      implicit none
c
      integer iy, i, ic
c
c reset variables after harvest date and data is written 
c to io.f and crop diagnostic output files
c have to make sure there is no memory in the system
c of upper canopy or lower canopy vegetation when crops
c replace natural vegetation!
c
               call const(plai,      npoi*npft, 0.01)
               call const(thrlai,    npoi*npft, 0.0)
               call const(peaklai,   npoi*npft, 0.0)
               call const(lai_decl_per_day, npoi*npft, 0.0)
               call const(grnfraccrop, npoi*npft, 0.0)
               call const(cbiol,     npoi*npft, 0.0)
               call const(cbios,     npoi*npft, 0.0)    
               call const(cbior,     npoi*npft, 0.0)   
               call const(cbiow,     npoi*npft, 0.0) 
               call const(cbiog,     npoi*npft, 0.0)  
               call const(hui,       npoi*npft, 0.0) 
               call const(aleaf,     npoi*npft, 0.0) 
               call const(aroot,     npoi*npft, 0.0)
               call const(astem,     npoi*npft, 0.0)      
               call const(arepr,     npoi*npft, 0.0)  
               call const(awood,     npoi*npft, 0.0)
               call const(aybprod,   npoi*npft, 0.0) 
               call const(ayrprod,   npoi*npft, 0.0)  
               call const(ayabprod,  npoi*npft, 0.0) 
               call const(aylprod,   npoi*npft, 0.0) 
               call const(harvidx,   npoi*npft, 0.0) 
               call const(leafout,   npoi*npft, 0.0) 
               call const(htmx,      npoi*1,    0.0) 
               call const(cumlvs,    npoi*npft, 0.0) 
               call const(plaimx,    npoi*npft, 0.0) 
               call const(dpgf,      npoi*npft, 0.0)
               call const(dpmature,  npoi*npft, 0.0)
               call const(biomass,   npoi*npft, 0.0) 
               call const(totnuptake,npoi*npft, 0.0) 
               call const(tnplant,   npoi*npft, 0.0) 
               call const(totnfix,   npoi*npft, 0.0) 
               call const(idpp,      npoi*npft, 0.0) 
               call const(fixn,      npoi*npft, 0.0) 
               call const(gddplant,  npoi*npft, 0.0) 
               call const(crmclim,   npoi*npft, 0.0) 
               call const(crmact,    npoi*npft, 0.0) 
               call const(crmplant,    npoi*npft, 0.0) 
               call const(emerge_day,  npoi*npft, 9999.0)
               call const(emerge_gdd,  npoi*npft, -1.0)
               call const(emerge_gddtsoi,  npoi*npft, -1.0)
               call const(grainday,  npoi*npft, 9999.0)
               call const(graingdd,  npoi*npft, -1.0)
               call const(gddtsoi,   npoi*npft, 0.0) 
               call const(fertinput, npoi*npft, 0.0) 
c
c grid cell variables - not dependent on pft
c
               call const(sai,     npoi*2, 0.0)
               call const(fu,      npoi,   0.0)
               call const(lai,     npoi*2, 0.0)
               call const(zbot,    npoi*2, 0.0) 
               call const(ztop,    npoi*2, 0.0)  
               call const(totbiou, npoi,   0.0) 
               call const(totbiol, npoi,   0.0) 
               call const(totlaiu, npoi,   0.0) 
               call const(totlail, npoi,   0.0) 
               call const(vf,      npoi,   0.0) 
c
c
      return
      end
c
c ----------------------------------------------------------------------------
      subroutine rotation(irestart, irstyear, iyear)
c ----------------------------------------------------------------------------
c
c subroutine used to plant different crop types in different years
c to simulate typical crop rotations in the United States 
c could be modified in future if rotations would include natural
c vegetation types such as grasslands or biofuel (switchgrass) crops
c
c currently three main types of rotations are used
c if iwheat eq:
c 2: maize/soybean rotation (alternating)
c 3: maize/soybean/spring wheat
c 4: soybean/winter wheat/maize
c
c note: by doing continuous winter wheat, land is fallow 
c only from harvest (june-july) through planting (Sept-Nov).
c
c uses:
c
      use comgrid
      use compar
      use comsum
      use comveg
      use comcrop
      use comnitr
      use comwork, only: cdummy
      use combcs,  only: cmask
c
      implicit none
c
c local variables
c
      real
     >  xfrac
c
      integer
     >  iyear,
     >  irestart,
     >  irstyear,
     >  iyrrot,
     >  i,
     >  idiv,
     >  idx
c
c begin grid
c in this case, irotation also is the number of crops in a specified rotation
c
c assumes that natural vegetation existence arrays are set to 0.0
c in climate.f (existence) 
c
      if (irestart .eq. 1) then
          iyrrot = irstyear
      else 
          iyrrot = 1950  ! base year to start crop rotations
      endif
c
c look at the fraction remainer to determine which crop should
c be planted
c
      !Initialization to avoid mod(*,0) for idiv (Y.Li)
      idiv = 1
      if (irotation .eq. 1) idiv = 3
      if (irotation .eq. 2) idiv = 2
      if (irotation .eq. 3) idiv = 3
      if (irotation .eq. 4) idiv = 3
      if (irotation .eq. 5) idiv = 2
c
      xfrac = mod ((iyear - iyrrot), idiv)  


      do 100 i = 1, npoi
        if(cmask(i) == 0 ) cycle
c
c 2:
c two-crop rotation (standard soybean/corn)
c alternate between even / odd years
c
        if (irotation .eq. 2) then
          if (xfrac .eq. 0) then
            isoybean = 1
            imaize   = 0
            iwheat   = 0
            exist(i,13) = 1.0
            exist(i,14) = 0.0
            exist(i,15) = 0.0
c
          else
            isoybean = 0
            imaize   = 1
            iwheat   = 0
            exist(i,13) = 0.0
            exist(i,14) = 1.0
            exist(i,15) = 0.0
          endif
c 3:
c rotation with 3 crops (corn, corn, soybean) 
c
        elseif (irotation .eq. 3) then
          if (xfrac .eq. 0) then 
            isoybean = 0
            imaize   = 1
            iwheat   = 0
            exist(i,13) = 0.0
            exist(i,14) = 1.0
            exist(i,15) = 0.0
c
          elseif (xfrac .eq. 1.0) then
            isoybean = 1
            imaize   = 0
            iwheat   = 0
            exist(i,13) = 1.0
            exist(i,14) = 0.0
            exist(i,15) = 0.0
c     
          else
            isoybean = 0
            imaize   = 1
            iwheat   = 0
            exist(i,13) = 0.0
            exist(i,14) = 1.0
            exist(i,15) = 0.0
          endif
c
c 4:
c 3 crop rotation with winter wheat and soybean planted in same year
c winter wheat harvested in year 2
c maize grown in year 3
c 
        elseif (irotation .eq. 4) then 
c
c soybean planted/harvested  
c winter wheat planted
c
          if (xfrac .eq. 0.0) then 
            isoybean = 1
            imaize   = 0
            iwheat   = 2
            exist(i,13) = 1.0
            exist(i,14) = 0.0
            exist(i,15) = 1.0
c
c winter wheat harvested in year 2
c
          elseif (xfrac .eq. 1.0) then
            isoybean = 0
            imaize   = 0
            iwheat   = 2
            exist(i,13) = 0.0
            exist(i,14) = 0.0
            exist(i,15) = 1.0
c     
c maize planted/harvested in year 3
c
          else
            isoybean = 0
            imaize   = 1
            iwheat   = 0
            exist(i,13) = 0.0
            exist(i,14) = 1.0
            exist(i,15) = 0.0
          endif
        
        else if(irotation .eq. 5) then
c two-crop rotation (standard soy/corn)
c alternate between even / odd years
c
          if (xfrac .eq. 0) then
            isoybean = 0
            imaize   = 1
            iwheat   = 0
            exist(i,13) = 0.0
            exist(i,14) = 1.0
            exist(i,15) = 0.0
c
          else
            isoybean = 1
            imaize   = 0
            iwheat   = 0
            exist(i,13) = 1.0
            exist(i,14) = 0.0
            exist(i,15) = 0.0
          endif
        endif
 100    continue
c
c
c---------------------------note-------------------------------------------
c  when iyr_pft = 1, set exist(npoi,idx) from read-in NetCDF file
c
c  Y.Li 2014-06-18
c
c-------------------------End note-----------------------------------------
c
      if(iyr_pft.eq.1) then

        !read-in NetCDF file to cdummy
        call read_NetCDF('pft',iyear)

        !initialize and assign the plant type to exist()
        !Important: need to determine whether isoybean, imaize and iwheat need
        !to be modified, or extended to depend on grid point
        isoybean=0;imaize=0;iwheat=0
        exist(:,13:15) = 0
        do i =1, npoi
          if(cmask(i) == 0 ) cycle

          idx = int(cdummy(i))

          if(idx == 13)     then
            isoybean= 1
            icropsum(i) = 1  !simulate the crop
          elseif(idx == 14) then
            imaize  = 1
            icropsum(i) = 1  !simulate the crop
          elseif(idx == 15) then
            iwheat  = 1
            icropsum(i) = 1  !simulate the crop
          else
            !Y.Li (temporarily disable this warning, this warning might or might not needed)
            !write(6,*) '(rotation) Warning: read-in pft not in the range (13~15) at grid point: ',i
            !write(6,*) 'read-in type index: ',idx
          end if

          if(idx>=1.and.idx<=15) then !current crop range
            exist(i,idx) = 1.0
          elseif(idx /= 0) then
            write(6,*) 'Error! out of range read-in pft value at grid: ',idx,i
          end if

        end do !on i =1, npoi

      end if !on if(iyr_pft.eq.1)

      return
      end
c
c ----------------------------------------------------------------------------
      subroutine planting(irestart, irstyear, iyear0,iyear,imonth, iday, jday,ffact)
c ----------------------------------------------------------------------------
c
c subroutine to determine planting dates for crops based on observations
c and/or climatic constraints
c
c uses:
c
      use comgrid
      use compar
      use comatm
      use comsoi
      use comsum
      use comveg
      use comsno
      use comnitr
      use comcrop
      use combcs
      use comwork, only: cdummy
      use management_para

c
      implicit none
c
c local variables
c
      real
     >  ffact
c
      integer
     >  imonth,
     >  iday,
     >  jday,
     >  iyear,
     >  iyear0,
     >  irestart,
     >  irstyear,
     >  i, j, k, n

      integer:: days_sincePlanting, mngt_idx
c 
c WJS 01.27.10 Code that was here to set cfert* variables from hist.fert.1945.1996
c has been moved to read_cfert.f, which is called from main program
c

c        
c in order to only allow a crop to be planted once each year
c initialize cropplant = 0., but hold it = 1 through the end of
c the year
c
c initialize other variables that are calculated for crops
c on an annual basis in cropresidue subroutine
c
      if (iday .eq. 1 .and. imonth .eq. 1) then 
c
        do i = 1, npoi
          if(cmask(i) == 0 ) cycle
c 
          do j = scpft, ecpft   ! crop plant functional types only
c
c make sure variables aren't changed at beginning of the year
c for a crop the is currently planted (e.g. winter wheat)
c if else if statement were removed in the next sequence, then
c winter wheat grown continuously would amount to a wheat/fallow
c rotation because wheat would only be planted every other year 
c 
            if (croppresent(i,j) .eq. 0)  then
              cropplant(i,j) = 0.0
              pdate(i,j)     = 0.0
              idop(i,j)      = 999
c
c for continuous, annual winter wheat type crop
c
            else if (croppresent(i,j) .eq. 1 .and.
     >               j .eq. 15            .and.
     >               iwheat .eq. 2)        then
              cropplant(i,j) = 0.0
            endif
c
            harvdate(i,j)  = 999
            hdate(i,j)     = 0.0
            matdate(i,j)   = 0.0
            matgdd(i,j)    = 0.0
            matgrnfraccrop(i,j) = 0.0
            emerge_day(i,j) = 9999.
            emerge_gdd(i,j) = -1.0
            emerge_gddtsoi(i,j) = -1.0
            grainday(i,j)  = 9999.
            graingdd(i,j)  = -1.0
            harvidx(i,j)   = 0.0
            croplaimx(i,j) = 0.0
            grainn(i,j)    = 0.0
            cropyld(i,j)   = 0.0
            dmleaf(i,j)    = 0.0
            dmstem(i,j)    = 0.0
            dmroot(i,j)    = 0.0
            dmresidue(i,j) = 0.0
            dmyield(i,j)   = 0.0
            dmcrop(i,j)    = 0.0
            cropfixn(i,j)  = 0.0
            residuen(i,j)  = 0.0
            cropn(i,j)     = 0.0
            nconcl(i,j)    = 0.0
            nconcs(i,j)    = 0.0
            nconcr(i,j)    = 0.0
            nconcg(i,j)    = 0.0
            cntops(i,j)    = 40.0
            cnroot(i,j)    = 60.0
            fertinput(i,j) = 0.0
            apar(i,j)      = 0.0
            fkill_gdd_missed(i,j) = 0.0
            dev_fraction_attained(i,j) = 0.0

          end do
        end do
      end if
c
c Note that the following call is outside the do loops, because it
c contains its own do loops.
c
      if (management_prescribed) then
         call planting_prescribed(jday)
      else
         call planting_prognostic(iyear0, iyear, imonth, iday, jday)
      end if
c
c
c---------------------------note-------------------------------------------
c  when iyr_fert = 1, set variable <fertnitro> from read-in NetCDF file,
c  the original mechanism below this method (if (iday .eq. 1 .and. imonth.eq.1)...)
c  will be skipped all the time during the simulation 
c  
c possibly should be integrated into the original fertilizer procedure (when this function is in practical use)
c
c  Y.Li 2014-06-23
c
c-------------------------End note-----------------------------------------
c

      if(iyr_fert.eq.1) then

        !--temporary, need to confirm with the triggering mechanism
        if(jday .eq. 1) then

          !read-in NetCDF file to cdummy
          call read_NetCDF('fert',iyear)

          do i = 1, npoi
            if(cmask(i) == 0 ) cycle

            do j = scpft, ecpft     ! crop plant functional types only

              !--set all the fertnitro to 0 at jday, then only determined by the read-in value
              fertnitro(i,j) = 0.0

              if(exist(i,j) .eq. 1) then

                !assign read-in value, apply the fertilization rate by ffact
                !fertnitro(i,j) = cdummy(i) * ffact 
                fertnitro(i,j) = cdummy(i)*1e-04 * ffact 

              end if !on if(exist(i,j))

c assign annual fertilizer input in kg/ha; 
c
              fertinput(i,j)  = fertnitro(i,j)*1e+04

            end do !on i = 1, npoi
          end do !on j = scpft, ecpft

        end if !on if (jday)

        !--finish fertilizing, goto the end of the subroutine      
        goto 9999

      end if !on if(iyr_pft.eq.1)
c
c

      do i = 1, npoi
        if(cmask(i) == 0 ) cycle
c 
        do j = scpft, ecpft     ! crop plant functional types only
c
c add fertilizer nitrogen input for each crop planted (kg m-2)
c on the planting date
c either input here, or read from a gridded dataset
c is treated a single, broadcast pulse at planting to the top soil layer
c this fertilizer is assumed to be ammonium nitrate, in a 50/50 ratio
c of NH4/NO3
c
c define amount of fertilizer added when crop is planted
c use historical changes for years between 1945-1996 (Alexander et al.,) 
c only add fertilizer on day of planting - use croppresent funtion to
c make sure that in successive years, idop from previous year doesn't
c get applied here.
c
c also, don't cycle through all crop types - only the one that
c is planted this year...otherwise it will zero out the fertilizer
c values for pfts 13, 14 if going through all pfts through 15 
c
c
! Y.Li 2016-04-20: extend jday >=idop to enable daily fertilizer since planting date
! the structure already re-organized to easily integrate <iyr_fert> part here

           if (jday .eq. idop(i,j) .and. croppresent(i,j) .eq. 1
     >         .and. exist(i,j) .eq. 1.) then 
c

              if (iyear .lt. fert_year1) then ! Before we have any data: use a constant value
                fertnitro(i,13) = 0.0000       ! soybeans - kg_n m-2 y-1
                fertnitro(i,14) = 0.0000       ! maize    - kg_n m-2 y-1
                fertnitro(i,15) = 0.0000       ! wheat    - kg_n m-2 y-1

                !fertnitro(i,13) = 0.0002       ! soybeans - kg_n m-2 y-1
                !fertnitro(i,14) = 0.0009       ! maize    - kg_n m-2 y-1
                !fertnitro(i,15) = 0.0002       ! wheat    - kg_n m-2 y-1

              else if (iyear .gt. fert_lastyear) then  !For years past when we have any data: use the last year of data (we will multiply this by ffact below)

                fertnitro(i,13) = fertsoy(i,fert_nyears)   * 1e-04 ! soybeans - kg_n m-2 y-1
                fertnitro(i,14) = fertmaize(i,fert_nyears) * 1e-04 ! maize    - kg_n m-2 y-1
                fertnitro(i,15) = fertwheat(i,fert_nyears) * 1e-04 ! wheat    - kg_n m-2 y-1
c
c$$$c modified for new updated fertilizer data from Simon 4.27.07
c$$$c
c$$$                  fertnitro(i,13) = fert05soy(i)   * 1e-04 ! soybeans - kg_n m-2 y-1
c$$$                  fertnitro(i,14) = fert05maize(i) * 1e-04 ! maize    - kg_n m-2 y-1
c$$$c                 fertnitro(i,14) = fert05sorg(i)  * 1e-04 ! sorghum  - kg_n m-2 y-1
c$$$                  fertnitro(i,15) = fert05wheat(i) * 1e-04 ! wheat    - kg_n m-2 y-1
c
              else  ! fert_year1 <= iyear <= fert_lastyear: years for which we have data
                 
                 if(flg_mngt == 1 .and. trim(fertilizer_mngt_mode) == 'years_ave') then

                   fertnitro(i,13) = sum(fertsoy  (i,fertMngt_ys_idx:fertMngt_ye_idx))/real(fertMngt_nyears) * 1e-04   ! soybeans - kg_n m-2 y-1
                   fertnitro(i,14) = sum(fertmaize(i,fertMngt_ys_idx:fertMngt_ye_idx))/real(fertMngt_nyears) * 1e-04   ! maize    - kg_n m-2 y-1
                   fertnitro(i,15) = sum(fertwheat(i,fertMngt_ys_idx:fertMngt_ye_idx))/real(fertMngt_nyears) * 1e-04   ! wheat    - kg_n m-2 y-1

                 else
                   fertnitro(i,13) = fertsoy(i,iyear+1-fert_year1)   * 1e-04       ! soybeans - kg_n m-2 y-1
                   fertnitro(i,14) = fertmaize(i,iyear+1-fert_year1) * 1e-04       ! maize    - kg_n m-2 y-1
                   fertnitro(i,15) = fertwheat(i,iyear+1-fert_year1) * 1e-04       ! wheat    - kg_n m-2 y-1

                 end if  !on if flg_mngt & fertilizer_mngt_mode
 
              endif !on if (iyear .lt. fert_year1)
c
c constant 2000 level
c set unfertilized cases for soybeans 4.24.07
c WJS 01.27.10: If you want to use this commented-out code, you should change the hard-coded 51
c to something like fert_nyears, as is done above
c
c             fertnitro(i,13) = fertsoy(i,51)   * 1e-04       ! soybeans - kg_n m-2 y-1
c             fertnitro(i,13) = 0.0                           ! soybeans - kg_n m-2 y-1
c             fertnitro(i,14) = fertmaize(i,51) * 1e-04       ! maize    - kg_n m-2 y-1
c             fertnitro(i,15) = fertwheat(i,51) * 1e-04       ! wheat    - kg_n m-2 y-1
c
c take care of new fertilizer data that doesn't exist for Canada
c use the old values from the calculated fertilizer use over the US (historical)

              if (fertnitro(i,13) .gt. 2e+03) then
                call fill_fert(i, 13, cfertsoy, iyear)
              end if
              if (fertnitro(i,14) .gt. 2e+03) then
                call fill_fert(i, 14, cfertmaize, iyear)
              end if
              if (fertnitro(i,15) .gt. 2e+03) then
                call fill_fert(i, 15, cfertwheat, iyear)
              end if

c
c
c constant mid 1990s levels for grid cells outside of US
c (This code can replace the calls to fill_fert, above)
c
c               if (fertnitro(i,13) .gt. 1e+03) fertnitro(i,13) = cfertsoy(1996) 
c               if (fertnitro(i,14) .gt. 1e+03) fertnitro(i,14) = cfertmaize(1996) 
c               if (fertnitro(i,15) .gt. 1e+03) fertnitro(i,15) = cfertwheat(1996) 


c
c For years > fert_lastyear, we multiply the fertilization rate by ffact:
              if (iyear .gt. fert_lastyear) then
                fertnitro(i,13) = fertnitro(i,13) * ffact
                fertnitro(i,14) = fertnitro(i,14) * ffact
                fertnitro(i,15) = fertnitro(i,15) * ffact
              end if

c
c
c or constant fertilizer application each year
c
c                fertnitro(i,13) = 0.0025      ! soybeans - kg_n m-2 y-1
c                fertnitro(i,14) = 0.0180      ! maize    - kg_n m-2 y-1
c                fertnitro(i,15) = 0.0080      ! wheat    - kg_n m-2 y-1
c
c assign annual fertilizer input in kg/ha

                fertinput(i,j)  = fertnitro(i,j) * 1.e+04
c
             !store the one-year fertilizer when daily distribution is requested
             if(trim(fertilizer_mngt_mode) == 'dat') then

               mngt_fertnitro0_13(i) = fertnitro(i,13)
               mngt_fertnitro0_14(i) = fertnitro(i,14)
               mngt_fertnitro0_15(i) = fertnitro(i,15)

               !at planting date, assign the value specified; fertMngt_pct = [0,1]
               days_sincePlanting = jday - idop(i,j)
               mngt_idx = days_sincePlanting + 1  !index starts from 1; days_sincePlanting starts from 0
               fertnitro(i,13) = fertMngt_pct(mngt_idx)*mngt_fertnitro0_13(i)
               fertnitro(i,14) = fertMngt_pct(mngt_idx)*mngt_fertnitro0_14(i)
               fertnitro(i,15) = fertMngt_pct(mngt_idx)*mngt_fertnitro0_15(i)

             end if !on if fertilizer_mngt_mode

           else !on if (jday .eq. idop(i,j)...

             if(flg_mngt == 1 .and. trim(fertilizer_mngt_mode) == 'dat') then
               if(jday > idop(i,j) .and. croppresent(i,j) .eq. 1 .and. exist(i,j) .eq. 1.) then

                 days_sincePlanting = jday - idop(i,j)
                 mngt_idx = days_sincePlanting + 1  !index starts from 1; days_sincePlanting starts from 0
                 fertnitro(i,13) = fertMngt_pct(mngt_idx)*mngt_fertnitro0_13(i)
                 fertnitro(i,14) = fertMngt_pct(mngt_idx)*mngt_fertnitro0_14(i)
                 fertnitro(i,15) = fertMngt_pct(mngt_idx)*mngt_fertnitro0_15(i)

               end if
             else  !no daily distribution

               if (exist(i,j) .eq. 1.0) then   !default: 0 for other days
                 fertnitro(i,j) = 0.0       
               end if

             end if  !on if flg_mngt & fertilizer_mngt_mode

           endif   ! planting date
c
         end do !on j = scpft, ecpft
       end do !on i = 1, npoi
c
c end of loop
c
c

9999  continue

c return to main program
c 
      return
      end
c
c ---------------------------------------------------------------------------------
      subroutine planting_prognostic(iyear0, iyear, imonth, iday, jday)
c ---------------------------------------------------------------------------------
c
c Subroutine to determine planting dates prognostically, based on
c climatic constraints. This is only meant to be called by the planting
c subroutine. This also sets gddmaturity. But note that, at least as of
c 5-6-10, this can only be used with prognostic cultivars: it does not
c include the code that is required to set cultivars based on
c prescribed, read-in maps (see planting_prescribed for how cultivars
c are set when using prescribed cultivars).
c
c WJS 05.06.10: The code here was moved from the planting() subroutine.
c Note that there is still some commented-out code here that refers to
c prescribed planting dates & hybrids (read in from files): I kept the
c code & comments as is, including that commented-out code that is out
c of place in this "prognostic" subroutine.
c
c uses:
c
      use comgrid
      use comcrop
      use comveg
      use comsum
      use combcs,only: cmask
      use compft,only: tauleaf,tauroot,tauwood0
c
      implicit none
c
c subroutine arguments
c
      integer
     >  iyear0,
     >  iyear,
     >  imonth,
     >  iday,
     >  jday

c local variables
c
      real   
     >  ptemp(npft), 
     >  pmmin(npft), 
     >  pdmin(npft), 
     >  pmmax(npft), 
     >  pdmax(npft), 
     >  pmintemp(npft), 
     >  sumhy(npoi,npft),
     >  sumdp(npoi,npft)

      real
     >  yc

      integer
     >  i,
     >  j,
     >  iy

      integer:: coverCrop_idx !cover crop index: 11 (c4 grasses) or 12 (c3 grasses)
      real::    tauwood_tmp   !local tauwood for c3 biomass (cover crops)

c
c minimum temperature required for crop planting and vegetative growth
c from EPIC model parameterizations
c
c the planting range in the US for maize is typically from April 10 - May 10
c the most active planting date in  US for soybean is typically May 15 - June 20 
c spring wheat planting dates are typically early April through mid-May - in line with maize 
c winter wheat is from Sept. 1 through early Nov., and is typically planted within 10-14
c days of the first likely frost event.
c
      ptemp(13)   = 284.16          ! soybean minimum 10 day average temperature for planting (K)
      ptemp(14)   = 283.16          ! maize minimum   10 day average temperature for planting (K)
      pmintemp(13)= 279.16          ! soybean minimum 10 day average min temp for planting (K) 
c     pmintemp(14)= 277.16          ! maize minimum 10 day average min temp for planting (K) 
      pmintemp(14)= 279.16          ! maize minimum 10 day average min temp for planting (K) 
      pmmin(13)   = 5               ! soybean earliest month to plant (month)
      pmmin(14)   = 4               ! maize earliest month to plant (month)
      pdmin(13)   = 1               ! soybean earliest day in earliest month to plant (day)
      pdmin(14)   = 1               ! maize earliest day in earliest month to plant (day)
c
c spring wheat planting date requirements
c note: iwheat = 1 (spring), iwheat = 2 (winter) 
c
      if (iwheat .eq. 1) then
        ptemp(15)   = 280.16        ! spring wheat minimum 10 day average temperature for planting (K) 
        pmintemp(15)= 272.16        ! spring wheat minimum 10 day average min temp threshold for planting (K)
        pmmin(15)   = 4             ! spring wheat earliest month to plant (month)
        pdmin(15)   = 1             ! spring wheat earliest day in earliest month to plant (day)
c
c typically winter wheat is planted when average minimum temperature gets to about 40 F
c and is planted no later than November 
c
      else if (iwheat .eq. 2) then  ! winter wheat
        pmintemp(15) = 278.16       ! average 5-day minimum temperature needed to plant winter wheat (~ 40 F) 
        pmmin(15)    = 9            ! earliest month in which winter wheat is planted (month)
        pdmin(15)    = 1            ! earliest day (of earliest month) in which winter wheat is planted (day) 
        pmmax(15)    = 11           ! latest possible month to plant winter wheat (month)
        pdmax(15)    = 30           ! latest day in latest month to plant winter wheat (day) 
      endif

      do i = 1, npoi
        if(cmask(i) == 0 ) cycle
c 
       do j = scpft, ecpft       ! crop plant functional types only

         if (exist(i,j) .eq. 1. .and. croppresent(i,j) .ne. 1 .and.
     >     cropplant(i,j) .eq. 0.0) then 

c gdd needed for * chosen crop and a likely hybrid (for that region) *
c to reach full physiological maturity
c
c based on accumulated seasonal average growing degree days from April 1 - Sept 30 (inclusive)
c for corn and soybeans in the United States -
c decided upon by what the typical average growing season length is
c
c and the gdd needed to reach maturity in those regions
c
c first choice is used for spring wheat and/or soybeans and maize
c
c         idop(i,j)       = max(90, int(xinpdate(i,iyear-iyear0+1)))  ! has to be later or equal to April 1  
c
c ------
c CJK 8/1/06
c for constant, average planting date input from *.nc file
c
c         idop(i,j)       = int(xinavepdate(i))
c ------
c
          if  ((j .eq. 13 .or. j .eq. 14) .or.
     >         (j .eq. 15 .and. iwheat .eq. 1))  then 
c
             if (iyear .eq. iyear0) then         
c
c for constant planting date
c
c               if (jday .eq. idop(i,j)) then
c
                if (a10td(i)         .gt. ptemp(j)     .and.         ! 10-day average soil temperature
     >             a10tmin(i)       .gt. pmintemp(j)  .and. 
     >             imonth           .ge. pmmin(j)     .and.         ! impose earliest planting date 
     >             iday             .ge. pdmin(j))    then 
cc               jday             .eq. idop(i,j)    .and.               
cc    >          gdd8(i)          .ge. gddmin(j))    then         ! impose limit on growing season length needed
c
                  croplive(i,j)   = 1.0        ! initialize freeze kill function to 1.0 - crops living 
                  cropplant(i,j)  = 1.0        ! initialize freeze kill function to 1.0 - crops living 
                  croppresent(i,j) = 1
                  grnfraccrop(i,j) = 1.0
                  idop(i,j)         = jday     ! has to be later or equal to April 1  
                  if (j .eq. 13) soydop(i,iyear)   = jday
                  if (j .eq. 14) corndop(i,iyear)  = jday
                  if (j .eq. 15) whtdop(i,iyear)   = jday
c
                  gddmaturity(i,13) = max(950., min (gdd10(i), hybgdd(j))) 
                  gddmaturity(i,14) = max(950., min (gdd10(i)  * 0.90, hybgdd(j)))  ! WJS (4-22-10): Changed base from gdd8 to gdd10
                  gddmaturity(i,15) = max(950., min (gdd0c(i), hybgdd(j))) 
c
c                 gddmaturity(i,14) = hybgdd(14) 
c -----
c CJK 8/1/06
c for constant, average hybrid input from *.nc file
c
c                 gddmaturity(i,14) = max(950.0, min(1800.,xinavehybrid(i)))
c -----
c
               endif 
c
             elseif (iyear .gt. iyear0 .and. iyear .lt. iyear0+5) then       ! after initial spinup for crop average planting dates
c ------
c CJK 8/1/06
c for constant, average planting date input from *.nc file
c
c                idop(i,j)       = int(xinavepdate(i))
c ------
c CJK 8/1/06
c for constant planting date
c
c                 if (jday .eq. idop(i,j)) then
c
c
              
c
               if (a10td(i)         .gt. ptemp(j)     .and.              ! 10-day average soil temperature
     >             a10tmin(i)       .gt. pmintemp(j)  .and. 
     >             imonth           .ge. pmmin(j)     .and.              ! impose earliest planting date 
     >             iday             .ge. pdmin(j))    then 
                     croplive(i,j)   = 1.0        ! initialize freeze kill function to 1.0 - crops living 
                     cropplant(i,j)  = 1.0        ! initialize freeze kill function to 1.0 - crops living 
                     croppresent(i,j) = 1
                     grnfraccrop(i,j) = 1.0
                     idop(i,j)         = jday     ! has to be later or equal to April 1  
                     if (j .eq. 13) soydop(i,iyear)   = jday
                     if (j .eq. 14) corndop(i,iyear)  = jday
                     if (j .eq. 15) whtdop(i,iyear)   = jday
c
                     gddmaturity(i,13) = max(950., min (gddcorn(i,iyear-1) *0.80, hybgdd(j)))  ! assign hybrid based on last year 
                     gddmaturity(i,14) = max(950., min (gddcorn(i,iyear-1) *0.90, hybgdd(j)))  ! assign hybrid based on last year 
                     gddmaturity(i,15) = max(950., min (gddcorn(i,iyear-1) *1.2,  hybgdd(j)))  ! assign hybrid based on last year 
c
c                 gddmaturity(i,14) = hybgdd(14) 
c -----
c CJK 8/1/06
c for constant, average hybrid input from *.nc file
c
c                 gddmaturity(i,14) = max(950.0, min(1800.,xinavehybrid(i)))
c -----
c
c
               endif
c
             else 
c
c CJK 8/1/06
c for constant planting dates and optimum hybrids with no adjustment
c
c              idop(i,j)         = int(xinavepdate(i))
c-------
               sumdp(i,j) = 0
               sumhy(i,j) = 0
c
c insert here - do iy from iyear-5 to  iyear-1 --> previous 5 year mean of hybrids planted
c keep planting date flexible for that year's weather conditions 
c
               yc = 5.0
c              yc = 11.0
c              do iy = 1949, 1959          ! 11 year averaging spinup for crops
               do iy = iyear-5,iyear-1     ! hybrid based on previous 5 year average - farm management 
c
                 if     (j .eq. 13) then 
                   sumhy(i,j) = sumhy(i,j) + gddcorn(i,iy) * 0.8
                   sumdp(i,j) = sumdp(i,j) + float(soydop(i,iy))
                 elseif (j .eq. 14) then 
                   sumhy(i,j) = sumhy(i,j) + gddcorn(i,iy) * 0.9
                   sumdp(i,j) = sumdp(i,j) + float(corndop(i,iy))
                 elseif (j .eq. 15) then
                   sumhy(i,j) = sumhy(i,j) + gddcorn(i,iy) * 1.2
                   sumdp(i,j) = sumdp(i,j) + float(whtdop(i,iy))
                 endif
               enddo
c
c               xinavehybrid(i) = sumhy/yc
c
                avehybrid(i,j)    = sumhy(i,j) / yc
                iavepdate(i,j)    = int(sumdp(i,j)/yc)
c
c              if (jday .eq. iavepdate(i,j)) then
c
c CJK 8/1/06
c              if (jday .eq. idop(i,j)) then
c
c

               if (a10td(i)         .gt. ptemp(j)     .and.              ! 10-day average soil temperature
     >             a10tmin(i)       .gt. pmintemp(j)  .and. 
     >             imonth           .ge. pmmin(j)     .and.         ! impose earliest planting date 
     >             iday             .ge. pdmin(j))    then 
                     croplive(i,j)   = 1.0        ! initialize freeze kill function to 1.0 - crops living 
                     cropplant(i,j)  = 1.0        ! initialize freeze kill function to 1.0 - crops living 
                     croppresent(i,j) = 1
                     grnfraccrop(i,j) = 1.0
c                    idop(i,j)       = iavepdate(i,j) 
                     idop(i,j)       = jday  
                     gddmaturity(i,j) = max(950.,min(avehybrid(i,j), hybgdd(j)))   ! using constant hybrid for location over entire period 

c                    gddmaturity(i,14) = hybgdd(14) 
c CJK 8/1/06
c                    gddmaturity(i,14) = max(950.0, min(1800., xinavehybrid(i)))
                endif

             endif   ! planting/year option
c
c              gddmaturity(i,13) = min (gdd10(i), hybgdd(j)) 
c              gddmaturity(i,14) = max(950., min (gdd10(i) * 0.85, hybgdd(j))) 
c              gddmaturity(i,15) = min (gdd0c(i), hybgdd(j))  ! for spring wheat  
c
c plant everywhere if not planted by middle of june  
c based on approximate planting dates designated by the USDA
c
c             else if (imonth .ge. 6 .and. iday .ge. 15) then
c                  croplive(i,j)   = 1.0        ! initialize freeze kill function to 1.0 - crops living 
c                  cropplant(i,j)  = 1.0        ! initialize freeze kill function to 1.0 - crops living 
c                  croppresent(i,j) = 1
c                  grnfraccrop(i,j) = 1.0
c                  idop(i,j)       = jday
c                  gddmaturity(i,13) = min (gdd10(i), hybgdd(j)) 
c                  gddmaturity(i,14) = max(950., min (gdd10(i) * 0.85, hybgdd(j))) 
cc                 gddmaturity(i,14) = min (hybgdd(j) , hybgdd(j)) 
c                  gddmaturity(i,15) = min (gdd0c(i), hybgdd(j))  ! for spring wheat  
c
c winter wheat : use gdd0c as a limit to plant winter wheat
c
           elseif (j .eq. 15 .and. iwheat .eq. 2) then   ! plant winter wheat
c
c add check to only plant winter wheat after other crops (soybean, maize)
c have been harvested 
c
c *** remember order of planting is crucial - in terms of which crops you want
c to be grown in what order ***    
c
c in this case, corn or soybeans are assumed to be planted before
c wheat would be in any particular year that both pfts are allowed
c to grow in the same grid cell (e.g., double-cropping)  
c
            if (iyear .eq. iyear0) then
c
             if (a5tmin(i)        .le. pmintemp(j) .and.
     >           imonth           .ge. pmmin(j)    .and.
     >           iday             .ge. pdmin(j)    .and.
     >           (harvdate(i,13)  .ne. 999         .or.
     >            harvdate(i,14)  .ne. 999         .or.
     >            irotation       .eq. 0)          .and. 
     >            gdd0c(i)        .ge. gddmin(j))         then
c
                   croplive(i,j)     = 1.0     
                   cropplant(i,j)    = 1.0     
                   croppresent(i,j) = 1
                   grnfraccrop(i,j) = 1.0
                   idop(i,j)         = jday
                   whtdop(i,iyear)   = jday
c                  gddmaturity(i,15) = hybgdd(j) 
                   gddmaturity(i,15) = max(950., min (gdd0c(i) * 0.90, hybgdd(j))) 
c
c plant winter wheat at latest possible date 
c and after all other crops were harvested for that year
c
             elseif (imonth  .ge. pmmax(j)     .and.
     >               iday    .ge. pdmax(j)     .and.
     >              (harvdate(i,13)  .ne. 999  .or.
     >               harvdate(i,14)  .ne. 999  .or. 
     >               irotation       .eq. 0)   .and. 
     >               gdd0c(i).ge. gddmin(j)) then  
c
                       croplive(i,j)     = 1.0     
                       cropplant(i,j)    = 1.0     
                       croppresent(i,j) = 1
                       grnfraccrop(i,j) = 1.0
                       idop(i,j)         = jday
c                      gddmaturity(i,15) = hybgdd(j) 
                       gddmaturity(i,15) = max(950., min (gdd0c(i), hybgdd(j))) 
             endif
c
           elseif (iyear .gt. iyear0 .and. iyear .lt. iyear0+5) then       ! after initial spinup for crop average planting dates
c
             if (a5tmin(i)        .le. pmintemp(j) .and.
     >           imonth           .ge. pmmin(j)    .and.
     >           iday             .ge. pdmin(j)    .and.
     >           (harvdate(i,13)  .ne. 999         .or.
     >            harvdate(i,14)  .ne. 999         .or.
     >            irotation       .eq. 0)          .and. 
     >            gdd0c(i)        .ge. gddmin(j))         then
c
                   croplive(i,j)     = 1.0     
                   cropplant(i,j)    = 1.0     
                   croppresent(i,j) = 1
                   grnfraccrop(i,j) = 1.0
                   idop(i,j)         = jday
                   whtdop(i,iyear)   = jday
                   gddmaturity(i,15) = max(950., min (gddcorn(i,iyear-1) * 1.20, hybgdd(j)))  ! assign hybrid based on last year 
c
c plant winter wheat at latest possible date 
c and after all other crops were harvested for that year
c
             elseif (imonth  .ge. pmmax(j)     .and.
     >               iday    .ge. pdmax(j)     .and.
     >              (harvdate(i,13)  .ne. 999  .or.
     >               harvdate(i,14)  .ne. 999  .or. 
     >               irotation       .eq. 0)   .and. 
     >               gdd0c(i).ge. gddmin(j)) then  
c
                       croplive(i,j)     = 1.0     
                       cropplant(i,j)    = 1.0     
                       croppresent(i,j) = 1
                       grnfraccrop(i,j) = 1.0
                       idop(i,j)         = jday
                       gddmaturity(i,15) = max(950., min (gddcorn(i,iyear-1) * 1.20, hybgdd(j)))  ! assign hybrid based on last year 
             endif
c
           else
               sumdp(i,j) = 0
               sumhy(i,j) = 0
c
c insert here - do iy from iyear-5 to  iyear-1 --> previous 5 year mean of hybrids planted
c keep planting date flexible for that year's weather conditions 
c
               yc = 5.0
c              yc = 11.0
c              do iy = 1949, 1959          ! 11 year averaging spinup for crops
               do iy = iyear-5,iyear-1     ! hybrid based on previous 5 year average - farm management 
                   sumhy(i,j) = sumhy(i,j) + gddcorn(i,iy) * 1.2
                   sumdp(i,j) = sumdp(i,j) + float(whtdop(i,iy))
               enddo
c
                avehybrid(i,j)    = sumhy(i,j) / yc
                iavepdate(i,j)    = int(sumdp(i,j)/yc)
c
c              if (jday .eq. iavepdate(i,j)) then
c
                if (a5tmin(i)        .le. pmintemp(j) .and.
     >           imonth           .ge. pmmin(j)    .and.
     >           iday             .ge. pdmin(j)    .and.
     >           (harvdate(i,13)  .ne. 999         .or.
     >            harvdate(i,14)  .ne. 999         .or.
     >            irotation       .eq. 0)          .and. 
     >            gdd0c(i)        .ge. gddmin(j))         then
c
                   croplive(i,j)     = 1.0     
                   cropplant(i,j)    = 1.0     
                   croppresent(i,j) = 1
                   grnfraccrop(i,j) = 1.0
c                  idop(i,j)         = iavepdate(i,j) 
                   idop(i,j)         = jday
                   gddmaturity(i,15) = max(950., min (avehybrid(i,j), hybgdd(j)))  ! assign hybrid based on last year 
c
c plant winter wheat at latest possible date 
c and after all other crops were harvested for that year
c
             elseif (imonth  .ge. pmmax(j)     .and.
     >               iday    .ge. pdmax(j)     .and.
     >              (harvdate(i,13)  .ne. 999  .or.
     >               harvdate(i,14)  .ne. 999  .or. 
     >               irotation       .eq. 0)   .and. 
     >               gdd0c(i).ge. gddmin(j)) then  
c
                       croplive(i,j)     = 1.0     
                       cropplant(i,j)    = 1.0     
                       croppresent(i,j) = 1
                       grnfraccrop(i,j) = 1.0
                       idop(i,j)         = jday
c                      idop(i,j)         = iavepdate(i,j) 
                       gddmaturity(i,15) = max(950., min (avehybrid(i,j), hybgdd(j)))  ! assign hybrid based on last year 
             endif
           endif  ! (year check)
          endif

           !remove cover crops at the planting date (Y.Li)
           if(flg_coverCrops /= 0 .and. (j==13 .or. j==14) ) then  !to plant soybean or corn
             if(idop(i,j) == jday) then !crop is being planting

               if(flg_coverCrops == 1) then  !c4 grasses
                 coverCrop_idx = 11
               elseif(flg_coverCrops == 2) then  !c3 grasses
                 coverCrop_idx = 12
               else
                 write(6,*) '(planting_prognostic)Error! unsupported flg_coverCrops: ',flg_coverCrops
                 write(6,*) 'stopping the program...'
                 stop 
               end if

               exist(i,coverCrop_idx) = 0   !kill c3/c4 grass

               !move biomass to litter pool
               tauwood_tmp = tauwood0(coverCrop_idx)/2.0    !from vegetation.f, pft cannot exist, then tauwood = taufin

               falll(i) = falll(i) + cbiol(i,coverCrop_idx) / tauleaf(coverCrop_idx)
               fallr(i) = fallr(i) + cbior(i,coverCrop_idx) / tauroot(coverCrop_idx)
               fallw(i) = fallw(i) + cbiow(i,coverCrop_idx) / tauwood_tmp

               cbiol(i,coverCrop_idx) = 0.0
               cbior(i,coverCrop_idx) = 0.0
               cbiow(i,coverCrop_idx) = 0.0

             end if !on if idop
           end if !on if flg_coverCrops

         endif  ! crop existence

        end do  ! j = scpft, ecpft


        
      end do  ! i = 1, npoi



      return
      end

c ---------------------------------------------------------------------------------
      subroutine planting_prescribed(jday)
c ---------------------------------------------------------------------------------
c
c Subroutine to determine planting dates based on prescribed, read-in
c variables. This is only meant to be called by the planting subroutine.
c This also sets gddmaturity. But note that, at least as of 5-6-10, this
c can only be used with prescribed cultivars: it does not include the
c logic that is required for prognostic cultivars (see
c planting_prognostic for how this logic works when using prognostic
c cultivars).
c
c uses:
c
      use comgrid
      use compar
      use comcrop
      use comveg
      use combcs,only: cmask
      use compft,only: tauleaf,tauroot,tauwood0
c
      implicit none
c
c Arguments
c 
      integer jday  ! current day of year
c
c Local variables
c
      integer i, j

      integer:: coverCrop_idx !cover crop index: 11 (c4 grasses) or 12 (c3 grasses)
      real:: tauwood_tmp      !local tauwood for c3 biomass (cover crops)
c
      do i = 1, npoi
        if(cmask(i) == 0 ) cycle

        do j = scpft, ecpft     ! loop over crop plant functional types only

          if (exist(i,j) .eq. 1. .and. croppresent(i,j) .ne. 1 .and.
     >      cropplant(i,j) .eq. 0.0) then  
        
c Note that there is a single xinpdate variable that applies to all crop
c types (and similarly for xinhybrid); in practice this means that you
c need to rename files depending on the crop type used in a given run,
c so that the proper input file is read for that crop. This also means
c that crop rotations cannot be used in conjunction with
c management_prescribed=TRUE.
c
c Note also that xinpdate is used directly, unlike in Chris's old code
c where there was an imposed minimum of 90. Chris had a comment
c associated with idop - "has to be later or equal to April 1" - but I
c am not sure why that is. It is possible that that minimum should be
c applied here to avoid breaking anything.
            if (jday .eq. nint(xinpdate(i))) then ! time to plant!
              croplive(i,j) = 1.0
              cropplant(i,j) = 1.0
              croppresent(i,j) = 1
              grnfraccrop(i,j) = 1.0
              idop(i,j) = jday
              
c Set gddmaturity using the prescribed, read-in cultivar variable:
              gddmaturity(i,j) = xinhybrid(i)

              !remove cover crops at the planting date (Y.Li)
              if(flg_coverCrops /= 0 .and. (j==13 .or. j==14) ) then  !to plant soybean or corn
                if(idop(i,j) == jday) then !crop is being planting

                  if(flg_coverCrops == 1) then  !c4 grasses
                    coverCrop_idx = 11
                  elseif(flg_coverCrops == 2) then  !c3 grasses
                    coverCrop_idx = 12
                  else
                    write(6,*) '(planting_prognostic)Error! unsupported flg_coverCrops: ',flg_coverCrops
                    write(6,*) 'stopping the program...'
                    stop 
                  end if

                  exist(i,coverCrop_idx) = 0   !kill c3/c4 grass

                  !move biomass to litter pool
                  tauwood_tmp = tauwood0(coverCrop_idx)/2.0    !from vegetation.f, pft cannot exist, then tauwood = taufin

                  falll(i) = falll(i) + cbiol(i,coverCrop_idx) / tauleaf(coverCrop_idx)
                  fallr(i) = fallr(i) + cbior(i,coverCrop_idx) / tauroot(coverCrop_idx)
                  fallw(i) = fallw(i) + cbiow(i,coverCrop_idx) / tauwood_tmp

                  cbiol(i,coverCrop_idx) = 0.0
                  cbior(i,coverCrop_idx) = 0.0
                  cbiow(i,coverCrop_idx) = 0.0

                end if !on if idop
              end if !on if flg_coverCrops

            end if
          end if  ! crop existence
        end do
      end do

      return
      end
      

c ---------------------------------------------------------------------------------
      subroutine fill_fert(i, pft, cfert, iyear)
c ---------------------------------------------------------------------------------
c
c Set fertnitro(i, pft) to the appropriate year's data from cfert.
c This should be called for points i that have missing fertilizer data.
c
      use comgrid
      use compar
      use comcrop
      use comnitr
c
      implicit none
c
c Subroutine arguments
      integer 
     >  i,     ! current grid cell
     >  pft,   ! pft for which to fill fertnitro
     >  iyear  ! current year

      real cfert(*)  ! the cfert* array to use to fill fertnitro(i, pft) (e.g., cfertsoy, cfertmaize, cfertwheat)

      if (iyear .lt. cfert_year1) then  ! years before we have data: use a constant value
        fertnitro(i,pft) = 0.  
      else if (iyear .gt. cfert_lastyear) then  ! years after the dataset ends: use last year of data
        fertnitro(i,pft) = cfert(cfert_lastyear)  
      else  ! cfert_year1 <= iyear <= cfert_lastyear: years for which we have data
        fertnitro(i,pft) = cfert(iyear)
      end if

      return
      end
c
c ---------------------------------------------------------------------------------
      subroutine vernalization(iyear,imonth, iday, jday)
c ---------------------------------------------------------------------------------
c
c last updated by C. Kucharik 2.19.2010
c
c * * * only call subroutine for winter wheat * * *
c
c subroutine calculates vernalization and photoperiod effects on gdd accumulation
c in winter wheat varieties. Thermal time accumulation is reduced in 1st period until
c plant is fully vernalized. During this time of emergence to spikelet formation,  
c photoperiod can also have a drastic effect on plant development.
c
c
c uses:
c
      use comgrid
      use compar
      use comatm
      use comsoi
      use comsum
      use comveg
      use comsno
      use comnitr
      use comcrop
      use combcs,only: cmask
c
      implicit none
c
c local variables
c
        real  p1d, p1v,
     >        tcrown,
     >        vd,vd1, vd2,
     >        tkil, tbase,
     >        hti 
c
        integer  i,
     >           iyear,
     >           imonth,
     >           iday,
     >           jday
 
c photoperiod factor calculation      
c genetic constant - can be modified
c
        p1d = 0.004  ! average for genotypes from Ritchey, 1991.
c                    ! Modeling plant and soil systems: Wheat phasic development
        p1v = 0.003  ! average for genotypes from Ritchey, 1991.  
c
        do 100 i = 1, npoi
          if(cmask(i) == 0 ) cycle
c
c only go through routine if winter wheat has been planted, is living,
c and the factor is not already 1.0
c 
         if (croplive(i,15) .eq. 1.0 .and. vf(i) .ne. 1.0) then 

c
c modified by C. Kucharik 2.19.2010
c
c pgs 52-53 "Wheat Phasic Devlopment" by J.T. Ritchie, in the book
c Modeling plant and soil systems; Agronomy No. 31; edited by John Hanks
c and J.T. Ritchie
c
c The vernalizaton routine is taken directly from the CERES-Wheat model code
c
c       jj   = latindex(i)
c       xlat = latscale(jj) 
c       s1 = sin(xlat * 0.01745)
c       c1 = cos(xlat * 0.01745)
c       dec  = 0.4093 * sin(0.0172 * (jday - 82.2)) 
c       dlv  = ((-s1 * sin(dec) - 0.1047)/(c1 * cos(dec)))
c       if (dlv .lt. -0.87) dlv = -0.87
c       hrlt = 7.639 * acos(dlv)  
c       df   = 1 - p1d * (20.0 - hrlt)**2
c
c daylength calculation in IBIS is in minutes - convert to hours
c
          df(i)   = 1 - p1d * (20.0 - daylength(i)/60.0)**2
c
c for all equations - temperatures must be in degrees (C)
c calculate temperature of crown of crop (e.g., 3 cm soil temperature)
c snow depth in centimeters
c
          if (td(i) .lt. 273.16) then
             tcrown = 2.0 + (td(i) - 273.16) * (0.4 + 0.0018 *
     >                (min(adsnod(i)*100., 15.0) - 15.0)**2)
          else
             tcrown = td(i) - 273.16 
          endif
c 
c vernalization factor calculation
c if vf(i) = 1.  then plant is fully vernalized - and thermal time
c accumulation in phase 1 will be unaffected 
c cumulative vd will need to be reset to 0 at planting of crop
c
          if (vf(i) .ne. 1.0) vd = 0.
          if (vf(i) .ne. 1.0 .and. tmax(i) .gt. 273.16) then
             if (tmin(i) .le. 288.16) then
               vd1      = 1.4 - 0.0778 * (tcrown)
               vd2      = 0.5 + 13.44/((tmax(i) - tmin(i) + 3.0)**2) * tcrown
               vd       = max(0.0, min(1., vd1, vd2))
               cumvd(i) = cumvd(i) + vd
             endif
c
             if (cumvd(i) .lt. 10 .and. tmax(i) .gt. 303.16) then
               cumvd(i) = cumvd(i) - 0.5 * (tmax(i) - 303.16)   
             endif 
             cumvd(i) = max(0.0, cumvd(i))    ! must be greater than 0 
c
             vf(i) = 1.0 - p1v * (50.0 - cumvd(i))
             vf(i) = max(0.0, min(vf(i), 1.0))       ! must be between 0 - 1
          endif
c             
c calculate cold hardening of plant
c determines for winter wheat varieties whether the plant has completed
c a period of cold hardening to protect it from freezing temperatures.  If  it has
c not, then exposure could result in death or killing of plants.
c
c there are two distinct phases of hardening 
c
c note:  tcrown temperatures are in deg C, others for tmax, tmin are given in deg K.
c
       tbase = 0.0
c
         if (tmin(i) .le. 270.16 .or. hdidx(i) .ne. 0) then
           hti =  1.0
           if (hdidx(i) .ge. hti) then   ! done with phase 1
             if(tcrown .le. tbase + 0.0) then
               hdidx(i) = hdidx(i) + 0.083
               if (hdidx(i) .gt. hti * 2.0) hdidx(i) = hti * 2.0  
             endif
c
             if (tmax(i) .ge. tbase + 283.16) then
               hdidx(i) = hdidx(i) - 0.02 * (tmax(i) - 283.16)
               if (hdidx(i) .gt. hti) hdidx(i) = hdidx(i) - 0.02 * (tmax(i) - 283.16) 
               hdidx(i) = max(0.0, hdidx(i))
             endif
c
           elseif (tcrown .ge. tbase - 1.0) then
             if (tcrown .le. tbase + 8.0) then
               hdidx(i) = hdidx(i) + 0.1 - ((tcrown - tbase + 3.5)**2 / 506.0)
               if (hdidx(i) .ge. hti .and. tcrown .le. tbase + 0.0) then
                 hdidx(i) = hdidx(i) + 0.083
                 if (hdidx(i) .gt. hti * 2.0) hdidx(i) = hti * 2.0
               endif
             endif 
c
             if (tmax(i) .ge. tbase + 283.16) then
               hdidx(i) = hdidx(i) - 0.02 * (tmax(i) - 283.16)
               if (hdidx(i) .gt. hti) hdidx(i) = hdidx(i) - 0.02 * (tmax(i) - 283.16) 
               hdidx(i) = max(0.0, hdidx(i))
             endif
           endif
c      
c calculate what the wheat killing temperature 
c there is a linear inverse relationship between 
c hardening of the plant and the killing temperature or
c threshold that the plant can withstand 
c when plant is fully-hardened (hdidx = 2), the killing threshold is -18 C
c
c will have to develop some type of relationship that reduces LAI and
c biomass pools in response to cold damaged crop 
c
           if (tmin(i) .le. 267.16) then
             tkil = (tbase - 6.0) - 6.0 * hdidx(i)         
             if (tkil .ge. tcrown) then
               if ((0.95 - 0.02 * (tcrown - tkil)**2) .ge. 0.02) then
                 write (*,*)  'crop damaged by cold temperatures'
               else
                 croplive(i,15) = 0.0   ! kill winter wheat
c WJS (05.16.10): Should anything else be done here? For example, should harvdate be set? Should croppresent be set to 0?
                 write (*,*)  '95% of crop killed by cold temperatures'
               endif
             endif
           endif
          endif
c
        else if (croplive(i,13) .eq. 1 .or. croplive(i,14) .eq. 1) then
          vf(i) = 1.0
c
        endif  ! croplive, vf check
c
 100   continue 
c
       return
       end
c
c ---------------------------------------------------------------------
      subroutine phenocrop(iyear, iyear0,imonth, iday, jday)
c ---------------------------------------------------------------------
c 
c subroutine which determine the phenological development of each crop
c type including allocation changes of photosynthate to 
c leaf, stem, root, and reproductive organs (pod, seed, grain) 
c leaf area index increase and decline
c subsequent carbon in biomass pools (leaf,stem,root,wood,reproductive)  
c
c Conversion factors for carbon to dry matter are from Penning DeVries
c et al., 1983 and Penning DeVries et al., 1989
c 
c values are fraction of carbon in dry matter
c
c leaves: 0.459
c stem  : 0.494
c root  : 0.467
c
c maize cob and grain : 0.491
c soybean pod and seed: 0.527 
c
c all of the phenology changes are based on the total number of gdd needed
c to change to the next phase - based on fractions of the total gdd typical
c for  that region based on the April 1 - Sept 30 window of development  
c
c Phase 1: Planting to leaf emergence
c Phase 2: Leaf emergence to beginning of grain fill (general LAI accumulation)
c Phase 3: Grain fill to physiological maturity and subsequently harvest (LAI decline)  
c
c uses:
c
      use comgrid
      use compar
      use comatm
      use comsoi
      use comsum
      use comveg
      use comsno
      use comnitr
      use comcrop
      use compft
      use combcs,only: cmask
c
      implicit none
c
c
      real huileaf(npft),         ! heat unit index needed to attain leaf emergence after planting 
     >     huigrain(npft)         ! heat unit index needed to reach vegetative maturity 
*     >     laicons(npft),        ! constant used in lai senescence equation
*     >     allconsl(npft),       ! constant used in dynamic allocation equations for leaves (decline in allocation)
*     >     allconss(npft),       ! constant used in dynamic allocation equations for stem (decline in allocation) 
*     >     laimx(npft),          ! maximum lai designated for crops from EPIC and EPICphase  models
*     >     tkill(npft),          ! minimum daily temperature threshold used to kill crops
*     >     mxmat(npft),          ! maximum length of growing season to harvest/physiological maturity 
*     >     mxdgfi(npft),         ! maximum number of days needed to reach grain fill stage
*     >     mxgddgf(npft)         ! maximum gdd past initiation of grain fill to reach phys maturity/harvest
c
      real 
     >     laidecl(npoi,npft)     ! decline in leaf area for crop
*     >     tauroot(npft),        ! time constant for root turnover 
*     >     arooti(npft),         ! initial allocation to crop fine roots 
*     >     arootf(npft),         ! end of growing season final allocation to crop fine roots
*     >     aleaff(npft),         ! end of growing season final allocation to crop leaf area
*     >     astemf(npft),         ! end of growing season final allocation to crop stem biomass
*     >     fleafi(npft),         ! initial fraction of aboveground carbon allocation (stem/leaf) allocated to leaf 
*     >     fleaf(npft),          ! fraction of aboveground carbon allocation (stem/leaf) allocated to leaf 
*     >     declfact(npft)        ! factor helping to control LAI decline at end of season
c
      real pc,
     >     ddays,
     >     ddfac,
     >     tthreshold,
*    >     xminlai, 
     >     dtt,
     >     explf,
     >     tix,
     >     xn,
     >     plag,
     >     abiol,
     >     crmcorn,
     >     aspecla,
     >     grainfillfrac,
     >     lai_at_harv,
     >     grnfracdecline1,
     >     grnfracdecline2,
     >     grnfracdecline2_start,
     >     grnfracdecline_postmaturity
c
      integer
     >     iyear,
     >     iyear0,
     >     imonth, 
     >     iday,
     >     jday, 
     >     nplants,
     >     my_daystoharv,     ! number of days between maturity and harvest in this grid cell, for the current pft
     >     i, j, k, l, n  

      logical
     >  maturity_reached_today,
     >  harvest_reached_today

c 
c Set some parameters
c
c LAI at harvest: After maturity, LAI decline will be calculated to achieve this LAI at the time of harvest
      parameter(lai_at_harv = 1.0)
c
c Daily decline in grnfraccrop during first part of the grain-fill period and during the second part of the grain-fill period:
c WJS 05.26.10: These values of 0.0025 and 0.014 are based on Valentinuz & Tollenaar 2004, "Vertical Profile of Leaf Senescence
c during the Grain-Filling Period in Older and Newer Maize Hybrids" - specifically, using the more recent hybrids from that paper
      parameter(grnfracdecline1 = 0.0025)
      parameter(grnfracdecline2 = 0.014)
c
c Fraction of the grain-fill period (in a GDD sense) at which we switch from using grnfracdecline1 to grnfracdecline2. Note that
c Valentinuz & Tollenaar (2004) divided the grain-fill period into halves based on number of days, whereas here we use GDD. A value
c of 0.5 means that we start using grnfracdecline2 once we have accumulated half (or more) of the required GDD between the start of
c grain-fill and maturity.
      parameter (grnfracdecline2_start = 0.5)
c
c Daily decline in grnfraccrop during the post-maturity period. 
c WJS 05.26.10: This value of 0.04 was chosen so that, if green fraction is ~ 50% at maturity (which will be the case for a 60-day
c grain-fill period, using grnfracdecline1 = 0.0025, grnfracdecline2 = 0.014 and grnfracdecline2_start = 0.5), then it will take ~ 2
c weeks for the leaves to turn fully brown, which Chris Kucharik suggested as a reasonable amount of time.
      parameter (grnfracdecline_postmaturity = 0.04)
c
c phenology for additional leaf drop - if drought related or temperature related at
c end of growing season     
c
       ddays      = 7.0
       ddfac      = 1.0 / ddays
       tthreshold = 273.16
c
c number of corn plants per square meter
c this is only important if using the leaf expansion equations
c of Ritchie, that is temperature dependent.  Our standard procedure
c here however is to use the allocation of C to leaf (aleaf) and
c specific leaf area (specla) to accumulate LAI during the season
c
       nplants    = 7 
c
      call vernalization(iyear,imonth,iday,jday)
c
c begin global grid
c
      do 100 i = 1, npoi
       if(cmask(i) == 0 ) cycle
c
       if (icropsum(i) .gt. 0.0) then
c
        aplantn(i) = 0.0
c
c only want to look at top 1.5 meters for plant available inorganic n
c
        do 115  k = 1, nsoilay 
c
          aplantn(i) = aplantn(i) + smsoil(i,k) + smsoln(i,k)
c
 115    continue

c time constant for root longevity (in days)
c
*           tauroot(13) = 200  ! 0.5% day-1 turnover rate   
*           tauroot(14) = 200  ! 0.5% day-1 turnover rate  
*           tauroot(15) = 200  ! 0.5% day-1 turnover rate 
c
c constants to control rate of lai decline after peak
c
*           laicons(13) = 2.0      
*           laicons(14) = 5.0 
*           laicons(15) = 3.5 
c
c constants to control rate of leaf and stem allocation decline after shift to grain 
c lower number gives more rapid decline in allocation - typical range between 2-5 
c
*           allconsl(13) = 2.0      
*           allconss(13) = 5.0 
c
*           allconsl(14) = 5.0      
*           allconss(14) = 2.0 
c
*           allconsl(15) = 3.0      
*           allconss(15) = 1.0 
c
c maximum potential lai for crops according to EPIC logic
c and EPIC-phase parameters 
c
*           laimx(13)   = 6.0      ! (e.g. soybeans)
*           laimx(14)   = 5.0      ! (e.g. maize) 
*           laimx(15)   = 7.0      ! (e.g. wheat) 
c
c crop phenology (gdd thresholds) controlled by gdd needed for maturity (physiological) which
c is based on the average gdd accumulation and hybrids in United States from April 1 - Sept 30  
c
c Phase 1:
c ==========
c threshold for attaining leaf emergence (based on fraction of gdd(i) -- climatological average)
c Hayhoe and Dwyer, 1990, Can. J. Soil Sci 70:493-497
c Carlson and Gage, 1989, Agric. For. Met., 45: 313-324 
c J.T. Ritchie, 1991: Modeling Plant and Soil systems
c
           huileaf(13)  = lfemerg(13)  * gddmaturity(i,13) 
           huileaf(14)  = lfemerg(14)  * gddmaturity(i,14) 
           huileaf(15)  = lfemerg(15)  * gddmaturity(i,15)  ! typically between 3-7% in wheat 
c
c Phase 2:
c ==========
c from leaf emergence to beginning of grain-fill period 
c this hypothetically occurs at the end of tassling, not the beginning 
c tassel initiation typically begins at 0.5-0.55 * gddmaturity 
c
c calculate linear relationship between huigrain fraction and relative
c maturity rating for maize
c
           if (management_prescribed) then
c Set huigrain using the prescribed, read-in xingrnfill variable; note
c that there is a single xingrnfill variable that applies to all crop
c types; in practice this means that you need to rename files depending
c on the crop type used in a given run, so that the proper input file is
c read for that crop. This also means that crop rotations cannot be used
c in conjunction with management_prescribed=TRUE.
c
c Note that for prescribed management, the CRM-based correction (as is
c done in the 'else' clause) doesn't need to be done: This is implicitly
c taken into account in creating the prescribed management files.
             huigrain(13) = xingrnfill(i) * gddmaturity(i,13)
             huigrain(14) = xingrnfill(i) * gddmaturity(i,14)
             huigrain(15) = xingrnfill(i) * gddmaturity(i,15)
           else
c Set huigrain using fixed parameters across space (and do an adjustment
c for corn using crm rating). This isn't really any more "prognostic"
c than the approach used if management_prescribed=true, so a better way
c to think about the difference (i.e., rather than prescribed vs.
c prognostic) is that with management_prescribed we use grid cell by
c grid cell prescribed management, whereas with
c management_prescribed=false we set the huigrain variable (and other
c variables, in the planting subroutine) without resorting to grid cell
c by grid cell prescriptions.
             crmcorn      = max(73., min((gddmaturity(i,14)+ 53.683)/13.882,135.))
             huigrain(14) = -0.002  * (crmcorn - 73.) + grnfill(14)
             huigrain(14) = min(max(huigrain(14),grnfill(14) - 0.1), grnfill(14)) 
             huigrain(14) = huigrain(14)   * gddmaturity(i,14) ! from Cabelguenne et al. 1999
c     
             huigrain(13) = grnfill(13)    * gddmaturity(i,13) ! from Cabelguenne et al. 1999 
c            huigrain(14) = grnfill(14)    * gddmaturity(i,14)  ! from Cabelguenne et al. 1999
             huigrain(15) = grnwht(iwheat_index) * gddmaturity(i,15)             
           end if
c
c set initial and final root npp allocation for crops 
c arooti is initial allocation (leaf emergence);
c arootf is final   allocation (crop maturity - harvest);
c
*            arooti(13) = 0.50    ! soybean (from Lal et al. 1999) 
*            arooti(14) = 0.45    ! maize 
c            arooti(15) = 0.40    ! wheat
*            arooti(15) = 0.30    ! wheat
c     
*            arootf(13) = 0.20    ! soybean (from Lal et al. 1999)
*            arootf(14) = 0.05    ! maize 
c            arootf(15) = 0.05    ! wheat
*            arootf(15) = 0.00    ! wheat
c
c set final leaf and stem allocation
c aleaff is ending allocation to leaves after beginning of reproductive stage
c astemf is ending allocation to stem after beginning of reproductive stage
c
*            aleaff(13) = 0.00    ! soybean
*            aleaff(14) = 0.00    ! maize 
*            aleaff(15) = 0.00    ! wheat
c
*            astemf(13) = 0.30    ! soybean
*            astemf(14) = 0.00    ! maize 
*            astemf(15) = 0.05    ! wheat
c
c lai decline factor
c
*            declfact(13)   = 1.05
*            declfact(14)   = 1.05
*            declfact(15)   = 1.05
c
c
c set fraction of aboveground allocation that goes to leaf vs. stem
c if lai has not reached a maximum value allowed (before reproductive phenology stage) 
c
*            fleafi(13)   = 0.575             ! soybean 
*            fleafi(14)   = 0.575             ! maize 
             fleafi(15)   = fleafiwht(iwheat_index)  ! wheat          
c
*            xminlai = 0.01
c
c freeze kill temperature threshold (K)
c -5 celsius 
c
*            tkill(13)  = 268.16
*            tkill(14)  = 268.16
*            tkill(15)  = 268.16
c
c maximum number of days before automatic switch to grain fill
c initiation.  
c
*            mxdgfi(13) = 100.0
*            mxdgfi(14) = 110.0
*            mxdgfi(15) = 170.0
c
c max allowable degree days past grain fill initiation
c
*            mxgddgf(13)  = 850.0
*            mxgddgf(14)  = 925.0
             mxgddgf(15)  = mgddgf(iwheat_index)  ! wheat
            
c
c max days from planting to physiological maturity
c will be different for spring versus winter wheat
c 
*            mxmat(13)    = 150.0 
*            mxmat(14)    = 165.0
             mxmat(15)    = mxmatwht(iwheat_index)   ! wheat
             
c
c crop phenology
c daily routine  
c
c allocation rules for crops based on maturity and linear decrease of amount allocated
c to roots over course of the growing season
c
         do 80 j = scpft, ecpft 
c
           if (croppresent(i,j) .eq. 1) then
c
c calculate fraction allocated to leaf (from J. Norman allocation curve)
c bfact and fleafi are set in params.crp
c
               fleaf(j) = fleafi(j) * (exp(-bfact(j)) - exp(-bfact(j) *
     >                               gddplant(i,j)/huigrain(j))) / (exp(-bfact(j))-1) 
c
c calculate accumulated growing degree days since planting (gddplant) 
c determine if growing degree days calculated from top layer soil temperature
c are enough for leaf emergence to occur 
c
             hui(i,j)       = gddplant(i,j)        
             leafout(i,j)   = gddtsoi(i,j)  
c
             if (management_prescribed) then
c use values specified by grid cell
c
c Note that there is a single xindaystoharv variable that applies to all
c crop types; in practice this means that you need to rename files
c depending on the crop type used in a given run, so that the proper
c input file is read for that crop. This also means that crop rotations
c cannot be used in conjunction with management_prescribed=TRUE.
               my_daystoharv = nint(xindaystoharv(i))
             else
c use a constant value for each pft
               my_daystoharv = nint(daystoharv(j))
             end if
c
             laidecl(i,j) = 0.0   
c
c calculate days past planting
c
c             if (iwheat .eq. 2 .and. croppresent(i,j) .eq. 1) then
c               idpp(i,j)  = jday + (ndaypy - idop(i,j))
c             else
c               idpp(i,j)  = jday - idop(i,j)
c             endif
c
                idpp(i,j) = idpp(i,j) + 1
c
c crop phenology from leaf emergence to start of leaf decline   
c determine allocation coefficients to help calculate increase in lai  
c and fine root biomass
c
             if (leafout(i,j) .ge. huileaf(j) .and.
     >           hui(i,j) .lt. huigrain(j)) then  
c
c Phase 1 completed:
c ==================
c if hui is less than the number of gdd needed for filling of grain
c leaf emergence also has to have taken place for lai changes to occur and carbon
c assimilation 
c
c calculate day of year that leaf emergence occurs; also, gdd at emergence
c note that emerge_day is initialized to 9999 at beginning of each growing season
c
               if (emerge_day(i,j) .gt. 366) then  ! leaf emergence has not occurred here yet
                 emerge_day(i,j) = real(jday)
                 emerge_gdd(i,j) = gddplant(i,j)
                 emerge_gddtsoi(i,j) = gddtsoi(i,j)
               end if
c
c allocation rules for crops based on maturity and linear decrease of amount allocated
c to roots over course of the growing season
c
               awood(i,j) = 0.0
c           
c check to see if lai is at maximum allowable 
c
               if (peaklai(i,j) .eq. 1) then 
                 aleaf(i,j) = 0.0
                 arepr(i,j) = 0.0
c CJK 9-23-04    astem(i,j) = 0.0
c                aroot(i,j) = 1. - arepr(i,j) - aleaf(i,j) - astem(i,j)
c **********
c
                 aroot(i,j) = min(1.0, (arooti(j) - (arooti(j) - arootf(j))
     >                        * min(1.0,hui(i,j)/gddmaturity(i,j)))) 
                 aroot(i,j) = max(0.0, aroot(i,j))
                 astem(i,j) = 1.0 - aroot(i,j) - aleaf(i,j) - arepr(i,j)
c
c                 
            
               else
                 aroot(i,j) = min(1.0, (arooti(j) - (arooti(j) - arootf(j))
     >                        * min(1.0,hui(i,j)/gddmaturity(i,j)))) 
                 aroot(i,j) = max(0.0, aroot(i,j))
                 aleaf(i,j) = max(0.0,(1.0 - aroot(i,j)) * fleaf(j)) 
                 astem(i,j) = 1.0 - aroot(i,j) - aleaf(i,j)
                 arepr(i,j) = 0.0
c
               endif
c
c calculate actual lai increase based on npp and allocation rules in ibis 
c only increment lai if it hasn't reached maximum allowable value 
c
               if (peaklai(i,j) .eq. 0) then
                 tlai(i,j) = plai(i,j) + (specla(j) * aleaf(i,j)
     >                        * max(0.0, adnpp(i,j)))

                 if (tlai(i,j) .ge. laimx(j)) then
                    aleaf(i,j) = min(1.0 - aroot(i,j), laimx(j) -
     >                        plai(i,j)) / (specla(j) * adnpp(i,j)) 
                    aleaf(i,j) = max(0.0, aleaf(i,j))
                    astem(i,j) = 1.0 - aroot(i,j) - aleaf(i,j) 
c
c CJK other possible source of over allocation
c above astem could be set to 0.0
c
                    peaklai(i,j) = 1
c
                 endif
c
c original phenology - apply to soybeans
c if maize - employ leaf expansion routine based on a modified
c ceres gdd approach
c
c                if     (j .eq. 13) then
c
                   plai(i,j)    = plai(i,j) + (specla(j) * aleaf(i,j)
     >                          * max(0.0, adnpp(i,j)))
c
c               elseif (j .eq. 14) then 
c
c 
c                 anpplai      = specla(j) * aleaf(i,j) * max(0., adnpp(i,j))
c
c==============================================================================
c addition of ceres-maize logic for leaf expansion and relationship to lai
c as a function of daily gdd accumulation
c==============================================================================
c
c calculate daily thermal time
c
c                   dtt = td(i) - baset(j)
c
c                   if (cumlvs(i,j) .lt. 5) then
c                     pc = 1.0
c                    explf = 73.5 
c                    explf = 63.0 
c                   else 
c                     pc = 1.0
c                    explf = 48.3 
c                     explf = 40.0
c                   endif
c
c                   cumlvs(i,j) = cumlvs(i,j) + max(0.0, dtt)/(explf * pc) 
c                   tix = max(0.0, dtt)/(explf * pc)
c
c calculate the number of the newest emerging leaf
c
c                   xn  = cumlvs(i,j) + 1.0
c
c ceres-maize and other models impose a stronger drop-off in
c leaf area accumulations for higher leaf numbers, but this
c requires the knowledge of the TLNO - or total number of
c leaves for the plant - which is something we don't or won't have
c in this model version
c
c                   if (xn .lt. 4) then
c                     plag = 5.55 * xn * tix 
c                   else if (xn .ge. 4 .and. xn .le. 13) then 
c                     plag = 0.261 * xn * xn * xn * tix 
c                   else
c                     plag = 611.5 * tix
c                   endif
c
c convert leaf area expanded (cm2) per plant to m2 basis
c
c                   plag = max(0.0, plag * nplants * 0.0001)
c
c                   abiol     = max(0.001, aleaf(i,j) * max(0.0, adnpp(i,j)))
c
c calculate specific leaf area just for added leaf area and biomass this time step
c
c                  aspecla   = plag / abiol
c
c impose limit on specific leaf area - limiting conditions for carbon assimilated
c
c                   if (aspecla .gt. 30.0) aspecla = 30.0
c
c                  plai(i,j) = plai(i,j) + abiol * aspecla 
c
c                  endif
c
c==============================================================================
c
               endif
c
c hold ending allocation values to stem and leaf for use by equations after shift to
c reproductive phenology stage begins
c
             astemi(i,j) = astem(i,j)
             aleafi(i,j) = aleaf(i,j)
c
c shift allocation either when enough gdd are accumulated or maximum number
c of days has elapsed since planting
c
c            else if ((hui(i,j) .ge. huigrain(j) .or. idpp(i,j) .ge. mxdgfi(j)) 
             else if (hui(i,j) .ge. huigrain(j) 
     >                 .and. croplive(i,j) .eq. 1.) then  ! only enter this block if we haven't reached maturity yet
c
              dpgf(i,j) = max(0.0, gddplant(i,j) - huigrain(j))             
c
c WJS (06.04.10): It would be more clear if grainday and graingdd were
c set inside the same conditional, as is done for emerge_day and
c emerge_gdd above, since grainday and graingdd end up getting set at
c the same time anyway - i.e., the first day of the grainfill period.
c
c calculate day of year that grain fill begins
c grainday is initialized to 9999 at beginning of each growing season
c
              grainday(i,j) = min(grainday(i,j), real(jday))
c
c calculate gdd at the time grain fill begins
c graingdd is initialized to a value < 0 at beginning of each growing season
c
              if (graingdd(i,j) .lt. 0) then
                graingdd(i,j) = hui(i,j)
              end if
c
c Phase 2 completed:
c
               awood(i,j) = 0.0
               aroot(i,j) = min(1.0, (arooti(j) - (arooti(j) - arootf(j))
     >                    * min(1.0, hui(i,j)/gddmaturity(i,j)))) 
               aroot(i,j) = max(0.0, aroot(i,j))
c
               if (thrlai(i,j) .lt. 0.0001) thrlai(i,j) = plai(i,j)

c
c lai plateau -- reached threshold -- according to heat unit index  (accumulated GDD logic)        
c set lai and hui threshold  (to be used in leaf phenology (senescence) calculations)
c shift aboveground allocation to grain/fruit
c
c
               templai(i,j)   = plai(i,j) 
c
c lai decline based on thermal time past grain fill  
c add check in to make sure huigrain is <= gddplant in this equation.  If it
c isn't because the days past planting is used to shift grain allocation, set
c huigrain value to the gddplant value when the shift took place so LAI starts
c the proper decline.
c 
c lai decline based thermal time accumulation past grain fill 
c add check in to make sure huigrain is <= gddplant in this equation.  If it
c isn't because the days past planting is used to shift grain allocation, set
c huigrain value to the gddplant value when the shift took place so LAI starts
c the proper decline.
c 
                plai(i,j) = max(thrlai(i,j) * (1.0-min(max((gddplant(i,j)-
     >                      huigrain(j)),0.0)/(0.55*gddmaturity(i,j)), 1.0)
     >                      **laicons(j)), xminlai)
c
c calculate decrease in lai for purpose of updating aboveground biomass pools
c
                 laidecl(i,j)   = max(0.0, templai(i,j) - plai(i,j)) 
                 if (astemi(i,j) .gt. astemf(j)) then
c
                   astem(i,j) = max(astem(i,j) * (1.0-min((hui(i,j)-
     >                       huigrain(j))/((gddmaturity(i,j)*declfact(j))
     >                      -huigrain(j)),1.0)**
     >                       allconss(j)),astemf(j)) 
                   astem(i,j) = max(0.0,astem(i,j)) 
                 endif
c
c CJK 9-23-04
c ***********
c
                 astem(i,j) = 0.0  ! CJK 9-23-04
                 aroot(i,j) = 0.0  ! CJK 9-23-04 
c
                 if (aleafi(i,j) .gt. aleaff(j)) then
c
                    aleaf(i,j) = max(aleaf(i,j) * (1.0-min((hui(i,j)-
     >                       huigrain(j))/((gddmaturity(i,j)*declfact(j))
     >                      -huigrain(j)),1.0)**
     >                       allconsl(j)),aleaff(j)) 
                    aleaf(i,j) = max(0.0,aleaf(i,j)) 
c
                 endif
                 arepr(i,j)   = 1.0 - aroot(i,j) - astem(i,j) - aleaf(i,j) - awood(i,j)
c
c Handle decline in green fraction during the grain-fill period
c
c fraction of the grainfill period that we have completed
                 grainfillfrac = (hui(i,j) - huigrain(j)) / (gddmaturity(i,j) - huigrain(j))

                 if (grainfillfrac .lt. grnfracdecline2_start) then
                   grnfraccrop(i,j) = max(0.0, grnfraccrop(i,j) - grnfracdecline1)
                 else
                   grnfraccrop(i,j) = max(0.0, grnfraccrop(i,j) - grnfracdecline2)
                 end if

               else if (croplive(i,j) .eq. 0) then
c The crop is present, but no longer living. This will be true between maturity and harvest                 

                 dpmature(i,j) = dpmature(i,j) + 1

c We expect NPP to be <= 0 after maturity; print a warning if this is not the case
                 if (adnpp(i,j) .gt. 0) then
                   print *,
     >               'WARNING: NPP .ne. 0 post-maturity for point ',
     >               i, ', pft ', j
                 end if

c The allocation fractions shouldn't matter, since NPP should be <= 0.
c However, I am ensuring that the sum of the allocation fractions is 1,
c just in case they do matter. I am arbitrarily setting astem(i,j) to 1,
c since allocation to the stem seems like it would have the smallest
c impact on the model's results.
                 awood(i,j) = 0.0
                 aleaf(i,j) = 0.0
                 arepr(i,j) = 0.0
                 aroot(i,j) = 0.0
                 astem(i,j) = 1.0

c Handle LAI decline between maturity and harvest
                 templai(i,j) = plai(i,j)
                 plai(i,j) = max(plai(i,j) - lai_decl_per_day(i,j), xminlai)
                 laidecl(i,j) = max(0.0, templai(i,j) - plai(i,j))

c Handle decline in green fraction between maturity and harvest
                 grnfraccrop(i,j) = max(0.0, grnfraccrop(i,j) - grnfracdecline_postmaturity)


               else
c WJS (06.03.10): I am adding this 'else' clause to catch situations
c where we are not in any of the above phenology periods. I believe this
c will just be the case before leaf emergence.
c
c WJS (06.03.10): The allocation fractions shouldn't matter very much,
c since NPP should be near zero (it might not be exactly 0, since plai
c is initialized to 0.01). However, I am ensuring that the sum of the
c allocation fractions is 1. I am arbitrarily setting astem(i,j) to 1,
c since allocation to the stem seems like it would have the smallest
c impact on the model's results, and it kind of makes sense in terms of
c what might happen during the pre-emergence period.
                 awood(i,j) = 0.0
                 aleaf(i,j) = 0.0
                 arepr(i,j) = 0.0
                 aroot(i,j) = 0.0
                 astem(i,j) = 1.0

               endif
c
c keep track of total biomass production for the entire year, and the
c aboveground value to calculate harvest index
c
               aybprod(i,j) = aybprod(i,j) + 
     >                       aleaf(i,j) * max(0.0,adnpp(i,j)) +            
     >                       astem(i,j) * max(0.0,adnpp(i,j)) +            
     >                       aroot(i,j) * max(0.0,adnpp(i,j)) +            
     >                       awood(i,j) * max(0.0,adnpp(i,j)) +            
     >                       arepr(i,j) * max(0.0,adnpp(i,j))             
c
               ayabprod(i,j) = ayabprod(i,j) + 
     >                       aleaf(i,j) * max(0.0,adnpp(i,j)) +            
     >                       astem(i,j) * max(0.0,adnpp(i,j)) +            
     >                       arepr(i,j) * max(0.0,adnpp(i,j))             
c
c keep track of annual total root production carbon
c
               ayrprod(i,j) = ayrprod(i,j) +
     >                        aroot(i,j) * max(0.0,adnpp(i,j))
c
c update carbon reservoirs using an analytical solution
c to the original carbon balance differential equation
c 
c keep track of total carbon allocated to leaves for litterfall
c calculation
c
               aylprod(i,j) = aylprod(i,j) + aleaf(i,j) *
     >                        max (0.0, adnpp(i,j))
c
               cbiol(i,j) = cbiol(i,j) -
     >                      (laidecl(i,j)/specla(j)) + 
     >                      aleaf(i,j) * max (0.0, adnpp(i,j)) 
c
               cbiog(i,j) = cbiog(i,j) + arepr(i,j) * max (0.0, adnpp(i,j))
               cbiog(i,j) = max(0.0, cbiog(i,j)) 
c
               cbios(i,j) = cbios(i,j) + astem(i,j) * max (0.0, adnpp(i,j))
               cbios(i,j) = max(0.0, cbios(i,j))
c
               cbior(i,j) = cbior(i,j) * exp(-1.0/tauroot(j)) + 
     >                      aroot(i,j) * tauroot(j) * max(0.0,adnpp(i,j)) *
     >                      (1.0 - exp(-1.0/tauroot(j)))

               cbior(i,j) = max(0.0, cbior(i,j)) 
c
               cbiow(i,j) = 0.0 
c
               biomass(i,j) = cbiol(i,j) + cbiog(i,j) + cbior(i,j)
     >                      + cbios(i,j) + cbiow(i,j)
c
c keep track of aboveground annual npp 
c
               ayanpp(i,j)  = (aleaf(i,j) + arepr(i,j) + astem(i,j)
     >                         + awood(i,j)) * adnpp(i,j) + ayanpp(i,j)
c
c---------------------------------------------------------------------------------
c check for climatic and phenological limits on maturity, growth, and harvest date

               maturity_reached_today = .FALSE.
               harvest_reached_today = .FALSE.
c
c check to see if minimum temperature has fallen below freeze
c kill threshold for 3 consecutive days and if lai is above a minimum, plant will
c be damaged/killed.  This function is more for spring freeze events
c or for early fall freeze events
c
c currently simulates too many grid cells that are killed by
c freezing temperatures 
c
c spring wheat is affected by this, winter wheat kill function
c is determined in crops.f - is a more elaborate function of
c cold hardening of the plant
c
               if (tmin(i) .le. tkill(j)) then
                 ccdays(i,j) = ccdays(i,j) + 1
               else
                 ccdays(i,j) = 0
               endif
c===============================================================
c removed on March 12 2002 - C. Kucharik
c until it can be a bit more refined, or used at a smaller scale.
c we really have no way of validating this routine
c too difficult to implement on 0.5 degree scale grid cells 
c
                if (ccdays(i,j) .ge. 3    .and. 
     >             hui(i,j)    .ge. 0.6*gddmaturity(i,j)   .and.
c     >             plai(i,j)   .ge. 0.25 .and. 
c     >             imonth      .ge. 7    .and.
c     >             iday        .ge. 1    .and.
     >             croplive(i,j) .eq. 1  .and.
     >             ((j .eq. 13 .or. j .eq. 14) .or.
     >               j .eq. 15 .and.   iwheat  .eq. 1)) then

c Freeze-kill is like reaching maturity prematurely                  
                     maturity_reached_today = .true.

c Record the fraction of the season (in a gdd sense) that the crop "missed" due to freeze kill
                     if (hui(i,j) .lt. gddmaturity(i,j)) then
                       fkill_gdd_missed(i,j) = (gddmaturity(i,j) - hui(i,j))/gddmaturity(i,j)
                     end if
               endif
c
c ================================================================
c if accumulated gdd past grain fill initiation exceeds limit
c or number of days past planting reaches a maximum, the crop has 
c reached physiological maturity
c
c crop could either be live or dead at this stage - these limits
c could lead to reaching physiological maturity or determining
c a harvest date for a crop that was killed by an early frost (see above)
c
c orig CJK Jan 27 03 if (dpgf(i,j) .ge. mxgddgf(j) .or.
c    >              idpp(i,j) .ge. mxmat(j)) then
c
                if (croplive(i,j) .eq. 1 .and.
     >              (hui(i,j) .ge. gddmaturity(i,j) .or.
     >              idpp(i,j) .ge. mxmat(j))) then
                  maturity_reached_today = .true.
                endif

c
c ================================================================
c if it's Dec. 31, then we force maturity for spring crops, in order to
c avoid problems that arise because a bunch of initialization happens on
c Jan. 1
c
                if (croplive(i,j) .eq. 1 .and.
     >            (j .eq. 13 .or. j .eq. 14 .or. 
     >            (j .eq. 15 .and. iwheat .eq. 1)) .and. 
     >            (imonth .eq. 12 .and. iday .eq. 31)) then
                  maturity_reached_today = .true.

                  print *, '*** Forcing Dec. 31 maturity for ', i, j
                end if

c ================================================================
c if one of the above checks led us to reach maturity today, handle that
c
                if (maturity_reached_today) then
                  croplive(i,j)     = 0.0
                  matdate(i,j) = jday 
                  matgdd(i,j) = hui(i,j)
                  matgrnfraccrop(i,j) = grnfraccrop(i,j)
                  dev_fraction_attained(i,j) = hui(i,j) / gddmaturity(i,j)
                  dpmature(i,j) = 0.0

c calculate LAI decline per day required to achieve LAI=lai_at_harv at harvest
                  if (plai(i,j) .gt. lai_at_harv .and. my_daystoharv .ge. 1) then
                    lai_decl_per_day(i,j) = (plai(i,j) - lai_at_harv)/my_daystoharv
                  else  ! either LAI < lai_at_harv already, or harvest will happen immediately
                    lai_decl_per_day(i,j) = 0.
                  end if
                end if

c ================================================================
c check whether we reach harvest today
c
                if (croplive(i,j) .eq. 0 .and. 
     >            nint(dpmature(i,j)) .eq. my_daystoharv) then

                  harvest_reached_today = .true.
                end if

c
c ================================================================
c if it's Dec. 31, then we force harvest for spring crops, in order to
c avoid problems that arise because a bunch of initialization happens on
c Jan. 1
c 
c Note that we only do this if croplive = 0 - i.e., we never harvest a
c live crop. However, a check done above should ensure that a crop is
c always forced to reach maturity by Dec. 31.
c
                if (croplive(i,j) .eq. 0 .and.
     >            (j .eq. 13 .or. j .eq. 14 .or. 
     >            (j .eq. 15 .and. iwheat .eq. 1)) .and. 
     >            (imonth .eq. 12 .and. iday .eq. 31)) then
                  harvest_reached_today = .true.

                  print *, '*** Forcing Dec. 31 harvest for ', i, j
                end if

c ================================================================
c handle reaching harvest today
c
                if (harvest_reached_today) then
                  croppresent(i,j)  = 0.0
                  harvdate(i,j)     = jday
                  grnfraccrop(i,j)  = 0.0 ! turn all vegetation to brown
                  plai(i,j)         = 0.25 ! simulates remaining stubble/mulch

                  !plant c4/c3 grasses as cover crops after corn or soybean harvest (Y.Li)
                  if(flg_coverCrops /= 0 .and. (j==13 .or. j==14) ) then
                    exist(i,:)  = 0

                    if(flg_coverCrops == 1) then        !c4 grasses
                      exist(i,11) = 1
                    else if (flg_coverCrops == 2) then  !c3 grasses
                      exist(i,12) = 1    
                    end if

                  end if  !on if flg_coverCrops
                  
                end if

c---------------------------------------------------------------------------------
           endif ! croppresent
          
c
 80     continue
c
       endif  ! crop existence (cropsum)
 100    continue
c
c call to update crop lai, fractions, etc.
c
        call cropupdate
c
        call cropresidue(iyear,iyear0,jday) 
c
        if (iday .eq. 31 .and. imonth .eq. 12) then 
          call cropoutput(iday,imonth,iyear,jday) 
        endif 
c
c return to main program
c 
      return
      end
c
c ---------------------------------------------------------------------
      subroutine cropupdate
c ---------------------------------------------------------------------
c
c uses:
c
      use comgrid
      use compar
      use comatm
      use comsoi
      use comsum
      use comveg
      use comsno
      use comcrop
      use comnitr
      use compft
      use combcs,only: cmask
c
      implicit none
c
*     real  laimx(npft)
c
      real  
*     >      xminlai,
*     >      ztopmxsoy,   ! maximum height of soybean canopy
*     >      ztopmxwht,   ! maximum height of wheat canopy
*     >      ztopmxmze,   ! maximum height of maize canopy
     >      avglail
c
      integer
     > i, j
c 
*      xminlai    = 0.010
c
c maximum crop height values (from EPIC and CERES models) 
c
*      ztopmxsoy   = 0.75  ! maximum height (m) of soybeans
*      ztopmxwht   = 1.2   ! maximum height (m) of wheat 
*      ztopmxmze   = 2.5   ! maximum height (m) of maize
c     
c maximum allowable lai
c 
*      laimx(13)   = 6.0      ! maximum lai of soybeans
*      laimx(14)   = 5.0      ! maximum lai of maize 
*      laimx(15)   = 7.0      ! maximum lai of wheat 
c
c begin global grid
c
      do 100 i = 1, npoi
        if(cmask(i) == 0 ) cycle

       if (icropsum(i) .gt. 0.0) then
c
c maintain minimum value of leaf carbon in crops that exist
c
        do 90 j = scpft, ecpft 
         if (exist(i,j) .eq. 1.0 .and. croppresent(i,j) .eq. 1) then
c 
          cbiol(i,j)   = max((exist(i,j) * xminlai/specla(j)),
     >                   cbiol(i,j))
c         plai(i,j)    = cbiol(i,j) * specla(j) 
          biomass(i,j) = cbiol(i,j) + cbiog(i,j) + cbior(i,j)
     >                   + cbiow(i,j)
c
c check for maximum plant leaf area index
c
          if (plai(i,j) .gt. plaimx(i,j)) plaimx(i,j) = plai(i,j)
         endif
c
 90   continue
c
c crop canopy single sided leaf area index (area-weighted)
c
          avglail = plai(i,13) + plai(i,14) + plai(i,15)           
c
c crop canopy fractions
c
        frac(i,13)  = plai(i,13)  /
     >               max (avglail, epsilon)
c
        frac(i,14)  = plai(i,14)  /
     >               max (avglail, epsilon)
c
        frac(i,15)  = plai(i,15)  /
     >               max (avglail, epsilon)
c
c calculate total crop leaf are index
c
        totlail(i) = plai(i,13) + plai(i,14) + plai(i,15)
c
c
        fl(i) = totlail(i) / 1.0
c
c constrain fractional cover
c
c       fl(i) = max(0.25, min(0.975, fl(i)))
        fl(i) = max(0.025, min(0.975, fl(i)))
c
c calculate the crop canopy leaf area index using the fractional vegetation cover
c
        lai(i,1) = avglail / fl(i)
c
c C. Kucharik  04.02.01
c calculate greenness fraction of crop canopy
c if plant optical properties were ever changed - due to browning of
c vegetation in the future on a daily basis, do it here...this will effect leaf nir and vis
c transmittance and reflectance in radiation.f
c greenfracl(i) is lower canopy total green fraction - which is also being calculated
c for natural vegetation in vegetation.f 
c
       greenfracl(i) = 0.0
       do 80 j = scpft, ecpft
        greenfracl(i) = greenfracl(i) + frac(i,j) * grnfraccrop(i,j)  
 80    continue
c
c calculate total crop canopy biomass
c
        totbiol(i) = biomass(i,13) + biomass(i,14) + biomass(i,15) 
c
c calculate crop bottom/top height
c ztop is calculated in phenology subroutine based leaf area index 
c is the weighted mean of plai of each crop pft in the grid cell
c will only be important if we allow more than one crop pft to exist
c within the same grid cell
c
        zbot(i,1) = 0.02
c
        if(lai(i,1) /= 0) then !Y.Li
          ztop(i,1) = plai(i,13)/lai(i,1) * ztopmxsoy *
     >                (min(plai(i,13)/(laimx(13)-1.0),1.0))**2 +
     >                plai(i,14)/lai(i,1) * ztopmxmze *
     >                (min(plai(i,14)/(laimx(14)-1.5),1.0))**2 +                              
     >                plai(i,15)/lai(i,1) * ztopmxwht *
     >                (min(plai(i,15)/(laimx(15)-1.0),1.0))**2                              
        end if !on if(lai(i,1)/=0)
c
        htmx(i,1) = max(htmx(i,1),ztop(i,1))
        ztop(i,1) = max(0.05, max(htmx(i,1),ztop(i,1)))
c
c calculate stem area index for crops 
c
        sai(i,1) = 0.20 * plai(i,13) + 0.10 * plai(i,14) + 0.20 * plai(i,15) 
c
c calculate annual aboveground npp total
c
        ayanpptot(i) = ayanpp(i,13) + ayanpp(i,14) + ayanpp(i,15)
c
c end of loop
c
        endif  ! existence
 100  continue
c
c return to main program
c 
      return
      end
c
c ---------------------------------------------------------------------
      subroutine cropresidue(iyear,iyear0,jday)
c ---------------------------------------------------------------------
c
c routine that calculates the amount of residue that is input back to
c the soil at the end of growing season from crop at harvest
c
c also keeps track of effects of residue and plant nitrogen uptake on
c the nitrogen budget
c
c
c uses:
c
      use comgrid
      use compar
      use comatm
      use comsoi
      use comsum
      use comveg
      use comcrop
      use comnitr
      use combcs,only: cmask
c
c
      real     
*     >         maxhi(npft),        ! Maximum harvest index allowed for crops
*     >         fyield(npft),       ! Adjustment factor for fraction of cbiog that
                                   ! is just grain - measured in the field 
*     >         cgrain(npft),       ! carbon fraction in grain
*     >         convfact(npft),     ! conversion factor to go to bu/acre 
     >         dumg,
     >         rdm,
     >         fdml,
     >         fdms,
     >         fnresidue
c
      integer
     >         i, j,
     >         iyear,
     >         iyear0,
     >         jday
c
c carbon in grain and conversion factors 
c
*      cgrain(13)   = 0.45
*      cgrain(14)   = 0.45
*      cgrain(15)   = 0.45
c
*      convfact(13) = 150.0
*      convfact(14) = 159.46
*      convfact(15) = 150.0
c
c set maximum allowable harvest indicies
c
*      maxhi(13)  = 0.38 ! soybean (from Cabelguenne et al. 1999
c
*      maxhi(14)  = 0.60 ! maize  (from EPIC and EPICphase  models) 
c
*      maxhi(15)  = 0.50 ! wheat  (from EPIC, EPIC phase, and CERES-wheat)  
c
c fraction of total alloction in reproductive stage (cbiog) that
c is actual harvested grain yield
c
*      fyield(13) = 0.85 ! from Lal et al. 1999 
*                        ! Agr. For. Met. 93:  53-70.
*      fyield(14) = 1.00
*      fyield(15) = 0.85 ! major storage organ is the ear, in which 85% is the grain
*                        ! Penning deVries, 1989 
c
c begin global grid
c
      do 100 i = 1, npoi
        if(cmask(i) == 0 ) cycle
c
        if (icropsum(i) .gt. 0.0) then
c
c zero out litter fall rates
c
        falll(i) = 0.0
        fallr(i) = 0.0
        fallw(i) = 0.0      
c
        do 90 j = scpft, ecpft
c
c calculate CRM values from Pioneer regression relationships 
c Jan 28 03 CJK
c
          crmclim(i,j)    = max(73., min((gddmaturity(i,j)+53.683)/13.882,135.))
c         crmact(i,14)    = max(73., min((gdd10this(i)+53.683)/13.882,135.))
c
c crmact should designate what CRM rating resulted from GDD accumulation between
c killing frost/freeze events (-2.2 C)  
c
          crmact(i,14)    = max(73., min((gddfzcorn(i)+53.683)/13.882,135.))
c
c gddplant gets reinitialized to 0.0 at harvest date, so save value here
c Note that we only stop updating crmplant when croplive is no longer 1
c (i.e., we don't include the gdd that are accumulated between maturity and harvest)
c
          if (gddplant(i,j) .gt. 0.0 .and. croplive(i,j) .eq. 1.) then
            crmplant(i,j)   = max(73., min((gddplant(i,j)+53.683)/13.882,135.))
          endif
c
c only write out values at harvest date, and re-initialize crop variables
c at this time - this allows for the same crop (e.g., wheat) to be grown
c across two consecutive calendar years   
c
          if (exist(i,j) .eq. 1.0 .and. harvdate(i,j) .eq. jday) then
c
            hdate(i,j) = harvdate(i,j)
            pdate(i,j) = idop(i,j)
c
            ayabprod(i,j) = max(cbiol(i,j), ayabprod(i,j))
c
c added so division by zero doesn't take place in those places where
c production was zero 
c
c adjust actual grain yield value
c
            dumg         = cbiog(i,j) 
            cbiog(i,j)   = cbiog(i,j) * fyield(j) 
c 
c add excess (i.e. pod of soybeans) to stem storage pool of plant 
c
            cbios(i,j)   = max(0.0, cbios(i,j) + (dumg - cbiog(i,j)))
c
c impose upper limit on harvest index as per JMN suggestion 
c if harvest index is > allowable, put excess to stem which will add to more
c litterfall and adjust the grain accordingly
c might have to revisit logic with respect to what actually get harvested, versus
c what is left in the field after harvest as litter input to the soil 
c 
            harvidx(i,j)   =  cbiog(i,j) / ayabprod(i,j)
c
            croplaimx(i,j) = plaimx(i,j)
c
            if (harvidx(i,j) .gt. maxhi(j)) then
c
              harvidx(i,j) = maxhi(j) 
              dumg         = cbiog(i,j) 
              cbiog(i,j)   = maxhi(j) * ayabprod(i,j)
c 
c add excess to stem of plant
c
              cbios(i,j)   = cbios(i,j) + (dumg - cbiog(i,j))
c
            endif
c
c
c calculate n in grain (kg/ha)
c
            grainn(i,j)   = (cbiog(i,j) / cfrac(j)) * fngrain(i,j) * 1e+04
c
c calculate crop yield in bu/acre
c
            cropyld(i,j) = cbiog(i,j) * convfact(j) / cgrain(j) 
c
c yield in Mg/ha dry matter
c
            dmyield(i,j) = cbiog(i,j) / cgrain(j) * 10.0
c
c calculate stem and leaf dry matter (Mg/ha)  
c
            dmleaf(i,j)    = aylprod(i,j) / cfrac(j) * 10.0
            dmstem(i,j)    = cbios(i,j)   / cfrac(j) * 10.0
            dmroot(i,j)    = cbior(i,j)   / cfrac(j) * 10.0
            dmresidue(i,j) = dmleaf(i,j)  + dmstem(i,j)
c
c calculate dry matter production of crop total 
c
            dmcrop(i,j)    = dmresidue(i,j) + dmroot(i,j) +
     >                       dmyield(i,j)
c
c calculate aboveground residue dry matter total
c
            rdm  = dmleaf(i,j) + dmstem(i,j) 
c
c calculate fractions for leaf and stem
c
            if(rdm /= 0.0) then !Y.Li
              fdml = dmleaf(i,j) / rdm 
              fdms = dmstem(i,j) / rdm
            else
              fdml = 0.0
              fdms = 0.0
            end if !on if(rdm /= 0.0)

            !fdml = dmleaf(i,j) / rdm 
            !fdms = dmstem(i,j) / rdm
c
            fnresidue      = fnleaf(i,j) * fdml + fnstem(i,j) * fdms  
c
c calculate amount of N in aboveground residue (kg/ha) 
c
            residuen(i,j)  = (fnleaf(i,j) * dmleaf(i,j) 
     >                     + fnstem(i,j) * dmstem(i,j)) * 1e+04
c
c assign leaf, stem, root, and grain nitrogen concentrations
c to new variables (in percent)
c
            nconcl(i,j)    = fnleaf(i,j)   * 100.0
            nconcs(i,j)    = fnstem(i,j)   * 100.0
            nconcr(i,j)    = fnroot(i,j)   * 100.0
            nconcg(i,j)    = fngrain(i,j)  * 100.0
c   
c assign total nitrogen plant uptake to new variable for
c purposes of outputting data at end of year (kg/ha)
c
            cropn(i,j)     = totnuptake(i,j) * 1.e+04  
c
c assign total nitrogen fixation for crop to new variable for
c purposes of outputting this data at end of calendar year
c units of kg/ha
c
            cropfixn(i,j)  = totnfix(i,j) * 1.e+04
c
c carbon nitrogen ratio of plant residue goes to biogeochem.f
c and fine roots
c 
            if (fnresidue .gt. 0.0) then
              cntops(i,j) = min(cfrac(j) / fnresidue, 200.0) 
            else
              cntops(i,j) = 60.0
            endif
c
            if (fnroot(i,j) .gt. 0.0) then
              cnroot(i,j) = min(cfrac(j) / fnroot(i,j), 200.0) 
            else
              cnroot(i,j) = 80.0
            endif
c  
c assume that stem and leaf are both included in leaf litterfall value
c these annual total values (falll, fallr, fallw) cannot be changed 
c on a daily basis because biogeochem.f uses the annual total, split
c between each day of the year equally 
c carbon returned as residue
c
            falll(i) = falll(i) +  cbios(i,j) + aylprod(i,j)
            fallr(i) = fallr(i) +  ayrprod(i,j) 
            fallw(i) = fallw(i) +  cbiow(i,j)
c
c re-initialize crop variables
c for current crop at harvest date 
c
            plai(i,j)        = 0.01
            thrlai(i,j)      = 0.0
            peaklai(i,j)     = 0.0
            lai_decl_per_day(i,j) = 0.0
            ccdays(i,j)      = 0.0
            cbiol(i,j)       = 0.0
            cbior(i,j)       = 0.0
            cbios(i,j)       = 0.0
            cbiog(i,j)       = 0.0
            cbiow(i,j)       = 0.0
            hui(i,j)         = 0.0
            aybprod(i,j)     = 0.0
            ayrprod(i,j)     = 0.0
            ayabprod(i,j)    = 0.0
            aylprod(i,j)     = 0.0
            leafout(i,j)     = 0.0
            htmx(i,1)        = 0.0
            cumlvs(i,j)      = 0.0
            plaimx(i,j)      = 0.0
            dpgf(i,j)        = 0.0
            dpmature(i,j)    = 0.0
            biomass(i,j)     = 0.0
            totnuptake(i,j)  = 0.0
            tnplant(i,j)     = 0.0
            totnfix(i,j)     = 0.0
            idpp(i,j)        = 0.0
            gddplant(i,j)    = 0.0
            gddtsoi(i,j)     = 0.0
            sai(i,1)         = 0.0
            fu(i)            = 0.0
            lai(i,1)         = 0.0
            zbot(i,1)        = 0.0
            ztop(i,1)        = 0.0
            totbiol(i)       = 0.0
            totlail(i)       = 0.0  
            vf(i)            = 0.0  ! vernalization factor for winter wheat
            arepr(i,j)         = 0.0
           endif
 90     continue
c
       endif ! crop existence
c  
 100  continue       
c
c call to map out vegetation classes
c subroutine is within vegetation.f
c
      call vegmap
c
c return to main program
c
      return
      end
c
c ---------------------------------------------------------------------
      subroutine irrigation(iday, imonth)
c ---------------------------------------------------------------------
c
c routine that calculates the amount of water that is applied to a
c managed ecosystem on a daily basis (mm/day)
c
c based on average daily water content in the soil
c
c this amount will be evenly applied in timestep increments (start and stop time)
c based on the duration of the event - similar to the way precipitation
c events are handled 
c
c uses:
c
      use comgrid
      use compar
      use comatm
      use comsoi
      use comsum
      use comveg
      use comcrop
      use comnitr
      use combcs,only: cmask
c
      implicit none
c
c local variables
c
      integer
     >  iday,
     >  imonth,
     >  i, j, k
c
      real
     >  awcmax,
     >  awc 

c
      do 100 i = 1, npoi
        if(cmask(i) == 0 ) cycle
c
        if (irrigate(i) .eq. 1) then
c 
          if (iday .eq. 1 .and. imonth .eq. 1) then
            totirrig(i) = 0.0
          endif
c 
          awcmax = 0.0
          awc    = 0.0 
c 
          do 110 k = 1, 5       ! top 5 soil layers = 100 cm 
c 
c calculate the maximum available water in top 1 m (root zone) at field capacity
c 
            awcmax = awcmax + max(0.0, (sfield(i,k) - swilt(i,k))) *
     >        hsoi(k) * poros(i,k) * 100
c 
c calculate actual amount of water available to plant - current water content
c based on daily average water and ice content in soil layers down to a meter 
c 
            awc    = awc    + max(0.0, (adwisoilay(i,k)     +
     >        (1. - adwisoilay(i,k))        *
     >        adwsoilay(i,k)) - swilt(i,k)) *
     >        hsoi(k) * poros(i,k)          * 100
c 
 110      continue    
c 
c irrigation will occur if :
c 
c * the crop has been planted
c * the minimum daily temperature is > 5 C
c * the 5-day running mean temperature is > 10 C 
c * the actual soil water content in the top 1 m of soil
c is less than or equal to 50% of the maximum value at field capacity,
c 
          if (awc .le. 0.50 * awcmax 
     >      .and. (tmin(i) .gt. 278.16 .and. a5td(i) .gt. 283.16)
     >      .and. (icropsum(i) .gt. 0.)) then
c .and. (croplive(i,j) .gt. 0)  ! have to add in that irrigation only occurs when crop is planted! 
c 
c irrigation (applied above the canopy)  is used to make up the difference to field capacity
c convert awc values from cm to mm - this is a per day value 
c so it is consistent with the values used in weather.f (diurnal) for precipitation 
c 
c set upper limit based on literature search typical of what a farmer could
c apply in a typical day of irrigation and what rates of application are
c 
            xirrig(i) = min(150.0, max(0.0, awcmax - awc) * 10.0)
c 
          else 
c 
            xirrig(i) = 0.0
c 
          endif
          
        else  ! not irrigated
          
          xirrig(i) = 0.0
c 
        end if  ! if (irrigate(i))
c
 100  continue
c
c return to main program
c
      return
      end
c
c ---------------------------------------------------------------------
      subroutine nitrostress(istep,iday,imonth) 
c ---------------------------------------------------------------------
c
c subroutine calculates the effects of the amount of available inorganic
c nitrogen on carbon assimilation in crops 
c
c strictly speaking, stressn* is multiplied to the vmax parameters 
c used in the photosynthesis calculations
c
c uses:
c
      use comgrid
      use compar
      use com1d
      use comatm
      use comsoi
      use comsum
      use comveg
      use comcrop
      use comnitr
      use combcs,only: cmask
c
      implicit none
c
c local variables
c
      real 
*    >      smax,           ! maximum value stressn can have
     >      tnsupply,       ! total nitrogen supply from soil to plant 
     >      gfn,            ! function controlling nitrogen uptake
*    >      cnmax,          ! maximum allowable residue c/n ratio
     >      wsupply,        ! supply of water to plant from transpiration stream (kg/h20/m2/s)
*    >      alphacc,
*    >      gnmin,
     >      awc,
     >      sumnval,
     >      fnmin,
     >      wc,
     >      dmplant,
     >      fdmleaf,
     >      fdmstem,
     >      fdmroot,
     >      fdmgrain,
     >      dmres,
     >      flres,
     >      fsres,
     >      tngrain,
     >      tnleaf,
     >      tnstem,
     >      tnroot, 
     >      fnsmax,
     >      fnrmax,
     >      fnpmax, 
     >      f1, f2, f3
c
c biological fixation
c
      real
     >      gs,
     >      fxg,
     >      rdepth,
     >      fc,
     >      wp,
     >      sm, 
     >      fxw,
     >      rd,
     >      fxn,
     >      fxr,
*    >      availn,
     >      fnresidue
c
      real  plantn(npoi)    ! plant available nitrogen total
c
*      real  fnlfmx(npft),
*     >      fngrmx(npft),
*     >      sratio(npft),
*     >      rratio(npft),
*     >      fnopt(npft)
c
      integer iday,
     >        imonth,
     >        istep,
     >        iter,
     >        iter2,
     >        i, j, k
c
*        cfrac(13)    = 0.50  ! fraction of dry matter that is carbon   
*        cfrac(14)    = 0.50    
*        cfrac(15)    = 0.45    
c
*        fnlfmx(13)   = 0.025 ! average leaf nitrogen content  
*        fnlfmx(14)   = 0.013 ! average leaf nitrogen content
c
*        if (iwheat .eq. 1) fnlfmx(15)   = 0.008 ! average leaf nitrogen content  
*        if (iwheat .eq. 2) fnlfmx(15)   = 0.010 ! average leaf nitrogen content 
c
         fnlfmx(15) = fnlfmxw(iwheat_index) 
c
c        fnlfmx(15)   = 0.060 ! average leaf nitrogen content  
c
*        fngrmx(13)   = 0.040 ! maximum allowed grain nitrogen fraction
*        fngrmx(14)   = 0.017
c
*        if (iwheat .eq. 1) fngrmx(15)   = 0.0300 ! maximum allowed grain nitrogen fraction
*        if (iwheat .eq. 2) fngrmx(15)   = 0.0225 ! maximum allowed grain nitrogen fraction
c
         fngrmx(15) = fngrmxw(iwheat_index)
c
*        sratio(13)   = 0.40  ! typical ratio of stem to leaf nitrogen content
*        sratio(14)   = 0.05  
*        sratio(15)   = 0.40  ! typical ratio of stem to leaf nitrogen content
c
*        rratio(13)   = 0.75  ! typical ratio of root to leaf nitrogen content
*        rratio(14)   = 0.75  ! typical ratio of root to leaf nitrogen content
*        rratio(15)   = 1.00  ! typical ratio of root to leaf nitrogen content
c  
c
*        fnopt(13) = 0.0075   ! minimum leaf nitrogen concentration - stress onset - soybean
c        fnopt(14) = 0.03125  ! leaf nitrogen concentration - stress onset - maize 
*        fnopt(14) = 0.02850  ! leaf nitrogen concentration - stress onset - maize 
c
*        if (iwheat .eq. 1) fnopt(15) = 0.0175   ! stress onset for leaf nitrogen - spring wheat 
*        if (iwheat .eq. 2) fnopt(15) = 0.0100   ! stress onset for leaf nitrogen - winter wheat 
c
         fnopt(15) = fnoptw(iwheat_index)
c
      do 100 i = 1, npoi
        if(cmask(i) == 0 ) cycle
c
*        alphac         = 0.0005   ! minimum uptake rate applied to nitrogen (mm/timestep) 
*        gnmin         = 0.001    ! minimum nitrogen fraction allowed in grain  
*        smax          = 1.05     ! maximum nitrogen stress factor value     
*        availn        = 1.0      ! scaling variable to adjust for plant capability in capturing
                                 ! nitrogen from pools - above 1.0 means plant can take
                                 ! some up in excess of what is in the transpiration stream 
        fngrain(i,13) = 0.035  
        fngrain(i,14) = 0.013
        fngrain(i,15) = 0.020  
c
        tnsupply  = 0.0
        awc       = 0.0
        sumnval   = 0.0
c
c initialize layer nitrogen uptake
c
        do 180 k = 1, nsoilay
c
          tnuptake(i,k) = 0.0
          anuptake(i,k) = 0.0
c
 180    continue
c
        do 200 j = scpft, ecpft 
c
          if (exist(i,j) .eq. 1.0) then
            if (croplive(i,j) .eq. 1.0) then 
c
*             cnmax       = 95 
              fnmin       = cfrac(j) / cnmax
              stressn(i,j)= 1.0
              gfn         = 1.0
c
c calculate the total nitrogen supply rate (kg/m2/day) for each soil layer based on
c 1) the total daily water uptake for crops (lower canopy) (mm/day - stats.f)
c 2) total available nitrogen pool to roots (soil solution) (kg-no3/m2)
c 3) and available water content (mm)
c
c NOTE:  at this time, logic in IBIS cannot be used to determine
c the uptake of nitrogen for each specific pft (mixed in each grid
c cell because upsoil is for the entire lower canopy...will have
c to weight it for now on the lai of that pft [frac(i,j)] 
c
c since it is being called each timestep, use instantaneous values
c of soil ice and moisture
c  
        do 210 k = 1, nsoilay 
c
c calculate water content in each layer - based on EPIC parameterizations
c that look at actual water content in mm and not available water 
c
              wc          =  max(0.0, (wisoi(i,k)   +
     >                       (1.0 - wisoi(i,k))      *
     >                       wsoi(i,k))) * hsoi(k)  * 
     >                       poros(i,k)             * 1000
c
c
c alphac is minimum uptake rate to account for nitrogen usage
c even when transpiration is small (plant still able to take up
c nitrogen)
c
c allow plant to take up nitrogen in excess of the
c transpiration stream at low rates early in the
c season 
c
c supply of nitrogen to crops is from roots corresponding to [l]ower
c canopy in model - since this routine is being called each timestep
c use the value of upsoil(i,k) from canopy.f rather than the value from
c stats.f which is the daily average (adupsoil)   
c upsoil is in units of mm/m2/s of transpiration
c
             wsupply      =  upsoil(i,k) * dtime 
c
c make sure that water content is not zero - 
c set to small limit
c
             wc = max(1.0, wc)
c
c the total nitrogen uptake from the layer comes from the total n
c in the layer both in solution and in soil - leachable n is only
c that portion that is in the solution
c
c value of tnuptake for this layer is used in subroutine in solute
c leaching algorithm as a net sink of nitrogen to the layer 
c
c only allow uptake in layers that have roots
c make sure nitrogen uptake only occurs while crop is active
c the minimum rate will only be applied when the plant is not
c experiencing moisture stress
c
              if (froot(k,1) .gt. 0.005 .and. tnpptot(i) .gt. 0.0) then
                tnuptake(i,k) = max(alphac * stressl(i,k), wsupply) * availn *
     >                          (smsoil(i,k) + smsoln(i,k)) 
              else
                tnuptake(i,k) = 0.0
              endif
c
 210    continue
c
c            endif
           
        if (aybprod(i,j) .gt. 0.0  .and.
     >      aylprod(i,j) .gt. 0.0) then
c
c for the purpose of dealing with total nitrogen uptake,
c we have to use year to date total carbon production
c in these equations because some root and leaf biomass
c has been adjusted due to phenology in the model
c
            dmplant   =  aybprod(i,j)  / cfrac(j)
            fdmleaf   =  (aylprod(i,j) / cfrac(j)) / dmplant
            fdmstem   =  (cbios(i,j)   / cfrac(j)) / dmplant
            fdmroot   =  (ayrprod(i,j) / cfrac(j)) / dmplant
            fdmgrain  =  (cbiog(i,j)   / cfrac(j)) / dmplant
c
            dmres     =  (aylprod(i,j) + cbios(i,j)) / cfrac(j) 
            flres     =  (aylprod(i,j) / cfrac(j)) / dmres
            fsres     =  (cbios(i,j)   / cfrac(j)) / dmres  
c
            fnplant(i,j)   =  max(0.0, totnuptake(i,j) / dmplant) 
c
c maintain minimum nitrogen concentration in leaf and stem (potential residue)
c
            iter  = 0
            iter2 = 0
c
c           write(*,*) fnleaf(i,j), fnresidue, cbiog(i,j) 
 400        fnleaf(i,j) = (fnplant(i,j) - fngrain(i,j) * fdmgrain) /
     >                    (fdmleaf + sratio(j) * fdmstem +
     >                     rratio(j) * fdmroot)
c
            fnleaf(i,j) = max(0.0, fnleaf(i,j))
            fnstem(i,j) = fnleaf(i,j) * sratio(j)
            fnroot(i,j) = fnleaf(i,j) * rratio(j)
c
            fnresidue   = fnleaf(i,j) * flres + fnstem(i,j) * fsres 
            if (fnresidue .gt. fnmin) iter2 = 1
            if (fnresidue .le. fnmin) iter  = 1 

c            if (iter .eq. 1 .and. fngrain(i,j) .gt. gnmin
c     >                      .and. fdmgrain     .gt. 0.0
c     >                      .and. iter2        .eq. 0) then 
c                    fngrain(i,j) = max(gnmin, fngrain(i,j) * 0.99)
c                    write(*,*) 'taking from grain', fngrain(i,j)
c                    goto 400
c         
            if (iter2 .eq. 1 .and. fngrain(i,j) .lt. fngrmx(j) 
     >                      .and. fdmgrain     .gt. 0.0
     >                      .and. iter         .eq. 0) then 
                    fngrain(i,j) = min(fngrmx(j), fngrain(i,j) * 1.01)
                    goto 400
            endif
c
c ------------------------------------------------------------------------------
c calculate nitrogen content in various pools
c
            tngrain = fngrain(i,j) * fdmgrain * dmplant
            tnleaf  = fnleaf(i,j)  * fdmleaf  * dmplant      
            tnstem  = fnstem(i,j)  * fdmstem  * dmplant      
            tnroot  = fnroot(i,j)  * fdmroot  * dmplant      
c
            tnplant(i,j) = tngrain + tnleaf + tnstem + tnroot 
c
            fnsmax   = sratio(j)  * fnlfmx(j)
            fnrmax   = rratio(j)  * fnlfmx(j) 
c
            fnpmax   = fnlfmx(j) * fdmleaf + fnsmax * fdmstem +
     >                 fnrmax * fdmroot + fngrmx(j) * fdmgrain               
c
c
c calculate function controlling rate of nitrogen uptake
c based on plant nitrogen concentration and maximum value
c there is a chance early in growing season that fnplant
c could be higher than fnpmax - thus maximize the below
c fraction to be .le. 1.0 
c
             gfn = 1. - min(1.0, fnplant(i,j) / fnpmax) ** 1.0 
c
c calculate the annual running total of the actual nitrogen
c uptake for all pfts in lower canopy (crops)
c
c adjust nitrogen uptake by plant by gfn factor - equal in each layer
c
             do 310 k = 1, nsoilay 
               anuptake(i,k) = tnuptake(i,k) * gfn
               tnsupply      = tnsupply + anuptake(i,k)
 310         continue 
c
             totnuptake(i,j) = totnuptake(i,j) + tnsupply              
c
             totnuptake(i,j) = max(0.0, min((1.0 - cfrac(j)) * dmplant, totnuptake(i,j))) 
c
c calculate stress parameter using rectangular hyperbola which
c relates leaf nitrogen concentration to vmax        
c
c rectangular hyperbola 
c f1 and f2  control the shape of the stress response function 
c which spans from 0 to 1.0 for leaf n concentrations from 0 to 4 percent 
c
c s-shaped curve nitrogen limitation effect from epic model
c
              f1    = 8.5 
              f2    = 11.0 
c
c ratio of leaf nitrogen concentration to the optimal maximum
c for corn/maize
c
             f3 = 2 * (fnleaf(i,j) / fnopt(j))
c
             stressn(i,j) = min (smax, (f3 / (f3 + exp(f1 - f2 * f3))))
c
             stressn(i,j) = max (0.10, stressn(i,j)) 
c
c----------------------------------------------------------------------
c biological fixation of nitrogen through symbiosis in soybeans
c----------------------------------------------------------------------
c
c Key reference:
c M. Cabelguenne et al., Agricultural systems 60: 175-196, 1999.
c this module is taken from the new epicphase model
c
c the amount of daily n-fixation by the plant is based on a fraction
c of the total daily n-uptake.  It is controlled by three main factors:
c
c * growth stage of the crop (0-1) *
c * soil moisture            (0-1) * 
c * nitrogen in rooting zone (0-1) * 
c
c the growth stage factor (fxp) inhibits fixation in young plants
c and old plants, and peaks between 30-55% of the crop cycle
c
c the soil water content factor (fxw) reduces n-fixation when the 
c water content in the top 0.3 m is less than 85% of field capacity
c
c the soil nitrogen (plant available) factor (fxn) reduces n-fixation
c when the nitrogen amount in the root zone is greater than 100 kg/ha  
c
c  
             if (j .eq. 13) then
c
c calculate growth stage and factor (fraction of total average gdd)
c
               gs = hui(i,j) / gddmaturity(i,j)
c
               if     (gs .le. 0.15 .or.  gs .ge. 0.75) then
                 fxg = 0.0
               elseif (gs .gt. 0.15 .and. gs .le. 0.30) then
                 fxg = 6.67 * gs - 1.0
               elseif (gs .gt. 0.30 .and. gs .le. 0.55) then
                 fxg = 1.0
               else
                 fxg = 3.75 - 5.0 * gs   
               endif
c
c calculate effect of soil moisture in top 25-30 cm
c
c              rdepth = 1. / (hsoi(1) + hsoi(2) + hsoi(3))  
               rdepth = 1. / (hsoi(1) + hsoi(2))  
c
               fc = 0.0
               wp = 0.0
               sm = 0.0
               plantn(i) = 0.0
c
               do 220 k = 1, 2
                 fc = fc + sfield(i,k) * hsoi(k)
                 wp = wp + swilt(i,k)  * hsoi(k)
                 sm = sm + wsoi(i,k)   * hsoi(k)
c
c calculate available plant nitrogen total in these layers
c
                 plantn(i) = plantn(i) + smsoil(i,k) + smsoln(i,k)
c
 220           continue 
c
               fc = fc * rdepth
               wp = wp * rdepth
               sm = sm * rdepth
               sm = min(sm, 0.85 * (fc - wp) + wp)
c
               fxw = (sm - wp) / (0.85 * (fc - wp)) 
c
c calculate effect of plant available nitrogen pool
c 
               rd = 1.0   ! rooting depth in meters

               if     (plantn(i) .gt. 0.0300) then
                 fxn = 0.0
               elseif (plantn(i) .le. 0.0100) then
                 fxn = 1.0
               else
                 fxn = 1.5 - 0.005 * (plantn(i) * 10000) / rd
               endif
c
c equation for fxn has to be in kg/ha nitrogen and meters for rd 
c
c CJK replaced 2/1/2006  
c the 1.5 at the end of the following equation was to increase
c the annual n-fixation rates by 50% because of simulated errors
c and a low bias compared to observations
c
               fxr  = min(1.0, fxw, fxn) * fxg 
c
               fixn(i,j) = fxr * tnsupply 
c              fixn(i,j) = 0.0 
c 
c fixn is thus calculated each timestep
c
c update plant available nitrogen pool for fixation from soybean
c
             else  ! non-soybean crop - no fixation
c
               fixn(i,j) = 0.0
c
             endif    ! soybeans only fix nitrogen
c
             totnfix(i,j) = totnfix(i,j) + fixn(i,j)
c
             do 240 k = 1, nsoilay 
c
c critical: what quantity do we add the nitrogen fixation
c to?  previous? depends on what order subroutines are called
c 
               smsoil(i,k) = smsoil(i,k) + fixn(i,j) * froot(k,1)
c
 240         continue
c
           endif   ! production gt 0
         endif    ! crop plant  
        endif     ! crop existence 
c
 200  continue    ! crop pft 
c
 100  continue    ! grid cell 
c
c return to program 
c
      return
      end
c ---------------------------------------------------------------------
      subroutine cropoutput(iday,imonth,iyear,jday) 
c ---------------------------------------------------------------------
c
c output crop variables for a single grid cell - diagnostic output 
c
c uses:
c
      use comgrid
      use compar
      use comatm
      use comsoi
      use comsum
      use comveg
      use comcrop
      use comnitr
      use combcs,only: cmask
c
      implicit none
c
c local variables
c
      integer
     >         iday,
     >         jday,
     >         imonth,
     >         iyear,
     >         i,j,k,l,n
c
      real
     >         tnsoi,
     >         nuptake,
     >         hi,
     >         grain,
     >         fyld,
     >         yld

c define soil layer for which you want output for
c calculate total immobile/mobile N pool totals to
c specified soil depth (isoilay input)
c
      do 100 i = 1, npoi 
        if(cmask(i) == 0 ) cycle
c 
                tsnimm(i) = 0.
                tsnmob(i) = 0.
                tnsoi     = 0.
c
c total soil inorganic nitrogen pools available to plant for uptake
c
                do 330 k = 1, isoilay
                  tsnimm(i) = tsnimm(i) + smsoil(i,k)     ! immobile pool
                  tsnmob(i) = tsnmob(i) + smsoln(i,k)     ! mobile pool (e.g., nitrate) 
                  tnsoi = tnsoi + smsoil(i,k) + smsoln(i,k)  
 330            continue
c
c output data for model studies
c cropn is in kg/ha, totnvegn is not 
c
                  nuptake = 0.0
                  do 335 j = scpft, ecpft
                    nuptake = nuptake + cropn(i,j)
 335              continue
                    nuptake = nuptake + totnvegn(i) * 1.e+04
c
c write out variables 
c
                  do 380 n = scpft, ecpft
c
                    if (exist(i,n) .eq. 1.0 .and. cropplant(i,n) .eq. 1.
     >                  .and. i .eq. 1 .and. harvdate(i,n) .ne. 999) then 
c
c grid cell 1 corresponds to that which is closest to lat, lon for ARL
c
                       open(16, file = 'crop.output.dat',status = 'unknown')
                       write(16,340) 
     >                             iyear,
     >                             cropyld(i,n),                 
c
c     >                             iyear,
c     >                             n,
c     >                             tnsoi*1.e04,
     >                             fertinput(i,n),
c    >                             ydeposn(i)*1.e04,
     >                             aynmintot(i)*1.e04,
c    >                             ayimmtot(i)*1.e04,
     >                             nuptake, 
c     >                             yno3leach(i),
c     >                             concn(i),
c    >                             snbalance(i),
c     >                             cropyld(i,n),
c    >                             dmyield(i,n),
     >                             harvidx(i,n),
c     >                             dmleaf(i,n),
c     >                             dmstem(i,n),
c     >                             dmroot(i,n),
c     >                             dmyield(i,n),
c     >                             dmcrop(i,n),
c     >                             cntops(i,n),
c    >                             cnroot(i,n),
     >                             nconcl(i,n),
     >                             nconcs(i,n),
     >                             nconcr(i,n),
     >                             nconcg(i,n),
     >                             grainn(i,n),
c     >                             drntot(i),
c    >                             ayprcp(i),
     >                             croplaimx(i,n),
     >                             cropfixn(i,n),
c    >                             totcsoi(i),
     >                             idop(i,n),
     >                             harvdate(i,n) 
                     endif
 340                format(i6,5f8.2,4f8.4,3f8.2,i5,i5)
c 340                format(i6,f8.2)
c
 380    continue  ! loop through all pfts
c
 100    continue
c
      return
      end

      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! purpose : initialize pft when iyr_pft == 1
      !           
      !
      ! note    : use first year pft data for initialization
      !
      !           Last updated 2015-08-03 by Yuwei Li
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !
      subroutine initial_pft(irestart,iyear0,irstyear)
        use comcrop, only: isoybean, imaize, iwheat
        use comgrid, only: npoi
        use comwork, only: cdummy

        implicit NONE
        !--input
        integer:: irestart,iyear0,irstyear
        !--local
        integer:: first_year,i,idx

        if(irestart == 1) then
          first_year = irstyear
        elseif(irestart == 0) then
          first_year = iyear0
        else
          write(6,*) "(initial_pft)Error! read-in irestart: ",irestart
        end if

        write(6,*) "iyr_pft == 1, using first year data for initialization..."
        call read_NetCDF('pft',first_year)

        isoybean=0;imaize=0;iwheat=0
        do i =1, npoi
          idx = int(cdummy(i))

          if(idx == 13)     then
            isoybean= 1
          elseif(idx == 14) then
            imaize  = 1
          elseif(idx == 15) then
            iwheat  = 1
          end if

        end do !on i =1, npoi

        return
      end subroutine initial_pft
