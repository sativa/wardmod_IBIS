c io.f  last update 11.06.03 C. Kucharik
c       for new climate datasets
c
c    #     ####
c    #    #    #
c    #    #    #
c    #    #    #
c    #    #    #
c    #     ####
c This file contains these subroutines
c wrestart
c wdaily
c wmonthly
c wyearly
c readit
c read_fert_data
c read_one_fert_file
c read_management_header
c read_management_data
c read_one_management_variable
c set_irrigation_map
c restart
c coldstart
c rdanom
c inird
c rdday
c diaginit
c wdiag
c
c See end of file for information on how to add new code to read a
c file or write a file.
c
c ---------------------------------------------------------------------
      subroutine wrestart (nday, iyear, iyear0)
c ---------------------------------------------------------------------
c
c this subroutine writes the restart values of:
c
c  fsnocov = fractional snow cover
c  tsno    = temperature of snow
c  hsno    = snow depth
c  tsoi    = soil temperature
c  wisoi   = soil ice content
c  wsoi    = soil moisture content
c  smsoil  = immobile inorganic nitrogen in soil
c  smsoln  = mobile inorganic nitrogen in soil
c  cbiol   = carbon in leaf biomass pool
c  cbiow   = carbon in woody biomass pool
c  cbior   = carbon in fine root biomass pool
c  cnroot  = c/n ratio of crop roots
c  cntops  = c/n ratio of tops of crop plants
c  sapfrac = sapwood fraction
c  decompl = litter decomposition factor
c  decomps = soil decomposition factor
c  clitlm  = leaf metabolic litter
c  clitls  = leaf structural litter
c  clitll  = leaf lignin litter
c  clitrm  = root metabolic litter
c  clitrs  = root structural litter
c  clitrl  = root lignin litter
c  clitwm  = woody metabolic litter
c  clitws  = woody structural litter
c  clitwl  = woody lignin litter
c  falll   = annual leaf litterfall
c  fallr   = annual fine root turnover 
c  fallw   = annual wood litterfall
c  totcmic = total microbial carbon
c  csoislop= slow soil carbon, protected humus
c  csoislon= slow soil carbon, nonprotected humus
c  csoipas = passive soil carbon
c  aplantn = available plant nitrogen
c  gdd0    = growing degree days 0
c  gdd0c   = growing degree days 0 for wheat between April 1 and Sept 30
c  gdd5    = growing degree days 5
c  gdd8    = growing degree days 8 
c  gdd10   = growing degree days 10 
c  tc      = coldest monthly temperature
c  tw      = warmest monthly temperature
c  wipud   = ice content of puddles per soil area
c  wpud    = liquid content of puddles per soil area
c  agddu   = annual accumulated growing degree days for bud burst, upper canopy
c  agddl   = annual accumulated growing degree days for bud burst, lower canopy
c  tempu   = cold-phenology trigger for trees
c  templs  = cold-phenology trigger for shrubs
c  greenfracl3 = fraction of green vegetation in C3 grasses
c  greenfracl4 = fraction of green vegetation in C4 grasses
c  Tavgann = average annual air temperature (purely from climatology)
c  PPTavgann = average annual precipitation (purely from climatology)
c  a10td    = 10-day avg daily temp
c  a10ts    = 10-day avg daily soil (1) temp
c  a10ancub = 10-day average canopy photosynthesis rate - broadleaf
c  a10ancuc = 10-day average canopy photosynthesis rate - conifers
c  a10ancls = 10-day average canopy photosynthesis rate - shrubs
c  a10ancl4 = 10-day average canopy photosynthesis rate - c4 grasses
c  a10ancl3 = 10-day average canopy photosynthesis rate - c3 grasses
c  a10scalparamu = 10-day average canopy scaling parameter - upper canopy
c  a10scalparaml = 10-day average canopy scaling parameter - lower canopy
c  a10daylightu = 10-day average daylight - upper canopy
c  a10daylightl = 10-day average daylight - lower canopy
c  a11soiltd = 11-day average surface soil temperature
c  a3tdmin = 3-day average daily minimum air temperature
c (NOTE: a10ancuc is not used at this point, so its wrestart entry 
c is commented out)
c ---------------------------------------------------------------------
c
c instantaneous output for restarts
c
      use comgrid
      use compar
      use comsoi
      use comsno
      use comsum
      use comveg
      use comwork
      use comcrop
      use comnitr
c
      implicit none
c
c Arguments
c
      integer nday,         ! number of days run since iyear0
     >        iyear,        ! this calendar year
     >        iyear0        ! initial year

c
c local variables
c
      integer lf,           ! number of characters in directory name
     >        n, k,         ! loop indices
     >        idies,        ! file indice (?) for netcdf
     >        istat,        ! error flag for netcdf      
     >        nyears        ! years of run iyear0
c
      integer istart(4),
     >        icount(4)         ! for writing restart vars
c
      character*21  tunits
      character*200 fdir         ! used to construct odd/even file names
      character*10  cdate        ! date to use in history attribute in files
      character*10  tdate        ! character date for time step
      character*80  dimnames(4)  ! names of dimensions for restart vars
      character*80  pftdef(npft) ! plant functional type defs (not used now)
      character*80  filen        ! file name
c
      real slayers(nsnolay),    ! index for snow layers
     >     depthsoi(nsoilay),   ! soil layer depths
     >     pindex(npft),        ! index for pfts
     >     ftime(1),            ! floating point time value; the dummy (1) is added for more robust compiling (Y.Li)
     >     tweight(1)           ! time weight (# days/sample); dummy (1) added

! dummy variables for a better compiling (more robust) (Y.Li)
      real:: dummy_vals3rd(1)   !for  *.nc 
c
c External
c
      integer lenchr,           ! Function: Find length of character string
     > NF_PUT_ATT_TEXT,         ! netcdf function
     > NF_GLOBAL                ! netcdf function
c
c ---------------------------------------------------------------------
c
      data istart / 1,1,1,1 /, 
     >     icount / nlon,nlat,1,1 /
      icount(1) = nlonsub
      icount(2) = nlatsub
c
      tweight = 1.

c
c check to see if iyear is odd or even, construct appropriate file name root
c
      if (mod(iyear,2) .eq. 0) then
         fdir = trim(myprocDir)//'restart/even'
      else
         fdir = trim(myprocDir)//'restart/odd'
      end if
      lf = lenchr(fdir)
c
c tdate is december of this year, max year = 999
c
      tdate='DEC000'//char(0)//char(0)//char(0)//char(0)
      nyears = iyear - iyear0 + 1

      !for restart run when restart files may not exist, yet nyears > 2 with the above computation and no file will then be created
      if(irestart_glo == 1) then   
        nyears = iyear - iyrlast_glo
      end if

      if (nyears .lt. 10) then
         write(tdate(6:6),'(i1)') nyears
       else if (nyears .lt. 100) then
         write(tdate(5:6),'(i2)') nyears
      else
         write(tdate(4:6),'(i3)') nyears
      end if
c
c initialize snow layer indicies, pft names, etc
c
      if (nyears .le. 2) then
c         call date(cdate)
         ftime = nday
c
c time units is days since Dec 31 of the year before iyear0
c
         tunits = 'days since 0000-12-31'
         write(tunits(12:15),'(i4)') iyear0-1
c
         do 5 n = 1, nsnolay
            slayers(n) = float(n)
 5       continue
c
         depthsoi(1) = hsoi(1)
         do 10 n = 2, nsoilay
            depthsoi(n) = depthsoi(n-1)+hsoi(n)
 10      continue
c
c define pft by index
c
         do 12 n = 1, npft
            pindex(n) = n
 12      continue
c
c and by character label
c
         pftdef(1) = 'trbrevtr - tropical broadleaf evergreen trees'
     >    //char(0)
         pftdef(2) = 
     >    'trbrdetr - tropical broadleaf drought-deciduous trees'
     >    //char(0)
         pftdef(3) = 
     >    'wtbrevtr - warm-temperate broadleaf evergreen trees'
     >    //char(0)
         pftdef(4) = 'tecoevtr - temperate conifer evergreen trees'
     >    //char(0)
         pftdef(5) = 
     >    'tebrdetr - temperate broadleaf cold-deciduous trees'//char(0)
         pftdef(6) = 'bocoevtr - boreal conifer evergreen trees'
     >    //char(0)
         pftdef(7) = 'bocodetr - boreal conifer cold-deciduous trees'
     >    //char(0)
         pftdef(8) = 
     >    'bobrdetr - boreal broadleaf cold-deciduous trees'//char(0)
         pftdef(9) = 'evsh - evergreen shrubs'//char(0)
         pftdef(10) = 'desh - deciduous shrubs'//char(0)
         pftdef(11) = 'c4gr - warm (c4) grasses'//char(0)
         pftdef(12) = 'c3gr - cool (c3) grasses'//char(0)
         pftdef(13) = 'c3 crop - soybean'//char(0)
         pftdef(14) = 'c4 crop - maize'//char(0)
         pftdef(15) = 'c3 crop - wheat'//char(0)
c
         dimnames(1) = 'longitude'
         dimnames(2) = 'latitude'
c
c dimnames(3) is set for each variable seperately
c
         dimnames(4) = 'time'
      end if
cc
cc dummy variable example, 3-d - copy & modify for new variable.
cc If new variable is 4-d, see tsno to see how 4-d differs from 3-d.
cc
c      filen = fdir(1:lf)//'dummyv.nc'
c      if (nyears .le. 2) then
c         call inifile(idies,filen,
c     >    'restart file for dummyv',
c     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'dummyv',
c     >    'instantaneous dummyv','dummyvs-units',3,dimnames,
c     >    OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (dummyv, cdummy)
c      icount(3) = 1
c      call writevar(filen,'dummyv',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wrestart, dummyv'
c         stop 1
c      end if
c
c fractional snow cover
c
      filen = fdir(1:lf)//'fsnocov.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for fractional snow cover',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'fsnocov',
     >    'instantaneous fractional snow cover','fraction',3,dimnames,
     >    OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (fi, cdummy)
      icount(3) = 1
      call writevar(filen,'fsnocov',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, fsnocov'
         stop 1
      end if
c
c temperature of snow layers
c
      filen = fdir(1:lf)//'tsno.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for snow temperature',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'snowlayer','snow layers top to bottom','',nsnolay,slayers,
     >    'down',tunits,'gregorian',istat)
         dimnames(3) = 'snowlayer'
         call inivar(idies,'tsno',
     >    'instantaneous snow cover temperature','degK',4,dimnames,
     >    OCEAN,istat)
         call endini(idies,istat)
      end if
      do 15 k = 1, nsnolay
         call vec2arr (tsno(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 15   continue
      icount(3) = nsnolay
      call writevar(filen,'tsno',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, tsno'
         stop 1
      end if
c
c thickness of snow layers
c
      filen = fdir(1:lf)//'hsno.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for snow layer thickness',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'snowlayer','snow layers top to bottom','',nsnolay,slayers,
     >    'down',tunits,'gregorian',istat)
         dimnames(3) = 'snowlayer'
         call inivar(idies,'hsno','instantaneous snow layer thickness',
     >    'meters',4,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      do 20 k = 1, nsnolay
         call vec2arr (hsno(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 20   continue
      icount(3) = nsnolay
      call writevar(filen,'hsno',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, hsno'
         stop 1
      end if
c
c temperature of soil layers
c
      filen = fdir(1:lf)//'tsoi.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for soil temperature',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'soillayer','depth of soil layer bottom','meter',nsoilay,
     >    depthsoi,'down',tunits,'gregorian',istat)
         dimnames(3) = 'soillayer'
         call inivar(idies,'tsoi','instantaneous soil temperature',
     >    'degK',4,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      do 25 k = 1, nsoilay
         call vec2arr (tsoi(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 25   continue
      icount(3) = nsoilay
      call writevar(filen,'tsoi',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, tsoi'
         stop 1
      end if
c
c ice content of soil
c
      filen = fdir(1:lf)//'wisoi.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for soil ice content',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'soillayer','depth of soil layer bottom','meter',nsoilay,
     >    depthsoi,'down',tunits,'gregorian',istat)
         dimnames(3) = 'soillayer'
         call inivar(idies,'wisoi',
     >    'instantaneous fraction of soil pore space containing ice',
     >    'fraction',4,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      do 30 k = 1, nsoilay
         call vec2arr (wisoi(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 30   continue
      icount(3) = nsoilay
      call writevar(filen,'wisoi',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, wisoi'
         stop 1
      end if
c
c water content of soil
c
      filen = fdir(1:lf)//'wsoi.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for soil water content',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'soillayer','depth of soil layer bottom','meter',nsoilay,
     >    depthsoi,'down',tunits,'gregorian',istat)
         dimnames(3) = 'soillayer'
         call inivar(idies,'wsoi',
     >    'instantaneous fraction of soil pore space containing water',
     >    'fraction',4,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      do 35 k = 1, nsoilay
         call vec2arr (wsoi(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 35   continue
      icount(3) = nsoilay
      call writevar(filen,'wsoi',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, wsoi'
         stop 1
      end if
c
c immobile inorganic nitrogen in soil - available to plants, but not leachable 
c
      filen = fdir(1:lf)//'smsoil.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for immobile soil inorganic nitrogen',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'soillayer','depth of soil layer bottom','meter',nsoilay,
     >    depthsoi,'down',tunits,'gregorian',istat)
         dimnames(3) = 'soillayer'
         call inivar(idies,'smsoil',
     >    'instantaneous immobile inorganic soil nitrogen',
     >    'kg/m2',4,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      do 36 k = 1, nsoilay
         call vec2arr (smsoil(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 36   continue
      icount(3) = nsoilay
      call writevar(filen,'smsoil',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, smsoil'
         stop 1
      end if
c
c mobile inorganic nitrogen in soil - available to plants, and potentially leachable 
c
      filen = fdir(1:lf)//'smsoln.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for leachable soil inorganic nitrogen',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'soillayer','depth of soil layer bottom','meter',nsoilay,
     >    depthsoi,'down',tunits,'gregorian',istat)
         dimnames(3) = 'soillayer'
         call inivar(idies,'smsoln',
     >    'instantaneous mobile inorganic soil nitrogen',
     >    'kg/m2',4,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      do 37 k = 1, nsoilay
         call vec2arr (smsoln(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 37   continue
      icount(3) = nsoilay
      call writevar(filen,'smsoln',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, smsoln'
         stop 1
      end if
c
c carbon in leaf biomass pool
c
      filen = fdir(1:lf)//'cbiol.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for carbon in leaf biomass pool',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'pft','plant functional type','_',npft,pindex,'',
     >    tunits,'gregorian',istat)
c add global attribute to define pfts with text, use netcdf low-level command
         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',
     >    npft*80,pftdef)
         dimnames(3) = 'pft'
         call inivar(idies,'cbiol',
     >    'instantaneous carbon in leaf biomass pool','kg/m^2',
     >    4,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      do 40 k = 1, npft
         call vec2arr (cbiol(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 40   continue
      icount(3) = npft
      call writevar(filen,'cbiol',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, cbiol'
         stop 1
      end if
c
c carbon in wood
c
      filen = fdir(1:lf)//'cbiow.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for carbon in wood biomass pool',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'pft','plant functional type','_',npft,pindex,'',
     >    tunits,'gregorian',istat)
c add global attribute to define pfts with text, use netcdf low-level command
         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',
     >    npft*80,pftdef)
         dimnames(3) = 'pft'
         call inivar(idies,'cbiow',
     >    'instantaneous carbon in wood biomass pool','kg/m^2',
     >    4,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      do 45 k = 1, npft
         call vec2arr (cbiow(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 45   continue
      icount(3) = npft
      call writevar(filen,'cbiow',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, cbiow'
         stop 1
      end if
c
c carbon in root
c
      filen = fdir(1:lf)//'cbior.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for carbon in root biomass pool',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'pft','plant functional type','_',npft,pindex,'',
     >    tunits,'gregorian',istat)
c add global attribute to define pfts with text, use netcdf low-level command
         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',
     >    npft*80,pftdef)
         dimnames(3) = 'pft'
         call inivar(idies,'cbior',
     >    'instantaneous carbon in root biomass pool','kg/m^2',
     >    4,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      do 50 k = 1, npft
         call vec2arr (cbior(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 50   continue
      icount(3) = npft
      call writevar(filen,'cbior',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, cbior'
         stop 1
      end if
c
c c/n ratio of crop roots 
c
      filen = fdir(1:lf)//'cnroot.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for cn ratio of crop roots',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'pft','plant functional type','_',npft,pindex,'',
     >    tunits,'gregorian',istat)
c add global attribute to define pfts with text, use netcdf low-level command
         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',
     >    npft*80,pftdef)
         dimnames(3) = 'pft'
         call inivar(idies,'cnroot',
     >    'end of year cn ratio of crop roots','dimensionless',
     >    4,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      do 52 k = 1, npft
         call vec2arr (cnroot(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 52   continue
      icount(3) = npft
      call writevar(filen,'cnroot',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, cnroot'
         stop 1
      end if
c
c c/n ratio of tops of crop plants 
c
      filen = fdir(1:lf)//'cntops.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for cn ratio of crop tops',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'pft','plant functional type','_',npft,pindex,'',
     >    tunits,'gregorian',istat)
c add global attribute to define pfts with text, use netcdf low-level command
         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',
     >    npft*80,pftdef)
         dimnames(3) = 'pft'
         call inivar(idies,'cntops',
     >    'end of year cn ratio of crop tops','dimensionless',
     >    4,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      do 54 k = 1, npft
         call vec2arr (cntops(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 54   continue
      icount(3) = npft
      call writevar(filen,'cntops',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, cntops'
         stop 1
      end if
c
c sapwood fraction
c
      filen = fdir(1:lf)//'sapfrac.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for sapwood fraction',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'sapfrac','instantaneous sapwood fraction',
     >    'fraction',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (sapfrac, cdummy)
      icount(3) = 1
      call writevar(filen,'sapfrac',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, sapfrac'
         stop 1
      end if
c
c litter decomposition factor 
c
      filen = fdir(1:lf)//'decompl.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for litter decomposition factor',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'decompl',
     >    'litter decomposition factor','dimensionless',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (decompl, cdummy)
      icount(3) = 1
      call writevar(filen,'decompl',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, decompl'
         stop 1
      end if
c
c soil decomposition factor 
c
      filen = fdir(1:lf)//'decomps.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for soil decomposition factor',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'decomps',
     >    'soil decomposition factor','dimensionless',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (decomps, cdummy)
      icount(3) = 1
      call writevar(filen,'decomps',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, decomps'
         stop 1
      end if
c
c leaf metabolic litter
c
      filen = fdir(1:lf)//'clitlm.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for leaf metabolic litter',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'clitlm',
     >    'instantaneous leaf metabolic litter carbon','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (clitlm, cdummy)
      icount(3) = 1
      call writevar(filen,'clitlm',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, clitlm'
         stop 1
      end if
c
c leaf structural carbon litter
c
      filen = fdir(1:lf)//'clitls.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for leaf structural litter',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'clitls',
     >    'instantaneous leaf structural litter carbon','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (clitls, cdummy)
      icount(3) = 1
      call writevar(filen,'clitls',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, clitls'
         stop 1
      end if
c
c leaf lignin carbon litter
c
      filen = fdir(1:lf)//'clitll.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for leaf lignin litter',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'clitll',
     >    'instantaneous leaf lignin litter carbon','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (clitll, cdummy)
      icount(3) = 1
      call writevar(filen,'clitll',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, clitll'
         stop 1
      end if
c
c root metabolic litter
c
      filen = fdir(1:lf)//'clitrm.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for root metabolic litter',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'clitrm',
     >    'instantaneous root metabolic litter carbon','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (clitrm, cdummy)
      icount(3) = 1
      call writevar(filen,'clitrm',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, clitrm'
         stop 1
      end if
c
c root structural litter
c
      filen = fdir(1:lf)//'clitrs.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for root structural litter',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'clitrs',
     >    'instantaneous root structural litter carbon','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (clitrs, cdummy)
      icount(3) = 1
      call writevar(filen,'clitrs',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, clitrs'
         stop 1
      end if
c
c root lignin litter
c
      filen = fdir(1:lf)//'clitrl.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for root lignin litter',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'clitrl',
     >    'instantaneous root lignin litter carbon','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (clitrl, cdummy)
      icount(3) = 1
      call writevar(filen,'clitrl',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, clitrl'
         stop 1
      end if
c
c woody metabolic litter
c
      filen = fdir(1:lf)//'clitwm.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for woody metabolic litter',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'clitwm',
     >    'instantaneous woody metabolic litter carbon','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (clitwm, cdummy)
      icount(3) = 1
      call writevar(filen,'clitwm',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, clitwm'
         stop 1
      end if
c
c woody structural litter
c
      filen = fdir(1:lf)//'clitws.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for woody structural litter',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'clitws',
     >    'instantaneous woody structural litter carbon','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (clitws, cdummy)
      icount(3) = 1
      call writevar(filen,'clitws',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, clitws'
         stop 1
      end if
c
c woody lignin litter
c
      filen = fdir(1:lf)//'clitwl.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for woody lignin litter',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'clitwl',
     >    'instantaneous woody lignin litter carbon','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (clitwl, cdummy)
      icount(3) = 1
      call writevar(filen,'clitwl',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, clitwl'
         stop 1
      end if
c
c annual leaf litterfall
c
      filen = fdir(1:lf)//'falll.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for annual leaf litterfall',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'falll',
     >    'annual leaf litterfall carbon','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (falll, cdummy)
      icount(3) = 1
      call writevar(filen,'falll',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, falll'
         stop 1
      end if
c
c annual fine root turnover 
c
      filen = fdir(1:lf)//'fallr.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for annual fine root turnover',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'fallr',
     >    'annual fine root turnover carbon','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (fallr, cdummy)
      icount(3) = 1
      call writevar(filen,'fallr',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, fallr'
         stop 1
      end if
c
c annual wood turnover 
c
      filen = fdir(1:lf)//'fallw.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for annual woody turnover',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'fallw',
     >    'annual wood turnover carbon','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (fallw, cdummy)
      icount(3) = 1
      call writevar(filen,'fallw',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, fallw'
         stop 1
      end if
c
c total microbial carbon
c
      filen = fdir(1:lf)//'totcmic.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for total microbial carbon',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'totcmic',
     >    'instantaneous total microbial carbon','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (totcmic, cdummy)
      icount(3) = 1
      call writevar(filen,'totcmic',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, totcmic'
         stop 1
      end if
c
c slow soil carbon, protected humus
c
      filen = fdir(1:lf)//'csoislop.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for slow soil carbon, protected humus',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'csoislop',
     >    'instantaneous slow soil carbon protected humus',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (csoislop, cdummy)
      icount(3) = 1
      call writevar(filen,'csoislop',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, csoislop'
         stop 1
      end if
c
c slow soil carbon, nonprotected humus
c
      filen = fdir(1:lf)//'csoislon.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for slow soil carbon, nonprotected humus',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'csoislon',
     >    'instantaneous slow soil carbon nonprotected humus',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (csoislon, cdummy)
      icount(3) = 1
      call writevar(filen,'csoislon',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, csoislon'
         stop 1
      end if
c
c passive soil carbon
c
      filen = fdir(1:lf)//'csoipas.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for passive soil carbon',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'csoipas',
     >    'instantaneous passive soil carbon','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (csoipas, cdummy)
      icount(3) = 1
      call writevar(filen,'csoipas',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, csoipas'
         stop 1
      end if
c
c available plant nitrogen 
c
      filen = fdir(1:lf)//'aplantn.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for available plant nitrogen soil pool',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'aplantn',
     >    'instantaneous plant available soil nitrogen ','kg/m^2',
     >    3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (aplantn, cdummy)
      icount(3) = 1
      call writevar(filen,'aplantn',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, aplantn'
         stop 1
      end if
c
c
c growing degree days
c
      filen = fdir(1:lf)//'gdd0.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for growing degree days above 0 deg_C',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'gdd0',
     >    'instantaneous growing degree days above 0 deg_C',
     >    'days degC',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (gdd0, cdummy)
      icount(3) = 1
      call writevar(filen,'gdd0',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, gdd0'
         stop 1
      end if
c
c growing degree days - wheat
c
      filen = fdir(1:lf)//'gdd0c.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for wheat - April 1 - Sept 30 GDD above 0 C',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'gdd0c',
     >    'instantaneous growing degree days above 0 deg_C',
     >    'days degC',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (gdd0c, cdummy)
      icount(3) = 1
      call writevar(filen,'gdd0c',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, gdd0c'
         stop 1
      end if
c
      filen = fdir(1:lf)//'gdd5.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for growing degree days above 5 deg_C',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'gdd5',
     >    'instantaneous growing degree days above 5 deg_C',
     >    'days degC',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (gdd5, cdummy)
      icount(3) = 1
      call writevar(filen,'gdd5',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, gdd5'
         stop 1
      end if
c
      filen = fdir(1:lf)//'gdd8.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for growing degree days above 8 deg_C',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'gdd8',
     >    'instantaneous growing degree days above 8 deg_C',
     >    'days degC',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (gdd8, cdummy)
      icount(3) = 1
      call writevar(filen,'gdd8',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, gdd8'
         stop 1
      end if
c
      filen = fdir(1:lf)//'gdd10.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for growing degree days above 10 deg_C',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'gdd10',
     >    'instantaneous growing degree days above 10 deg_C',
     >    'days degC',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (gdd10, cdummy)
      icount(3) = 1
      call writevar(filen,'gdd10',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, gdd10'
         stop 1
      end if
c
c coldest monthly temperature
c
      filen = fdir(1:lf)//'tc.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for coldest monthly temperature',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'tc',
     >    'instantaneous coldest monthly temperature',
     >    'degC',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (tc, cdummy)
      icount(3) = 1
      call writevar(filen,'tc',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, tc'
         stop 1
      end if
c
c warmest monthly temperature
c
      filen = fdir(1:lf)//'tw.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for warmest monthly temperature',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'tw',
     >    'instantaneous warmest monthly temperature',
     >    'degC',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (tw, cdummy)
      icount(3) = 1
      call writevar(filen,'tw',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, tw'
         stop 1
      end if
c
c ice content of puddles
c
      filen = fdir(1:lf)//'wipud.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for ice content of puddles',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'wipud',
     >    'instantaneous ice content of puddles',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (wipud, cdummy)
      icount(3) = 1
      call writevar(filen,'wipud',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, wipud'
         stop 1
      end if
c
c liquid content of puddles
c
      filen = fdir(1:lf)//'wpud.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for liquid content of puddles',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'wpud',
     >    'instantaneous liquid water content of puddles',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (wpud, cdummy)
      icount(3) = 1
      call writevar(filen,'wpud',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, wpud'
         stop 1
      end if
c
c annual accumulated growing degree days for bud burst, upper canopy
c
      filen = fdir(1:lf)//'agddu.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for upper canopy growing degree days',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'agddu',
     >    'instantaneous growing degree days uc',
     >    'days degC',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (agddu, cdummy)
      icount(3) = 1
      call writevar(filen,'agddu',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, agddu'
         stop 1
      end if
c
c annual accumulated growing degree days for bud burst, lower canopy
c
      filen = fdir(1:lf)//'agddl.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for lower canopy growing degree days',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'agddl',
     >    'instantaneous growing degree days lc',
     >    'days degC',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (agddl, cdummy)
      icount(3) = 1
      call writevar(filen,'agddl',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, agddl'
         stop 1
      end if
c
c cold-phenology trigger for trees
c
      filen = fdir(1:lf)//'tempu.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for cold phenology trigger for trees',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'tempu',
     >    'cold phenology trigger for trees',
     >    '_',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (tempu, cdummy)
      icount(3) = 1
      call writevar(filen,'tempu',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, tempu'
         stop 1
      end if
c
c cold-phenology trigger for shrubs
c
      filen = fdir(1:lf)//'templs.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for cold phenology trigger for shrubs',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'templs',
     >    'cold phenology trigger for shrubs',
     >    '_',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (templs, cdummy)
      icount(3) = 1
      call writevar(filen,'templs',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, templs'
         stop 1
      end if
c
c fraction of green vegetation in C3 grasses
c
      filen = fdir(1:lf)//'greenfracl3.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for fraction of green vegetation in C3 grasses',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'greenfracl3',
     >    'fraction of green vegetation in C3 grasses',
     >    '_',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (greenfracl3, cdummy)
      icount(3) = 1
      call writevar(filen,'greenfracl3',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, greenfracl3'
         stop 1
      end if
c
c fraction of green vegetation in C4 grasses
c
      filen = fdir(1:lf)//'greenfracl4.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for fraction of green vegetation in C4 grasses',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'greenfracl4',
     >    'fraction of green vegetation in C4 grasses',
     >    '_',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (greenfracl4, cdummy)
      icount(3) = 1
      call writevar(filen,'greenfracl4',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, greenfracl4'
         stop 1
      end if
c
c average annual air temperature (purely from climatology)
c
      filen = fdir(1:lf)//'Tavgann.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for average annual air temperature',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'Tavgann',
     >    'average annual air temperature',
     >    'deg C',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (Tavgann, cdummy)
      icount(3) = 1
      call writevar(filen,'Tavgann',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, Tavgann'
         stop 1
      end if
c
c average annual precipitation (purely from climatology)
c
      filen = fdir(1:lf)//'PPTavgann.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for average annual precipitation',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'PPTavgann',
     >    'average annual precipitation',
     >    'mm/day',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (PPTavgann, cdummy)
      icount(3) = 1
      call writevar(filen,'PPTavgann',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, PPTavgann'
         stop 1
      end if
c
c 10-day average daily air temperature
c
      filen = fdir(1:lf)//'a10td.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for 10-day average daily air T',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'a10td',
     >    '10-day average daily air T',
     >    'K',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (a10td, cdummy)
      icount(3) = 1
      call writevar(filen,'a10td',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, a10td'
         stop 1
      end if
c
c 10-day average soil temperature average of first layer
c
      filen = fdir(1:lf)//'a10ts.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for 10-day average daily soil T',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'a10ts',
     >    '10-day average daily air T',
     >    'K',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (a10ts, cdummy)
      icount(3) = 1
      call writevar(filen,'a10ts',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, a10ts'
         stop 1
      end if
c
c 10-day average canopy photosynthesis rate - broadleaf
c
      filen = fdir(1:lf)//'a10ancub.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for 10-day average canopy photosynth. rate,'
     >    //' broadleaf',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'a10ancub',
     >    '10-day average canopy photosynth. rate, broadleaf',
     >    'mol_co2 m-2 s-1',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (a10ancub, cdummy)
      icount(3) = 1
      call writevar(filen,'a10ancub',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, a10ancub'
         stop 1
      end if
ccc
ccc 10-day average canopy photosynthesis rate - conifer
ccc
cc      filen = fdir(1:lf)//'a10ancuc.nc'
cc      if (nyears .le. 2) then
cc         call inifile(idies,filen,
cc     >    'restart file for 10-day average canopy photosynth. rate,'
cc     >    //' conifer',
cc     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
cc     >    latscale,'','none',
cc     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
cc         dimnames(3) = 'time'
cc         call inivar(idies,'a10ancuc',
cc     >    '10-day average canopy photosynth. rate, conifer',
cc     >    'mol_co2 m-2 s-1',3,dimnames,OCEAN,istat)
cc         call endini(idies,istat)
cc      end if
cc      call vec2arr (a10ancuc, cdummy)
cc      icount(3) = 1
cc      call writevar(filen,'a10ancuc',istart,icount,cdummy,ftime,
cc     > tweight,tdate,istat)
cc      if (istat .ne. 0) then
cc         write(*,*) 'ERROR in wrestart, a10ancuc'
cc         stop 1
cc      end if
c
c 10-day average canopy photosynthesis rate - shrubs
c
      filen = fdir(1:lf)//'a10ancls.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for 10-day average canopy photosynth. rate,'
     >    //' shrubs',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'a10ancls',
     >    '10-day average canopy photosynth. rate, shrubs',
     >    'mol_co2 m-2 s-1',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (a10ancls, cdummy)
      icount(3) = 1
      call writevar(filen,'a10ancls',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, a10ancls'
         stop 1
      end if
c
c 10-day average canopy photosynthesis rate - c4 grasses
c
      filen = fdir(1:lf)//'a10ancl4.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for 10-day average canopy photosynth. rate,'
     >    //' c4 grasses',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'a10ancl4',
     >    '10-day average canopy photosynth. rate, c4 grasses',
     >    'mol_co2 m-2 s-1',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (a10ancl4, cdummy)
      icount(3) = 1
      call writevar(filen,'a10ancl4',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, a10ancl4'
         stop 1
      end if
c
c 10-day average canopy photosynthesis rate - c3 grasses
c
      filen = fdir(1:lf)//'a10ancl3.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for 10-day average canopy photosynth. rate,'
     >    //' c3 grasses',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'a10ancl3',
     >    '10-day average canopy photosynth. rate, c3 grasses',
     >    'mol_co2 m-2 s-1',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (a10ancl3, cdummy)
      icount(3) = 1
      call writevar(filen,'a10ancl3',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, a10ancl3'
         stop 1
      end if
c
c 10-day average canopy scaling parameter - upper canopy
c
      filen = fdir(1:lf)//'a10scalparamu.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for 10-day average canopy scaling parameter,'
     >    //' upper canopy',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'a10scalparamu',
     >    '10-day average canopy scaling parameter, upper canopy',
     >    '_',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (a10scalparamu, cdummy)
      icount(3) = 1
      call writevar(filen,'a10scalparamu',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, a10scalparamu'
         stop 1
      end if
c
c 10-day average canopy scaling parameter - lower canopy
c
      filen = fdir(1:lf)//'a10scalparaml.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for 10-day average canopy scaling parameter,'
     >    //' lower canopy',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'a10scalparaml',
     >    '10-day average canopy scaling parameter, lower canopy',
     >    '_',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (a10scalparaml, cdummy)
      icount(3) = 1
      call writevar(filen,'a10scalparaml',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, a10scalparaml'
         stop 1
      end if
c
c 10-day average daylight - upper canopy
c
      filen = fdir(1:lf)//'a10daylightu.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for 10-day average daylight,'
     >    //' upper canopy',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'a10daylightu',
     >    '10-day average daylight, upper canopy',
     >    'W m-2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (a10daylightu, cdummy)
      icount(3) = 1
      call writevar(filen,'a10daylightu',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, a10daylightu'
         stop 1
      end if
c
c 10-day average daylight - lower canopy
c
      filen = fdir(1:lf)//'a10daylightl.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for 10-day average daylight,'
     >    //' lower canopy',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'a10daylightl',
     >    '10-day average daylight, lower canopy',
     >    'W m-2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (a10daylightl, cdummy)
      icount(3) = 1
      call writevar(filen,'a10daylightl',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, a10daylightl'
         stop 1
      end if
c
c 11-day average surface soil temperature
c
      filen = fdir(1:lf)//'a11soiltd.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for 11-day average surface soil temperature,'
     >    //' surface soil',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'a11soiltd',
     >    '11-day average surface soil temperature',
     >    'deg K',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (a11soiltd, cdummy)
      icount(3) = 1
      call writevar(filen,'a11soiltd',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, a11soiltd'
         stop 1
      end if
c
c 3-day average daily minimum air temperature
c
      filen = fdir(1:lf)//'a3tdmin.nc'
      if (nyears .le. 2) then
         call inifile(idies,filen,
     >    'restart file for 3-day average daily minimum air temperature,'
     >    //' air temperature',
     >    'ibis wrestart',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'a3tdmin',
     >    '3-day average daily minimum air temperature',
     >    'deg K',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (a3tdmin, cdummy)
      icount(3) = 1
      call writevar(filen,'a3tdmin',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wrestart, a3tdmin'
         stop 1
      end if
c
      return
c
      end
c
c
c ---------------------------------------------------------------------
      subroutine wdaily (nday, iyear, iyear0, ndayr, jday, irstyear, irestart)
c ---------------------------------------------------------------------
c
c writes out daily files
c
c ---------------------------------------------------------------------
c
      use netcdf_wrappers
      use comgrid
      use compar
      use comsum
      use comveg
      use comwork
      use comsoi
      use comatm
      use comcrop
      use comnitr
c
      implicit none
c
c Arguments
c
      integer nday,   ! number of days run since iyear0
     >        ndayr,  ! 
     >        jday,
     >        irstyear,
     >        iyear,  ! this calendar year
     >        iyear0, ! very first year ever for this run sequence
     >        irestart  ! 0: normal mode; 1: restart mode
c
      integer NF_PUT_ATT_TEXT,   ! netcdf function
     >        NF_GLOBAL          ! netcdf function
c
c local variables
c
      integer  mstep,         ! this time step for netcdf file
     >        i,j,k,l,m,      ! loop indice
     >        idies,          ! netcdf file indice
     >        istat           ! netcdf error flag

      integer istart(4), icount(4) ! for writing vars
      real    pindex(npft)         ! index used for pfts and canopies
      real    dindex(nsoilay)      ! index used for soil layers
c
      character*10 cdate        ! date to use in history attribute in files
      character*10 tdate        ! character date for time step
      character*13 canopies(2)  ! canopy definitions
      character*21 tunits       ! time units
      character*80 dimnames(4)  ! names of dimensions for vars
      character*80 pftdef(npft) ! plant functional type defs (not currently used)
      character*80 filen        ! file name
c
      real ftime(1),              ! real form of nday  dummy (1)
     >     tweight(1)             ! number of days in daily average = 1. dummy (1)

      real:: dummy_vals3rd(1)     !dummy for *.nc 
c
c ---------------------------------------------------------------------
c
      data istart / 1,1,1,1 /,
     >     icount / nlon,nlat,1,1 /
      icount(1) = nlonsub
      icount(2) = nlatsub
c
c     current time value, step, time weight
c
c The following has been edited to allow a restart with no initial
c daily output file (ie. output file is created upon restart)
c
      ftime = nday
      if (irestart .eq. 0) then
        mstep = nday
      else
        mstep = nday - ndayr
      end if
      tweight = 1.
c
c tdate is julian day (3 char) followed by iyear (4 char)
c
      tdate='0000000'//char(0)//char(0)//char(0)
      if (nday .lt. 10) then
         write(tdate(3:3),'(i1)') nday
      else if (nday .lt. 100) then
         write(tdate(2:3),'(i2)') nday
      else
         write(tdate(1:3),'(i3)') nday
      end if
      if (iyear .ge. 1000) then
         write(tdate(4:7),'(i4)') iyear
      else if (iyear .lt. 10) then
         write(tdate(7:7),'(i1)') iyear
      else if (iyear .lt. 100) then
         write(tdate(6:7),'(i2)') iyear
      else
         write(tdate(5:7),'(i3)') iyear
      end if
c
c first time only
c
      if (mstep.eq.1) then
c
c initialize snow layer indices, pft names, etc
c
c         call date(cdate)
c
c time units is days since Dec 31 of the year before iyear0
c
         tunits = 'days since 0000-12-31'
         write(tunits(12:15),'(i4)') iyear0-1
c
         dimnames(1) = 'longitude'
         dimnames(2) = 'latitude'
c        dimnames(3) is set for each variable seperately
c        dimnames(4) not used in this subr for now, but define for future use
         dimnames(4) = 'time'

c
c define plant functional types, canopies with indices
c
         do i = 1, npft
           pindex(i) = i
         enddo
c and with characters
         pftdef(1) = 'trbrevtr - tropical broadleaf evergreen trees'
     >    //char(0)
         pftdef(2) = 
     >    'trbrdetr - tropical broadleaf drought-deciduous trees'
     >    //char(0)
         pftdef(3) = 
     >    'wtbrevtr - warm-temperate broadleaf evergreen trees'
     >    //char(0)
         pftdef(4) = 'tecoevtr - temperate conifer evergreen trees'
     >    //char(0)
         pftdef(5) = 
     >    'tebrdetr - temperate broadleaf cold-deciduous trees'//char(0)
         pftdef(6) = 'bocoevtr - boreal conifer evergreen trees'
     >    //char(0)
         pftdef(7) = 
     >    'bocodetr - boreal conifer cold-deciduous trees'//char(0)
         pftdef(8) = 
     >    'bobrdetr - boreal broadleaf cold-deciduous trees'//char(0)
         pftdef(9) = 'evsh - evergreen shrubs'//char(0)
         pftdef(10) = 'desh - deciduous shrubs'//char(0)
         pftdef(11) = 'c4gr - warm (c4) grasses'//char(0)
         pftdef(12) = 'c3gr - cool (c3) grasses'//char(0)
         pftdef(13) = 'c3 crop - soybean'//char(0)
         pftdef(14) = 'c4 crop - corn'//char(0)
         pftdef(15) = 'c3 crop - wheat'//char(0)

         canopies(1) = 'lower canopy'//char(0)
         canopies(2) = 'upper canopy'//char(0)
c
      end if
cc
cc dummy variable example, 3-d - copy & modify for new variable.
cc
c      filen = trim(myprocDir)//'output/daily/dummyv.nc'
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'daily average dummyv',
c     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'dummyv','average dummyv',
c     >    'dummyvs-units',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (addummyv, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'dummyv',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wdaily, dummyv'
c         stop 1
c      end if
c
c daily windspeed
cc
c      filen = trim(myprocDir)//'output/daily/ud.nc'
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'daily average windspeed',
c     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'ud','average windspeed',
c     >    'm/s',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (adud, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'ud',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wdaily, ud'
c         stop 1
c      end if
c
c
c rainfall
c
      filen = trim(myprocDir)//'output/daily/rain.nc'
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'daily average rainfall',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'rain','average rainfall',
     >    'mm/day',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (adrain, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'rain',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, rain'
         stop 1
      end if
c
c cloudiness
c
c      filen = trim(myprocDir)//'output/daily/cloud.nc'
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'daily average cloudiness',
c     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,
c     >    latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'cloud','average cloudiness',
c     >    '%',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (cloud, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'cloud',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wdaily, cloud'
c         stop 1
c      end if
c
c rh
c
c      filen = trim(myprocDir)//'output/daily/rh.nc'
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'daily average relative humidity',
c     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,
c     >    latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'rh','average rh',
c     >    '%',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (adrh, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'rh',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wdaily, rh'
c         stop 1
c      end if
c
c snowfall
c
c      filen = trim(myprocDir)//'output/daily/snow.nc'
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'daily average snowfall',
c     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,
c     >    latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'snow','average snowfall',
c     >    'mm/day',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (adsnow, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'snow',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wdaily, snow'
c         stop 1
c      end if
cc
cc aet
cc
c      filen = trim(myprocDir)//'output/daily/aet.nc'
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'daily average aet',
c     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'aet','average aet',
c     >    'mm/day',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (adaet, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'aet',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wdaily, aet'
c         stop 1
c      end if
cc
cc trunoff
cc
      filen = trim(myprocDir)//'output/daily/trunoff.nc'
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'daily average total runoff',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'trunoff','average total runoff',
     >    'mm/day',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (adtrunoff, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'trunoff',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, trunoff'
         stop 1
      end if
c
c srunoff
c
      filen = trim(myprocDir)//'output/daily/srunoff.nc'
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'daily average surface runoff',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'srunoff','average surface runoff',
     >    'mm/day',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (adsrunoff, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'srunoff',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, srunoff'
         stop 1
      end if
c
c drainage
c
      filen = trim(myprocDir)//'output/daily/drainage.nc'
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'daily average drainage',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'drainage','average drainage',
     >    'mm/day',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (addrainage, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'drainage',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, drainage'
         stop 1
      end if
cc
cc wsoi
cc
      filen = trim(myprocDir)//'output/daily/wsoi.nc'
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'daily average soil moisture',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'wsoi','average soil moisture',
     >    'fraction',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (adwsoi, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'wsoi',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, wsoi'
         stop 1
      end if

c wisoi
c
      filen = trim(myprocDir)//'output/daily/wisoi.nc'
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'daily average soil ice',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'wisoi','average soil ice',
     >    'fraction',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (adwisoi, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'wisoi',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, wisoi'
         stop 1
      end if
c
cc snod
cc
      filen = trim(myprocDir)//'output/daily/snod.nc'
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'daily average snow depth',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'snod','average snow depth',
     >    'meters',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (adsnod, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'snod',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, snod'
         stop 1
      end if
cc
cc snof
cc
      filen = trim(myprocDir)//'output/daily/snof.nc'
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'daily average snow fraction',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'snof','average snow fraction',
     >    'fraction',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (adsnof, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'snof',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, snof'
         stop 1
      end if
cc
cc co2ratio
cc
      filen = trim(myprocDir)//'output/daily/co2ratio.nc'
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'daily average ratio of root to total soil co2 flux',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'co2ratio','average co2 ratio',
     >    'fraction',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (adco2ratio, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'co2ratio',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, co2ratio'
         stop 1
      end if
cc
cc co2mic
cc
      filen = trim(myprocDir)//'output/daily/co2mic.nc'
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'daily flux of carbon due to soil microbe co2 flux',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'co2mic','soil microbe carbon flux',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (adco2mic, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'co2mic',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, co2mic'
         stop 1
      end if
c
c nitrate concentrations in solution 
c
      do 54 l = 1, nsoilay
        dindex(l) = l  
 54   continue
    
      filen = trim(myprocDir)//'output/daily/nconc.nc'
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'solute nitrate concentration',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'level','soil depth','cm',nsoilay,dindex,'',
     >    tunits,'gregorian',istat)
cc add global attribute to define pfts with text, use netcdf low-level command
c         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'depth-definition',
c     >    nsoilay*80,soildef)
         dimnames(3) = 'level'
         call inivar(idies,'csoln','nitrate concentration each layer',
     >    'layer',4,dimnames,OCEAN,istat)
      end if
      do 12 m = 1, nsoilay 
         call vec2arr (csoln(1,m), cdummy((m-1)*nlonsub*nlatsub + 1))
 12   continue
      istart(3) = 1
      istart(4) = mstep
      icount(3) = nsoilay 
      call writevar(filen,'csoln',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, csoln'
         stop 1
      end if

c crop daily lai 
c 'lev' replaced 'pft' so we could look at output with GrADS
c
c
      filen = trim(myprocDir)//'output/daily/plai.nc'

      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'daily pft lai',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'level','plant functional type','pft',npft,pindex,'',
     >    tunits,'gregorian',istat)
c add global attribute to define pfts with text, use netcdf low-level command
         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',
     >    npft*80,pftdef)
         dimnames(3) = 'level'
         call inivar(idies,'plai','plai for each pft',
     >    'm2/m2',4,dimnames,OCEAN,istat)
         call endini(idies,istat)
       end if
      do 220 k = 1, npft
         call vec2arr (plai(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 220    continue
      istart(3) = 1
      istart(4) = mstep
      icount(3) = npft
      call writevar(filen,'plai',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, plai'
         stop 1
      end if
c
c biomass, by pft, and category 
c additions by CJK 3/5/99 for crop output variables
c 'lev' replaced 'pft' so we could look at output with GrADS
c
c      filen = trim(myprocDir)//'output/daily/biomass.nc'
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'daily biomass of carbon',
c     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,
c     >    'level','plant fuctional type','pft',npft,pindex,'',
c     >    tunits,'gregorian',istat)
c add global attribute to define pfts with text, use netcdf low-level command
c         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',
c     >    npft*80,pftdef)
c         dimnames(3) = 'level'
c         call inivar(idies,'biomass','biomass for each pft',
c     >    'kg/m^2',4,dimnames,OCEAN,istat)
c         dimnames(3) = 'level'
c         call inivar(idies,'cbiol','total biomass of leaves
c     >     for each  pft',
c     >    'kg/m^2',4,dimnames,OCEAN,istat)
c         dimnames(3) = 'level'
c         call inivar(idies,'cbior','total fine root biomass
c     >     for each  pft',
c     >    'kg/m^2',4,dimnames,OCEAN,istat)
c         dimnames(3) = 'level'
c         call inivar(idies,'cbiog','total grain biomass for
c     >     each crop pft',
c     >    'kg/m^2',4,dimnames,OCEAN,istat)
c         dimnames(3) = 'level'
c         call inivar(idies,'cbios','total stem biomass for
c     >     each crop pft',
c     >    'kg/m^2',4,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c
c      call write4dvar (biomass, npft, mstep, filen, 'biomass', ftime, tweight, tdate)
c      call write4dvar (cbiol, npft, mstep, filen, 'cbiol', ftime, tweight, tdate)
c      call write4dvar (cbior, npft, mstep, filen, 'cbior', ftime, tweight, tdate)
c      call write4dvar (cbiog, npft, mstep, filen, 'cbiog', ftime, tweight, tdate)
c      call write4dvar (cbios, npft, mstep, filen, 'cbios', ftime, tweight, tdate)
c
c
cc gdd8
cc
c      filen = trim(myprocDir)//'output/daily/gdd8.nc'
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'accumulation of gdd base 8 C',
c     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'gdd8','accumulation of gdd base 8 C',
c     >    'deg C',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (gdd8this, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'gdd8',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wdaily, gdd8'
c         stop 1
c      end if
cc
c phenology data 
c additions by CJK 12/22/99 for crop output variables
c
c      filen = trim(myprocDir)//'output/daily/phenology.nc'
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'phenological variables',
c     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,
c     >    'pft','plant fuctional type','_',npft,pindex,'',
c     >    tunits,'gregorian',istat)
c add global attribute to define pfts with text, use netcdf low-level command
c         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',
c     >    npft*80,pftdef)
c         dimnames(3) = 'pft'
c         call inivar(idies,'gddplant','accumulated gdd since
c     >     each crop pft planting ',
c     >    'gdd',4,dimnames,OCEAN,istat)
c         dimnames(3) = 'pft'
c         call inivar(idies,'hui','heat unit index',
c     >    'dimensionless',4,dimnames,OCEAN,istat)
c         dimnames(3) = 'pft'
c         call inivar(idies,'croplive','switch to tell if
c     >     crop planting has taken place',
c     >    'binary',4,dimnames,OCEAN,istat)
c      endif
c
c      call write4dvar (gddplant, npft, mstep, filen, 'gddplant', ftime, tweight, tdate)
c      call write4dvar (hui, npft, mstep, filen, 'hui', ftime, tweight, tdate)
c      call write4dvar (croplive, npft, mstep, filen, 'croplive', ftime, tweight, tdate)
c
c 
c bottom and top of vegetation canopies
c added CJK 3/11/99 to analyze crops/grasses
c 'lev' replaced 'canopy' so we could look at output with GrADS
c
c      filen = trim(myprocDir)//'output/daily/zcanopy.nc'
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'daily height of vegetation canopies',
c     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,
c     >    'level','lev','canopy',2,pindex,'up',
c     >    tunits,'gregorian',istat)
c add global attribute to define canopies with text, use netcdf low-level com
c         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'canopy_def',
c     >    26+5,'1='//canopies(1)//' 2='//canopies(2))
c         dimnames(3) = 'level'
c         call inivar(idies,'zbot',
c     >    'bottom heights of lower and upper canopies',
c     >    'meters',4,dimnames,OCEAN,istat)
c         call inivar(idies,'ztop',
c     >    'top heights of lower and upper canopies',
c     >    'meters',4,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      do 20 k = 1, 2
c         call vec2arr (zbot(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
c 20   continue
c      istart(3) = 1
c      istart(4) = mstep
c      icount(3) = 2
c      call writevar(filen,'zbot',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wdaily, zbot'
c         stop 1
c      end if
c      do 25 k = 1, 2
c         call vec2arr (ztop(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
c 25   continue
c      istart(3) = 1
c      istart(4) = mstep
c      icount(3) = 2
c      call writevar(filen,'ztop',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wdaily, ztop'
c         stop 1
c      end if
c
cc templs 
cc
c      filen = trim(myprocDir)//'output/daily/templs.nc'
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'index based on gdd for growth/sensecence',
c     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'templs','index based on gdd for growth/sensescence',
c     >    'dimensionless',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (templs, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'templs',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wdaily, templs'
c         stop 1
c      end if
c
c bottom and top height of lower and upper canopies
c
c      filen = trim(myprocDir)//'output/daily/zcanopy.nc'
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'daily height of vegetation canopies',
c     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,
c     >    'canopy','canopy','_',2,pindex,'up',
c     >    tunits,'gregorian',istat)
c add global attribute to define canopies with text, use netcdf low-level com
c         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'canopy_def',
c     >    26+5,'1='//canopies(1)//' 2='//canopies(2))
c         dimnames(3) = 'canopy'
c         call inivar(idies,'zbot',
c     >    'bottom heights of lower and upper canopies',
c     >    'meters',4,dimnames,OCEAN,istat)
c         call inivar(idies,'ztop',
c     >    'top heights of lower and upper canopies',
c     >    'meters',4,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      do 20 k = 1, 2
c         call vec2arr (zbot(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
c 20   continue
c      istart(3) = 1
c      istart(4) = mstep
c      icount(3) = 2
c      call writevar(filen,'zbot',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wdaily, zbot'
c         stop 1
c      end if
c      do 25 k = 1, 2
c         call vec2arr (ztop(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
c 25   continue
c      istart(3) = 1
c      istart(4) = mstep
c      icount(3) = 2
c      call writevar(filen,'ztop',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wdaily, ztop'
c         stop 1
c      end if
c
c upper and lower canopy daily lai
c
      filen = trim(myprocDir)//'output/daily/laicanopy.nc'
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'daily lai-leaf area index of upper and lower vegetation canopies',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    '','','',1,dummy_vals3rd,'',
     >    tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'laiu',
     >    'daily lai of upper canopy',
     >    'm2/m2',3,dimnames,OCEAN,istat)
       call inivar(idies,'lail',
     >    'daily lai of lower canopy',
     >    'm2/m2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if

      call vec2arr (lai(1,1), cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'laiu',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
        write(*,*) 'ERROR in wdaily, laiu'
         stop 1
      end if
      call vec2arr (lai(1,2), cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'lail',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, lail'
         stop 1
      end if
c
c lower canopy daily green & brown lai
c
      filen = trim(myprocDir)//'output/daily/lailower.nc'
      if (mstep .eq. 1) then
        call inifile(idies,filen,
     >    'daily green & brown lai-leaf area index of lower vegetation canopy',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    '','','',1,dummy_vals3rd,'',
     >    tunits,'gregorian',istat)
        dimnames(3) = 'time'
        call inivar(idies,'glail',
     >    'daily green lai of lower canopy',
     >    'm2/m2',3,dimnames,OCEAN,istat)
        call inivar(idies,'blail',
     >    'daily brown lai of lower canopy',
     >    'm2/m2',3,dimnames,OCEAN,istat)
        call endini(idies,istat)
      end if

      work(1:npoi) = fl(:)*lai(:,1)*greenfracl(:)  ! green lai of lower canopy
      call vec2arr (work(1:npoi), cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'glail',istart,icount,cdummy,ftime,
     >  tweight,tdate,istat)
      if (istat .ne. 0) then
        write(*,*) 'ERROR in wdaily, glail'
        stop 1
      end if

      work(1:npoi) = fl(:)*lai(:,1)*(1. - greenfracl(:))  ! brown lai of lower canopy
      call vec2arr (work(1:npoi), cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'blail',istart,icount,cdummy,ftime,
     >  tweight,tdate,istat)
      if (istat .ne. 0) then
        write(*,*) 'ERROR in wdaily, blail'
        stop 1
      end if
c
c crop daily lai 
c 'lev' replaced 'pft' so we could look at output with GrADS
c
c      filen = trim(myprocDir)//'output/daily/plai.nc'
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'daily pft lai',
c     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,
c     >    'level','plant functional type','pft',npft,pindex,'',
c     >    tunits,'gregorian',istat)
cc add global attribute to define pfts with text, use netcdf low-level command
c         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',
c     >    npft*80,pftdef)
c         dimnames(3) = 'level'
c         call inivar(idies,'plai','plai for each pft',
c     >    'm2/m2',4,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      do 222 k = 1, npft
c         call vec2arr (plai(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
c 222    continue
c      istart(3) = 1
c      istart(4) = mstep
c      icount(3) = npft
c      call writevar(filen,'plai',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wdaily, plai'
c         stop 1
c      end if
c
c frac - fraction of canopy occupied by pft  MM 4/29/10
c 'lev' replaced 'pft' so we could look at output with GrADS
c
c      filen = trim(myprocDir)//'output/daily/frac.nc'
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'fraction of canopy occupied by pft',
c     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,
c     >    'level','plant functional type','pft',npft,pindex,'',
c     >    tunits,'gregorian',istat)
ccc add global attribute to define pfts with text, use netcdf low-level command
c         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',
c     >    npft*80,pftdef)
c         dimnames(3) = 'level'
c         call inivar(idies,'frac','fraction of canopy occupied by pft',
c     >    'fraction',4,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      do 222 k = 1, npft
c         call vec2arr (frac(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
c 222    continue
c      istart(3) = 1
c      istart(4) = mstep
c      icount(3) = npft
c      call writevar(filen,'frac',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wdaily, frac'
c         stop 1
c      end if
c
c total plant nitrogen uptake and stress factor by pft
c CJK 2-24-00 for crop output variables
c 'lev' replaced 'pft' so we could look at output with GrADS
c
      filen = trim(myprocDir)//'output/daily/plantn.nc'
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'daily crop nitrogen variabiles',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'level','plant functional type','pft',npft,pindex,'',
     >    tunits,'gregorian',istat)
c add global attribute to define pfts with text, use netcdf low-level command
         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',
     >    npft*80,pftdef)
         dimnames(3) = 'level'
         call inivar(idies,'stressn','nitrogen stress factor for each pft',
     >    'dimensionless',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'totnuptake','total nitrogen uptake 
     >     for each  pft',
     >    'kg/m^2',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'fnleaf','nitrogen fraction in leaf 
     >     for each  pft',
     >    'fraction',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'fnstem','nitrogen fraction in stem 
     >     for each  pft',
     >    'fraction',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'fnroot','nitrogen fraction in root 
     >     for each  pft',
     >    'fraction',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'fngrain','nitrogen fraction in grain 
     >     for each  pft',
     >    'fraction',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'fnplant','nitrogen fraction in entire 
     >     plant for each  pft',
     >    'fraction',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'tnplant','nitrogen total in plant 
     >     for each  pft',
     >    'fraction',4,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if

      call write4dvar (stressn, npft, mstep, filen, 'stressn', ftime, tweight, tdate)
      call write4dvar (totnuptake, npft, mstep, filen, 'totnuptake', ftime, tweight, tdate)
      call write4dvar (fnleaf, npft, mstep, filen, 'fnleaf', ftime, tweight, tdate)
      call write4dvar (fnstem, npft, mstep, filen, 'fnstem', ftime, tweight, tdate)
      call write4dvar (fnroot, npft, mstep, filen, 'fnroot', ftime, tweight, tdate)
      call write4dvar (fngrain, npft, mstep, filen, 'fngrain', ftime, tweight, tdate)
      call write4dvar (fnplant, npft, mstep, filen, 'fnplant', ftime, tweight, tdate)
      call write4dvar (tnplant, npft, mstep, filen, 'tnplant', ftime, tweight, tdate)
c
c
c daily rate of nitrate leaching 
c 
      filen = trim(myprocDir)//'output/daily/leachr.nc'
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'rate of nitrogen leached from profile (kg N m-2 y-1)',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'dnileach','rate of nitrogen leaching',
     >    'kg n ha-1 y-1',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (dnileach, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'dnileach',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, dnileach'
         stop 1
      end if
c
c instantaneous plant available inorganic nitrogen 
c 
      filen = trim(myprocDir)//'output/daily/aplantn.nc'
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'inorganic plant available inorganic nitrogen (kg N m-2 )',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'aplantn','inorganic plant available nitrogen',
     >    'kg n m-2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (aplantn, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'aplantn',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, aplantn'
         stop 1
      end if
c
c cumulative nitrogen leaching 
c 
      filen = trim(myprocDir)//'output/daily/cumnleach.nc'
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'accumulated nitrogen leached from profile at specific depth (kg N m-2)',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'ftot','cumulative nitrogen leaching',
     >    'kg n ha-1 ',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (ftot, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'ftot',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, ftot'
         stop 1
      end if
c
c cumulative drainage 
c 
      filen = trim(myprocDir)//'output/daily/cumdrainage.nc'
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'accumulated drainage from profile at specific depth (mm h20)',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'drntot','cumulative drainage',
     >    'mm h20 ',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (drntot, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'drntot',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, drntot'
         stop 1
      end if
c
c daily plant transpiration  
c 
      filen = trim(myprocDir)//'output/daily/trans.nc'
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'daily rate of canopy transpiration (upper and lower) (mm/day)',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'adtrans','daily canopy transpiration',
     >    'mm/day',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (adtrans, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'adtrans',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, adtrans'
         stop 1
      end if
c
c 
c daily total evaporation  
c 
      filen = trim(myprocDir)//'output/daily/evap.nc'
      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'daily rate of evaporation (mm/day)',
     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'adevap','daily evaporation',
     >    'mm/day',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (adevap, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'adevap',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wdaily, adevap'
         stop 1
      end if
c
c 
c daily ratio of transpiration to total ET  
c 
c      filen = trim(myprocDir)//'output/daily/tratio.nc'
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'daily ratio of trans to total ET (dimensionless)',
c     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'adtratio','ratio of trans to total ET',
c     >    'dimensionless',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (adtratio, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'adtratio',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wdaily, adtratio'
c         stop 1
c      end if
cc
c adprobfire - daily probability of fire
c
c      filen = trim(myprocDir)//'output/daily/adprobfire.nc'
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'average daily probability of fire',
c     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,
c     >    latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'adprobfire','average daily probability of fire',
c     >    'unitless',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c       end if
c      call vec2arr (adprobfire, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'adprobfire',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wdaily, adprobfire'        
c         stop 1
c      end if
c
c totsoi - total soil moisture (wisoi(i,1) + wsoi(i,1))
c
c      filen = trim(myprocDir)//'output/daily/totsoi.nc'
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'total soil moisture of first layer',
c     >    'ibis wdaily',cdate,nlonsub,lonscale,nlatsub,
c     >    latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'totsoi','total soil moisture of first layer',
c     >    'unitless',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c       end if
c      call vec2arr (totsoi, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'totsoi',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wdaily, totsoi'
c         stop 1
c      end if
      
      return
      end

c
c ---------------------------------------------------------------------
      subroutine wmonthly (nday, imonth, iyear, iyear0, irstyear, jday, irestart)
c ---------------------------------------------------------------------
c
c writes out monthly files
c
      use comgrid
      use compar
      use comsum
      use comwork
      use comcrop
c
      implicit none
c
c Arguments
c
      integer nday,    ! number of days run since iyear0
     >        jday,
     >        imonth,  ! this month
     >        iyear,   ! this calendar year
     >        irstyear,   ! this calendar year
     >        iyear0,  ! very first year ever for this run
     >        irestart ! 0: normal mode; 1: restart mode
c
c Local variables
c
      integer mstep,   ! time step in netcdf file
     >        idies,   ! netcdf file indice
     >        istat    ! netcdf error flag
c
      integer istart(4), icount(4) ! for writing vars
c
      character*10 cdate       ! date to use in history attribute in files
      character*10 tdate       ! character date for time step
      character*21 tunits      ! units for time
      character*80 dimnames(4) ! names of dimensions for vars
      character*80 filen       ! file name
      character*3  chmon(12)   ! month abbreviations
      character*4   num,num2
c
      real ftime(1),              ! real form of nday dummy (1)
     >     tweight(1)             ! number of days in monthly average, dummy(1)

      real:: dummy_vals3rd(1)     !dummy for *.nc 

      character*2 num_lay         !char for soil layer number
      integer:: ilayer
c
c ---------------------------------------------------------------------
c
      data chmon / 'JAN','FEB','MAR','APR','MAY','JUN',
     >             'JUL','AUG','SEP','OCT','NOV','DEC' /
      data istart / 1,1,1,1 /,
     >     icount / nlon,nlat,1,1 /
      icount(1) = nlonsub
      icount(2) = nlatsub
c
c     current time value, step, time weight
c
      ftime = nday
      if (irestart .eq. 0) then
        mstep = 12*(iyear-iyear0) + imonth
      else
        mstep = 12*(iyear-irstyear) + imonth
      end if
      tweight = ndaypm(imonth)

      if(imonthout_glo == 2) then
        mstep = imonth 
        write(num,'(i4.4)')  iyear
        write(num2,'(i4.4)') iyear-1
      end if

c
c tdate is this month (3 char), followed by this year (4 char)
c
      tdate=chmon(imonth)//'0000'//char(0)//char(0)//char(0)
      if (iyear .ge. 1000) then
         write(tdate(4:7),'(i4)') iyear
      else if (iyear .lt. 10) then
         write(tdate(7:7),'(i1)') iyear
      else if (iyear .lt. 100) then
         write(tdate(6:7),'(i2)') iyear
      else
         write(tdate(5:7),'(i3)') iyear
      end if
c
c first time only
c
      if (mstep.eq.1) then
c
c initialize snow layer indices, pft names, etc
c
c         call date(cdate)
c
c time units is days since Dec 31 of the year before iyear0
c
         tunits = 'days since 0000-12-31'
         if(imonthout_glo == 2) then
           write(tunits(12:15),'(i4)') iyear-1
         else
           write(tunits(12:15),'(i4)') iyear0-1
         end if
c
         dimnames(1) = 'longitude'
         dimnames(2) = 'latitude'
c        dimnames(3) is set for each variable seperately
c        dimnames(4) not used in this subr for now, but define for future use
         dimnames(4) = 'time'
      end if !on if mstep == 1
cc
cc dummy variable example, 3-d - copy & modify for new variable.
cc
c      filen = trim(myprocDir)//'output/monthly/dummyv.nc'
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'monthly average dummyv',
c     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'dummyv','average dummyv',
c     >    'dummyvs-units',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (amdummyv, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'dummyv',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wmonthly, dummyv'
c         stop 1
c      end if
c
c temperature
c
c       if(imonthout_glo == 2) then
c         filen = trim(myprocDir)//'output/monthly/'//num//'_monthly_temp.nc'
c       else
c        filen = trim(myprocDir)//'output/monthly/temp.nc'
c       end if
c
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'monthly average air temperature',
c     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
c     >    latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'temp','average air temperature',
c     >    'C',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (amtemp, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'temp',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wmonthly, temp'
c         stop 1
c      end if
c
c rainfall
c
      if(imonthout_glo == 2) then
        filen = trim(myprocDir)//'output/monthly/'//num//'_monthly_rain.nc'
      else
        filen = trim(myprocDir)//'output/monthly/rain.nc'
      end if

      if (mstep .eq. 1) then
        call inifile(idies,filen,
     >    'monthly average rainfall',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
        dimnames(3) = 'time'
        call inivar(idies,'rain','average rainfall',
     >    'mm/day',3,dimnames,OCEAN,istat)
        call endini(idies,istat)
      end if

      call vec2arr (amrain, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'rain',istart,icount,cdummy,ftime,
     >  tweight,tdate,istat)
      if (istat .ne. 0) then
        write(*,*) 'ERROR in wmonthly, rain'
        stop 1
      end if
c
c cloudiness
c
c       if(imonthout_glo == 2) then
c         filen = trim(myprocDir)//'output/monthly/'//num//'_monthly_cloud.nc'
c       else
c         filen = trim(myprocDir)//'output/monthly/cloud.nc'
c       end if
c
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'monthly average cloudiness',
c     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
c     >    latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'cloud','average cloudiness',
c     >    '%',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (amcloud, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'cloud',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wmonthly, cloud'
c         stop 1
c      end if
c
c rh
c
c      if(imonthout_glo == 2) then 
c        filen = trim(myprocDir)//'output/monthly/'//num//'_monthly_rh.nc'
c      else
c        filen = trim(myprocDir)//'output/monthly/rh.nc'
c      end if
c
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'monthly average relative humidity',
c     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
c     >    latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'rh','average rh',
c     >    '%',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (amrh, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'rh',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wmonthly, rh'
c         stop 1
c      end if
c
c snowfall
c
c      if(imonthout_glo == 2) then
c        filen = trim(myprocDir)//'output/monthly/'//num//'_monthly_snow.nc'
c      else
c        filen = trim(myprocDir)//'output/monthly/snow.nc'
c      end if
c
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'monthly average snowfall',
c     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
c     >    latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'snow','average snowfall',
c     >    'mm/day',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (amsnow, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'snow',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wmonthly, snow'
c         stop 1
c      end if
c
c specific humidity
c
c      if(imonthout_glo == 2) then
c        filen = trim(myprocDir)//'output/monthly/'//num//'_monthly_qa.nc'
c      else
c        filen = trim(myprocDir)//'output/monthly/qa.nc'
c      end if
c
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'monthly average specific humidity',
c     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
c     >    latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'qa','average specific humidity',
c     >    'kg-h2o/kg-air',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (amqa, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'qa',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wmonthly, qa'
c         stop 1
c      end if
c
c evapotranspiration
c
      if(imonthout_glo == 2) then
        filen = trim(myprocDir)//'output/monthly/'//num//'_monthly_aet.nc'
      else
        filen = trim(myprocDir)//'output/monthly/aet.nc'
      end if

      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'monthly average aet',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'aet','average evapotranspiration',
     >    'mm/day',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'trans','average transpiration',
     >    'mm/day',3,dimnames,OCEAN,istat)
         !Y.Li add evap for THMB
         call inivar(idies,'evap','average evaporation',
     >    'mm/day',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amaet, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'aet',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, aet'
         stop 1
      end if
      call vec2arr (amtrans, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'trans',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, trans'
         stop 1
      end if
      call vec2arr ((amaet-amtrans), cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'evap',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, evap'
         stop 1
      end if
c
c transETratio
c
c      if(imonthout_glo == 2) then
c        filen = trim(myprocDir)//'output/monthly/'//num//'_monthly_transETratio.nc'
c      else
c        filen = trim(myprocDir)//'output/monthly/transETratio.nc'
c      end if
c
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'monthly average transpiration:ET ratio',
c     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
c     >    latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'transET','average trans:ET ratio',
c     >    'dimensionless',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (amtratio, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'transET',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wmonthly, transET'
c         stop 1
c      end if
c
c trunoff, srunoff, drainage
c
      if(imonthout_glo == 2) then
        filen = trim(myprocDir)//'output/monthly/'//num//'_monthly_runoff.nc'
      else
        filen = trim(myprocDir)//'output/monthly/runoff.nc'
      end if

      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'monthly average total runoff',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'trunoff','average total runoff',
     >    'mm/day',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'srunoff','average surface runoff',
     >    'mm/day',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'drainage','average drainage',
     >    'mm/day',3,dimnames,OCEAN,istat)

         do ilayer = 1,nsoilay
           write(num_lay,'(i2)') ilayer
           call inivar(idies,'drainage_layer'//adjustl(num_lay),'average drainage for layer '//adjustl(num_lay),
     >    'mm/day',3,dimnames,OCEAN,istat)

         end do !on do ilayer = 1,nsoilay

         dimnames(3) = 'time'
         call inivar(idies,'totnleach','average total n leaching',
     >    'kg/ha/day',3,dimnames,OCEAN,istat)

         do ilayer = 1,nsoilay
           write(num_lay,'(i2)') ilayer
           call inivar(idies,'totnleach_layer'//adjustl(num_lay),'average total n leaching for layer '//adjustl(num_lay),
     >      'kg/ha/day',3,dimnames,OCEAN,istat)
         end do !on do ilayer = 1,nsoilay

         dimnames(3) = 'time'
         call inivar(idies,'no3leach','average no3 leaching',
     >    'kg/ha/day',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amtrunoff, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'trunoff',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, trunoff'
         stop 1
      end if
      call vec2arr (amsrunoff, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'srunoff',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, srunoff'
         stop 1
      end if
      call vec2arr (amdrainage, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'drainage',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, drainage'
         stop 1
      end if

      !drainage for each layer
      do ilayer = 1, nsoilay

        write(num_lay,'(i2)') ilayer
        call vec2arr (amdrainage_layer(:,ilayer), cdummy)
        istart(3) = mstep
        icount(3) = 1
        call writevar(filen,'drainage_layer'//adjustl(num_lay),istart,icount,cdummy,ftime,
     >                tweight,tdate,istat)
        if (istat .ne. 0) then
           write(*,*) 'ERROR in wmonthly, drainage_layer'//adjustl(num_lay)
           stop 1
        end if

      end do !on ilayer

      call vec2arr (amtotnleach(:,9), cdummy)  !the original one at layer 9; keep this for compatability
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'totnleach',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, totnleach'
         stop 1
      end if

      !totnleach_layer* for each layer
      do ilayer = 1, nsoilay

        write(num_lay,'(i2)') ilayer
        call vec2arr (amtotnleach(:,ilayer), cdummy)  !the original one at layer 9; keep this for compatability
        istart(3) = mstep
        icount(3) = 1
        call writevar(filen,'totnleach_layer'//adjustl(num_lay),istart,icount,cdummy,ftime,
     >                tweight,tdate,istat)
        if (istat .ne. 0) then
           write(*,*) 'ERROR in wmonthly, totnleach_layer'//adjustl(num_lay)
           stop 1
        end if

      end do !on do ilayer

      call vec2arr (amno3leach, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'no3leach',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, no3leach'
         stop 1
      end if
c
c soil temperature
c
      if(imonthout_glo == 2) then
        filen = trim(myprocDir)//'output/monthly/'//num//'_monthly_tsoi.nc'
      else
        filen = trim(myprocDir)//'output/monthly/tsoi.nc'
      end if

      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'monthly average soil temperature',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'tsoi','average soil temperature',
     >    'degC',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amtsoi, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'tsoi',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, tsoi'
         stop 1
      end if
c
c soil moisture, ice, volumetric water content, plant available water
c
      if(imonthout_glo == 2) then
        filen = trim(myprocDir)//'output/monthly/'//num//'_monthly_wsoi.nc'
      else
        filen = trim(myprocDir)//'output/monthly/wsoi.nc'
      end if

      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'monthly average soil moisture',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'wsoi','average soil moisture',
     >    'fraction',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'wisoi','average soil ice',
     >    'fraction',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'vwc','average volumetric water content',
     >    'fraction',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'awc',
     >    'average plant available water content',
     >    'cm',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amwsoi, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'wsoi',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, wsoi'
         stop 1
      end if
      call vec2arr (amwisoi, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'wisoi',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, wisoi'
         stop 1
      end if
      call vec2arr (amvwc, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'vwc',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, vwc'
         stop 1
      end if
      call vec2arr (amawc, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'awc',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, awc'
         stop 1
      end if
c
c snow depth
c
      if(imonthout_glo == 2) then
        filen = trim(myprocDir)//'output/monthly/'//num//'_monthly_snod.nc'
      else
        filen = trim(myprocDir)//'output/monthly/snod.nc'
      end if

      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'monthly average snow depth',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'snod','average snow depth',
     >    'meters',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amsnod, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'snod',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, snod'
         stop 1
      end if
c
c snow fraction
c
      if(imonthout_glo == 2) then
        filen = trim(myprocDir)//'output/monthly/'//num//'_monthly_snof.nc'
      else
        filen = trim(myprocDir)//'output/monthly/snof.nc'
      end if

      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'monthly average snow fraction',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'snof','average snow fraction',
     >    'm^2/m^3',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amsnof, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'snof',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, snof'
         stop 1
      end if
c
c solar radiation
c
c      if(imonthout_glo == 2) then
c        filen = trim(myprocDir)//'output/monthly/'//num//'_monthly_solar.nc'
c      else
c        filen = trim(myprocDir)//'output/monthly/solar.nc'
c      end if
c
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'monthly average incident solar radiation',
c     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
c     >    latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'solar','average incident solar radiation',
c     >    'W/m^2',3,dimnames,OCEAN,istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'solarnet','average net solar radiation',
c     >    'W/m^2',3,dimnames,OCEAN,istat)         
c         call endini(idies,istat)
c      end if
c      call vec2arr (amsolar, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'solar',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wmonthly, solar'
c         stop 1
c      end if
c      call vec2arr (amsolarnet, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'solarnet',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wmonthly, solarnet'
c         stop 1
c      end if
c
c net radiation
c
c      if(imonthout_glo == 2) then
c        filen = trim(myprocDir)//'output/monthly/'//num//'_monthly_netrad.nc'
c      else
c        filen = trim(myprocDir)//'output/monthly/netrad.nc'
c      end if
c
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'monthly average net radiation',
c     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
c     >    latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'netrad','average net radiation',
c     >    'W/m^2',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (amnetrad, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'netrad',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wmonthly, netrad'
c         stop 1
c      end if
c
c albedo
c
c      if(imonthout_glo == 2) then
c        filen = trim(myprocDir)//'output/monthly/'//num//'_monthly_albedo.nc'
c      else 
c        filen = trim(myprocDir)//'output/monthly/albedo.nc'
c      end if
c
c       if (mstep .eq. 1) then
c        call inifile(idies,filen,
c     >    'monthly average albedo',
c     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'albedo','average albedo',
c     >    'fraction',3,dimnames,OCEAN,istat)
c          call endini(idies,istat)
c       end if
c       call vec2arr (amalbedo, cdummy)
c       istart(3) = mstep
c       icount(3) = 1
c       call writevar(filen,'albedo',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c       if (istat .ne. 0) then
c        write(*,*) 'ERROR in wmonthly, albedo'
c        stop 1
c       end if
c
c downward and upward infrared radiation
c
c      if(imonthout_glo == 2) then
c        filen = trim(myprocDir)//'output/monthly/'//num//'_monthly_ir.nc'
c      else
c        filen = trim(myprocDir)//'output/monthly/ir.nc'
c      end if
c
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'monthly average infrared radiation',
c     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
c     >    latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'irdown','average downward IR',
c     >    'W/m^2',3,dimnames,OCEAN,istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'irup','average upward IR',
c     >    'W/m^2',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (amirdown, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'irdown',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wmonthly, irdown'
c         stop 1
c      end if
c      call vec2arr (amirup, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'irup',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wmonthly, irup'
c         stop 1
c      end if
c
c sensible heat flux
c
      if(imonthout_glo == 2) then
        filen = trim(myprocDir)//'output/monthly/'//num//'_monthly_sens.nc'
      else
        filen = trim(myprocDir)//'output/monthly/sens.nc'
      end if

      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'monthly average sensible heat flux',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'sens','average sensible heat flux',
     >    'W/m^2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amsens, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'sens',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, sens'
         stop 1
      end if
c
c latent heat flux
c
      if(imonthout_glo == 2) then
        filen = trim(myprocDir)//'output/monthly/'//num//'_monthly_latent.nc'
      else
        filen = trim(myprocDir)//'output/monthly/latent.nc'
      end if

      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'monthly average latent heat flux',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'latent','average latent heat flux',
     >    'W/m^2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amlatent, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'latent',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, latent'
         stop 1
      end if
c
c leaf area index upper and lower
c
      if(imonthout_glo == 2) then
        filen = trim(myprocDir)//'output/monthly/'//num//'_monthly_lai.nc'
      else
        filen = trim(myprocDir)//'output/monthly/lai.nc'
      end if

      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >   'monthly average leaf area index for upper and lower canopies',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'laiu','average lai upper canopy',
     >    'fraction',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'lail','average lai lower canopy',
     >    'fraction',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'glail','average green lai lower canopy',
     >    'fraction',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'blail','average brown lai lower canopy',
     >    'fraction',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amlaiu, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'laiu',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, laiu'
         stop 1
      end if
      call vec2arr (amlail, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'lail',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, lail'
         stop 1
      end if
      call vec2arr (amglail, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'glail',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, glail'
         stop 1
      end if
      call vec2arr (amblail, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'blail',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, blail'
         stop 1
      end if
c
c total net primary productivity
c
      if(imonthout_glo == 2) then
        filen = trim(myprocDir)//'output/monthly/'//num//'_monthly_npptot.nc'
      else
        filen = trim(myprocDir)//'output/monthly/npptot.nc'
      end if

      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'monthly total net primary productivity of carbon',
     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'npptot','total npp of carbon',
     >    'kg m-2 month-1',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (amnpptot, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'npptot',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, npptot'
         stop 1
      end if
c
c co2ratio
c
c      if(imonthout_glo == 2) then
c        filen = trim(myprocDir)//'output/monthly/'//num//'_monthly_co2ratio.nc'
c      else
c        filen = trim(myprocDir)//'output/monthly/co2ratio.nc'
c      end if
c
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'monthly ratio of root respiration to total soil respiration',
c     >    'ibis wmonthly',cdate,nlonsub,lonscale,nlatsub,
c     >    latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'co2ratio',
c     >    'ratio of root to soil respiration',
c     >    'fraction',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (amco2ratio, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'co2ratio',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wmonthly, co2ratio'
c         stop 1
c      end if
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine wyearly (nday,iyear,iyear0, irstyear, irestart)
c ---------------------------------------------------------------------
c
c writes out yearly files
c
      use netcdf_wrappers
      use comgrid
      use compar
      use comsum
      use comveg
      use comwork
      use comsoi
      use comatm
      use comcrop
      use comnitr
c
      implicit none
c
c Arguments
c
      integer nday,   ! number of days run since iyear0
     >        iyear,  ! this calendar year
     >        irstyear,  ! 
     >        iyear0, ! very first year ever for this run
     >        irestart  ! 0: normal mode; 1: restart mode

c
c local variables
c
      integer mstep,  ! time for this year step
     >        idies,  ! netcdf file indice
     >        istat,  ! netcdf error flag
     >        i,k,l,m ! loop indices        
c
      integer istart(4), icount(4) ! for writing vars
      real    pindex(npft)         ! index used for pfts and canopies
      real    dindex(nsoilay)      ! index used for soil layers
c
      character*10 cdate        ! date to use in history attribute in files
      character*10 tdate        ! character date for time step
      character*13 canopies(2)  ! canopy definitions
      character*21 tunits       ! time units
      character*80 dimnames(4)  ! names of dimensions for vars
      character*80 pftdef(npft) ! plant functional type defs (not currently used)
      character*80 filen        ! file name
      character*4   num,num2
c
      real ftime(1),            ! real form of nday
     >     tweight(1)           ! number of days in yearly average

      real:: dummy_vals3rd(1)     !dummy for *.nc 
c
c External
c 
      integer NF_PUT_ATT_TEXT,   ! netcdf function
     >        NF_GLOBAL          ! netcdf function
c
c ---------------------------------------------------------------------
c
      data istart / 1,1,1,1 /,
     >     icount / nlon,nlat,1,1 /
      icount(1) = nlonsub
      icount(2) = nlatsub
c
c current time value and step: make ftime Jan 1 of this year 
c instead of Dec 31
c
      ftime = nday - ndaypy + 1
      tweight = ndaypy
      if (irestart .eq. 0) then
        mstep = iyear - iyear0 + 1
      else
        mstep = iyear - irstyear + 1
      end if

      if(iyearout_glo == 2) then
        mstep = 1 
        write(num,'(i4.4)')  iyear
        write(num2,'(i4.4)') iyear-1
      end if
c
c tdate is ANN (3 char) followed by this year (4 char)
c
      tdate='ANN0000'//char(0)//char(0)//char(0)
      if (iyear .ge. 1000) then
         write(tdate(4:7),'(i4)') iyear
      else if (iyear .lt. 10) then
         write(tdate(7:7),'(i1)') iyear
      else if (iyear .lt. 100) then
         write(tdate(6:7),'(i2)') iyear
      else
         write(tdate(5:7),'(i3)') iyear
      end if
c
c first time only
c
      if (mstep .eq. 1) then
c
c         call date(cdate)
c
c time units is days since Dec 31 of the year before iyear0
c
         tunits = 'days since 0000-12-31'
         if(iyearout_glo == 2) then
           write(tunits(12:15),'(i4)') iyear-1
         else
           write(tunits(12:15),'(i4)') iyear0-1
         end if
c
c dimension names
c
         dimnames(1) = 'longitude'
         dimnames(2) = 'latitude'
c        dimnames(3) is set for each variable seperately
         dimnames(4) = 'time'
c
c define plant functional types, canopies with indices
c
         do i = 1, npft
           pindex(i) = i
         enddo
c and with characters
         pftdef(1) = 'trbrevtr - tropical broadleaf evergreen trees'
     >    //char(0)
         pftdef(2) = 
     >    'trbrdetr - tropical broadleaf drought-deciduous trees'
     >    //char(0)
         pftdef(3) = 
     >    'wtbrevtr - warm-temperate broadleaf evergreen trees'
     >    //char(0)
         pftdef(4) = 'tecoevtr - temperate conifer evergreen trees'
     >    //char(0)
         pftdef(5) = 
     >    'tebrdetr - temperate broadleaf cold-deciduous trees'//char(0)
         pftdef(6) = 'bocoevtr - boreal conifer evergreen trees'
     >    //char(0)
         pftdef(7) = 
     >    'bocodetr - boreal conifer cold-deciduous trees'//char(0)
         pftdef(8) = 
     >    'bobrdetr - boreal broadleaf cold-deciduous trees'//char(0)
         pftdef(9) = 'evsh - evergreen shrubs'//char(0)
         pftdef(10) = 'desh - deciduous shrubs'//char(0)
         pftdef(11) = 'c4gr - warm (c4) grasses'//char(0)
         pftdef(12) = 'c3gr - cool (c3) grasses'//char(0)
         pftdef(13) = 'c3 crop - soybean'//char(0)
         pftdef(14) = 'c4 crop - corn'//char(0)
         pftdef(15) = 'c3 crop - wheat'//char(0)

         canopies(1) = 'lower canopy'//char(0)
         canopies(2) = 'upper canopy'//char(0)
c
      end if
cc
cc dummy variable example, 4-d var whose 3rd dim is a character dim (pft)
cc copy and modify for a new variable
cc
c      filen = trim(myprocDir)//'output/yearly/dummyv.nc'
c      if (mstep .eq. 1) then
c         call inifilec(idies,filen,
c     >    'annual dummyv',
c     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,
c     >    'pft','plant fuctional type','_',npft,80,pftdef,
c     >    tunits,'gregorian',istat)
c         dimnames(3) = 'pft'
c         call inivar(idies,'dummyv','dummyv for each pft',
c     >    'dummyvs-units',4,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      do 5 k = 1, npft
c         call vec2arr (aydummyv(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
c 5    continue
c      istart(3) = 1
c      istart(4) = mstep
c      icount(3) = npft
c      call writevar(filen,'dummyv',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, dummyv'
c         stop 1
c      end if
cc
cc dummy variable example, 3-d
cc
c      filen = trim(myprocDir)//'output/yearly/dummyv.nc'
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'dummyv','annual dummyv',cdate,nlonsub,lonscale,
c     >    nlatsub,latscale,'','none','none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'dummyv','annual dummyv',
c     >    'dummyvs-units',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (aydummyv, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'dummyv',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, dummyv'
c         stop 1
c      end if
c
c net primary productivity, by pft and total
c
      if(iyearout_glo == 2) then
        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_npp.nc'
      else
        filen = trim(myprocDir)//'output/yearly/npp.nc'
      end if

      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual net primary productivity of carbon',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'pft','plant fuctional type','_',npft,pindex,'',
     >    tunits,'gregorian',istat)
c add global attribute to define pfts with text, use netcdf low-level command
         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',
     >    npft*80,pftdef)
         dimnames(3) = 'pft'
         call inivar(idies,'npp','npp of carbon for each pft',
     >    'kg m-2 year-1',4,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'npptot','total npp',
     >    'kg m-2 year-1',3,dimnames,OCEAN,istat)
         call inivar(idies,'anpptot','total above-ground npp',
     >    'kg m-2 year-1',3,dimnames,OCEAN,istat)
         dimnames(3) = 'pft'
         call inivar(idies,'gpp','gpp of carbon for each pft',
     >    'kg m-2 year-1',4,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'gpptot','total gpp',
     >    'kg m-2 year-1',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      do 5 k = 1, npft
         call vec2arr (aynpp(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 5    continue
      istart(3) = 1
      istart(4) = mstep
      icount(3) = npft
      call writevar(filen,'npp',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, npp'
         stop 1
      end if
      call vec2arr (aynpptot, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'npptot',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, npptot'
         stop 1
      end if
      call vec2arr (ayanpptot, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'anpptot',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, anpptot'
         stop 1
      end if
      do k = 1, npft
         call vec2arr (aygpp(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
      end do
      istart(3) = 1
      istart(4) = mstep
      icount(3) = npft
      call writevar(filen,'gpp',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, gpp'
         stop 1
      end if
      call vec2arr (aygpptot, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'gpptot',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, gpptot'
         stop 1
      end if
c
c evapotranspiration
c
      if(iyearout_glo == 2) then
        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_aet.nc'
      else
        filen = trim(myprocDir)//'output/yearly/aet.nc'
      end if

      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual average evapotranspiration',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'aet','average evapotranspiration',
     >    'mm/year',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'trans','average transpiration',
     >    'mm/year',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (ayaet, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'aet',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, aet'
         stop 1
      end if
      call vec2arr (aytrans, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'trans',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, trans'
         stop 1
      end if
c
c trunoff, srunoff, drainage, rratio, tratio
c
      if(iyearout_glo == 2) then
        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_runoff.nc'
      else
        filen = trim(myprocDir)//'output/yearly/runoff.nc'
      end if

      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual runoff',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'trunoff','total runoff',
     >    'mm/year',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'srunoff','surface runoff',
     >    'mm/year',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'totirrig','total irrigation water applied',
     >    'mm/year',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'ayprcp','total annual precipitation',
     >    'mm/year',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'drainage','drainage',
     >    'mm/year',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'rratio','average runoff ratio',
     >    'fraction',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'tratio','average transpiration ratio',
     >    'fraction',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (aytrunoff, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'trunoff',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, trunoff'
         stop 1
      end if
      call vec2arr (aysrunoff, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'srunoff',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, srunoff'
         stop 1
      end if
      call vec2arr (totirrig, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'totirrig',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, totirrig'
         stop 1
      end if
      call vec2arr (ayprcp, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'ayprcp',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, ayprcp'
         stop 1
      end if
      call vec2arr (aydrainage, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'drainage',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, drainage'
         stop 1
      end if
      call vec2arr (ayrratio, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'rratio',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, rratio'
         stop 1
      end if
      call vec2arr (aytratio, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'tratio',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, tratio'
         stop 1
      end if
c
c gdd information - climate and current year 
c
c      if(iyearout_glo == 2) then
c        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_gdd.nc'
c      else
c        filen = trim(myprocDir)//'output/yearly/gdd.nc'
c      end if
c
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'annual gdd info',
c     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,
c     >    latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'gdd8clim','climate gdd8',
c     >    'degrees C',3,dimnames,OCEAN,istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'gdd10clim','climate gdd10',
c     >    'degrees C',3,dimnames,OCEAN,istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'gdd0clim','climate gdd0',
c     >    'degrees C',3,dimnames,OCEAN,istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'currgdd8','current year gdd8',
c     >    'degrees C',3,dimnames,OCEAN,istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'currgdd10','current year gdd10',
c     >    'degrees C',3,dimnames,OCEAN,istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'currgdd0','current year gdd0',
c     >    'degrees C',3,dimnames,OCEAN,istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'gddfzcorn','current year gdd between freeze events',
c     >    'degrees C',3,dimnames,OCEAN,istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'gddfzsoy','current year gdd between freeze events',
c     >    'degrees C',3,dimnames,OCEAN,istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'gsdays','number of days between freeze events',
c     >    'days',3,dimnames,OCEAN,istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'iniday','last day in spring for freeze',
c     >    'day',3,dimnames,OCEAN,istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'endday','first day in fall for freeze',
c     >    'day',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (gdd8, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'gdd8clim',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, gdd8clim'
c         stop 1
c      end if
c      call vec2arr (gdd10, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'gdd10clim',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, gdd10clim'
c         stop 1
c      end if
c      call vec2arr (gdd0c, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'gdd0clim',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, gdd0clim'
c         stop 1
c      end if
c      call vec2arr (gdd8this, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'currgdd8',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, currgdd8'
c         stop 1
c      end if
c      call vec2arr (gdd10this, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'currgdd10',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, currgdd10'
c         stop 1
c      end if
c      call vec2arr (gdd0cthis, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'currgdd0',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, currgdd0'
c         stop 1
c      end if
c      call vec2arr (gddfzcorn, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'gddfzcorn',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, gddfzcorn'
c         stop 1
c      end if
c      call vec2arr (gddfzsoy, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'gddfzsoy',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, gddfzsoy'
c         stop 1
c      end if
c      call vec2arr (gsdays, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'gsdays',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, gsdays'
c         stop 1
c      end if
c      call vec2arr (iniday, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'iniday',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, iniday'
c         stop 1
c      end if
c      call vec2arr (endday, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'endday',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, endday'
c         stop 1
c      end if
c
c soil moisture, soil ice, volumetric water content, plant available water
c
      if(iyearout_glo == 2) then
        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_wsoi.nc'
      else
        filen = trim(myprocDir)//'output/yearly/wsoi.nc'
      end if

      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual average soil moisture',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'wsoi','average soil moisture',
     >    'fraction',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'wisoi','average soil ice',
     >    'fraction',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'vwc','average volumetric water content',
     >    'fraction',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'awc','average volumetric water content',
     >    'cm',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (aywsoi, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'wsoi',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, wsoi'
         stop 1
      end if
      call vec2arr (aywisoi, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'wisoi',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, wisoi'
         stop 1
      end if
      call vec2arr (ayvwc, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'vwc',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, vwc'
         stop 1
      endif
      call vec2arr (ayawc, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'awc',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, awc'
         stop 1
      endif
c
c soil temperature
c
      if(iyearout_glo == 2) then
        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_tsoi.nc'
      else
        filen = trim(myprocDir)//'output/yearly/tsoi.nc'
      end if

      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual average soil temperature',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'tsoi','average soil temperature',
     >    'degrees C',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (aytsoi, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'tsoi',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, tsoi'
         stop 1
      endif
cc

c extreme min air temp over year
c
c      if(iyearout_glo == 2) then
c        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_tcmin.nc'
c      else
c        filen = trim(myprocDir)//'output/yearly/tcmin.nc'
c      end if
c
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'minimum daily temp over year',
c     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,
c     >    latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'tcmin','minimum daily temp over year',
c     >    'degrees C',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (tcmin, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'tcmin',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c        write(*,*) 'ERROR in wyearly, tcmin'
c         stop 1
c      endif

cc solar radiation
cc
c      if(iyearout_glo == 2) then
c        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_solar.nc'
c      else
c        filen = trim(myprocDir)//'output/yearly/solar.nc'
c      end if
c
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'annual average solar incident radiation',
c     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'solar','average solar radiation',
c     >    'W/m^2',3,dimnames,OCEAN,istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'solarnet','average net solar radiation',
c     >    'W/m^2',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (aysolar, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'solar',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, solar'
c         stop 1
c      end if
c      call vec2arr (aysolarnet, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'solarnet',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, solarnet'
c         stop 1
c      end if
cc
cc albedo
cc
c      if(iyearout_glo == 2) then
c        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_albedo.nc'
c      else
c        filen = trim(myprocDir)//'output/yearly/albedo.nc'
c      end if
c
c     if (mstep .eq. 1) then
c        call inifile(idies,filen,
c    >    'annual average albedo',
c    >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c    >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c        dimnames(3) = 'time'
c        call inivar(idies,'albedo','average albedo',
c    >    'fraction',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c     end if
c     call vec2arr (ayalbedo, cdummy)
c     istart(3) = mstep
c     icount(3) = 1
c     call writevar(filen,'albedo',istart,icount,cdummy,ftime,
c    > tweight,tdate,istat)
c     if (istat .ne. 0) then
c        write(*,*) 'ERROR in wyearly, albedo'
c        stop 1
c     end if
cc
cc upward and downward infrared radiation
cc
c      if(iyearout_glo == 2) then
c        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_ir.nc'
c      else
c        filen = trim(myprocDir)//'output/yearly/ir.nc'
c      end if
c
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'annual average solar infrared radiation',
c     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'irdown','average downward ir',
c     >    'W/m^2',3,dimnames,OCEAN,istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'irup','average upward ir',
c     >    'W/m^2',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (ayirdown, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'irdown',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, irdown'
c         stop 1
c      end if
c      call vec2arr (ayirup, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'irup',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, irup'
c         stop 1
c      end if
c
c sensible heat flux
c
c      if(iyearout_glo == 2) then
c        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_sens.nc'
c      else
c        filen = trim(myprocDir)//'output/yearly/sens.nc'
c      end if
c
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'annual average sensible heat flux',
c     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,
c     >    latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'sens','average sensible heat flux',
c     >    'W/m^2',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (aysens, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'sens',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, sens'
c         stop 1
c      end if
c
c firelength - length of fire season
c
c      if(iyearout_glo == 2) then
c        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_firelength.nc'
c      else
c        filen = trim(myprocDir)//'output/yearly/firelength.nc'
c      end if
c
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'length of fire season',
c     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,
c     >    latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'firelength','length of fire season',
c     >    'Unitless',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (firelength, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'firelength',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c        write(*,*) 'ERROR in wyearly, firelength'
c         stop 1
c      end if
c
c firefrac - fraction of grid cell lost to fire
c
c      if(iyearout_glo == 2) then
c        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_firefrac.nc'
c      else
c        filen = trim(myprocDir)//'output/yearly/firefrac.nc'
c      end if
c
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'fraction of grid cell lost to fire',
c     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,
c     >    latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'firefrac','fraction of grid cell lost to fire',
c     >    'Unitless',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (firefrac, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'firefrac',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, firefrac'
c         stop 1
c      end if

cc
cc latent heat flux
cc
c      if(iyearout_glo == 2) then
c        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_latent.nc'
c      else
c        filen = trim(myprocDir)//'output/yearly/latent.nc'
c      end if
c
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'annual average latent heat flux',
c     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'latent','average latent heat flux',
c     >    'W/m^2',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (aylatent, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'latent',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, latent'
c         stop 1
c      end if
cc
cc net radiation
cc
c      if(iyearout_glo == 2) then
c        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_netrad.nc'
c      else
c        filen = trim(myprocDir)//'output/yearly/netrad.nc'
c      end if
c
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'annual average net radiation',
c     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'netrad','average net radiation',
c     >    'W/m^2',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (aynetrad, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'netrad',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, netrad'
c         stop 1
c      end if
c
c lai, by pft, total upper canopy, total lower canopy
c
      if(iyearout_glo == 2) then
        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_plai.nc'
      else
        filen = trim(myprocDir)//'output/yearly/plai.nc'
      end if

      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual leaf area index',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'pft','plant fuctional type','_',npft,pindex,'',
     >    tunits,'gregorian',istat)
c add global attribute to define pfts with text, use netcdf low-level command
         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',
     >    npft*80,pftdef)
         dimnames(3) = 'pft'
         call inivar(idies,'plai','leaf area index for each pft',
     >    'fraction',4,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'totlaiu','total lai for upper canopy',
     >    'fraction',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'totlail','total lai for lower canopy',
     >    'fraction',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      do 10 k = 1, npft
         call vec2arr (plai(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 10   continue
      istart(3) = 1
      istart(4) = mstep
      icount(3) = npft
      call writevar(filen,'plai',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, plai'
         stop 1
      end if
      call vec2arr (totlaiu, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'totlaiu',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, totlaiu'
         stop 1
      end if
      call vec2arr (totlail, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'totlail',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, totlail'
         stop 1
      end if
c
c biomass, by pft, upper canopy, lower canopy
c additions by CJK 3/5/99 for crop output variables
c this includes looking at total biomass production (not remaining standing)
c aboveground and total, along with biomass in grain and harvest index as a
c function of plant functional type 
c 'lev' replaced 'pft' so we could look at output with GrADS
c
      if(iyearout_glo == 2) then
        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_biomass.nc'
      else
        filen = trim(myprocDir)//'output/yearly/biomass.nc'
      end if

      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual biomass of carbon',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'level','plant fuctional type','pft',npft,pindex,'',
     >    tunits,'gregorian',istat)
c add global attribute to define pfts with text, use netcdf low-level command
         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',
     >    npft*80,pftdef)
         dimnames(3) = 'level'
         call inivar(idies,'biomass','biomass for each pft',
     >    'kg/m^2',4,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'totbiou','total biomass for upper canopy',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'totbiol','total biomass for lower canopy',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
          call endini(idies,istat)
      end if
      do 15 k = 1, npft
         call vec2arr (biomass(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
 15   continue
      istart(3) = 1
      istart(4) = mstep
      icount(3) = npft
      call writevar(filen,'biomass',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, biomass'
         stop 1
      end if
c 
      call vec2arr (totbiou, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'totbiou',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, totbiou'
         stop 1
      end if
      call vec2arr (totbiol, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'totbiol',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, totbiol'
        stop 1
      end if
c
c crops output
c 'lev' replaced 'pft' so we could look at output with GrADS
c
      if(iyearout_glo == 2) then
        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_crops.nc'
      else
        filen = trim(myprocDir)//'output/yearly/crops.nc'
      end if

      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual crop growth variables',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,
     >    'level','crop type','pft',npft,pindex,'',
     >    tunits,'gregorian',istat)
c add global attribute to define pfts with text, use netcdf low-level command
         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',
     >    npft*80,pftdef)
         dimnames(3) = 'level'
         call inivar(idies,'matdate','crop maturity date',
     >    'day of year',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'harvdate','crop harvest date',
     >    'day of year',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'plantdate','crop planting date',
     >    'day of year',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'harvidx','harvest index',
     >    'fraction',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'croplaimx','maximum crop lai',
     >    'm^2/m^2',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'grainn','nitrogen removed by grain',
     >    'kg/ha',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'cropyld','crop yield',
     >    'bu/ac',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'dmyield','crop yield dry matter',
     >    'Mg/ha',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'dmleaf','leaf dry matter',
     >    'Mg/ha',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'dmstem','stem dry matter',
     >    'Mg/ha',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'dmroot','root dry matter',
     >    'Mg/ha',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'dmresidue','aboveground residue dry matter',
     >    'Mg/ha',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'dmcrop','total crop dry matter',
     >    'Mg/ha',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'residuen','total nitrogen in residue',
     >    'kg/ha',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'nconcl','leaf nitrogen concentration',
     >    'percent',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'nconcs','stem nitrogen concentration',
     >    'percent',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'nconcr','root nitrogen concentration',
     >    'percent',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'nconcg','grain nitrogen concentration',
     >    'percent',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'cropn','total crop nitrogen uptake',
     >    'kg/ha',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'cropfixn','nitrogen fixation',
     >    'kg/ha',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'cntops','cn ratio of aboveground residue ',
     >    'dimensionless',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'cnroot','cn ratio of  fineroots ',
     >    'dimensionless',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'fertilizer','fertilizer applied ',
     >    'kg/ha',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'emerge_day','day of year of leaf emergence',
     >    'day',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'emerge_gdd','GDD since planting at the time of leaf emergence, using air temp.',
     >    'degrees C',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'emerge_gddtsoi','soil GDD accumulated since planting at the time of leaf emergence',
     >    'degrees C',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'grainday','grainday',
     >    'day',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'graingdd','GDD since planting at the time of grainfill initiation',
     >    'degrees C',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'matgdd','GDD since planting at the time of maturity',
     >    'degrees C',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'matgrnfraccrop','green fraction at the time of maturity',
     >    'fraction',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'hybrid','hybridgdd',
     >    'degrees C',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'crmclim','CRM rating for climatology',
     >    'days',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'crmact','best CRM rating for this year',
     >    'days',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'crmplant','CRM rating for this year based on planting date',
     >    'days',4,dimnames,OCEAN,istat)
c WJS (08.02.10) Commenting out APAR until absorb calculation is fixed in stats.f (see comment there and todo item)
c         dimnames(3) = 'level'
c         call inivar(idies,'apar','total APAR',
c     >    'MJ/m2',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'fkill','fraction of the season (in a gdd sense) that the crop missed due to freeze kill',
     >      'fraction',4,dimnames,OCEAN,istat)
         dimnames(3) = 'level'
         call inivar(idies,'devfraction',
     >      'fraction of development attained at the time the crop reaches maturity (hui/gddmaturity)',
     >      'fraction',4,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if

      call write4dvar (matdate, npft, mstep, filen, 'matdate', ftime, tweight, tdate)
      call write4dvar (hdate, npft, mstep, filen, 'harvdate', ftime, tweight, tdate)
      call write4dvar (pdate, npft, mstep, filen, 'plantdate', ftime, tweight, tdate)
      call write4dvar (harvidx, npft, mstep, filen, 'harvidx', ftime, tweight, tdate)
      call write4dvar (croplaimx, npft, mstep, filen, 'croplaimx', ftime, tweight, tdate)
      call write4dvar (grainn, npft, mstep, filen, 'grainn', ftime, tweight, tdate)
      call write4dvar (cropyld, npft, mstep, filen, 'cropyld', ftime, tweight, tdate)
      call write4dvar (dmyield, npft, mstep, filen, 'dmyield', ftime, tweight, tdate)
      call write4dvar (dmleaf, npft, mstep, filen, 'dmleaf', ftime, tweight, tdate)
      call write4dvar (dmstem, npft, mstep, filen, 'dmstem', ftime, tweight, tdate)
      call write4dvar (dmroot, npft, mstep, filen, 'dmroot', ftime, tweight, tdate)
      call write4dvar (dmresidue, npft, mstep, filen, 'dmresidue', ftime, tweight, tdate)
      call write4dvar (dmcrop, npft, mstep, filen, 'dmcrop', ftime, tweight, tdate)
      call write4dvar (residuen, npft, mstep, filen, 'residuen', ftime, tweight, tdate)
      call write4dvar (nconcl, npft, mstep, filen, 'nconcl', ftime, tweight, tdate)
      call write4dvar (nconcs, npft, mstep, filen, 'nconcs', ftime, tweight, tdate)
      call write4dvar (nconcr, npft, mstep, filen, 'nconcr', ftime, tweight, tdate)
      call write4dvar (nconcg, npft, mstep, filen, 'nconcg', ftime, tweight, tdate)
      call write4dvar (cropn, npft, mstep, filen, 'cropn', ftime, tweight, tdate)
      call write4dvar (cropfixn, npft, mstep, filen, 'cropfixn', ftime, tweight, tdate)
      call write4dvar (cntops, npft, mstep, filen, 'cntops', ftime, tweight, tdate)
      call write4dvar (cnroot, npft, mstep, filen, 'cnroot', ftime, tweight, tdate)
      call write4dvar (fertinput, npft, mstep, filen, 'fertilizer', ftime, tweight, tdate)
      call write4dvar (gddmaturity, npft, mstep, filen, 'hybrid', ftime, tweight, tdate)
      call write4dvar (crmclim, npft, mstep, filen, 'crmclim', ftime, tweight, tdate)
      call write4dvar (crmact, npft, mstep, filen, 'crmact', ftime, tweight, tdate)
      call write4dvar (crmplant, npft, mstep, filen, 'crmplant', ftime, tweight, tdate)
      call write4dvar (emerge_day, npft, mstep, filen, 'emerge_day', ftime, tweight, tdate)
      call write4dvar (emerge_gdd, npft, mstep, filen, 'emerge_gdd', ftime, tweight, tdate)
      call write4dvar (emerge_gddtsoi, npft, mstep, filen, 'emerge_gddtsoi', ftime, tweight, tdate)
      call write4dvar (grainday, npft, mstep, filen, 'grainday', ftime, tweight, tdate)
      call write4dvar (graingdd, npft, mstep, filen, 'graingdd', ftime, tweight, tdate)
      call write4dvar (matgdd, npft, mstep, filen, 'matgdd', ftime, tweight, tdate)
      call write4dvar (matgrnfraccrop, npft, mstep, filen, 'matgrnfraccrop', ftime, tweight, tdate)
c WJS (08.02.10) Commenting out APAR until absorb calculation is fixed in stats.f (see comment there and todo item)
c      call write4dvar (apar, npft, mstep, filen, 'apar', ftime, tweight, tdate, 1.e-06)
      call write4dvar (fkill_gdd_missed, npft, mstep, filen, 'fkill', ftime, tweight, tdate)
      call write4dvar (dev_fraction_attained, npft, mstep, filen, 'devfraction', ftime, tweight, tdate)
c 
c soil texture 
c
c      do 55 l = 1, nsoilay
c        dindex(l) = l  
c 55   continue
c    
c      if(iyearout_glo == 2) then
c        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_soitext.nc'
c      else
c        filen = trim(myprocDir)//'output/yearly/soitext.nc'
c      end if
c
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'soil texture',
c     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,
c     >    'level','soil depth','cm',nsoilay,dindex,'',
c     >    tunits,'gregorian',istat)
c add global attribute to define pfts with text, use netcdf low-level command
c         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'depth-definition',
c     >    nsoilay*80,soildef)
c         dimnames(3) = 'level'
c         call inivar(idies,'soisand','sand fraction for each layer',
c     >    'fraction',4,dimnames,OCEAN,istat)
c         dimnames(3) = 'level'
c         call inivar(idies,'soiclay','clay fraction for each layer',
c     >    'fraction',4,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c       end if
c
c      call write4dvar (soisand, nsoilay, mstep, filen, 'soisand', ftime, tweight, tdate)
c      call write4dvar (soiclay, nsoilay, mstep, filen, 'soiclay', ftime, tweight, tdate)
c
c
c soil carbon: rootbio, totalit, totrlit, totcsoi, totcmic
c
      if(iyearout_glo == 2) then
        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_csoi.nc'
      else
        filen = trim(myprocDir)//'output/yearly/csoi.nc'
      end if

      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual total soil carbon',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'rootbio','total live root biomass carbon',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'totalit','total above ground litter carbon',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'totrlit','total below ground litter carbon',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'totcsoi','total soil carbon w/o litter',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'totcmic','total microbial carbon',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (ayrootbio, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'rootbio',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, rootbio'
         stop 1
      end if
      call vec2arr (ayalit, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'totalit',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, totalit'
         stop 1
      end if
      call vec2arr (ayblit, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'totrlit',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, totrlit'
         stop 1
      end if
      call vec2arr (aycsoi, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'totcsoi',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, totcsoi'
         stop 1
      end if
      call vec2arr (aycmic, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'totcmic',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, totcmic'
         stop 1
      end if
c
c litterfall: falll, fallr, fallw
c
c      if(iyearout_glo == 2) then
c        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_litterfall.nc'
c      else
c        filen = trim(myprocDir)//'output/yearly/litterfall.nc'
c      end if
c
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'annual litterfall',
c     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,
c     >    latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'falll',
c     >    'total annual leaf litterfall',
c     >    'kg/m^2',3,dimnames,OCEAN,istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'fallr',
c     >    'total below ground annual root turnover',
c     >    'kg/m^2',3,dimnames,OCEAN,istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'fallw','total wood litterfall',
c     >    'kg/m^2',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (falll, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'falll',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, falll'
c         stop 1
c      end if
c      call vec2arr (fallr, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'fallr',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, fallr'
c         stop 1
c      end if
c      call vec2arr (fallw, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'fallw',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, fallw'
c         stop 1
c      end if
c
c soil nitrogen: totanlit, totrnlit, totnsoi, nmintot
c and residue carbon : nitrogen ratio for crops
c
c nitrogen leaching
c
c nitrogen uptake by natural vegetation
c
      if(iyearout_glo == 2) then
        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_nsoi.nc'
      else
        filen = trim(myprocDir)//'output/yearly/nsoi.nc'
      end if

      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual total soil nitrogen',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'totanlit',
     >    'total above ground litter nitrogen',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'totrnlit',
     >    'total below ground litter nitrogen',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'totnsoi','total soil organic nitrogen w/o litter',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'nmin','net nitrogen mineralization',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'nimmob','total nitrogen immobilized',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'ndepos','total nitrogen deposition',
     >    'kg/m2',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'nfixnat','nitrogen fixation by natural vegetation',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'totnleach','annual total inorganic nitrogen leaching',
     >    'kg/hectare',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'no3leach','annual total nitrate leaching',
     >    'kg/hectare',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'nbalance','annual inorganic nitrogen balance',
     >    'kg/hectare',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'drainage','annual total drainage',
     >    'mm/year',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'no3conc','flow-weighted mean nitrate concentration',
     >    'mg/liter',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'natvegnup','total annual inorganic nitrogen uptake by natural veg',
     >    'kg/m2/yr',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'insoiln','inorganic, immobile soil n to specified layer',
     >    'kg/ha/day',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'insolutn','inorganic, mobile soil n in solution to specified layer',
     >    'kg/ha/day',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (ayanlit, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'totanlit',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, totanlit'
         stop 1
      end if
      call vec2arr (aybnlit, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'totrnlit',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, totrnlit'
         stop 1
      end if
      call vec2arr (aynsoi, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'totnsoi',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, totnsoi'
         stop 1
      end if
      call vec2arr (aynmintot, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'nmin',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, nmin'
         stop 1
      end if
      call vec2arr (ayimmtot, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'nimmob',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, nimmob'
         stop 1
      end if
      call vec2arr (ydeposn, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'ndepos',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, ndepos'
         stop 1
      end if
      call vec2arr (yfixsoin, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'nfixnat',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, nfixnat'
         stop 1
      end if
      call vec2arr (ftot, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'totnleach',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, totnleach'
         stop 1
      end if
      call vec2arr (yno3leach, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'no3leach',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, no3leach'
         stop 1
      end if
      call vec2arr (snbalance, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'nbalance',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, nbalance'
         stop 1
      end if
      call vec2arr (drntot, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'drainage',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, drainage'
         stop 1
      end if
c
      call vec2arr (concn, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'no3conc',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, no3conc'
         stop 1
      end if
c
      call vec2arr (totnvegn, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'natvegnup',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, natvegnup'
         stop 1
      end if

      call vec2arr (tsnimm, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'insoiln',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, insoiln'
         stop 1
      end if
      call vec2arr (tsnmob, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'insolutn',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wmonthly, insolutn'
         stop 1
      end if
cc
cc total litter
cc
      if(iyearout_glo == 2) then
        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_totlit.nc'
      else
        filen = trim(myprocDir)//'output/yearly/totlit.nc'
      end if

      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual total litter carbon',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'totlit','total litter carbon',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (totlit, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'totlit',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, totlit'
         stop 1
      end if
cc
cc total wood litter
cc
      if(iyearout_glo == 2) then
        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_clitw.nc'
      else
        filen = trim(myprocDir)//'output/yearly/clitw.nc'
      end if

      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual total wood litter carbon',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'clitw','total wood litter carbon',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr ( (clitwm+clitws+clitwl), cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'clitw',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, clitw'
         stop 1
      end if
cc
cc total litterfall
cc
      if(iyearout_glo == 2) then
        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_totfall.nc'
      else
        filen = trim(myprocDir)//'output/yearly/totfall.nc'
      end if

      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual total litterfall carbon',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'totfall','total litterfall',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (totfall, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'totfall',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, totfall'
         stop 1
      end if
c
cc total soil carbon in slow pool
cc
      if(iyearout_glo == 2) then
        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_csoislo.nc'
      else
        filen = trim(myprocDir)//'output/yearly/csoislo.nc'
      end if

      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual total soil carbon in slow pool',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'csoislo','total soil carbon in slow pool',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (csoislo, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'csoislo',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, csoislo'
         stop 1
      end if
cc
cc total soil carbon in passive pool
cc
      if(iyearout_glo == 2) then
        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_csoipas.nc'
      else
        filen = trim(myprocDir)//'output/yearly/csoipas.nc'
      end if

      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual total soil carbon in pasive pool',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'csoipas','total soil carbon in passive pool',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (csoipas, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'csoipas',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, csoipas'
         stop 1
      end if
c
c co2 carbon exchange: net ecosystem, microbial resp, root resp, soil resp
c
      if(iyearout_glo == 2) then
        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_co2fluxes.nc'
      else
        filen = trim(myprocDir)//'output/yearly/co2fluxes.nc'
      end if

      if (mstep .eq. 1) then
         call inifile(idies,filen,
     >    'annual total carbon from exchange of co2',
     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,
     >    latscale,'','none',
     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
         dimnames(3) = 'time'
         call inivar(idies,'neetot',
     >    'total net ecosystem echange carbon',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'co2mic','total microbe respiration carbon',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'co2root','total root respiration carbon',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         dimnames(3) = 'time'
         call inivar(idies,'co2soi','total soil respiration carbon',
     >    'kg/m^2',3,dimnames,OCEAN,istat)
         call endini(idies,istat)
      end if
      call vec2arr (ayneetot, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'neetot',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, neetot'
         stop 1
      end if
      call vec2arr (ayco2mic, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'co2mic',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, co2mic'
         stop 1
      end if
      call vec2arr (ayco2root, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'co2root',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, co2root'
         stop 1
      end if
      call vec2arr (ayco2soi, cdummy)
      istart(3) = mstep
      icount(3) = 1
      call writevar(filen,'co2soi',istart,icount,cdummy,ftime,
     > tweight,tdate,istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in wyearly, co2soi'
         stop 1
      end if
c
c vegetation type
c
c      if(iyearout_glo == 2) then
c        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_vegtype0.nc'
c      else
c        filen = trim(myprocDir)//'output/yearly/vegtype0.nc'
c      end if
c
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'annual vegetation type - ibis classification',
c     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,
c     >    latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'vegtype0','vegetation type',
c     >    '_',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (vegtype0, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'vegtype0',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, vegtype0'
c         stop 1
c      end if
c
c fractional cover of upper and lower canopies
c
c      if(iyearout_glo == 2) then
c        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_fcover.nc'
c      else
c        filen = trim(myprocDir)//'output/yearly/fcover.nc'
c      end if
c
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'annual fractional cover of canopies','ibis wyearly',
c     >    cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'fu','fractional cover of upper canopy',
c     >    'fraction',3,dimnames,OCEAN,istat)
c         call inivar(idies,'fl','fractional cover of lower canopy',
c     >    'fraction',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (fu, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'fu',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, fu'
c         stop 1
c      end if
c      call vec2arr (fl, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'fl',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, fl'
c         stop 1
c      end if
cc
cc sapwood fraction
cc
c      if(iyearout_glo == 2) then
c        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_sapfrac.nc'
c      else
c        filen = trim(myprocDir)//'output/yearly/sapfrac.nc'
c      end if
c
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'annual sapwood fraction',
c     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,'','none',
c     >    'none',1,dummy_vals3rd,'',tunits,'gregorian',istat)
c         dimnames(3) = 'time'
c         call inivar(idies,'sapfrac','sapwood fraction',
c     >    'fraction',3,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      call vec2arr (sapfrac, cdummy)
c      istart(3) = mstep
c      icount(3) = 1
c      call writevar(filen,'sapfrac',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, sapfrac'
c         stop 1
c      end if
c
c bottom and top of vegetation canopies
c
c      if(iyearout_glo == 2) then
c        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_zcanopy.nc'
c      else
c        filen = trim(myprocDir)//'output/yearly/zcanopy.nc'
c      end if
c
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'annual height of vegetation canopies',
c     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,
c     >    'canopy','canopy','_',2,pindex,'up',
c     >    tunits,'gregorian',istat)
cc add global attribute to define canopies with text, use netcdf low-level com
c         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'canopy_def',
c     >    26+5,'1='//canopies(1)//' 2='//canopies(2))
c         dimnames(3) = 'canopy'
c         call inivar(idies,'zbot',
c     >    'bottom heights of lower and upper canopies',
c     >    'meters',4,dimnames,OCEAN,istat)
c         call inivar(idies,'ztop',
c     >    'top heights of lower and upper canopies',
c     >    'meters',4,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      do 20 k = 1, 2
c         call vec2arr (zbot(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
c 20   continue
c      istart(3) = 1
c      istart(4) = mstep
c      icount(3) = 2
c      call writevar(filen,'zbot',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, zbot'
c         stop 1
c      end if
c      do 25 k = 1, 2
c         call vec2arr (ztop(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
c 25   continue
c      istart(3) = 1
c      istart(4) = mstep
c      icount(3) = 2
c      call writevar(filen,'ztop',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, ztop'
c         stop 1
c      end if
c
c existence of pfts
c
c      if(iyearout_glo == 2) then
c        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_exist.nc'
c      else
c        filen = trim(myprocDir)//'output/yearly/exist.nc'
c      end if
c
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'annual existence of each plant functional type',
c     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,
c     >    'pft','plant fuctional type','_',npft,pindex,'',
c     >    tunits,'gregorian',istat)
cc add global attribute to define pfts with text, use netcdf low-level command
c         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',
c     >    npft*80,pftdef)
c         dimnames(3) = 'pft'
c         call inivar(idies,'exist','existence for each pft',
c     >    '_',4,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      do 30 k = 1, npft
c         call vec2arr (exist(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
c 30   continue
c      istart(3) = 1
c      istart(4) = mstep
c      icount(3) = npft
c      call writevar(filen,'exist',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, exist'
c         stop 1
c      end if
c
c
c fraction of canopy as pft
c
c      if(iyearout_glo == 2) then
c        filen = trim(myprocDir)//'output/yearly/'//num//'_yearly_frac.nc'
c      else
c        filen = trim(myprocDir)//'output/yearly/frac.nc'
c      end if
c
c      if (mstep .eq. 1) then
c         call inifile(idies,filen,
c     >    'fraction of canopy occupied by each pft',
c     >    'ibis wyearly',cdate,nlonsub,lonscale,nlatsub,latscale,
c     >    'pft','plant fuctional type','_',npft,pindex,'',
c     >    tunits,'gregorian',istat)
cc add global attribute to define pfts with text, use netcdf low-level command
c         istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',
c     >    npft*80,pftdef)
c         dimnames(3) = 'pft'
c         call inivar(idies,'frac','fraction of canopy',
c     >    '_',4,dimnames,OCEAN,istat)
c         call endini(idies,istat)
c      end if
c      do 35 k = 1, npft
c         call vec2arr (frac(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
c 35   continue
c      istart(3) = 1
c      istart(4) = mstep
c      icount(3) = npft
c      call writevar(filen,'frac',istart,icount,cdummy,ftime,
c     > tweight,tdate,istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in wyearly, frac'
c         stop 1
c      end if
c
      return
      end
c
c ---------------------------------------------------------------------
      subroutine readit (isimveg,snorth,ssouth,swest,seast,iwest,jnorth)
c ---------------------------------------------------------------------
c
c reads in initialization files and initializes some fields
c
      use clim_file_utils
      use comgrid
      use compar
      use cominput
      use combcs
      use comsoi
      use comveg
      use comwork
      use comcrop
      use comnitr
      use, INTRINSIC :: IEEE_ARITHMETIC,only: IEEE_IS_NAN
c
      implicit none
c
c Arguments
c
      integer isimveg,  ! dynamic vegetation (1) or static (0)
     >        iwest,    !
     >        jnorth    !       
c
      real snorth, ssouth, swest, seast
c
c Local variables
c

      integer  istat,    ! netcdf error flag
     > i, j, ntime,      ! loop indices
     > jj, ii, 
     > ndim,             ! number of dimensions
     > nlpoints,         ! number of land points
     > ldim

      integer:: idx,n_cmask
c
      real xlat
c
      integer jsouth, ieast
c
      character*80 filen
      character*39 pathn
c
      integer istart(4), icount(4) ! for reading vars
      integer(kind=4):: t_s,t_e                !clock counter, to check the total elapsed time for the whole program
      integer(kind=4):: clock_rate, clock_max  !see system_clock() for reference

c
c ---------------------------------------------------------------------
c
      data istart / 1,1,1,1 /, icount / nlon,nlat,1,1 /

! skip reading land mask, latitudes and longitudes when flg_udg == 1
! (Y.Li)
c
c 2-d surface and vegetation arrays
c
      call system_clock(t_s,clock_rate,clock_max)

      if(flg_udg == 0 ) then
        icount(3) = 1
c       pathn = '/Volumes/wet/kucharik/IBIS_input/input/'
c
c land mask, latitudes, and longitudes
c
c       filen = pathn//'surta.nc'
        filen = 'input/surta.nc'
        aname = 'surta'
        call readvar(filen,aname,'level',istart,icount,
     >   work2d,lonscale,latscale,cdummy(1),cdummy(ndim4),istat)
        if (istat.lt.0) then
           write (*,9000)
           print *, 'while reading surta'
           stop 1
        end if
c
c
        if (abs(abs(lonscale(1)-lonscale(2))-xres).gt.0.001 .or.
     >      abs(abs(latscale(1)-latscale(2))-yres).gt.0.001) then
           write (*,9000)
           write(*,*) 'resolution mismatch!'
           write(*,*) 'xres, yres in comgrid = ', xres, yres
           write(*,*) 'xres, yres from input/surta.nc = ',
     >     abs(lonscale(1)-lonscale(2)), abs(latscale(1)-latscale(2))
           stop 1
        end if
c
c subset the grid if not default whole grid
c
        if (snorth.lt.latscale(1) .or. ssouth.gt.latscale(nlat) .or.
     >       swest.gt.lonscale(1) .or.  seast.lt.lonscale(nlon)) then
          jnorth = 0
          if (snorth .lt. (latscale(nlat)+latscale(nlat-1))/2.) then
            jnorth = nlat
          else if (snorth .ge. (latscale(1)+latscale(2))/2.) then
            jnorth = 1
          else
            do 1 j = nlat-1,1,-1
              if (snorth .ge. (latscale(j)+latscale(j+1))/2.) jnorth = j
 1          continue
          end if
          jsouth = 0
          if (ssouth .lt. (latscale(nlat)+latscale(nlat-1))/2.) then
            jsouth = nlat
          else if (ssouth .ge. (latscale(1)+latscale(2))/2.) then
            jsouth = 1
          else
            do 2 j = nlat-1,1,-1
              if (ssouth .ge. (latscale(j)+latscale(j+1))/2.) jsouth = j
 2          continue
          end if
c
          iwest = 0
          if (swest .lt. (lonscale(1)+lonscale(2))/2.) then
            iwest = 1
          else if (swest .ge. (lonscale(nlon)+lonscale(nlon-1))/2.) then
            iwest = nlon
          else
            do 3 i = 2, nlon
              if(swest .ge. (lonscale(i)+lonscale(i-1))/2.) iwest=i
 3          continue
          end if
          ieast = 0
          if (seast .lt. (lonscale(1)+lonscale(2))/2.) then
            ieast = 1
          else if (seast .ge. (lonscale(nlon)+lonscale(nlon-1))/2.) then
            ieast = nlon
          else
            do 4 i = 2, nlon
              if(seast .ge. (lonscale(i)+lonscale(i-1))/2.) ieast=i
 4          continue
          end if
          nlonsub = ieast - iwest + 1
          nlatsub = jsouth - jnorth + 1
          istart(1) = iwest
          icount(1) = nlonsub
          istart(2) = jnorth
          icount(2) = nlatsub
        else
          iwest = 1
          ieast = nlon
          jnorth = 1
          jsouth = nlat
          nlonsub = nlon
          nlatsub = nlat
        end if
cc        print *, iwest, ieast, jnorth, jsouth
cc        print *, swest, seast, snorth, ssouth
c
c
c initialize lonindex, latindex for use in arr2vec, vec2arr, etc.
c and calculate the approximate the area of the land gridcells
c
        nlpoints = 0
c
c here, i/j refers to entire grid (1 to nlon/nlat), 
c ii/jj refers to subgrid (1 to nlonsub/nlatsub)
c
        do 5 j = jnorth, jsouth
c
          jj = j - jnorth + 1
          do 6 i = iwest, ieast
c
            ii = i - iwest + 1
            lmask(ii,jj) = nint(work2d(i,j))
c
            if (lmask(ii,jj).eq.1) then
c
              nlpoints = nlpoints + 1
              lonindex(nlpoints) = ii
              latindex(nlpoints) = jj
              xlat = latscale(j) * pi / 180.0
              garea(nlpoints) = yres * 111400.0 * xres * 111400.0 *
     >                          cos(xlat)
c
            end if
c
 6        continue
 5      continue
c
cc      print *, jnorth, jsouth, iwest, ieast
cc      print *, latscale(jnorth),latscale(jsouth),lonscale(iwest),lonscale(ieast)
cc      print *, ((lmask(i,j),i=1,nlonsub),j=1,nlatsub)
cc      print *, istart
cc      print *, icount
c
        do 7 j = jnorth, jsouth
          jj = j - jnorth + 1
          latscale(jj) = latscale(j)
 7      continue
c
        do 8 i = iwest, ieast
          ii = i - iwest + 1
          lonscale(ii) = lonscale(i)
 8      continue

      end if !on if(flg_udg == 0)
c
cc      print *, (lonscale(i),i=1,nlonsub)
cc      print *, (latscale(i),i=1,nlatsub)
cc      print *, lonindex
cc      print *, latindex
c
! 
! -- Overwrite the above sub-set grid information if user-defined grid
!    is used (flg_udg == 1) 
!
      if(flg_udg == 1 .or. flg_udg == 2) then

        nlonsub = udg_ie - udg_is + 1
        nlatsub = udg_je - udg_js + 1

        istart(1) = udg_is
        istart(2) = udg_js

        icount(1) = nlonsub
        icount(2) = nlatsub

        !--approximate lonscale and latscale for subgrid
        !  to lonscale(1~nlonsub) and latscale(1~nlatsub)
        idx = 1
        do i = udg_is, udg_ie
          lonscale(idx) = lon_llcorner + (i-1)*xres
          idx = idx + 1
        end do

        idx = 1
        do j = udg_js, udg_je
          if(flg_udg == 1) then
            latscale(idx) = lat_llcorner + (j-1)*yres
          elseif(flg_udg == 2) then  !from high latitude to low latitude
            latscale(idx) = lat_llcorner - (j-1)*yres
          end if

          idx = idx + 1
        end do

        idx = 1
        do j = 1, nlatsub
          do i = 1, nlonsub

            lonindex(idx) = i
            latindex(idx) = j
            xlat = latscale(j) * pi / 180.0
            garea(idx) = yres * 111400.0 * xres * 111400.0 *
     >                        cos(xlat)
            idx = idx + 1
          end do !on i
        end do !on j

        !some variables needed for compatability
        iwest  = udg_is
        jnorth = udg_js  !udg_js is jnorth for flg_udg=1 and flg_udg=2

        !cheat lmask and nlpoints
        lmask    = 1
        nlpoints = npoi
      end if  !on flg_udg == 1 or flg_udg == 2 
c
      if (nlpoints .ne. npoi) then
         write(*,9000)
         write(*,*) 'number of land points in input/surta.nc'
         write(*,*) 'does not match number of land points in comgrid'
         write(*,*) 'in surta =',nlpoints,' in comgrid =',npoi
         stop 1
      else
         write (*,9010) 
         write (*,9020) nlpoints
         write (*,9010) 
      end if

! summarize grid information
      write(6,'(A,i5,2(2x,f12.5))') ' total number, min, max (longtitude subgrid): ',nlonsub, lonscale(1),lonscale(nlonsub)
      write(6,'(A,i5,2(2x,f12.5))') ' total number, min, max (latitude   subgrid): ',nlatsub, latscale(1),latscale(nlatsub)
c
cc
cc dummy variable example, 4-d, but 3rd dim ('level') = 1
cc copy and chanve for a new variable
cc
c      filen = 'input/dummyv.nc'
c      aname = 'dummyv'
c      call readvar(filen,aname,'level',istart,icount,
c     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
c      if (istat.lt.0) then
c         write(*,9000)
c         print *, 'while reading dummyv'
c         stop 1
c      end if
c      call arr2vec (cdummy, xindummy)

!appropriate start and count should be constructed from above (udg or non-udg)
c
c cmask
c
      !initialize and read computational mask, if request
      cmask = 1  !default, every cell is to compute

      if(flg_udg == 2) then   !get <cmask> from <surta.nc>
        write(6,*) 'flg_udg == 2, assigning <surta.nc> data to <cmask>...' 
        filen = 'input/surta.nc'
        aname = 'surta'
        call readvar(filen,aname,'level',istart,icount,cmask, 
     >               cdummy(1),cdummy(1),cdummy(1),cdummy(1),istat)

        if(istat.lt.0) then
          write(6,*) 'Error in reading <surta.nc> to <cmask> from file: ',trim(filen)
          stop 'Error in readit'
        end if

      end if !on if flg_udg == 2

      if(flg_cmask == 1) then
        filen = 'input/'//trim(fileName_cmask)
        aname = trim(varName_cmask)

        call readvar(filen,aname,'None',istart,icount,cmask, 
     >               cdummy(1),cdummy(1),cdummy(1),cdummy(1),istat)
        if(istat.lt.0) then
          write(6,*) 'Error in reading computational mask file: ',trim(filen)
          stop 'Error in readit'
        end if
        write(6,*) 'read-in cmask file completed...'

      end if !on if flg_cmask

      !quick summary
      n_cmask = 0;
      do i = 1,npoi
        if(cmask(i) == 1) n_cmask = n_cmask+1
      end do
      write(6,*) 'number of computational cell(cmask=1): ',n_cmask, '; npoi = ',npoi
c
c tminavgann 
c
      if (fn_tminavgann .eq. ' ') then
        write(*,9000)
        print *, 'fn_tminavgann not specified'
        stop 1
      end if
      call climate_full_fname(filen, fn_tminavgann, 0, .false.)
      aname = var_tminavgann
      call readvar(filen,aname,level_var,istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat.lt.0) then
         write(*,9000)
         print *, 'while reading tminavgann'
         stop 1
      end if
      call arr2vec (cdummy, tminavgann)

c
c average crop planting date - held constant
c
c     filen = pathn//'pdate.ave.nc'
c      filen = 'input/pdate.ave.nc'
c      aname = 'pdate'
c      call readvar(filen,aname,'level',istart,icount,
c     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
c      if (istat.lt.0) then
c         write(*,9000)
c         print *, 'while reading pdate'
c         stop 1
c      end if
c      call arr2vec (cdummy, xinavepdate)
c
c average corn hybrid- held constant
c
c     filen = pathn//'hybrid.ave.nc'
c      filen = 'input/hybrid.ave.nc'
c      aname = 'hybrid'
c      call readvar(filen,aname,'level',istart,icount,
c     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
c      if (istat.lt.0) then
c         write(*,9000)
c         print *, 'while reading hybrid'
c         stop 1
c      end if
c      call arr2vec (cdummy, xinavehybrid)
c
c topography
c
c     filen = pathn//'topo.nc'
      filen = 'input/topo.nc'
      aname = 'topo'
      call readvar(filen,aname,'level',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat.lt.0) then
         write(*,9000)
         print *, 'while reading topo'
         stop 1
      end if
      call arr2vec (cdummy, xintopo)

c
c fixed vegetation map
c
      if (isimveg .le. 1) then
c        filen = pathn//'vegtype.nc'
         filen = 'input/vegtype.nc'
         aname = 'vegtype'
         call readvar(filen,aname,'level',istart,icount,cdummy,
     >    work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
         if (istat.lt.0) then
            write(*,9000)
            print *, 'while reading vegtype'
            stop 1
         end if
         call arr2vec (cdummy, xinveg)
      end if
cc 2-d soil array
cc
c     filen = 'nput/soil.nc'
c     aname = 'soil'
c     call readvar(filen,aname,'level',istart,icount,
c    > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
c     if (istat.lt.0) then
c        write(*,9000)
c        print *, 'while reading soil'
c        stop 1
c     end if
c     call arr2vec (cdummy, soita)
c
c delta t
c
c     filen = pathn//'deltat.nc'
      filen = 'input/deltat.nc'
      aname = 'deltat'
      call readvar(filen,aname,'level',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat.lt.0) then
         write(*,9000)
         print *, 'while reading deltat'
         stop 1
      end if
      call arr2vec (cdummy, deltat)
c
c 3-d soil texture array
c
c icount(3) is the 11 layers used in soil_text.nc
c
      icount(3) = 11 
      icount(4) = 1
c     filen = pathn//'soil_text.nc'
      filen = 'input/soil_text.nc'
      aname = 'domtext'
      call readvar(filen,aname,'layer',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat.lt.0) then
         write(*,9000)
         print *, 'while reading soil_text'
         stop 1
      end if
      do 14 j = 1, nsoilay
        call arr2vec (cdummy((j-1)*nlonsub*nlatsub + 1), domtext(1,j))
 14   continue

      if(flg_sandlayers == 1) then
        domtext(:,10) = 1
        domtext(:,11) = 1
      end if
 
      !temporary treatment
      do j = 1,nsoilay
        do i = 1,npoi
        
          if(domtext(i,j) .eq. 0) then
            write(6,*) '(io)Warning! soil type (domtext) in soil_text.nc = 0 for (cell#,layer#) ', i,j
            write(6,*) 'permanent treatment for this is needed. For now,reset it as 1'
            
            domtext(i,j) = 1
          end if !on if(domtext .eq.0)
 
        end do !on i
      end do !on j
      !end tmporary treatment
c
c WJS 01.21.10: Commented out the following code, since we don't have
c the necessary datasets at 5' resolution, and 'clay' is only referred
c to by commented-out code in the model.
c
c WJS 01.21.10: Also, I believe there is a big in this code, in that
c the read-in sand data is put into the clay array!
c$$$      icount(3) = 6 
c$$$      icount(4) = 1
c$$$c     filen = pathn//'soita.sand.nc'
c$$$      filen = 'input/soita.sand.nc'
c$$$      aname = 'sandpct'
c$$$      call readvar(filen,aname,'layer',istart,icount,
c$$$     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
c$$$      if (istat.lt.0) then
c$$$         write(*,9000)
c$$$         print *, 'while reading soita.sand'
c$$$         stop 1
c$$$      end if
c$$$      do 12 j = 1, 6
c$$$        call arr2vec (cdummy((j-1)*nlonsub*nlatsub + 1), clay(1,j))
c$$$ 12   continue
c$$$c
c$$$      icount(3) = 6 
c$$$      icount(4) = 1
c$$$c     filen = pathn//'soita.clay.nc'
c$$$      filen = 'input/soita.clay.nc'
c$$$      aname = 'claypct'
c$$$      call readvar(filen,aname,'layer',istart,icount,
c$$$     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
c$$$      if (istat.lt.0) then
c$$$         write(*,9000)
c$$$         print *, 'while reading soita.clay'
c$$$         stop 1
c$$$      end if
c$$$      do 13 j = 1, nsoilay
c$$$        call arr2vec (cdummy((j-1)*nlonsub*nlatsub + 1), clay(1,j))
c$$$ 13   continue
c
c 3-d climate arrays
c
      icount(3) = 1
      icount(4) = 12
c
      if (fn_wetd_mon_clim .eq. ' ') then
c     It is okay if we don't have any wetd data, but then we can't use the weather generator.
c     (When we try to use the weather generator, we first check whether xinwet is uninitialized.)
        call const(clmwet, npoi*12, un_init) 
      else       
        call climate_full_fname(filen, fn_wetd_mon_clim, 0, .false.)
        aname = var_wetd_mon_clim
        call readvar(filen,aname,level_var,istart,icount,
     >       cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
        if (istat.lt.0) then
          write(*,9000)
          print *, 'while reading wetd'
          stop 1
        end if
        do 15 ntime = 1,12
          call arr2vec (cdummy((ntime-1)*nlonsub*nlatsub + 1),
     >         clmwet(1,ntime))
 15     continue
      end if
c
      if (read_tmin_tmax) then
c     Read tmin and tmax; calculate t and trng

        if (fn_tmin_mon_clim .eq. ' ') then
          write(*,9000)
          print *, 'fn_tmin_mon_clim must be specified for read_tmin_tmax = true'
          stop 1
        end if
        call climate_full_fname(filen, fn_tmin_mon_clim, 0, .false.)
        aname = var_tmin_mon_clim
        call readvar(filen,aname,level_var,istart,icount,
     >       cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
        if (istat.lt.0) then
          write(*,9000)
          print *, 'while reading tmin'
          stop 1
        end if
        do ntime = 1,12
          call arr2vec (cdummy((ntime-1)*nlonsub*nlatsub + 1),
     >         clmtmin(1,ntime))
        end do        

        if (fn_tmax_mon_clim .eq. ' ') then
          write(*,9000)
          print *, 'fn_tmax_mon_clim must be specified for read_tmin_tmax = true'
          stop 1
        end if
        call climate_full_fname(filen, fn_tmax_mon_clim, 0, .false.)
        aname = var_tmax_mon_clim
        call readvar(filen,aname,level_var,istart,icount,
     >       cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
        if (istat.lt.0) then
          write(*,9000)
          print *, 'while reading tmax'
          stop 1
        end if
        do ntime = 1,12
          call arr2vec (cdummy((ntime-1)*nlonsub*nlatsub + 1),
     >         clmtmax(1,ntime))
        end do 

        if(flg_udg == 2) then  !initialize
          clmt = 0.0
        end if

        do ntime = 1, 12
          do i = 1, npoi
            !if(cmask(i)==0 .or. IEEE_IS_NAN(clmtmax(i, ntime)) .or. IEEE_IS_NAN(clmtmin(i, ntime))) cycle
            if(cmask(i)==0) cycle
            clmtrng(i, ntime) = clmtmax(i, ntime) - clmtmin(i, ntime)
            clmt(i, ntime) = tmax_wt * clmtmax(i, ntime) + tmin_wt * clmtmin(i, ntime)
          end do
        end do

      else  ! .not. read_tmin_tmax
c     Read t and trng; calculate tmin and tmax

        if (fn_temp_mon_clim .eq. ' ') then
          write(*,9000)
          print *, 'fn_temp_mon_clim must be specified for read_tmin_tmax = false'
          stop 1
        end if
        call climate_full_fname(filen, fn_temp_mon_clim, 0, .false.)
        aname = var_temp_mon_clim
        call readvar(filen,aname,level_var,istart,icount,
     >       cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
        if (istat.lt.0) then
          write(*,9000)
          print *, 'while reading temp'
          stop 1
        end if
        do 20 ntime = 1,12
          call arr2vec (cdummy((ntime-1)*nlonsub*nlatsub + 1),
     >         clmt(1,ntime))
 20     continue
c     
        if (fn_dtr_mon_clim .eq. ' ') then
          write(*,9000)
          print *, 'fn_dtr_mon_clim must be specified for read_tmin_tmax = false'
          stop 1
        end if
        call climate_full_fname(filen, fn_dtr_mon_clim, 0, .false.)
        aname = var_dtr_mon_clim
        call readvar(filen,aname,level_var,istart,icount,
     >       cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
        if (istat.lt.0) then
          write(*,9000)
          print *, 'while reading trange'
          stop 1
        end if
        do 25 ntime = 1,12
          call arr2vec (cdummy((ntime-1)*nlonsub*nlatsub + 1),
     >         clmtrng(1,ntime))
 25     continue

        do ntime = 1, 12
          do i = 1, npoi
            if(cmask(i)==0) cycle
            clmtmin(i, ntime) = clmt(i, ntime) - tmax_wt * clmtrng(i, ntime)
            clmtmax(i, ntime) = clmt(i, ntime) + tmin_wt * clmtrng(i, ntime)
          end do
        end do

      end if  ! read_tmin_tmax

c
      if (fn_prec_mon_clim .eq. ' ') then
        write(*,9000)
        print *, 'fn_prec_mon_clim not specified'
        stop 1
      end if
      call climate_full_fname(filen, fn_prec_mon_clim, 0, .false.)
      aname = var_prec_mon_clim
      call readvar(filen,aname,level_var,istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat.lt.0) then
         write(*,9000)
         print *, 'while reading prec'
         stop 1
      end if
      do 30 ntime = 1,12
         call arr2vec (cdummy((ntime-1)*nlonsub*nlatsub + 1),
     >    clmprec(1,ntime))
 30   continue
c
      if (fn_wspd_mon_clim .eq. ' ') then
        write(*,9000)
        print *, 'fn_wspd_mon_clim not specified'
        stop 1
      end if      
      call climate_full_fname(filen, fn_wspd_mon_clim, 0, .false.)
      aname = var_wspd_mon_clim
      call readvar(filen,aname,level_var,istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat.lt.0) then
         write(*,9000)
         print *, 'while reading wspd'
         stop 1
      end if
      do 35 ntime = 1,12
         call arr2vec (cdummy((ntime-1)*nlonsub*nlatsub + 1),
     >    clmw(1,ntime))
 35   continue
c
      if (read_radiation) then
c     Read rads, set cloud to un_init

        if (fn_rads_mon_clim .eq. ' ') then
          write(*,9000)
          print *, 'fn_rads_mon_clim must be specified for read_radiation = true'
          stop 1
        end if
        call climate_full_fname(filen, fn_rads_mon_clim, 0, .false.)
        aname = var_rads_mon_clim
        call readvar(filen,aname,level_var,istart,icount,
     >       cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
        if (istat.lt.0) then
          write(*,9000)
          print *, 'while reading rads'
          stop 1
        end if
        do ntime = 1,12
          call arr2vec (cdummy((ntime-1)*nlonsub*nlatsub + 1),
     >         clmrads(1,ntime))
        end do    
            
        call const(clmcld, npoi*12, un_init)

      else ! .not. read_radiation
c     Read cloud, set rads to un_init

        if (fn_cloud_mon_clim .eq. ' ') then
          write(*,9000)
          print *, 'fn_cloud_mon_clim must be specified for read_radiation = false'
          stop 1
        end if
        call climate_full_fname(filen, fn_cloud_mon_clim, 0, .false.)
        aname = var_cloud_mon_clim
        call readvar(filen,aname,level_var,istart,icount,
     >       cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
        if (istat.lt.0) then
          write(*,9000)
          print *, 'while reading cld'
          stop 1
        end if
        do 40 ntime = 1,12
          call arr2vec (cdummy((ntime-1)*nlonsub*nlatsub + 1),
     >         clmcld(1,ntime))
 40     continue

        call const(clmrads, npoi*12, un_init)

      end if  ! read_radiation
c
      if (fn_rh_mon_clim .eq. ' ') then
        write(*,9000)
        print *, 'fn_rh_mon_clim not specified'
        stop 1
      end if
      call climate_full_fname(filen, fn_rh_mon_clim, 0, .false.)
      aname = var_rh_mon_clim
      call readvar(filen,aname,level_var,istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat.lt.0) then
         write(*,9000)
         print *, 'while reading rh'
         stop 1
      end if
      do 45 ntime = 1,12
         call arr2vec (cdummy((ntime-1)*nlonsub*nlatsub + 1),
     >    clmq(1,ntime))
 45   continue
c
 
      call read_fert_data('input', iwest, jnorth)

c WJS 1.11.10 Commented out the following code specific to the ethanol-nitrate runs
c$$$c CJK 4.23.07 - modified for ethanol-nitrate runs
c$$$c Simon provided a new N-fertilizer dataset that had average
c$$$c rates for 2001-2005, so only 1 year of inputs were present
c$$$c
c$$$c corn 2005 average
c$$$c
c$$$      icount(4) = 1
c$$$c
c$$$      filen = 'input/frate_2005.corn.nc'
c$$$      aname = 'frate'
c$$$      call readvar(filen,aname,'level',istart,icount,
c$$$     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
c$$$      if (istat.lt.0) then
c$$$         write(*,9000)
c$$$         print *, 'while reading frate 2005 corn'
c$$$         stop 1
c$$$      end if
c$$$      call arr2vec (cdummy, fert05maize)
c$$$c
c$$$c soybean  2005 average
c$$$c
c$$$      filen = 'input/frate_2005.soy.nc'
c$$$      aname = 'frate'
c$$$      call readvar(filen,aname,'level',istart,icount,
c$$$     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
c$$$      if (istat.lt.0) then
c$$$         write(*,9000)
c$$$         print *, 'while reading frate 2005 soybean'
c$$$         stop 1
c$$$      end if
c$$$      call arr2vec (cdummy, fert05soy)
c$$$c
c$$$c wheat  2005 average
c$$$c
c$$$      filen = 'input/frate_2005.wheat.nc'
c$$$      aname = 'frate'
c$$$      call readvar(filen,aname,'level',istart,icount,
c$$$     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
c$$$      if (istat.lt.0) then
c$$$         write(*,9000)
c$$$         print *, 'while reading frate 2005 wheat'
c$$$         stop 1
c$$$      end if
c$$$      call arr2vec (cdummy, fert05wheat)
c$$$c
c$$$c sorghum  2005 average
c$$$c
c$$$      filen = 'input/frate_2005.sorg.nc'
c$$$      aname = 'frate'
c$$$      call readvar(filen,aname,'level',istart,icount,
c$$$     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
c$$$      if (istat.lt.0) then
c$$$         write(*,9000)
c$$$         print *, 'while reading frate 2005 sorghum'
c$$$         stop 1
c$$$      end if
c$$$      call arr2vec (cdummy, fert05sorg)
c$$$c
c time-dependent nitrogen deposition fertilizer dataset 
c begins in 1940, through 2000 
c
c WJS 01.21.10: Commented out the following code, since we don't have
c the necessary dataset at 5' resolution, and ndepfact is only referred
c to by commented-out code in the model
c
c WJS 02.03.10: IMPORTANT NOTE: If you use this code, you will need to
c change how the file is read to read one year at a time (as is done in
c subroutine read_one_fert_file), now that cdummy has been made smaller.
c
c$$$      icount(3) = 1
c$$$      icount(4) = 60 
c$$$c     filen = pathn//'ndep.1940.2000.nc'
c$$$      filen = 'input/ndep.1940.2000.nc'
c$$$      aname = 'ndep'
c$$$      call readvar(filen,aname,'level',istart,icount,
c$$$     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
c$$$      if (istat.lt.0) then
c$$$         write(*,9000)
c$$$         print *, 'while reading ndep '
c$$$         stop 1
c$$$      end if
c$$$      do 53  ntime = 1,60
c$$$         call arr2vec (cdummy((ntime-1)*nlonsub*nlatsub + 1),
c$$$     >    ndepfact(1,ntime))
c$$$ 53   continue
c
c WJS 05.06.10: Read in header file for multi-year planting date & cultivar-related maps:
c
      if (management_prescribed) then
        call read_management_header('input')
c
c WJS 05.06.10: Read in spinup planting date & cultivar-related variables. These variables will remained unchanged until the
c simulation year becomes >= management_year1, at which point the variables will be overwritten with the current year's data from
c the multi-year files (through the read_management_data subroutine).
c
        icount(3) = 1
        icount(4) = 1

        filen = 'input/management/planting_date_average.nc'
        aname = 'planting.date'
        call readvar(filen,aname,'level',istart,icount,
     >       cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
        if (istat.lt.0) then
          write(*,9000)
          print *, 'while reading planting.date'
          stop 1
        end if
        call arr2vec(cdummy, xinpdate)

        filen = 'input/management/cultivar_average.nc'
        aname = 'total.gdd'
        call readvar(filen,aname,'level',istart,icount,
     >       cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
        if (istat.lt.0) then
          write(*,9000)
          print *, 'while reading total.gdd'
          stop 1
        end if
        call arr2vec(cdummy, xinhybrid)

        filen = 'input/management/cultivar_average.nc'
        aname = 'grainfill'
        call readvar(filen,aname,'level',istart,icount,
     >       cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
        if (istat.lt.0) then
          write(*,9000)
          print *, 'while reading grainfill'
          stop 1
        end if
        call arr2vec(cdummy, xingrnfill)

        filen = 'input/management/cultivar_average.nc'
        aname = 'days.to.harvest'
        call readvar(filen,aname,'level',istart,icount,
     >       cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
        if (istat.lt.0) then
          write(*,9000)
          print *, 'while reading days.to.harvest'
          stop 1
        end if
        call arr2vec(cdummy, xindaystoharv)

      end if   ! if (management_prescribed)

c WJS 05.05.10: The following is OLD code inherited from Chris for reading in planting date & hybrid information; it has been
c superseded  by the above code.
c
c time-dependent planting date information for corn from
c control runs 
c i, time
c begins in 1940, through 1960
c
c WJS 02.03.10: IMPORTANT NOTE: If you use this code, you will need to
c change how the file is read to read one year at a time (as is done in
c subroutine read_one_fert_file), now that cdummy has been made smaller.
c
c       istart(1) = 1  
c       istart(2) = 1   
c       istart(3) = 14       ! maize pft
c       istart(4) = 1  
c
c      icount(1) = 53       ! nlons in output
c      icount(2) = 34       ! nlats in output
c       icount(1) = 2       ! nlons in output
c       icount(2) = 2       ! nlats in output
c       icount(3) = 1        ! number of pfts
c       icount(4) = 21       ! number of years to read 
c
c      filen = pathn//'crops.input.nc'
c      aname = 'plantdate'
c      call readvar(filen,aname,'level',istart,icount,
c     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
c      if (istat.lt.0) then
c          write(*,9000)
c         print *, 'while reading plantdate'
c         stop 1
c      end if
c      do 533  ntime = 1, icount(4)
c         call arr2vec (cdummy((ntime-1)*nlonsub*nlatsub + 1),
c     >    xinpdate(1,ntime))
c 533   continue
c
c read in gdd information - use gddfzcorn to determine typical hybrid  
c to plant for the whole period in holding constant 
c
c WJS 02.03.10: IMPORTANT NOTE: If you use this code, you may need to
c change how the file is read to read one year at a time (as is done in
c subroutine read_one_fert_file), now that cdummy has been made smaller.
c
c      istart(3) = 1        ! maize pft
c      icount(3) = 1
c      istart(4) = 6577 
c      icount(4) = 1 
c
c      filen = pathn//'gdd.input.nc'
c      aname = 'gddfzcorn'
c      call readvar(filen,aname,'level',istart,icount,
c     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
c      if (istat.lt.0) then
c         write(*,9000)
c         print *, 'while reading gddfzcorn '
c         stop 1
c      end if
c      do 540  ntime = 1,  icount(4) 
c         call arr2vec (cdummy((ntime-1)*nlonsub*nlatsub + 1),
c     >    xinhybrid(1,ntime))
c          write(*,*) xinhybrid(1,ntime)
c 540   continue

      call set_irrigation_map('input', iwest, jnorth)

c
c copy all climatology fields to clm+anom fields for spin up
c
      ndim = npoi*12
      call scopy (ndim, clmt, xint)
      call scopy (ndim, clmtrng, xintrng)
      call scopy (ndim, clmtmin, xintmin)
      call scopy (ndim, clmtmax, xintmax)
      call scopy (ndim, clmprec, xinprec)
      call scopy (ndim, clmcld, xincld)
      call scopy (ndim, clmrads, xinrads)
      call scopy (ndim, clmq, xinq)
      call scopy (ndim, clmw, xinwind)
      call scopy (ndim, clmwet, xinwet)
c
 9000 format (1x,'ERROR in subroutine readit')
 9010 format (1x,' ')
 9020 format (1x,'number of land points: ', i10)
c
c return to main program
c
      call system_clock(t_e,clock_rate,clock_max)
      write(6,'(A,g16.8,A)') '(readit)time to read init data',real(t_e-t_s)/(real(clock_rate)),' seconds'

      return
      end
c 
c ---------------------------------------------------------------------
      subroutine read_fert_data(pathn, iwest, jnorth)
c ---------------------------------------------------------------------
c
c Read in fertilizer maps and their associated header file
c
      use comgrid
      use compar
      use comcrop
c
      implicit none
c
c Arguments
      character*(*) pathn  ! path in which frate*.nc files and frate.hdr are located

      integer 
     >  iwest,             ! first longitude index to read
     >  jnorth             ! first latitude index to read
c
c Local variables
c
      logical verbose      ! if true, print out extra info
      integer funit        ! file unit assignment for input
      integer ios          ! I/O error status

      parameter (verbose = .FALSE.)
      parameter (funit = 101)

c ---------------------------------------------------------------------

c Read in frate.hdr file
c This file contains the following variables, one per line:
c  - fert_year1 (first year in frate* files)
c  - fert_nyears (number of years in frate* files)
      open(unit=funit,
     >  file=trim(pathn) // '/frate.hdr',
     >  status='old',
     >  iostat=ios)

      if (ios .ne. 0) then
        write(*,*) 'READ_FERT_DATA: Error opening file: ',
     >    trim(pathn) // '/frate.hdr'
        write(*,*) 'IOS = ', ios
        stop 1
      end if

      read (funit,*) fert_year1
      read (funit,*) fert_nyears

      close(funit)

      if (fert_nyears .le. 0) then
        write(*,*) 'ERROR: fert_nyears must be > 0'
        write(*,*) 'fert_nyears = ', fert_nyears
        write(*,*) 'in: ', trim(pathn) // '/frate.hdr'
        stop 1
      end if

      if (fert_nyears .gt. max_fert_years) then
        write(*,*) 'ERROR: fert_nyears > max_fert_years'
        write(*,*) 'fert_nyears = ', fert_nyears
        write(*,*) 'max_fert_years = ', max_fert_years
        write(*,*) 'Please increase max_fert_years and recompile'
        stop 1
      end if

      fert_lastyear = fert_year1 + fert_nyears - 1
      
      if (verbose) then
        write(*,*) 'fert_year1    = ', fert_year1
        write(*,*) 'fert_lastyear = ', fert_lastyear
      end if

c Read in the data files
      call read_one_fert_file(fertmaize, pathn, 'corn', iwest, jnorth)
      call read_one_fert_file(fertsoy, pathn, 'soy', iwest, jnorth)
      call read_one_fert_file(fertwheat, pathn, 'wheat', iwest, jnorth)

      return
      end


c ---------------------------------------------------------------------
      subroutine read_one_fert_file(dat, pathn, crop, iwest, jnorth)
c ---------------------------------------------------------------------
c
c Read in the fertilizer map for one crop
c
c WJS (5-4-10): The code here is very similar to the code in read_one_management_variable; these two subroutines could be combined
c into one, more general subroutine
c
      use comgrid
      use compar
      use comwork
      use comcrop
c
      implicit none
c
c Arguments
      real dat(npoi,*)     ! output array (e.g. fertmaize)
      character*(*) pathn  ! path in which frate*.nc files are located
      character*(*) crop   ! file name will be frate.<crop>.nc

      integer 
     >  iwest,             ! first longitude index to read
     >  jnorth             ! first latitude index to read

c
c Local variables
c
      integer istart(4),   ! for reading vars
     >        icount(4)    

      integer time_len     ! number of times in the file
      integer ntime        ! loop counter

      integer istat        ! netcdf error flag

      character*80 filen

      data istart / 1,1,1,1 /, icount / nlon,nlat,1,1 /
      istart(1) = iwest
      istart(2) = jnorth
      icount(1) = nlonsub
      icount(2) = nlatsub

      filen = trim(pathn) // '/frate.' // trim(crop) // '.nc'

c First make sure that the file contains the proper number of years
      call dimlen(filen, 'time', time_len, istat)
      if (istat .lt. 0) then
        write(*,*) 'ERROR IN READ_ONE_FERT_FILE'
        write(*,*) 'while getting length of time dimension from'
        write(*,*) trim(filen)
        stop 1
      end if

      if (time_len .ne. fert_nyears) then
        write(*,*) 'ERROR IN READ_ONE_FERT_FILE'
        write(*,*) 'time_len .ne. fert_nyears'
        write(*,*) 'time_len = ', time_len
        write(*,*) 'fert_nyears = ', fert_nyears
        write(*,*) 'file = ', trim(filen)
        stop 1
      end if

c Now read the data, one year at a time
c (reading all the years at once would overflow cdummy)
      aname = 'frate'
      do ntime = 1, fert_nyears
        istart(4) = ntime
        call readvar(filen,aname,'level',istart,icount,
     >    cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
        if (istat.lt.0) then
          write(*,*) 'ERROR IN READ_ONE_FERT_FILE'
          write(*,*) 'while reading ', trim(filen)
          stop 1
        end if
        call arr2vec (cdummy, dat(1,ntime))
      end do

      return
      end

c ---------------------------------------------------------------------
      subroutine read_management_header(pathn)
c ---------------------------------------------------------------------
c
c Read in header file for multi-year management data files (planting_date_and_cultivar.hdr)
c This file contains the following variables, one per line:
c  - management_year1 (first year in planting date & cultivar-related netcdf files)
c  - management_nyears (number of years in planting date & cultivar-related netcdf files)
c
c WJS (5-4-10): This code could be changed to instead read this information from the 'time' variable of one of the netcdf files, but
c that would be a little trickier; for now I am using this .hdr file partly for consistency with the frate.hdr file.
c
      use comcrop, only: management_year1, management_lastyear
c
      implicit none
c
c Arguments
      character*(*) pathn  ! path in which planting_date_and_cultivar.hdr is located
c
c Local variables
c
      logical verbose      ! if true, print out extra info
      character*80 filen   ! name of input file, relative to pathn
      integer funit        ! file unit assignment for input
      integer ios          ! I/O error status
      integer nyears

      parameter (verbose = .FALSE.)
      parameter (filen = 'management/planting_date_and_cultivar.hdr')
      parameter (funit = 101)  

      open(unit=funit,
     >  file=trim(pathn) // '/' // trim(filen),
     >  status='old',
     >  iostat=ios)

      if (ios .ne. 0) then
        write(*,*) 'READ_MANAGEMENT_HEADER: Error opening file: ',
     >    trim(pathn) // '/' // trim(filen)
        write(*,*) 'IOS = ', ios
        stop 1
      end if

      read (funit,*) management_year1
      read (funit,*) nyears

      close(funit)

      if (nyears .le. 0) then
        write(*,*) 'ERROR: nyears must be > 0'
        write(*,*) 'nyears = ', nyears
        write(*,*) 'in: ', trim(pathn) // '/' // trim(filen)
        stop 1
      end if

      management_lastyear = management_year1 + nyears - 1
      
      if (verbose) then
        write(*,*) 'management_year1    = ', management_year1
        write(*,*) 'management_lastyear = ', management_lastyear
      end if

      return
      end



c ---------------------------------------------------------------------
      subroutine read_management_data(pathn, iyear, iwest, jnorth, planting_input, cultivar_input)
c ---------------------------------------------------------------------
c
c Read the current year from the multi-year planting date & cultivar-related maps. However, if the current year is earlier than the
c first year in these multi-year maps, then don't read anything: Instead, we use the spin-up files that have already been read in
c (in readit). If the current year is later than the last year in these multi-year maps, then we halt the program, because we don't
c know how to handle this situation.
c
c PRE: The associated header file must already have been read in (through read_management_header), setting management_year1 and
c management_lastyear appropriately.
c
      use comcrop, only : management_year1, management_lastyear,
     >  xinpdate, xinhybrid, xingrnfill, xindaystoharv
      
c
      implicit none
c
c Arguments
      character*(*) pathn  ! path in which planting date & cultivar-related netcdf files are located

      integer 
     >  iyear,             ! current year of the simulation
     >  iwest,             ! first longitude index to read
     >  jnorth             ! first latitude index to read

c For planting_input and cultivar_input: File name is determined as follows: If planting_input='foo' and cultivar_input='bar', then
c the planting date data are read from a file named planting_date_foo.nc, and the cultivar data are read from a file named
c cultivar_bar.nc.
      character*(*)
     >  planting_input,  ! Suffix on the planting date input file name
     >  cultivar_input   ! Suffix on the cultivar input file name
c
c ---------------------------------------------------------------------

      if (iyear .lt. management_year1) then
        write(*,*) iyear, ': Using spin-up files for planting date & cultivar'
        return
      else if (iyear .gt. management_lastyear) then
        write(*,*) "ERROR IN READ_MANAGEMENT_DATA: iyear > management_lastyear currently unimplemented"
        write(*,*) 'iyear = ', iyear
        write(*,*) 'management_lastyear = ', management_lastyear
        stop 1
      else  ! management_year1 <= iyear <= management_lastyear
c Read in the data files for the current year
        call read_one_management_variable(xinpdate,
     >    trim(pathn)//'/management/planting_date_'//trim(planting_input)//'.nc', 'planting.date',
     >    iyear, iwest, jnorth)
        call read_one_management_variable(xinhybrid,
     >    trim(pathn)//'/management/cultivar_'//trim(cultivar_input)//'.nc', 'total.gdd',
     >    iyear, iwest, jnorth)
        call read_one_management_variable(xingrnfill,
     >    trim(pathn)//'/management/cultivar_'//trim(cultivar_input)//'.nc', 'grainfill',
     >    iyear, iwest, jnorth)
        call read_one_management_variable(xindaystoharv,
     >    trim(pathn)//'/management/cultivar_'//trim(cultivar_input)//'.nc', 'days.to.harvest',
     >    iyear, iwest, jnorth)
      end if

      return
      end


c ---------------------------------------------------------------------
      subroutine read_one_management_variable(dat, filen, varname, iyear, iwest, jnorth)
c ---------------------------------------------------------------------
c
c Read in one planting date or cultivar-related map, for the given year.
c
c PRE: management_year1 <= iyear <= management_lastyear
c
      use comgrid
      use compar
      use comwork
      use comcrop, only : management_year1
c
      implicit none
c
c Arguments
      real dat(npoi)         ! output array (e.g. xinpdate)
      character*(*) filen    ! file name containing the variable of interest (including full path)
      character*(*) varname  ! variable name in the given file

      integer 
     >  iyear,             ! calendar year to read (this is NOT the same as the time index in the file - e.g., iyear=1981, not 1)
     >  iwest,             ! first longitude index to read
     >  jnorth             ! first latitude index to read

c
c Local variables
c
      integer istart(4),   ! for reading vars
     >        icount(4)    

      integer time_len     ! number of times in the file

      integer istat        ! netcdf error flag

      data istart / 1,1,1,1 /, icount / nlon,nlat,1,1 /
      istart(1) = iwest
      istart(2) = jnorth
      icount(1) = nlonsub
      icount(2) = nlatsub

c Determine time index to read, and make sure that the file contains enough years
      istart(4) = iyear - management_year1 + 1
      call dimlen(filen, 'time', time_len, istat)
      if (istat .lt. 0) then
        write(*,*) 'ERROR IN READ_ONE_MANAGEMENT_VARIABLE'
        write(*,*) 'while getting length of time dimension from'
        write(*,*) trim(filen)
        stop 1
      end if
      if (time_len .lt. istart(4)) then
        write(*,*) 'ERROR IN READ_ONE_MANAGEMENT_VARIABLE'
        write(*,*) 'time_len .lt. istart(4)'
        write(*,*) 'time_len = ', time_len
        write(*,*) 'istart(4) = ', istart(4)
        write(*,*) 'file = ', trim(filen)
        stop 1
      end if

c Now read the data
      call readvar(filen,varname,'level',istart,icount,
     >  cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat.lt.0) then
        write(*,*) 'ERROR IN READ_ONE_MANAGEMENT_VARIABLE'
        write(*,*) 'while reading ', trim(filen)
        stop 1
      end if
      call arr2vec (cdummy, dat)

      return
      end

c ---------------------------------------------------------------------
      subroutine set_irrigation_map(pathn, iwest, jnorth)
c ---------------------------------------------------------------------
c
c Set spatial irrigation map, telling us whether each grid cell is
c irrigated. If irrigate_flag is 0 or 1, this map is constant everywhere
c (0 or 1, respectively). However, if irrigate_flag is 2, then we read
c this map from a file. In that case, the file should contain values of
c 0 (unirrigated) or 1 (irrigated), as floating point values (although
c we also round fractional values to the nearest integer).
c
      use comgrid
      use compar
      use comwork
      use comcrop, only : irrigate_flag, irrigate
c
      implicit none
c
c Arguments
      character*(*) pathn  ! path in which irrigation file is located

      integer 
     >  iwest,             ! first longitude index to read
     >  jnorth             ! first latitude index to read

c
c Local variables
c
      character *80 filen  ! name of input file, relative to pathn
      
      integer istart(4),   ! for reading vars
     >        icount(4)    

      integer istat        ! netcdf error flag

      parameter (filen = 'irrigated.nc')

      data istart / 1,1,1,1 /, icount / nlon,nlat,1,1 /
      istart(1) = iwest
      istart(2) = jnorth
      icount(1) = nlonsub
      icount(2) = nlatsub

      if (irrigate_flag .eq. 0) then  ! no irrigation anywhere
        irrigate(:) = 0
      else if (irrigate_flag .eq. 1) then  ! irrigation everywhere
        irrigate(:) = 1
      else if (irrigate_flag .eq. 2) then  ! spatially-varying irrigation
        call readvar((trim(pathn) // '/' // trim(filen)),'irrigated','level',istart,icount,
     >    cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
        if (istat.lt.0) then
          write(*,*) 'ERROR IN SET_IRRIGATION_MAP'
          write(*,*) 'while reading ', trim(filen)
          stop 1
        end if
        call arr2vec (cdummy, work(1:npoi))
        irrigate(:) = nint(work(1:npoi))
      else
        write(*,*) 'ERROR IN SET_IRRIGATION_MAP'
        write(*,*) 'Unknown value of irrigate_flag: ', irrigate_flag
        stop 1
      end if

      return
      end



c ---------------------------------------------------------------------
      subroutine restart (iyrlast)
c ---------------------------------------------------------------------
c
c reads in restart files, initializes some variables
c
c this subroutine reads the restart values of:
c
c  fsnocov = fractional snow cover
c  tsno    = temperature of snow
c  hsno    = snow depth
c  tsoi    = soil temperature
c  wisoi   = soil ice content
c  wsoi    = soil moisture content
c  smsoil = immobile soil inorganic nitrogen pool for plant uptake
c  smsoln = mobile soil inorganic nitrogen pool for plant uptake, and leachable 
c  cbiol   = carbon in leaf biomass pool
c  cbiow   = carbon in woody biomass pool
c  cbior   = carbon in fine root biomass pool
c  cntops  = c/n ratio of plant tops (crops)
c  cnroot  = c/n ratio of crop roots
c  sapfrac = sapwood fraction
c  decompl = litter decomposition factor
c  decomps = soil decomposition factor
c  clitlm  = leaf metabolic litter
c  clitls  = leaf structural litter
c  clitll  = leaf lignin litter
c  clitrm  = root metabolic litter
c  clitrs  = root structural litter
c  clitrl  = root lignin litter
c  clitwm  = woody metabolic litter
c  clitws  = woody structural litter
c  clitwl  = woody lignin litter
c  falll   = annual leaf litterfall 
c  fallr   = annual fine root turnover
c  fallw   = annual woody turnover
c  totcmic = total microbial carbon
c  csoislop= slow soil carbon, protected humus
c  csoislon= slow soil carbon, nonprotected humus
c  csoipas = passive soil carbon
c  gdd0    = growing degree days 0
c  gdd0c   = growing degree days 0 for wheat - accumulated April 1 - Sept 30
c  gdd5    = growing degree days 5
c  gdd8    = growing degree days 8 
c  gdd10   = growing degree days 10 
c  aplantn = available plant nitrogen
c  tc      = coldest monthly temperature
c  tw      = warmest monthly temperature
c  wipud   = ice content of puddles per soil area
c  wpud    = liquid content of puddles per soil area
c  agddu   = annual accumulated growing degree days for bud burst, upper canopy
c  agddl   = annual accumulated growing degree days for bud burst, lower canopy
c  tempu   = cold-phenology trigger for trees
c  templs  = cold-phenology trigger for shrubs
c  greenfracl3 = fraction of green vegetation in C3 grasses
c  greenfracl4 = fraction of green vegetation in C4 grasses
c  Tavgann = average annual air temperature
c  PPTavgann = average annual precipitation
c  a10td    = 10-day avg daily temp
c  a10ts    = 10-day avg daily soil (1) temp
c  a10ancub = 10-day average canopy photosynthesis rate - broadleaf
c  a10ancuc = 10-day average canopy photosynthesis rate - conifer
c  a10ancls = 10-day average canopy photosynthesis rate - shrubs
c  a10ancl4 = 10-day average canopy photosynthesis rate - c4 grasses
c  a10ancl3 = 10-day average canopy photosynthesis rate - c3 grasses
c  a10scalparamu = 10-day average canopy scaling parameter - upper canopy
c  a10scalparaml = 10-day average canopy scaling parameter - lower canopy
c  a10daylightu = 10-day average daylight - upper canopy
c  a10daylightl = 10-day average daylight - lower canopy
c  a11soiltd = 11-day average surface soil temperature
c  a3tdmin = 3-day daily minimum air temperature
c (NOTE: a10ancuc is not used at this point, so its restart entry 
c is commented out)
c
      use comgrid
      use compar
      use comsoi
      use comsno
      use combcs
      use comsum
      use comveg
      use comwork
      use comcrop
      use comnitr
c
      implicit none
c
c Arguments
c
      integer iyrlast
c
c Local variables
c
      integer istart(4), icount(4) ! for reading restart vars
c
      character*20 fdir  ! used to construct odd/even file names
      character*80 filen ! file name
c
      integer  lf,           ! number of characters in directory name
     >         i,            ! loop indices
     >         istat,        ! error flag for netcdf
     >         nlevel        ! loop indice on the soil layer
c
c External
c
      integer lenchr,       ! Function: Find length of character string
     >  NF_PUT_ATT_TEXT,    ! netcdf function
     >  NF_GLOBAL           ! netcdf function
c
c ---------------------------------------------------------------------
c
      data istart / 1,1,1,1 /, icount / nlon,nlat,1,1 /
      icount(1) = nlonsub
      icount(2) = nlatsub

      !assign appropriate istart for proc* to read the base files
      if(flg_mpi == 1) then
        istart(1) = udg_is-istart_all+1
        istart(2) = udg_js-jstart_all+1
      end if

c
c check to see if iyrlast is odd or even and read from the appropriate
c files for that year
c
      if (mod(iyrlast,2) .eq. 0) then
         !fdir = trim(myprocDir)//'restart/even'
         fdir = 'restart/even'
      else
         !fdir = trim(myprocDir)//'restart/odd'
         fdir = 'restart/odd'
      end if
      lf = lenchr(fdir)

      !Echo info
      write(6,'(\,A)') 'Reading restart data at path: '//trim(fdir)
cc
cc dummy variable example, 3-d - copy and change for ne variable
cc
c      icount(3) = 1
c      filen = fdir(1:lf)//'dummyv.nc'
c      call readvar(filen,'dummyv','',istart,icount,
c     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in restart, dummyvv'
c         stop 1
c      end if
c      call arr2vec (cdummy, dummyv)
cc
cc dummy variable example, 4-d, 3rd dim (snowlayer) = nsnolay
cc
c      icount(3) = nsnolay
c      filen = fdir(1:lf)//'dummyv.nc'
c      call readvar(filen,'dummyv','snowlayer',istart,icount,
c     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
c      if (istat .ne. 0) then
c         write(*,*) 'ERROR in restart, dummyv'
c         stop 1
c      end if
c      do 5 nlevel = 1,nsnolay
c         call arr2vec (cdummy((nlevel-1)*nlonsub*nlatsub + 1),
c     >    dummyv(1,nlevel))
c 5    continue
c
c fsnocov
c
      filen = fdir(1:lf)//'fsnocov.nc'
      icount(3) = 1
      call readvar(filen,'fsnocov','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, fsnocov'
         stop 1
      end if
      call arr2vec (cdummy, fi)
c
c nsnolay variables: tsno and hsno
c
      icount(3) = nsnolay
c
      filen = fdir(1:lf)//'tsno.nc'
      call readvar(filen,'tsno','snowlayer',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, tsno'
         stop 1
      end if
      do 5 nlevel = 1,nsnolay
         call arr2vec (cdummy((nlevel-1)*nlonsub*nlatsub + 1),
     >    tsno(1,nlevel))

         call restart_correction('tsno',tsno(1,nlevel),273.16,50.0)
 5    continue
c
      filen = fdir(1:lf)//'hsno.nc'
      call readvar(filen,'hsno','snowlayer',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, hsno'
         stop 1
      end if
      do 10 nlevel = 1,nsnolay
         call arr2vec (cdummy((nlevel-1)*nlonsub*nlatsub + 1),
     >    hsno(1,nlevel))
 10   continue
c
c nsoilay variables: tsoi, wisoi, wsoi
c and soil nitrogen
c
      icount(3) = nsoilay
c
      filen = fdir(1:lf)//'tsoi.nc'
      call readvar(filen,'tsoi','soillayer',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, tsoi'
         stop 1
      end if
      do 15 nlevel = 1,nsoilay
         call arr2vec (cdummy((nlevel-1)*nlonsub*nlatsub + 1),
     >    tsoi(1,nlevel))

         call restart_correction('tsoi',tsoi(1,nlevel),276.13,100.0)

 15   continue
c
      filen = fdir(1:lf)//'wisoi.nc'
      call readvar(filen,'wisoi','soillayer',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, wisoi'
         stop 1
      end if
      do 20 nlevel = 1,nsoilay
         call arr2vec (cdummy((nlevel-1)*nlonsub*nlatsub + 1),
     >    wisoi(1,nlevel))
 20   continue
c
      filen = fdir(1:lf)//'wsoi.nc'
      call readvar(filen,'wsoi','soillayer',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, wsoi'
         stop 1
      end if
      do 25 nlevel = 1,nsoilay
         call arr2vec (cdummy((nlevel-1)*nlonsub*nlatsub + 1),
     >    wsoi(1,nlevel))
 25   continue
c
c soil nitrogen variables
c
      filen = fdir(1:lf)//'smsoil.nc'
      call readvar(filen,'smsoil','soillayer',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, smsoil'
         stop 1
      end if
      do 26 nlevel = 1,nsoilay
         call arr2vec (cdummy((nlevel-1)*nlonsub*nlatsub + 1),
     >    smsoil(1,nlevel))
 26   continue
c
      filen = fdir(1:lf)//'smsoln.nc'
      call readvar(filen,'smsoln','soillayer',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, smsoln'
         stop 1
      end if
      do 27 nlevel = 1,nsoilay
         call arr2vec (cdummy((nlevel-1)*nlonsub*nlatsub + 1),
     >    smsoln(1,nlevel))
 27   continue
c
c npft variables
c
      icount(3) = npft
c
      filen = fdir(1:lf)//'cbiol.nc'
      call readvar(filen,'cbiol','pft',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, cbiol'
         stop 1
      end if
      do 30 nlevel = 1,npft
         call arr2vec (cdummy((nlevel-1)*nlonsub*nlatsub + 1),
     >    cbiol(1,nlevel))
 30   continue
c
      filen = fdir(1:lf)//'cbiow.nc'
      call readvar(filen,'cbiow','pft',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, cbiow'
         stop 1
      end if
      do 35 nlevel = 1,npft
         call arr2vec (cdummy((nlevel-1)*nlonsub*nlatsub + 1),
     >    cbiow(1,nlevel))
 35   continue
c
      filen = fdir(1:lf)//'cbior.nc'
      call readvar(filen,'cbior','pft',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, cbior'
         stop 1
      end if
      do 40 nlevel = 1,npft
         call arr2vec (cdummy((nlevel-1)*nlonsub*nlatsub + 1),
     >    cbior(1,nlevel))
 40   continue
c
      filen = fdir(1:lf)//'cntops.nc'
      call readvar(filen,'cntops','pft',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, cntops'
         stop 1
      end if
      do 42 nlevel = 1,npft
         call arr2vec (cdummy((nlevel-1)*nlonsub*nlatsub + 1),
     >    cntops(1,nlevel))
 42   continue
c
      filen = fdir(1:lf)//'cnroot.nc'
      call readvar(filen,'cnroot','pft',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, cnroot'
         stop 1
      end if
      do 45 nlevel = 1,npft
         call arr2vec (cdummy((nlevel-1)*nlonsub*nlatsub + 1),
     >    cnroot(1,nlevel))
 45   continue
c
c single level variables
c
      icount(3) = 1
c
      filen = fdir(1:lf)//'sapfrac.nc'
      call readvar(filen,'sapfrac','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, sapfrac'
         stop 1
      end if
      call arr2vec (cdummy, sapfrac)
c
      filen = fdir(1:lf)//'decompl.nc'
      call readvar(filen,'decompl','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, decompl'
         stop 1
      end if
      call arr2vec (cdummy, decompl)
c
      filen = fdir(1:lf)//'decomps.nc'
      call readvar(filen,'decomps','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, decomps'
         stop 1
      end if
      call arr2vec (cdummy, decomps)
c
      filen = fdir(1:lf)//'clitlm.nc'
      call readvar(filen,'clitlm','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, clitlm'
         stop 1
      end if
      call arr2vec (cdummy, clitlm)
c
      filen = fdir(1:lf)//'clitls.nc'
      call readvar(filen,'clitls','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, clitls'
         stop 1
      end if
      call arr2vec (cdummy, clitls)
c
      filen = fdir(1:lf)//'clitll.nc'
      call readvar(filen,'clitll','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, clitll'
         stop 1
      end if
      call arr2vec (cdummy, clitll)
c
      filen = fdir(1:lf)//'clitrm.nc'
      call readvar(filen,'clitrm','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, clitrm'
         stop 1
      end if
      call arr2vec (cdummy, clitrm)
c
      filen = fdir(1:lf)//'clitrs.nc'
      call readvar(filen,'clitrs','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, clitrs'
         stop 1
      end if
      call arr2vec (cdummy, clitrs)
c
      filen = fdir(1:lf)//'clitrl.nc'
      call readvar(filen,'clitrl','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, clitrl'
         stop 1
      end if
      call arr2vec (cdummy, clitrl)
c
      filen = fdir(1:lf)//'clitwm.nc'
      call readvar(filen,'clitwm','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, clitwm'
         stop 1
      end if
      call arr2vec (cdummy, clitwm)
c
      filen = fdir(1:lf)//'clitws.nc'
      call readvar(filen,'clitws','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, clitws'
         stop 1
      end if
      call arr2vec (cdummy, clitws)
c
      filen = fdir(1:lf)//'clitwl.nc'
      call readvar(filen,'clitwl','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, clitwl'
         stop 1
      end if
      call arr2vec (cdummy, clitwl)
c
      filen = fdir(1:lf)//'falll.nc'
      call readvar(filen,'falll','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, falll'
         stop 1
      end if
      call arr2vec (cdummy, falll)
c
      filen = fdir(1:lf)//'fallr.nc'
      call readvar(filen,'fallr','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, fallr'
         stop 1
      end if
      call arr2vec (cdummy, fallr)
c
      filen = fdir(1:lf)//'fallw.nc'
      call readvar(filen,'fallw','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, fallw'
         stop 1
      end if
      call arr2vec (cdummy, fallw)
c
      filen = fdir(1:lf)//'totcmic.nc'
      call readvar(filen,'totcmic','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, totcmic'
         stop 1
      end if
      call arr2vec (cdummy, totcmic)
c
      filen = fdir(1:lf)//'csoislop.nc'
      call readvar(filen,'csoislop','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, csoislop'
         stop 1
      end if
      call arr2vec (cdummy, csoislop)
c
      filen = fdir(1:lf)//'csoislon.nc'
      call readvar(filen,'csoislon','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, csoislon'
         stop 1
      end if
      call arr2vec (cdummy, csoislon)
c
      filen = fdir(1:lf)//'csoipas.nc'
      call readvar(filen,'csoipas','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, csoipas'
         stop 1
      end if
      call arr2vec (cdummy, csoipas)
c
      filen = fdir(1:lf)//'gdd0.nc'
      call readvar(filen,'gdd0','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, gdd0'
         stop 1
      end if
      call arr2vec (cdummy, gdd0)
c
      filen = fdir(1:lf)//'gdd0c.nc'
      call readvar(filen,'gdd0c','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, gdd0c'
         stop 1
      end if
      call arr2vec (cdummy, gdd0c)
c
      filen = fdir(1:lf)//'gdd5.nc'
      call readvar(filen,'gdd5','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, gdd5'
         stop 1
      end if
      call arr2vec (cdummy, gdd5)
c
      filen = fdir(1:lf)//'gdd8.nc'
      call readvar(filen,'gdd8','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, gdd8'
         stop 1
      end if
      call arr2vec (cdummy, gdd8)
c
      filen = fdir(1:lf)//'gdd10.nc'
      call readvar(filen,'gdd10','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, gdd10'
         stop 1
      end if
      call arr2vec (cdummy, gdd10)
c
      filen = fdir(1:lf)//'aplantn.nc'
      call readvar(filen,'aplantn','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, aplantn'
         stop 1
      end if
      call arr2vec (cdummy, aplantn)
c
      filen = fdir(1:lf)//'tc.nc'
      call readvar(filen,'tc','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, tc'
         stop 1
      end if
      call arr2vec (cdummy, tc)
c
      filen = fdir(1:lf)//'tw.nc'
      call readvar(filen,'tw','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, tw'
         stop 1
      end if
      call arr2vec (cdummy, tw)
c
      filen = fdir(1:lf)//'wipud.nc'
      call readvar(filen,'wipud','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, wipud'
         stop 1
      end if
      call arr2vec (cdummy, wipud)
c
      filen = fdir(1:lf)//'wpud.nc'
      call readvar(filen,'wpud','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, wpud'
         stop 1
      end if
      call arr2vec (cdummy, wpud)
c
      filen = fdir(1:lf)//'agddu.nc'
      call readvar(filen,'agddu','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, agddu'
         stop 1
      end if
      call arr2vec (cdummy, agddu)
c
      filen = fdir(1:lf)//'agddl.nc'
      call readvar(filen,'agddl','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, agddl'
         stop 1
      end if
      call arr2vec (cdummy, agddl)
c
      filen = fdir(1:lf)//'tempu.nc'
      call readvar(filen,'tempu','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, tempu'
         stop 1
      end if
      call arr2vec (cdummy, tempu)
c
      filen = fdir(1:lf)//'templs.nc'
      call readvar(filen,'templs','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, templs'
         stop 1
      end if
      call arr2vec (cdummy, templs)
c
      filen = fdir(1:lf)//'greenfracl3.nc'
      call readvar(filen,'greenfracl3','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, greenfracl3'
         stop 1
      end if
      call arr2vec (cdummy, greenfracl3)
c
      filen = fdir(1:lf)//'greenfracl4.nc'
      call readvar(filen,'greenfracl4','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, greenfracl4'
         stop 1
      end if
      call arr2vec (cdummy, greenfracl4)
c
      filen = fdir(1:lf)//'Tavgann.nc'
      call readvar(filen,'Tavgann','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, Tavgann'
         stop 1
      end if
      call arr2vec (cdummy, Tavgann)
c
      filen = fdir(1:lf)//'PPTavgann.nc'
      call readvar(filen,'PPTavgann','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, PPTavgann'
         stop 1
      end if
      call arr2vec (cdummy, PPTavgann)
c
      filen = fdir(1:lf)//'a10td.nc'
      call readvar(filen,'a10td','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, a10td'
         stop 1
      end if
      call arr2vec (cdummy, a10td)
c
      filen = fdir(1:lf)//'a10ts.nc'
      call readvar(filen,'a10ts','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, a10ts'
         stop 1
      end if
      call arr2vec (cdummy, a10ts)
c
      filen = fdir(1:lf)//'a10ancub.nc'
      call readvar(filen,'a10ancub','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, a10ancub'
         stop 1
      end if
      call arr2vec (cdummy, a10ancub)
ccc
cc      filen = fdir(1:lf)//'a10ancuc.nc'
cc      call readvar(filen,'a10ancuc','',istart,icount,
cc     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
cc      if (istat .ne. 0) then
cc         write(*,*) 'ERROR in restart, a10ancuc'
cc         stop 1
cc      end if
cc      call arr2vec (cdummy, a10ancuc)
c
      filen = fdir(1:lf)//'a10ancls.nc'
      call readvar(filen,'a10ancls','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, a10ancls'
         stop 1
      end if
      call arr2vec (cdummy, a10ancls)
c
      filen = fdir(1:lf)//'a10ancl4.nc'
      call readvar(filen,'a10ancl4','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, a10ancl4'
         stop 1
      end if
      call arr2vec (cdummy, a10ancl4)
c
      filen = fdir(1:lf)//'a10ancl3.nc'
      call readvar(filen,'a10ancl3','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, a10ancl3'
         stop 1
      end if
      call arr2vec (cdummy, a10ancl3)
c
      filen = fdir(1:lf)//'a10scalparamu.nc'
      call readvar(filen,'a10scalparamu','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, a10scalparamu'
         stop 1
      end if
      call arr2vec (cdummy, a10scalparamu)
c
      filen = fdir(1:lf)//'a10scalparaml.nc'
      call readvar(filen,'a10scalparaml','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, a10scalparaml'
         stop 1
      end if
      call arr2vec (cdummy, a10scalparaml)
c
      filen = fdir(1:lf)//'a10daylightu.nc'
      call readvar(filen,'a10daylightu','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, a10daylightu'
         stop 1
      end if
      call arr2vec (cdummy, a10daylightu)
c
      filen = fdir(1:lf)//'a10daylightl.nc'
      call readvar(filen,'a10daylightl','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, a10daylightl'
         stop 1
      end if
      call arr2vec (cdummy, a10daylightl)
c
      filen = fdir(1:lf)//'a11soiltd.nc'
      call readvar(filen,'a11soiltd','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, a11soiltd'
         stop 1
      end if
      call arr2vec (cdummy, a11soiltd)
c
      filen = fdir(1:lf)//'a3tdmin.nc'
      call readvar(filen,'a3tdmin','',istart,icount,
     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in restart, a3tdmin'
         stop 1
      end if
      call arr2vec (cdummy, a3tdmin)
c
c calculate tcmin
c
      do i = 1, npoi
         if(cmask(i)==0) cycle
         tcmin(i) = tc(i) + deltat(i)
      enddo
c
      call existence
c
      return
      end
c
c ---------------------------------------------------------------------
      subroutine coldstart
c ---------------------------------------------------------------------
c  
      use comgrid
      use compar
      use comsoi
      use comsno
c
      implicit none
c
c initialize some model variables for cold start conditions
c
      call const (fi, npoi, 0.0)
c
      call const (hsno, npoi*nsnolay, 0.0)
      call const (tsno, npoi*nsnolay, 273.16)
c
      call const (tsoi, npoi*nsoilay, 278.16)
      call const (wsoi, npoi*nsoilay, 0.50)
      call const (wisoi, npoi*nsoilay, 0.00)
c
c return to main program
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine rdanom(imonth,iyear,iy2,iwest,jnorth,next_month)
c ---------------------------------------------------------------------
c
c reads in monthly values for imonth (if next_month=false) or imonth+1
c (if next_month=true)
c The naming convention reflects its old behavior, which was to read in
c anomalies and add them to the climatological values
c
c this subroutine reads the monthly values of:
c
c temp - mean temperature (degC)
c trng - temperature range  (degC)
c tmin - min daily temperature (degC)
c tmax - max daily temperature (degC)
c prec - precipitation rate (mm/day)
c cld  - cloudiness (percent)
c rads - solar radiation (W m-2)
c rh   - relative humidity (percent)
c wspd - wind speed (m s-1)
c wetd - wet days per month
c
c  
      use clim_file_utils
      use comgrid
      use compar
      use cominput
      use combcs
      use comveg
      use comwork
c
      implicit none
c
c Arguments
c
      integer imonth, ! month
     >        iyear,  ! year
     >        iy2,    ! final year of the run
     >        iwest,  ! 1st lon index for subset
     >        jnorth  ! 1st lat index for subset

c If next_month is true, we read imonth+1 rather than imonth. This is
c used for the weather generator, for which we need the previous and
c next month's data in addition to the current month's data.
      logical next_month

c
c Local variables
c
      real anom(npoi)
c
      integer imon,   ! month to read: either current month or month + 1
     >        iyear_rd, ! year to read (handling possible shift into next year when add 1 to month)
     >        iyr,    ! number of years after begin of data file(istyrm)
     >        iyr_cal,! calendar year corresponding to iyr (similar to iyear_rd, but handling looping)
     >        istat,  ! error flag for netcdf
     >        i       ! loop indice on land points
c
      integer istart(4), icount(4) ! for reading rdanom vars
      character*80 filen
c
c ---------------------------------------------------------------------
c
      data istart / 1,1,1,1 /, icount / nlon,nlat,1,1 /
      istart(1) = iwest
      istart(2) = jnorth
      icount(1) = nlonsub
      icount(2) = nlatsub
c
c If timestep equals December of final year in anomaly file and also
c equals the final year of the run, and next_month is true, then set
c January anomalies to zero and return.
c
      if (iyear .eq. istyrm+nanom-1 .and. imonth .eq. 12 .and. next_month) then
        if (iyear .eq. iy2) then
c
          print *, 'WARNING: last month of run; no anomalies for January'
          print *, 'Using climatologies for month year ='
          print *, imonth+1,iyear+1
          do 4 i = 1, npoi
            if(cmask(i) ==0) cycle

            xint(i,1) = clmt(i,1)
            xintrng(i,1) = clmtrng(i,1)
            xintmin(i,1) = clmtmin(i,1)
            xintmax(i,1) = clmtmax(i,1)
            xinprec(i,1) = clmprec(i,1)
            xincld(i,1) = clmcld(i,1)
            xinrads(i,1) = clmrads(i,1)
            xinq(i,1) = clmq(i,1)
            xinwind(i,1) = clmw(i,1)
            xinwet(i,1) = clmwet(i,1)
 4        continue
          return
        end if
      end if
c
c determine which month to read
c
      iyear_rd = iyear
      imon = imonth
      if (next_month) then
        imon = imon + 1
        if (imon .eq. 13) then
          imon = 1
          iyear_rd = iyear_rd + 1
        end if
      end if
c
c Handle cycling through the years in the anomaly file, 
c if run continues past the end of the file
      iyr = year_to_read(iyear_rd, istyrm, overlap_start, overlap_end)
      if (iyr .lt. 0) then
        write(*,*) 'ERROR in rdanom: tried to read a year prior to the start of the files'
        write(*,*) 'iyear_rd = ', iyear_rd
        write(*,*) 'istyrm = ', istyrm
        stop 1
      end if
      istart(4) = iyr*12 + imon
      iyr_cal = iyr + istyrm

      if (iyr .eq. (overlap_start - istyrm) .and. iyear_rd .gt.
     >     overlap_start .and. imon .eq. 1) then
        print *, 'WARNING: Attempted to read past last month in anomaly file'
        print *, 'Looping back to the beginning of the file'
      end if
c
c      if (istart(4) .gt. 0) then
c        print *, 'rdanom reading month year step ='
c        print *, imon,iyr+istyrm,istart(4)
c      else
c        print *, 'WARNING, anomalies begin in year ',istyrm
c        print *, 'Not reading in anomalies for month year ='
c        print *, imon,iyr+istyrm
c        return
c      end if
c
c
c Note: For all variables, we only read the variable if the associated file name is not ' '.
c If the file name is ' ', we'll use the values that are already in the xin* variable,
c which are the values copied from the climatologies in readit. For this reason, it is important
c that all monthly xin* variables are set in readit!
c
c Also, note that, for variables that are only read if certain flags are set 
c (e.g., rads is only read if read_radiation is true), we do not explicitly set
c these variables to un_init here if they are not read. This is because they should
c already have been set to un_init in readit.
c
cc
cc     dummy variable example, 4-d, whose 3rd dim (level) = 1
cc
c      if (fn_dummyv_mon .ne. ' ') then
c        aname = var_dummyv_mon
c        call climate_full_fname(filen, fn_dummyv_mon, iyr_cal, .false.)  ! sets filen variable
c        call readvar(filen,aname,level_var,istart,icount,
c       > work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
c       > cdummy(3*nlonsub+1),istat)
c        if (istat .ne. 0) then
c           write(*,*) 'ERROR in rdanom, dummyv'
c           stop 1
c        end if
c        call arr2vec( work, anom )
c        do 10 i = 1, npoi
c           xindummyv(i,imon) = max (clmtrng(i,imon) + anom(i), 0.1)
c 10     continue
c      end if
c
c
c     mean temperature
c
      if (read_tmin_tmax) then
c     Read tmin and tmax; calculate t and trng
c Note: As of 01.28.10, the code in this 'if' clause is untested

        if (fn_tmin_mon .ne. ' ') then
          aname = var_tmin_mon
          call climate_full_fname(filen, fn_tmin_mon, iyr_cal, .false.)
          call readvar(filen,aname,level_var,istart,icount,
     >         work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     >         cdummy(3*nlonsub+1),istat)
          if (istat .ne. 0) then
            write(*,*) 'ERROR in rdanom, tmin'
            stop 1
          end if
          call arr2vec( work, anom )
          do i = 1, npoi
            if(cmask(i)==0) cycle
            
            xintmin(i,imon) = anom(i)
          end do
        end if        

        if (fn_tmax_mon .ne. ' ') then
          aname = var_tmax_mon
          call climate_full_fname(filen, fn_tmax_mon, iyr_cal, .false.)
          call readvar(filen,aname,level_var,istart,icount,
     >         work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     >         cdummy(3*nlonsub+1),istat)
          if (istat .ne. 0) then
            write(*,*) 'ERROR in rdanom, tmax'
            stop 1
          end if
          call arr2vec( work, anom )
          do i = 1, npoi
            if(cmask(i)==0) cycle
            xintmax(i,imon) = anom(i)
          end do
        end if   

        do i = 1, npoi
          if(cmask(i)==0) cycle
          xintrng(i, imon) = xintmax(i, imon) - xintmin(i, imon)
          xint(i, imon) = tmax_wt * xintmax(i, imon) + tmin_wt * xintmin(i, imon)
        end do

      else  ! .not. read_tmin_tmax
c     Read t and trng; calculate tmin and tmax

        if (fn_temp_mon .ne. ' ') then
          aname = var_temp_mon
          call climate_full_fname(filen, fn_temp_mon, iyr_cal, .false.)
          call readvar(filen,aname,level_var,istart,icount,
     >         work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     >         cdummy(3*nlonsub+1),istat)
          if (istat .ne. 0) then
            write(*,*) 'ERROR in rdanom, temp'
            stop 1
          end if
          call arr2vec( work, anom )
          do 5 i = 1, npoi
            if(cmask(i)==0) cycle
c     old CRU TS 1 dataset         xint(i,imon) = clmt(i,imon) + anom(i)
            xint(i,imon) = anom(i)
 5        continue
        end if 
c     c
c     c     temperature range
c     c
        if (fn_dtr_mon .ne. ' ') then
          aname = var_dtr_mon
          call climate_full_fname(filen, fn_dtr_mon, iyr_cal, .false.)
          call readvar(filen,aname,level_var,istart,icount,
     >         work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     >         cdummy(3*nlonsub+1),istat)
          if (istat .ne. 0) then
            write(*,*) 'ERROR in rdanom, trange'
            stop 1
          end if
          call arr2vec( work, anom )
          do 10 i = 1, npoi
            if(cmask(i)==0) cycle
c     xintrng(i,imon) = max (clmtrng(i,imon) + anom(i), 0.1)
            xintrng(i,imon) = max (anom(i), 0.1)
 10       continue
        end if

        do i = 1, npoi
          if(cmask(i)==0) cycle
          xintmin(i, imon) = xint(i, imon) - tmax_wt * xintrng(i, imon)
          xintmax(i, imon) = xint(i, imon) + tmin_wt * xintrng(i, imon)
        end do

      end if  ! read_tmin_tmax
c
c     precipitation rate
c
      if (fn_prec_mon .ne. ' ') then
        aname = var_prec_mon
        call climate_full_fname(filen, fn_prec_mon, iyr_cal, .false.)
        call readvar(filen,aname,level_var,istart,icount,
     >       work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     >       cdummy(3*nlonsub+1),istat)
        if (istat .ne. 0) then
          write(*,*) 'ERROR in rdanom, prec'
          stop 1
        end if
        call arr2vec( work, anom )
        do 15 i = 1, npoi
          if(cmask(i)==0) cycle
c     xinprec(i,imon) = clmprec(i,imon) + anom(i)
          xinprec(i,imon) = anom(i)
 15     continue
      end if
c
      if (read_radiation) then
c     Read rads
c Note: As of 01.28.10, the code in this 'if' clause is untested
        if (fn_rads_mon .ne. ' ') then
          aname = var_rads_mon
          call climate_full_fname(filen, fn_rads_mon, iyr_cal, .false.)
          call readvar(filen,aname,level_var,istart,icount,
     >         work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     >         cdummy(3*nlonsub+1),istat)
          if (istat .ne. 0) then
            write(*,*) 'ERROR in rdanom, rads'
            stop 1
          end if
          call arr2vec( work, anom )
          do i = 1, npoi
            if(cmask(i)==0) cycle
            xinrads(i,imon) = anom(i)
          end do
        end if

      else  ! .not. read_radiation
c     Read cloud
c
        if (fn_cloud_mon .ne. ' ') then
          aname = var_cloud_mon
          call climate_full_fname(filen, fn_cloud_mon, iyr_cal, .false.)
          call readvar(filen,aname,level_var,istart,icount,
     >         work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     >         cdummy(3*nlonsub+1),istat)
          if (istat .ne. 0) then
            write(*,*) 'ERROR in rdanom, cld'
            stop 1
          end if
          call arr2vec( work, anom )
          do 20 i = 1, npoi
            if(cmask(i)==0) cycle
c     xincld(i,imon) = clmcld(i,imon) + anom(i)
            xincld(i,imon) = anom(i)
 20       continue
        end if

      end if  ! read_radiation
c
c     relative humidity
c
      if (fn_rh_mon .ne. ' ') then
        aname = var_rh_mon
        call climate_full_fname(filen, fn_rh_mon, iyr_cal, .false.)
        call readvar(filen,aname,level_var,istart,icount,
     >       work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     >       cdummy(3*nlonsub+1),istat)
        if (istat .ne. 0) then
          write(*,*) 'ERROR in rdanom, rh'
          stop 1
        end if
        call arr2vec( work, anom )
        do 25 i = 1, npoi
          if(cmask(i)==0) cycle
c     xinq(i,imon) = clmq(i,imon) + anom(i)
          xinq(i,imon) = anom(i)
 25     continue
      end if
c
c     wind speed
c
      if (fn_wspd_mon .ne. ' ') then
        aname = var_wspd_mon
        call climate_full_fname(filen, fn_wspd_mon, iyr_cal, .false.)
        call readvar(filen,aname,level_var,istart,icount,
     >       work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     >       cdummy(3*nlonsub+1),istat)
        if (istat .ne. 0) then
          write(*,*) 'ERROR in rdanom, wspd'
          stop 1
        end if
        call arr2vec( work, anom )
        do 30 i = 1, npoi
          if(cmask(i)==0) cycle
c     xinwind(i,imon) = clmw(i,imon) + anom(i)
          xinwind(i,imon) = anom(i)
 30     continue
      end if
c
c     wet days
c
      if (fn_wetd_mon .ne. ' ') then
        aname = var_wetd_mon
        call climate_full_fname(filen, fn_wetd_mon, iyr_cal, .false.)
        call readvar(filen,aname,level_var,istart,icount,
     >       work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     >       cdummy(3*nlonsub+1),istat)
        if (istat .ne. 0) then
          write(*,*) 'ERROR in rdanom, wetd'
          stop 1
        end if
        call arr2vec( work, anom )
        do 35 i = 1, npoi
          if(cmask(i)==0) cycle
c          xinwet(i,imon) = clmwet(i,imon) + anom(i)
          xinwet(i,imon) = anom(i)
 35     continue
      end if
c
c
      return
      end
c
c
c
c ---------------------------------------------------------------------
      subroutine inird(file,istyr)
c ---------------------------------------------------------------------
c
c This subroutine obtains the year+1 of the year in the units attribute
c of a file.  Low-level netcdf commands are used.
c
      implicit none
c
      include 'netcdf.inc'
c
c Arguments
c
      character*(*) file
      integer istyr       
c
c Local Variables
c
      integer idies, istat, idtime, lf1
c
c      integer NF_OPEN,        ! netcdf function
c     >        NF_NOWRITE,     ! '
c     >        NF_NOERR,       ! '
c     >        NF_INQ_VARID,   ! '
c     >        NF_GET_ATT_TEXT ! '
c              
      character*80 units
c ---------------------------------------------------------------------
c
c open file
c
      istat = NF_OPEN(file,NF_NOWRITE,idies)
      if (istat .ne. NF_NOERR) then
         print *, 'Error in inird while trying to open file'
         print *, file
         print *, NF_STRERROR(istat)
         istyr = -1
         return
      end if
c
c get units attribute for time
c
      istat = NF_INQ_VARID(idies,'time',idtime)
      if (istat .ne. NF_NOERR) then
         print *, 'Error in inird while trying to get time id'
         print *, NF_STRERROR(istat)
         istyr = -1
         return
      end if
      units = ' '
      istat = NF_GET_ATT_TEXT(idies,idtime,'units',units)
      if (istat .ne. NF_NOERR) then
         print *, 'Error in inird while trying to get time units'
         print *, NF_STRERROR(istat)
         istyr = -1
         return
      end if
c
c put character units year into integer variable, add 1
c
      lf1 = index(units,'since') + 6
      read(units(lf1:lf1+3),'(i4)') istyr
c      read(units(12:15),'(i4)') istyr
      istyr = istyr + 1
      return
      end
c
c ---------------------------------------------------------------------
      subroutine rdday(jday, imonth, iyear, iwest, jnorth)
c ---------------------------------------------------------------------
c
c This subroutine reads in daily fields..
c
      use clim_file_utils
      use comgrid
      use compar
      use cominput
      use combcs
      use comwork
c
      implicit none
c
c Arguments
c
      integer jday,   ! day of the year
     >        imonth, ! month
     >        iyear,  ! year
     >        iwest,  ! 1st lon index for subset
     >        jnorth  ! 1st lat index for subset
c
c Local variables
c
      integer i,      ! loop indice on years after istyrd
     >        istat,  ! error flag for netcdf
     >        iyr,    ! index of year to read (0 = 1st year in file, 1 = 2nd, etc., handling looping around when we reach end of file)
     >        iyr_cal ! calendar year to read (handling looping around when we reach end of file)
c
      character*80 filen
      character*39 pathn
      integer istart(4), icount(4)

c     first_time: saved variable to track whether this is the first time
c     the daily files are read. If so, we do extra things that only have
c     to be done once, such as setting unread variables to un_init. 
      logical first_time  
c
c External
c
      logical is_leap        ! Function: is the given year a leap year?
c ---------------------------------------------------------------------
c
c     pathn = '/Volumes/wet/kucharik/IBIS_input/input/'
c
      save first_time
      data first_time / .TRUE. / ! this starts as true, but will be set to false after the daily files are read for the first time
      data istart / 1,1,1,1 /, icount / nlon,nlat,1,1 /
      istart(1) = iwest
      istart(2) = jnorth
      icount(1) = nlonsub
      icount(2) = nlatsub
c
c
c calculate counter for looping through anomaly data 
c
c WJS (4-16-10): I believe that this could lead to unexpected behavior
c (or even a crash of the code) in some cases, by having the model read
c a leap year in a non-leap year, or vice versa. For example, if climate
c data are present for the years 1949-1951, then when the model year is
c 1952, we will try to get that year's data from 1949, which will be bad
c when we try to access day 366. I haven't looked at this closely,
c though, to say definitively whether this is a problem.
c
      iyr = year_to_read(iyear, istyrd, overlap_start, overlap_end)
      if (iyr .lt. 0) then
        write(*,*) 'ERROR in rdday: tried to read a year prior to the start of the files'
        write(*,*) 'iyear = ', iyear
        write(*,*) 'istyrd = ', istyrd
        stop 1
      end if
      iyr_cal = iyr + istyrd

      if(first_time) then  !Y.Li
        write(6,*) 'Reading daily climate data from year: ',iyr_cal
      end if !on if(first_time)
c

      if (iyear .lt. istyrd) then
        print *, 'daily data begins in year ', istyrd
        print *, 'not reading in daily data'
        return
      end if
c
c count how many days since beginning of daily data
c daily data begin on Jan 1, istyrd
c
      if (one_file_per_year .or. (iyr_cal .eq. istyrd)) then
        istart(4) = jday
      else
        istart(4) = 0
        do 10 i = istyrd, iyr_cal-1
          istart(4) = istart(4) + 365
          if (is_leap(i)) then
            istart(4) = istart(4) + 1
          end if
 10     continue
        istart(4) = istart(4) + jday
      end if
c
c      if (istart(4) .gt. 0) then
c        print *, 'rdday reading day month year step ='
c        print *, jday, imonth, iyr_cal,istart(4)
c      else
c        print *, 'WARNING, anomalies begin in year ',istyrd
c        print *, 'Not reading in anomalies for day month year ='
c        print *, jday, imonth, iyr_cal
c        return
c      end if
c
c
c read daily precip
c
      if (fn_prec_daily .eq. ' ') then
        write(*,*) 'ERROR in rdday:'
        write(*,*) 'fn_prec_daily not specified'
        stop 1
      end if
      aname = var_prec_daily
      call climate_full_fname(filen, fn_prec_daily, iyr_cal, one_file_per_year)
      call readvar(filen,aname,level_var,istart,icount,
     > work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     > cdummy(3*nlonsub+1),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in rdday, prec'
         stop 1
      end if
      call arr2vec( work, xinprecd(1) )

      if (read_tmin_tmax) then
c     Read tmin and tmax; set t and trng to un_init if this is the first time we're reading daily files
c     (Note: unlike monthly variables, we do NOT calculate t and trng here, 
c      partly because of the difficulties that arise if these variables represent anomalies)
c
c read daily temp min 
c
        if (fn_tmin_daily .eq. ' ') then
          write(*,*) 'ERROR in rdday:'
          write(*,*) 'fn_tmin_daily must be specified for read_tmin_tmax = true'
          stop 1
        end if
        aname = var_tmin_daily
        call climate_full_fname(filen, fn_tmin_daily, iyr_cal, one_file_per_year)
        call readvar(filen,aname,level_var,istart,icount,
     >       work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     >       cdummy(3*nlonsub+1),istat)
        if (istat .ne. 0) then
          write(*,*) 'ERROR in rdday, tmin'
          stop 1
        end if
        call arr2vec( work, xintmind(1) )
c     
c     read daily temp max
c     
        if (fn_tmax_daily .eq. ' ') then
          write(*,*) 'ERROR in rdday:'
          write(*,*) 'fn_tmax_daily must be specified for read_tmin_tmax = true'
          stop 1
        end if
        aname = var_tmax_daily
        call climate_full_fname(filen, fn_tmax_daily, iyr_cal, one_file_per_year)
        call readvar(filen,aname,level_var,istart,icount,
     >       work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     >       cdummy(3*nlonsub+1),istat)
        if (istat .ne. 0) then
          write(*,*) 'ERROR in rdday, tmax'
          stop 1
        end if
        call arr2vec( work, xintmaxd(1) )
        
        if (first_time) then
          call const(xintd, npoi, un_init)
          call const(xintrngd, npoi, un_init)
        end if

      else  ! .not. read_tmin_tmax
c     Read t and trng; set tmin and tmax to un_init if this is the first time we're reading daily files
c     (Note: unlike monthly variables, we do NOT calculate tmin and tmax here, 
c      partly because of the difficulties that arise if these variables represent anomalies)
c
c read daily temp
c
        if (fn_temp_daily .eq. ' ') then
          write(*,*) 'ERROR in rdday:'
          write(*,*) 'fn_temp_daily must be specified for read_tmin_tmax = false'
          stop 1
        end if
        aname = var_temp_daily
        call climate_full_fname(filen, fn_temp_daily, iyr_cal, one_file_per_year)
        call readvar(filen,aname,level_var,istart,icount,
     >       work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     >       cdummy(3*nlonsub+1),istat)
        if (istat .ne. 0) then
          write(*,*) 'ERROR in rdday, temp'
          stop 1
        end if
        call arr2vec( work, xintd(1) )
c     
c     read daily trange
c     
        if (fn_dtr_daily .eq. ' ') then
          write(*,*) 'ERROR in rdday:'
          write(*,*) 'fn_dtr_daily must be specified for read_tmin_tmax = false'
          stop 1
        end if
        aname = var_dtr_daily
        call climate_full_fname(filen, fn_dtr_daily, iyr_cal, one_file_per_year)
        call readvar(filen,aname,level_var,istart,icount,
     >       work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     >       cdummy(3*nlonsub+1),istat)
        if (istat .ne. 0) then
          write(*,*) 'ERROR in rdday, trange'
          stop 1
        end if
        call arr2vec( work, xintrngd(1) )

        if (first_time) then
          call const(xintmind, npoi, un_init)
          call const(xintmaxd, npoi, un_init)
        end if

      end if  ! read_tmin_tmax
c     
      if (read_radiation) then
c     Read rads, set cloud to un_init if this is the first time we're reading daily files
c     
c read daily radiation
c     
        if (fn_rads_daily .eq. ' ') then
          write(*,*) 'ERROR in rdday:'
          write(*,*) 'fn_rads_daily must be specified for read_radiation = true'
          stop 1
        end if
        aname = var_rads_daily
        call climate_full_fname(filen, fn_rads_daily, iyr_cal, one_file_per_year)
        call readvar(filen,aname,level_var,istart,icount,
     >       work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     >       cdummy(3*nlonsub+1),istat)
        if (istat .ne. 0) then
          write(*,*) 'ERROR in rdday, rads'
          stop 1
        end if
        call arr2vec( work, xinradsd(1) )

        if (first_time) then
          call const(xincldd, npoi, un_init)
        end if

      else  ! .not. read_radiation
c     Read cloud, set rads to un_init if this is the first time we're reading daily files
c
c read daily cloudiness
c     
        if (fn_cloud_daily .eq. ' ') then
          write(*,*) 'ERROR in rdday:'
          write(*,*) 'fn_cloud_daily must be specified for read_radiation = false'
          stop 1
        end if
        aname = var_cloud_daily
        call climate_full_fname(filen, fn_cloud_daily, iyr_cal, one_file_per_year)
        call readvar(filen,aname,level_var,istart,icount,
     >       work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     >       cdummy(3*nlonsub+1),istat)
        if (istat .ne. 0) then
          write(*,*) 'ERROR in rdday, cld'
          stop 1
        end if
        call arr2vec( work, xincldd(1) )

        if (first_time) then
          call const(xinradsd, npoi, un_init)
        end if
      end if  ! read_radiation
c     
c read daily windspeed
c
      if (fn_wspd_daily .eq. ' ') then
        write(*,*) 'ERROR in rdday:'
        write(*,*) 'fn_wspd_daily not specified'
        stop 1
      end if
      aname = var_wspd_daily
      call climate_full_fname(filen, fn_wspd_daily, iyr_cal, one_file_per_year)
      call readvar(filen,aname,level_var,istart,icount,
     > work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     > cdummy(3*nlonsub+1),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in rdday, wspd'
         stop 1
      end if
      call arr2vec( work, xinwindd(1) )
c
c read daily humidity
c
      if (fn_rh_daily .eq. ' ') then
        write(*,*) 'ERROR in rdday:'
        write(*,*) 'fn_rh_daily not specified'
        stop 1
      end if
      aname = var_rh_daily
      call climate_full_fname(filen, fn_rh_daily, iyr_cal, one_file_per_year)
      call readvar(filen,aname,level_var,istart,icount,
     > work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     > cdummy(3*nlonsub+1),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in rdday, sphum'
         stop 1
      end if
      call arr2vec( work, xinqd(1) )
c
c
      first_time = .FALSE.

      return
      end
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! purpose : Read-in NetCDF files
!           Support type: soil, pft, fertilizer 
!
! note    : Read-in data communicated through variable <cdummy>
!
!           Last updated 2014-06-18 by Yuwei Li
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      subroutine read_NetCDF(ftype,year)

      use compar ,only: nlonsub, nlatsub, nsoilay,udg_is,udg_js
      use comwork,only: work,cdummy

      implicit NONE

      !--input
      integer ::   year       !dummy when ftype="soil"
      character*(*)  ftype      !option: soil, pft, fert
      
      !--local
      integer:: dimlon,dimlat,dim3rd,ierr
      integer:: lon_is,lat_js !start index in lon and lat direction
      integer,dimension(4):: istart,icount

      character*4   year_char
      character*80  varname
      character*120 fpath,fname
      character*120 filen,name3rd

      !Initialization
      write(year_char,100) year
      fpath = 'input/'

      lon_is = udg_is
      lat_js = udg_js

      dimlon = nlonsub
      dimlat = nlatsub

      !Seems like function <trim> not working properly, use ftype(1:3) to make sure the comparison works
      if( trim(ftype) /= 'soil'.and.trim(ftype(1:3)) /= 'pft'.and.trim(ftype) /= 'fert'.and.trim(ftype) /= 'test' ) then
        write(6,*) '(read_NetCDF)Error! Available options are: soil, pft, fert '
        write(6,*) 'current input type is: ',trim(ftype)
        goto 9999
      end if
 
      if( trim(ftype) == 'soil') then

        fname   = 'soil_layers.nc'
        varname = 'soiltype'
        name3rd = 'layer'

        dim3rd  = nsoilay

        write(6,*) 'Reading soils with layered data from a NetCDF file: ',trim(fname)

      elseif( trim(ftype(1:3)) == 'pft') then

        fname   = 'pft'//'_'//trim(year_char)//'.nc'
        varname = 'vegtype'
        name3rd = 'layer'

        dim3rd  = 1
       
        write(6,*) 'Reading land use on a year-by-year basis from NetCDF file: ', trim(fname)

      elseif( trim(ftype) == 'fert' ) then

        fname   = 'fert'//'_'//trim(year_char)//'.nc'
        varname = 'frate'
        name3rd = 'layer'
        dim3rd  = 1

        write(6,*) 'Reading fertilizer application on a year-by-year basis from NetCDF files: ',trim(fname)

      elseif (trim(ftype) == 'test' ) then  !TEST ONLY, REMOVE LATER

        write(6,*) 'Reading sample NetCDF file...'
        
        fname   = 'sample'//'_'//trim(year_char)//'.nc'
        varname = 'data'
        name3rd = 'layer'

        dim3rd  = 2


      end if  !on if trim(ftype)


      filen = trim(fpath)//trim(fname)

      istart  = (/lon_is, lat_js,      1, 1 /)
      icount  = (/dimlon, dimlat, dim3rd, 1 /)

      !variable <cdummy> used to store data, <work> used as junk
      cdummy = 1e20
      call readvar(filen,varname,name3rd,istart,icount,cdummy,work(1),work(1),work(1),work(1),ierr)
      if(ierr < 0) write(6,*) 'Error in reading NetCDF file: ', trim(filen)


9999  continue

  100 format(i4.4)
      return
      end subroutine read_NetCDF
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! purpose : Example to write NetCDF input file for testing
!
! note    : Last updated 2014-06-18 by Yuwei Li
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!

      subroutine write_netCDF()
      implicit NONE

      integer,parameter:: dim3rd = 1    !(for soil: 11)
      integer,parameter:: nlon = 1070, nlat = 705
      real,parameter:: lon_s = -96.5, lon_e = -90.2, lat_s = 40.3, lat_e = 43.5
    
      !variables for initilizations (file, var)
      integer:: idies,ierror,i,j,k
      integer:: n3rd,ndims
      
      real:: valmissing,delt
      real:: alons(nlon),alats(nlat)
      real,dimension(dim3rd):: vals3rd
    
      character*120 fpath,fname
      character*120 filen,title,source,history,name3rd,long3rd,units3rd,pos3rd,tunits,calendar

      character*80  varname,longname,units
      character(len=20),dimension(4)::dimnames
    
      !variables for writevar
      integer,dimension(4)  :: istart, icount
      real                  :: ftime,tweight
      real,dimension(nlon*nlat*dim3rd)     :: cdummy
      real,dimension(nlon,nlat,dim3rd) :: new_mydata
      character*10 cdate
    
      !--Initialization
      
        fpath = '/Shared/spaklab/wsc/agro-ibis/inputs/iowa/templates/'
        !fname = 'soil_layers.nc'
        !fname = 'fert_2011.nc'
        fname = 'pft_2011.nc'
        varname     = 'vegtype'                    !'soiltype'   !'frate'          !'vegtype'
        longname    = 'vegetation'                   !'soil type'  !'N-fertilizer'   !'vegetation'
        units       = ''                            !''           !'kg/ha'          !''
      
        filen    = trim(fpath)//trim(fname)         !name for new file
        title    = 'User-defined yearly nc file'    !title for file
        source   = 'source'                         !source for file
        history  = '2015-07-07'                     !date of creation, and possibly author
        name3rd  = 'level'                          !name of 3rd dimension
        long3rd  = 'abstract layer variable'        !long name for 3rd dimension variable
        units3rd = ''
        pos3rd   = ''                               !If 3rd dimension does not refer to a height or level, use ''.
        tunits   = 'days since 0000-12-31 00:00:00'
        calendar = 'gregorian'                      !type of calendar
      
        n3rd     = dim3rd                           !size of 3rd dimension

        !initialize longitude and latitude
        delt = (lon_e - lon_s)/real(nlon)
        alons(1) = lon_s
        do i = 2, nlon
          alons(i) = alons(i-1) + delt
        end do !on i=1,nlon

        delt = (lat_e - lat_s)/real(nlat)
        alats(1) = lat_s
        do i = 2, nlat
          alats(i) = alats(i-1) + delt
        end do !on i=1,nlat

        !layer
        do i=1,dim3rd
          vals3rd(i)  = real(i)
        end do
      
      !--Generate the *.nc file following procedure from the end of io.f 

      !initialize file
      call inifile(idies,filen,title,source,history,
     >              nlon,alons,nlat,alats,name3rd,long3rd,units3rd,n3rd,vals3rd,pos3rd,tunits,calendar,ierror)

      if(ierror < 0 ) stop 'Error in test_NetCDF'

      !initialize variables
      ndims       = 4
      dimnames(1) = 'longitude'
      dimnames(2) = 'latitude'
      dimnames(3) = 'level'
      dimnames(4) = 'time'
      valmissing  = 0.0
  
      call inivar(idies,varname,longname,units,ndims,dimnames,valmissing,ierror)
  
      !call endini after initializing ALL the variables (here only 1 variable is used)
      call endini(idies,ierror)

      !write dummy data into the file
      istart  = (/1,    1,    1      , 1 /)  !initialize, 3rd dim will change later
      icount  = (/nlon, nlat, dim3rd , 1 /)
      cdate   = 'ANN1980  '//char(0)
      ftime   = 1.0
      tweight = 365.0

      !!dummy data
      !do k=1,dim3rd
      !  new_mydata(:,:,k) = real(k)
      !  call vec2arr(new_mydata(:,:,k),cdummy((k-1)*nlon*nlat+1) )
      !end do

      !call writevar(filen,varname,istart,icount,cdummy,ftime,tweight,cdate,ierror)

      !if(ierror.ne.0) then
      !  stop '(test_NetCDF)Error in writing NetCDF file!'
      !else
      !  write(6,*)'** NetCDF file successfully written as**'
      !  write(6,*) trim(fpath)//trim(fname)
      !end if
       

      return
      end subroutine write_netCDF
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! purpose : Coupled Agro-IBIS and THMB through file communication
!
! note    : Read-in input from THMB
!
!         - A check will be conducted once for grid consistenty between
!           THMB and Agro-IBIS
!         - <varName> does NOT need to append $ as done in THMB
!
!                Last updated 2015-04-20 by Yuwei Li
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!

      subroutine agro_thmb_read_input(ndaypm,iyear,nyear,imonth,iday,
     >                                varArray,varName,fname)

      use compar, only: thmb_lakem,agro_thmb_path,stamp_eps,
     >                  checked_grid,grid_torlerance,nlonsub,nlatsub,
     >                  thmb_FnameID,thmb_Fname,commu_sleep,commu_MaxWait,first_day
      use comwork,only: lonscale,latscale,work,cdummy
      use comgrid,only: npoi

      implicit NONE

      !--Input
      integer:: iyear,nyear,imonth,iday
      integer:: ndaypm(12)
      real   :: varArray(npoi)       !variable data
      character*(*):: varName,fname  !variable and file to read

      !--local
      integer:: ilon,ilat
      integer:: dimlon,dimlat,dim3rd,ierr,iwait
      integer,dimension(4):: istart,icount
 
      real,dimension(3):: stamp, !time stamp for communication
     >                    key,   !key retrevied from communication file
     >                    check_file
      real:: grid_long(nlonsub),grid_latg(nlatsub)

      character*4   year_char
      character*200 fpath,filen
      character*80  name3rd
      logical::     file_ready, file_exist

      !--Initialization
  
      varArray = 1e20

      !special treament for stamp here, as to appropriate couple Agro-IBIS and THMB
      
      !default
      stamp(1)= real(iyear) 
      stamp(2)= real(imonth) 
      stamp(3)= real(iday - 1)  !input is previous day of THMB's results
      !special cases

      !case 1: cross a year
      if(imonth == 1 .and. iday == 1) then
        stamp(1) = real(iyear - 1)
        stamp(2) = real(12)
        stamp(3) = real(ndaypm(12))
      !case 2: cross a month
      elseif(iday == 1) then
        stamp(2) = real(imonth - 1)
        stamp(3) = real(ndaypm(imonth-1))
      end if

      check_file = -999.0

      dimlon  = nlonsub
      dimlat  = nlatsub

      name3rd = 'level'
      dim3rd  = 3

      file_ready = .false.
      iwait      = 0

      filen = trim(agro_thmb_path)//trim(fname)

      istart  = (/1,           1,      1, 1 /)
      icount  = (/dimlon, dimlat, dim3rd, 1 /)

      !--Begin the work
      !no read-in from THMB for the first day of the simulation,
      if( nyear==0 .and. imonth==1 .and. iday==1 ) then
        file_ready = .true.
        first_day  = 1
        write(6,*) 'First day of the simulation, no read-in from THMB'
      else
        first_day = 0
      end if

      !--read data
      do while(.not.(file_ready))

        !check file existance
        inquire(file=trim(agro_thmb_path)//trim(thmb_Fname),EXIST=file_exist)
        if(file_exist) then
          open(unit=thmb_FnameID,file=trim(agro_thmb_path)//trim(thmb_Fname),status='old',action='read',iostat=ierr)
          read(thmb_FnameID,*,end=998,err=998) check_file(1),check_file(2),check_file(3)
          close(thmb_FnameID)
        else
          write(6,*) 'Error!Communication file <thmb.communicate> not exist. Check the path or name'
          write(6,*) 'User specified path: ', trim(agro_thmb_path)
        end if

        !check stamp
        if(abs(stamp(1)-check_file(1))>stamp_eps .or.
     >     abs(stamp(2)-check_file(2))>stamp_eps .or.
     >     abs(stamp(3)-check_file(3))>stamp_eps ) then

           iwait = iwait + 1
           write(6,*) 'Read-in from Agro-IBIS not ready. Waiting...#', iwait
           call sleep(commu_sleep)

           if(iwait > commu_MaxWait) then
             write(6,*) 'Error!Exceed maximum waiting. MaxWait = ', commu_MaxWait
             write(6,*) 'Stopping the program...'
             stop
           end if
        else !file ready to read

           file_ready = .true.
           write(6,'(A,3(2x,f12.5))') 'File is ready to read from Agro-IBIS for time stamp ', stamp

           !variable <cdummy> used to store data, <work> used as junk
           cdummy = 1e20
           call readvar(filen,varname,name3rd,istart,icount,cdummy,grid_long,grid_latg,key,work(1),ierr)
           if(ierr < 0) write(6,*) 'Error in reading NetCDF file: ', trim(filen)

           !check grid consistency
           if(.not.(checked_grid)) then  !check grid
             write(6,'(A,f12.5)') 'checking grid consistency for Agro-IBIS and THMB coupling, with grid tolerance: ',
     >                             grid_torlerance

             do ilon = 1, nlonsub

               if( abs(lonscale(ilon) - grid_long(ilon)) > grid_torlerance ) then
                 write(6,*) 'Error! Grid not consistent for Agro-IBIS and THMB at ilon: ',ilon
                 write(6,*) '(lontitude) for Agro-IBIS: ',lonscale(ilon)
                 write(6,*) '(lontitude) for THMB: ',grid_long(ilon)
                 write(6,*) 'Stopping the program...'
                 stop 'Error in grid consistency between Agro-IBIS and THMB'
               end if

             end do

             do ilat = 1, nlatsub
               if( abs(latscale(ilat) - grid_latg(ilat)) > grid_torlerance) then
                 write(6,*) 'Error! Grid not consistent for Agro-IBIS and THMB at ilon: ',ilon
                 write(6,*) '(latitude) for Agro-IBIS: ',latscale(ilat)
                 write(6,*) '(latitude) for THMB: ',grid_latg(ilat)
                 write(6,*) 'Stopping the program...'
                 stop 'Error in grid consistency between Agro-IBIS and THMB'
               end if
             end do

             checked_grid = .true.
             write(6,*) 'Checked. Grid is consistent between Agro-IBIS and THMB'

           end if !on if(.not.(checked_grid))

           !check stamp
           if(abs(stamp(1)-key(1))>stamp_eps .or.
     >        abs(stamp(2)-key(2))>stamp_eps .or.
     >        abs(stamp(3)-key(3))>stamp_eps ) then

                write(6,*) 'Error in checking stamp for file: ',trim(filen)
                write(6,*) 'The Key retrieved is: ', key
                write(6,*) 'Current stamp: ', stamp

             else
               write(6,*) 'Sucessfully read-in file ',trim(filen)
               call arr2vec(cdummy,varArray)
           end if
           
 
        end if 

  998 continue

      end do !on do while

      !debug info
      !write(6,'(A,A,A,f15.8)') 'Maximum read-in value for variable ', trim(varName),' :', maxval(varArray)
  
   
      return
      end subroutine agro_thmb_read_input
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! purpose : Coupled Agro-IBIS and THMB through file communication
!
! note    : write output for THMB
!           aname:no appeneded $ for variable (unlike THMB)
! 
!                Last updated 2015-04-08 by Yuwei Li
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!

      subroutine agro_thmb_write_output(iyear,imonth,iday,
     >                                  varArray,varName,fname,flg_commu)

      use compar,  only: agro_thmb_path,nlonsub,nlatsub,agro_FnameID,agro_Fname
      use comwork, only: lonscale,latscale,cdummy
      use comgrid, only: npoi

      implicit NONE

      !--Input
      integer:: iyear,imonth,iday
      integer:: flg_commu  !flag for communication 1:write stamp to file <agro.communicate>; 0: don't write
                           !signal this flag when the last *.nc was written (when several *.nc files are needed)
      real,dimension(npoi):: varArray
      character*(*):: varName,fname

      !--Local
      integer,parameter:: dim3rd = 3   !for stamp (year, month, day)
      real,dimension(3):: stamp
    
      !variables for initilizations (file, var)
      integer:: idies,ierror,i,j,k
      integer:: n3rd,ndims
      
      real   :: valmissing
    
      character*200 fpath
      character*200 filen,title,source,history,name3rd,long3rd,units3rd,pos3rd,tunits,calendar

      character*80  longname,units
      character(len=20),dimension(4)::dimnames
    
      !variables for writevar
      integer,dimension(4)  :: istart, icount
      real                  :: ftime(1),tweight(1)  !dummy (1)
      character*10 cdate

      !--Initialization

      stamp(1) = real(iyear)
      stamp(2) = real(imonth)
      stamp(3) = real(iday)
      
      filen    = trim(agro_thmb_path)//trim(fname)    !name for new file
      title    = 'Agro-IBIS_Input_for_THMB'           !title for file
      source   = 'source'                             !source for file
      history  = 'history'                            !date of creation, and possibly author
      name3rd  = 'level'                              !name of 3rd dimension
      long3rd  = 'time stamp'                         !long name for 3rd dimension variable
      units3rd = ''
      pos3rd   = ''                                   !If 3rd dimension does not refer to a height or level, use ''.
      tunits   = 'days since 0000-12-31 00:00:00'
      calendar = 'gregorian'                          !type of calendar
                  
      n3rd     = dim3rd                               !size of 3rd dimension
      
      !--Generate the *.nc file following procedure from the end of io.f 

      !initialize file
      call inifile(idies,filen,title,source,history,
     >             nlonsub,lonscale,nlatsub,latscale,name3rd,long3rd,units3rd,n3rd,stamp,pos3rd,tunits,calendar,ierror)

      if(ierror < 0 ) stop 'Error in test_NetCDF'

      !initialize variables
      longname    = 'Input_data_For_THMB'
      units       = ''
      ndims       = 4
      dimnames(1) = 'longitude'
      dimnames(2) = 'latitude'
      dimnames(3) = 'level'
      dimnames(4) = 'time'
      valmissing  = 0.0
  
      call inivar(idies,varName,longname,units,ndims,dimnames,valmissing,ierror)
      call endini(idies,ierror) !call endini after initializing ALL the variables (here only 1 variable is used)

      istart  = (/1,       1,       1,    1 /)
      icount  = (/nlonsub, nlatsub, 1,    1 /)
      cdate   = 'ANN1980  '//char(0)
      ftime   = 1.0
      tweight = 365.0

      cdummy = 0.0  

      !call vec2arr(varArray,cdummy)
      call vec2arr_ver2(varArray,cdummy)
      call writevar(filen,varName,istart,icount,cdummy,ftime,tweight,cdate,ierror)

      if(ierror.ne.0) then
        stop '(agro_thmb_write_output)Error in writing NetCDF file!'
      else
        write(6,*) 'Finished writting nc file for ', trim(varName)
      end if


      if(flg_commu == 1) then
        open(unit=agro_FnameID,file=trim(agro_thmb_path)//trim(agro_Fname),status='replace',action='write',iostat=ierror)
        write(agro_FnameID,*) stamp(1), stamp(2), stamp(3)
        close(agro_FnameID)

        write(6,'(A,3(2x,f12.4))') 'Output files for Agro-IBIS completed for: ',stamp
      end if !on if(flg_commu)
   
      return
      end subroutine agro_thmb_write_output
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! purpose : write daily input for THMB
!
! note    : variables list
!
!           adsrunoff   : daily average surface runoff 
!           addrainage  : daily average drainage
!           adrain      : daily average rain
!           adevap      : daily average evaporation
!           adtotnleach : daily average nitrogen leach
!           
! 
!                Last updated 2015-12-07 by Yuwei Li
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!

      subroutine write_dailyOutput4THMB(jday,iyear,imonth,iday)

      use comwork, only: lonscale,latscale,cdummy,OCEAN
      use comgrid, only: npoi
      use comsum

      implicit NONE

      !--Input
      integer:: jday,iyear,imonth,iday

      !--Local
      integer,parameter:: dim3rd = 1   !just 1 for 3rd dim
    
      !variables for initilizations (file, var)
      integer:: idies,ierror,i,j,k,ilayer
      integer:: n3rd,ndims
      real(kind=4):: var3rd(1)
      
      real   :: valmissing
    
      character*200 varName,fname
      character*200 filen,title,source,history,name3rd,long3rd,units3rd,pos3rd,tunits,calendar
      character*4   num,num2
      character*2 num_lay         !char for soil layer number

      character*80  longname,units
      character(len=20),dimension(4)::dimnames
    
      !variables for writevar
      integer,dimension(4)  :: istart, icount
      real                  :: ftime(1),tweight(1)  !dummy (1)
      character*10 cdate

      !--Initialization
      write(num,'(i4.4)')  iyear
      write(num2,'(i4.4)') iyear-1
      
      filen    = trim(myprocDir)//'output/daily/'//num//'_daily4thmb.nc'
      title    = 'Agro-IBIS daily input for THMB'           !title for file
      source   = 'source'                             !source for file
      history  = 'history'                            !date of creation, and possibly author
      name3rd  = 'level'                              !name of 3rd dimension
      long3rd  = 'level'                              !long name for 3rd dimension variable
      var3rd   = 1.0                                  !dummy value 
      units3rd = ''
      pos3rd   = ''                                   !If 3rd dimension does not refer to a height or level, use ''.
      tunits   = 'days since '//num2//'-12-31 00:00:00'
      calendar = 'gregorian'                          !type of calendar
                  
      n3rd     = dim3rd                               !size of 3rd dimension
      
      !--Generate the *.nc file following procedure from the end of io.f 

      !initialize file
      if(imonth == 1 .and. iday == 1) then
        call inifile(idies,trim(filen),title,source,history, 
     >                nlonsub,lonscale,nlatsub,latscale,name3rd,long3rd,units3rd,n3rd,var3rd,pos3rd,tunits,calendar,ierror)

        if(ierror < 0 ) write(6,*) '(write_dailyOutput4THMB) Error! in create nc file'

        !initialize variables
        ndims       = 4
        dimnames(1) = 'longitude'
        dimnames(2) = 'latitude'
        dimnames(3) = 'level'
        dimnames(4) = 'time'
        valmissing  = 1e20

        !**adsrunoff
        varName     = 'adsrunoff'
        longname    = 'daily average surface runoff'
        units       = 'mm/day'
        call inivar(idies,trim(varName),longname,units,ndims,dimnames,valmissing,ierror)

        !**addrainage
        varname     = 'addrainage'
        longname    = 'daily average drainage'
        units       = 'mm/day'
        call inivar(idies,trim(varName),longname,units,ndims,dimnames,valmissing,ierror)
        
        do ilayer = 1, nsoilay

          write(num_lay,'(i2)') ilayer
          varname     = 'addrainage_layer'//adjustl(num_lay)
          longname    = 'daily average drainage for layer '//adjustl(num_lay)
          call inivar(idies,trim(varName),longname,units,ndims,dimnames,valmissing,ierror)

        end do !on ilayer

        !**adrain
        varName     = 'adrain'
        longname    = 'daily average rain'
        units       = 'mm/day'
        call inivar(idies,trim(varName),longname,units,ndims,dimnames,valmissing,ierror)

        !**adevap
        varName     = 'adevap'
        longname    = 'daily average evaporation'
        units       = 'mm/day'
        call inivar(idies,trim(varName),longname,units,ndims,dimnames,valmissing,ierror)

        !**adtotnleach
        varName     = 'adtotnleach'
        longname    = 'daily average nitrogen leaching'
        units       = 'kg/ha/day'
        call inivar(idies,trim(varName),longname,units,ndims,dimnames,valmissing,ierror)

        do ilayer = 1, nsoilay

          write(num_lay,'(i2)') ilayer
          varname     = 'adtotnleach_layer'//adjustl(num_lay)
          longname    = 'daily average nitrogen leaching for layer '//adjustl(num_lay)
          call inivar(idies,trim(varName),longname,units,ndims,dimnames,valmissing,ierror)

        end do !on ilayer

        call endini(idies,ierror) !call endini after initializing ALL the variables (here only 1 variable is used)
      end if !on if(imonth == 1 && iday == 1)


      !**write data
      istart  = (/1,       1,       1,    jday /)
      icount  = (/nlonsub, nlatsub, 1,       1 /)
      cdate   = 'ANN1980  '//char(0)
      ftime   = real(jday)
      tweight = 365.0


      !**adsrunoff
      varName = 'adsrunoff'
      cdummy  = 1e20  
      !call vec2arr(adsrunoff,cdummy)
      call vec2arr_ver2(adsrunoff,cdummy)
      call writevar(trim(filen),trim(varName),istart,icount,cdummy,ftime,tweight,cdate,ierror)

      !**addrainage
      varName = 'addrainage'
      cdummy  = 1e20  
      !call vec2arr(addrainage,cdummy)
      call vec2arr_ver2(addrainage,cdummy)
      call writevar(trim(filen),trim(varName),istart,icount,cdummy,ftime,tweight,cdate,ierror)

      do ilayer = 1, nsoilay

        write(num_lay,'(i2)') ilayer
        varname = 'addrainage_layer'//adjustl(num_lay)
        cdummy  = 1e20  
        call vec2arr_ver2(addrainage_layer(:,ilayer),cdummy)
        call writevar(trim(filen),trim(varName),istart,icount,cdummy,ftime,tweight,cdate,ierror)

      end do !on ilayer

      !**adrain
      varName = 'adrain'
      cdummy  = 1e20  
      !call vec2arr(adrain,cdummy)
      call vec2arr_ver2(adrain,cdummy)
      call writevar(trim(filen),trim(varName),istart,icount,cdummy,ftime,tweight,cdate,ierror)

      !**adevap
      varName = 'adevap'
      cdummy  = 1e20  
      !call vec2arr(adevap,cdummy)
      call vec2arr_ver2(adevap,cdummy)
      call writevar(trim(filen),trim(varName),istart,icount,cdummy,ftime,tweight,cdate,ierror)

      !**adtotnleach
      varName = 'adtotnleach'
      cdummy  = 1e20  
      !call vec2arr(adtotnleach,cdummy)
      call vec2arr_ver2(adtotnleach,cdummy)
      call writevar(trim(filen),trim(varName),istart,icount,cdummy,ftime,tweight,cdate,ierror)

      do ilayer = 1, nsoilay

        write(num_lay,'(i2)') ilayer
        varname = 'adtotnleach_layer'//adjustl(num_lay)
        cdummy  = 1e20  
        call vec2arr_ver2(adtotnleach_layer(:,ilayer),cdummy)
        call writevar(trim(filen),trim(varName),istart,icount,cdummy,ftime,tweight,cdate,ierror)

      end do !on ilayer
   
      return
      end subroutine write_dailyOutput4THMB

!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! purpose : Resolve NaN or unreasonable restart values
!
! note    :
!      varName  : for echo-use only to see which variable is corrected
!      vect     : 1D array for the read-in restart data
!      baseVal  : base value w.r.t the variable to compare
!      largeVal : if abs(difference) between the read-in data and <baseVal>
!                 larger than <largeVal>, value from previous point will
!                 be assigned to the current point
!           
! 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      subroutine restart_correction(varName,vect,baseVal,largeVal)
      !Correction for read-in restart variables
      ! if abs(vect) > largeNum, assign value from previous point to it 
        use comgrid, only: npoi
        use combcs,  only: cmask

        implicit NONE
        real:: vect(npoi),baseVal,largeVal
        character*(*):: varName
        integer i

        do i=1,npoi
          if(cmask(i)==0) cycle
          if( abs(vect(i)-baseVal)>largeVal .or. isnan(vect(i)) ) then
              write(6,*) 'Warning at <restart>, variable,i,val: ',trim(varName),i,vect(i)
              write(6,*) 'Use restart value from previous point#, val: ', i-1, vect(i-1)
              vect(i) = vect(i-1)
          end if
        end do !on i=1,npoi


      end subroutine restart_correction
c
c
c -------------------------------------------------------------------------
c HOW TO READ/WRITE NETCDF FILES IN IBIS
c -------------------------------------------------------------------------
c Reading/writing files in ibis is done through subroutines in
c ies-io.f.  The most important concept to understand is the
c relationship between the locations of points in an n-dimensional
c array and the values of istart and icount.  In FORTRAN, values in an
c array are stored so that the first dimension varies the fastest.  In
c C, the last dimension varies the fastest.  If you were to use ncdump
c on a netcdf file whose variable 'mydata' has 4 dimensions (latitude,
c longitude, level, time), the variable would be shown as follows:
c mydata(time, level, latitude, longitude)
c because ncdump was written in C and reflects C's ordering of dimensions.
c In this example time varies the the slowest and longitude the fastest.
c If you were to define this variable in a FORTRAN program so that
c time again varies the slowest and longitude the fastest, it would
c look like this:
c     real mydata(nlons, nlats, nlevels, ntimes)
c where nlons, nlats, nlevels, and ntimes are integer parameters of
c some specified size.  Looping through the data in the order of
c storage in memory would look like this:
c     do l = 1, ntimes
c       do k = 1, nlevels
c         do j = 1, nlats
c           do i = 1, nlons
c             mydata(i,j,k,l) = ...
c           enddo
c         enddo
c       enddo
c     enddo
c Since ies-io.f is FORTRAN code, keep in mind the FORTRAN
c representation will be flipped in comparison to what you see from
c ncdump.
c The netcdf interface reads and writes data using two integer vectors,
c istart, and icount.  Istart refers to the starting point for reading/writing
c along each dimension.  If you have a 4-d variable as above and want
c to write starting at the first point, istart would be a vector of
c length 4 and have the values (1,1,1,1). Icount refers to how far
c along each dimension data will be read/written.  If you wish to
c write the entire array in the sample FORTRAN code above, icount
c would be a vector of length 4 and have the values
c (nlons,nlats,nlevels,ntimes).
c Things get a little more complicated when you want to read/write
c only a portion of the variable.  In ibis we typically read/write a
c single time step, but possibly more than one level/pft and either an
c entire or partial lat/lon grid. Below are examples of istart and
c icount for various situations of reading/writing.
c 1) an entire lat/lon grid for only one (6th) time step and one (2nd) pft:
c istart is (1,1,2,6), icount is (nlons,nlats,1,1)
c
c 2) entire lat/lon arrays for all 9 pfts at one (6th) time step:
c istart is (1,1,1,6), icount is (nlons,nlats,npfts,1)
c
c 3) a single lat/lon point (at the ilon-th longitude and the ilat-th
c latitude) of a 3-d variable at all 100 time steps:
c istart is (ilon,ilat,1), icount is (1,1,100)
c Note that if istart and icount have been declared length 4, the 4th
c value in each is ignored when referencing a 3-d variable.
c
c 4) a subsection of a lat/lon grid, 20 points wide in longitude, 15
c points high in latitude, starting at point (ilon,ilat), at one (18th)
c level and 12 times, beginning at time step itime:
c istart is (ilon,ilat,18,itime) icount is (20,15,1,12)
c
c HOW TO ADD NEW CODE TO READ A FILE:
c To read a file, use subroutine readvar in ies-io.f. This subroutine
c assumes that the variable being read has dimensions called longitude,
c latitude, possibly time, and possibly another dimension, which is
c named in the call.  Only the bare essentials are returned.
c
c General call:
c     call readvar(filen,varname,name3d,istart,icount,values,
c    > alons,alats,vals3d,times,ierror)
c
c INPUT
c     filen - character*(*) - file name from which to read data
c     varname - character*(*) - name of variable from which to read
c     name3d - character*(*) - name of 3rd, nontime dimension variable.
c      Ignored if varname variable is only 3-d.
c     istart - integer(*) - starting points along each dimension of
c      the variable, for example the 1st point a 4-d variable is (1,1,1,1) 
c      and the 3rd level and 2nd time step of a 4d variable is (1,1,3,2).
c     icount - integer(*) - number of points along each dimension to read,
c      for example, to read in a single lat/lon grid from a 4-d variable,
c      icount would be (nlon,nlat,1,1)
c OUTPUT
c     values - real(*) - returned real values of the designated hyperslab
c     alons - real(*) - longitude vector for those points read in
c     alats - real(*) - latitude vector for those points read in
c     vals3d - real(*) - vector of values of the 3rd dimension (unchanged if
c      variable is 2- or 3-d, or if 3rd dimension has character values).
c     times - real(*) - time vector vector for those points read in (unchanged
c      if variable is 2-d).
c     ierror - integer - error code = 0 if no error, < 0 if error occurred
c
c Example: Read in an entire lat/lon grid for one (9th) pft (dimension 
c name is 'pft') and one (24th) time
c     parameter (nlons = 360, nlats=180)
c     real x(nlons,nlats), alons(nlons), alats(nlats), xjunk, time
c     integer istart(4), icount(4), ierr
c     istart(1) = 1
c     istart(2) = 1
c     istart(3) = 9
c     istart(4) = 24
c     icount(1) = nlons
c     icount(2) = nlats
c     icount(3) = 1
c     icount(4) = 1
c     call readvar('myfile.nc','myvar','pft',istart,icount,x,alons,alats,
c    > xjunk,time,ierr)
c     if (ierr .lt. 0) then
c       print *, 'Error occurred in readvar'
c     end if
c Note that above, I used a junk variable for the values of the 3rd
c dimension, because in ibis, the pft dimension is a character
c dimension.  If the 3rd dimension type is character, readvar does not
c return the value.
c
c Ibis-specific example: read in a variable which has 3rd dimension
c 'layer', using cdummy to store full array before extracting land
c points into the variable used directly in ibis.  The variable work is
c used to store returned info that won't be used later.
c Previously declared
c      istart(1) = 1        !begin at 1st longitude point
c      istart(2) = 1        !begin at 1st latitude point
c      istart(3) = 1        !begin at 1st layer
c      istart(4) = 1        !begin at 1st time
c      icount(1) = nlon     !read nlon points along longitude dimension
c      icount(2) = nlat     !read nlat points along latitude dimension
c Read in data in file input/soita.nc to ibis variable tex
c      icount(3) = nsoilay  !read nsoilay points along layer dimension
c      icount(4) = 1        !read one time step
c      filen = 'input/soita.nc'
c      aname = 'tex'
c      call readvar(filen,aname,'layer',istart,icount,
c     > cdummy,work(1),work(ndim4),work(2*ndim4),work(3*ndim4),istat)
c      if (istat.lt.0) then
c         write(*,9000)
c         print *, 'while reading soita'
c         stop 1
c      end if
c For each lat/lon grid, strip off ocean points. tex only has land points.
c      do 12 j = 1, nsoilay
c        call arr2vec (cdummy((j-1)*nlonsub*nlatsub + 1), tex(1,j))
c 12   continue
c
c HOW TO ADD NEW CODE TO WRITE TO A FILE:
c Writing to a file involves several steps.  First, if the file does
c not yet exist, you must create the file using inifile or inifilec.
c Then, you must create one or more variables using inivar.  Then after
c all initializing is done, declare the end of the initialization phase
c using endini.
c Finally, once the file and variable exist, you can write data into
c the variable using writevar.
c Creating and writing to files in seperate steps may seem
c complicated at first, but this 4-step method allows much flexibility.
c Not only do you have the choice of creating a variable without a 3rd
c dimension, with a real 3rd dimension, or with a character 3rd
c dimension, you may also have more than one variable within the same
c file.
c Note that when using ies-io.f routines if you have more
c than one variable in the file, and both have 4 dimensions, that 3rd
c dimension (level, pft, etc.) must be shared by both variables.  For
c example, you may have one variable whose 3rd dimension is pft and
c another variable which does not have a 3rd non-time dimension in the
c same file, but you may not have in the same file one variable whose
c 3rd dimension is pft and another variable whose 3rd dimension is level.
c
c 1) Initialize a file
c If you wish a file to contain only latitude, longitude, and time
c dimensions, or want it to also have a 3rd real dimension, use
c inifile.  If you want to the file to have a 3rd character dimension,
c use inifilec.
c General call for inifile:
c      call inifile(idies,filen,title,source,history,nlon,alons,
c     > nlat,alats,name3rd,long3rd,units3rd,n3rd,vals3rd,pos3rd,tunits,
c     > calendar,ierror)
c
c INPUT
c     filen - character*(*) - name for new file
c     title - character*(*) - title for file
c     source - character*(*) - source for file
c     history - character*(*) - date of creation, and possibly author
c     nlon - integer - number of point along longitude direction
c     alons - real(nlon) - values of longitudes
c     nlat - integer - number of point along latitude direction
c     alats - real(nlat) - values of latitudes
c     name3rd - character*(*) - name of 3rd dimension - use '' if you
c      don't want this dimension ('' is the empty string)
c     long3rd - character*(*) - long name for 3rd dimension variable
c      (ignored if nam3rd='')
c     n3rd - integer - size of 3rd dimension (ignored if name3rd='')
c     vals3rd - real(n3rd) - values along 3rd dimension (ignored if 
c      name3rd='')
c     pos3rd - character*(*) - if the 3rd dimension exists and refers to
c      a measured level or height, use 'up' if values increase with distance
c      above earth (as in height), or 'down' if values decrease with distance 
c      (as in pressure) (ignored if name3rd=''). If 3rd dimension does not
c      refer to a height or level, use ''.
c     tunits - character*(*) - units for time, must be in the form of days
c      since yyyy-mm-dd tt:tt:tt, where yyyy is a year or 0000, mm is a
c      month, dd is a day, and tt:tt:tt is (optional) time or 00:00:00.
c     calendar - character*(*) - type of calendar.  Choose from 'noleap',
c      'gregorian','n ka BP', etc.  Use iescal (in ies.f) if orbital 
c      parameters differ from modern.
c OUTPUT
c     idies - integer - id number of new file for later use
c     ierror - integer - error code, 0 = no error, < 0 = an error occured
c
c Example: Create a file which will hold fractional snow
c cover (only 3-d) and snow thickness for each layer (4-d).
c previously defined: 
c
c     parameter (nlons = 360, nlats = 180, nsnolay = 6)
c     real lonscale(nlons), latscale(nlats), slayers(nsnolay)
c     real snowc(nlons,nlats), snowh(nlons,nlats,nsnolay)
c     integer istat
c     character*80 cdate, tunits
c     alats = ...   ! create values for latitude somehow
c     alons = ...   ! create values for longitude somehow
c     slayers = ... ! create values for snow layers somehow
c     cdate = 'created on 6/29/97'
c     tunits = 'days since 1969-12-31'
c Now initialize a file with a 3rd real variable
c     filen = 'snow.nc'
c     call inifile(idies,filen,
c    > 'file for snow variables',
c    > 'C Molling, program.f v1.01',cdate,nlons,lonscale,nlats,latscale,
c    > 'snowlayer','snow layers top to bottom','',nsnolay,slayers,
c    > 'down',tunits,'gregorian',istat)
c     if (istat .lt. 0)then
c       print *, 'Error in inifile'
c     end if
c Note above, the empty string ('') is used because the 3rd dimension,
c snowlayer, does not have any units.  The 'positive' attribute for
c the 3rd dimension is 'down' because snowlayer 1 is at the top and the last
c layer is at the bottom (distance above the center of the earth is
c decreasing as snowlayer increases).  The calendar is gregorian
c because there are leap years included.  If all years are 365 days, you
c should use 'noleap'.  Units of time are days dince a date of the
c form yyyy-mm-dd.  Time units should use this form to be compatible with
c GrADS.  Other units should be compatible with Udunits (see
c www.unidata.ucar.edu).  The returned variable idies will be used in
c the subroutine that initializes a variable.
c
c General call for inifilec
c     call inifilec(idies,filen,title,source,history,nlon,alons,
c    > nlat,alats,name3rd,long3rd,units3rd,n3rd,len3rd,chars3rd,
c    > tunits,calendar,ierror)
c
c INPUT
c     filen - character*(*) - name for new file
c     title - character*(*) - title for file
c     source - character*(*) - source for file
c     history - character*(*) - date of creation, and possibly author
c     nlon - integer - number of point along longitude direction
c     alons - real(nlon) - values of longitudes
c     nlat - integer - number of point along latitude direction
c     alats - real(nlat) - values of latitudes
c     name3rd - character*(*) - name of 3rd dimension
c     long3rd - charqcter*(*) - long name for 3rd dimension variable
c     n3rd - integer - size of 3rd dimension
c     len3rd - integer length of chracter strings in vals3rd
c     chars3rd - character*len3rd(n3rd) - values along 3rd dimension
c     tunits - character*(*) - units for time, must be in the form of days
c      since yyyy-mm-dd tt:tt:tt, where yyyy is a year or 0000, mm is a
c      month, dd is a day, and tt:tt:tt is (optional) time or 00:00:00
c     calendar - charater*(*) - type of calendar.  Choose from 'noleap',
c      'gregorian','n ka BP', etc.  Use iescal if orbital parameters
c      differ from modern.
c OUTPUT
c     idies - integer - id number of new file for later use
c     ierror - integer - error code, 0 = no error, < 0 = an error occured
c
c Example: Create a file that has a character 3rd dimension.
c Defined previously - most variables as in previous example, plus...
c      parameter (npft = 9)
c      character*80 pftdef(npft)
c      pftdef(1) = 'boreal evergreens'
c      pftdef(2) = ...
c      etc...
c Initialize file with inifilec
c      filen = 'exist.nc'
c      call inifilec(idies,filen,
c     > 'annual existence of each plant functional type',
c     > 'ibis wyearly',cdate,nlon,lonscale,nlat,latscale,
c     > 'pft','plant fuctional type','',npft,80,pftdef,
c     > tunits,'gregorian',istat)
c The '80' above refers to the length of the character strings in pftdef.
c         dimnames(3) = 'pft'
c         call inivar(idies,'exist','existence for each pft',
c     >    '',4,dimnames,OCEAN,istat)
c End initialization phase
c         call endini(idies,istat)
c
c 2) Initialize a variable and end initialization phase
c After you initialize a file, you need to initialize a variable.
c Initializing reserves space for data to be written into the file, and
c saves information about the data (such as units, descriptive name, and
c a missing value).  You can use inivar to initialize any variable of 1
c or more dimensions, provided that those dimensions already exist in
c the file.  You may initialize more than one variable in a single
c file, as stated above.
c If you have a special value that denotes 'missing', such as using
c 1.e+36 to denote ocean grid cells, use this value for valmissing.  If
c you use 0. for valmissing, it is ignored.  Pick a value for valmissing
c which is well outside the valid range for the data.  For example,
c -99. would be a fine missing value for pressure, but not a good value
c for topography.
c The character array dimnames is used to store the
c names of the dimensions upon which the new variable will depend.
c For example, a 3-d variable would only depend on longitude, latitude,
c and time.  A 4-d variable would depend on those and also pft, or
c snowlater, or level, or some other dimension.  The dimension names
c must be in the same order of varying: fastest to slowest.  See the
c discussions for istart and icount above.  For example, the 4-d
c variable to go in the file snow.nc (above inifile example)  would have
c dimnames as follows
c      dimnames(1) = 'longitude'
c      dimnames(2) = 'latitude'
c      dimnames(3) = 'snowlayer'
c      dimnames(4) = 'time'
c
c General call for inivar:
c      call inivar(idies,varname,longname,units,ndims,dimnames,
c     > valmissing,ierror)
c
c INPUT
c     idies - integer - id number of a new file from a previous call
c      to inifile, inifilec, or iesopen
c     varname - charcter*(*) - name for new variable
c     longname - character*(*) - descriptive name for variable
c     units - character*(*) - units for variable, compatible with udunits
c     ndims - integer - number of dimensions for this variable
c     dimnames - character*(*)(ndims) - name of each dimension, in order
c     valmissing - real - missing value, ignored if valmissing=0.
c OUTPUT
c     ierror - integer - error code, 0 = no error, < 0 = error occurred
c
c Example: Define a 4-d and a 3-d variable to go into the file snow.nc.
c defined previously
c        character*80 dimnames(4)
c        real OCEAN
c        dimnames(1) = 'longitude'
c        dimnames(2) = 'latitude'
c        OCEAN = 9.e+20
c define 4-d variable
c        dimnames(3) = 'snowlayer'
c        dimnames(4) = 'time'
c        call inivar(idies,'hsno','snow layer thickness',
c     >   'meters',4,dimnames,OCEAN,istat)
c define 3-d variable in same file
c        dimnames(3) = 'time'
c        call inivar(idies,'snowf','snow cover fraction',
c     >   'fraction',3,dimnames,OCEAN,istat)
c end initialization phase
c         call endini(idies,istat)
c Notice that you need to change the value of dimnames(3) for the 3-d
c variable.  The value in dimnames(4) gets ignored.
c
c 3) End the initialization phase for the file
c When you are done initializing the file and the variables in the file, you
c must call endini.  This subroutine merely closes the file.  But by closing
c the file, this is a signal to the computer to write out all changes to this
c file, synchronizing the instructions in the program with what is written on
c the hard disk.  Since most computers buffer, that is save up, data until a
c there is a large chunk to write, it is essential that all buffered
c information for a netcdf file be written to disk before attempting to write
c data to the file.  The subroutine endini accomplishes this.  The subroutine
c endini should be called after initializing the last variable in the netcdf
c file (inivar) and before calling writevar.  See the example above.
c
c 4) Write to a variable
c After you have completed the initialization phase by initializing the file
c and initializing one or more variables in the file, you can write data 
c into the space reserved for each variable with subroutine writevar.
c You need to supply istart and icount, just as when you read a
c variable.  It is perfectly legal to write values in only a portion of the
c space reserved for the variable.  If nothing is written to a
c particular grid point, it has a special fill value supplied by netcdf
c when the variable was initialized. Notice that for ease of use, I
c have designed writevar to also write in the values of the time steps
c you are writing as well as the weighting for that time step.
c Time weighting refers to the number of days in that time sample.
c For example, some monthly means will be a mean value over 31 days,
c while others will be over only 30, or 28, or 29 days.  I assumed that
c most persons will calculate data and write it out one time step at a
c time, so it is nice to write out the time values and weight as you
c go along.  Another time-related item written is the character label
c for the time step.  Since after a few years, it's hard to keep track
c of what month is represented by a large number of days, I invented the
c date variable.  Date is a 10-character long string in which you can
c put a time label.  For example the time value 4745 may not be
c informative, but the data label 'JAN1913   ' is.  I usually use 9
c characters plus a null character at the end (char(0)).  Use strings
c like 'DJF001005 ' to denote time averages (Winter years 1 through 5).
c
c General call
c     call writevar(filen,varname,istart,icount,values,
c    > times,tweights,dates,ierror)
c
c INPUT
c     filen - character*(*) - file name to which to write data
c     varname - character*(*) - name of variable to which to write
c     istart - integer(*) - starting points along each dimension of
c      the variable, for example the 1st point a 4-d variable is (1,1,1,1) 
c      and the 3rd level and 2nd time step of a 4d variable is (1,1,3,2).
c     icount - integer(*) - number of points along each dimension to write,
c      for example, to write a single lat/lon grid to a 4-d variable,
c      icount would be (nlon,nlat,1,1).
c     values - real(*) - real values of the designated hyperslab
c     times - real(*) - time vector vector for those points written (ignored
c      if variable is 2-d).
c     tweights - real(*) - number of days in each time step in times (ignored
c      if variable is 2-d).
c     dates - character*10(*) - character labels for each time step in times,
c      last character should be a null (dates ignored if variable 2-d).
c OUTPUT
c     ierror - integer - error code = 0 if no error, < 0 if error occurred
c
c Example: write data to the 4-d variable initialized above
c previously defined
c     parameter (ndim3=nlons*nlats*nsnolay)
c     real cdummy(ndim3), hsno(npts,nsnolay)
c     real time, timewght
c     integer istart(4), icount(4)
c     character*10 cdate
c     time = 730.
c     timewght = 365.
c     istart(1) = 1
c     istart(2) = 1
c     istart(3) = 1
c     istart(4) = 2
c     icount(1) = nlons
c     icount(2) = nlats
c     icount(3) = nsnolay
c     icount(4) = 1
c     cdate = 'ANN1980  '//char(0)
c put land-only data in hsno into land and sea lat/lon grid in workspace
c     do 20 k = 1, nsnolay
c        call vec2arr (hsno(1,k), cdummy((k-1)*nlonsub*nlatsub + 1))
c 20  continue
c     call writevar('snow.nc','hsno',istart,icount,cdummy,time,
c    > timewght,cdate,istat)
c     if (istat .ne. 0) then
c        write(*,*) 'ERROR writing hsno'
c        stop 1
c     end if
c
c Some other thoughts:
c It's a good idea to end all strings stored in netcdf file with a
c null character.  Some less-robust C programs may fail when
c encountering non-null-ended strings, because in C, all strings end
c with a null character.
c Udunits provides several ways to express multiplication, powers,
c etc.  I like to use ^ to denote exponents (as in m^2 instead of just
c m2) because programs like Matlab see ^ as a LaTeX command to write the
c next character as a superscript.  Likewise I try to avoid using _,
c because this is the LaTeX command to make the next character a
c subscript.  So instead of using deg_C, I use degC, to prevent the C
c from becoming a subscript.
c The format of the netcdf files written by ies-io.f subroutines is
c meant to be compatible with the COARDS and CSM conventions (see
c www.cgd.ucar.edu:80/csm/experiments/output/format.html).
c Theoretically, you should be able to use NCO (see
c www.cgd.ucar.edu:80/cms/nco/index.html) and NCL (see
c ngwww.ucar.edu/ngdoc/ng4.1alpha/ug/ncl/ncloview.html) on the files.
c A few other packages that work with this format are GrADS
c (grads.iges.org/grads/head.html) and Ncview
c (meteora.ucsd.edu/~pierce/ncview_home_page.html).
c The routines in ies-io.f are just the beginning.  You can put many
c more things in a netcdf file by using routines in ies.f or by using
c the low-level netcdf commands directly.  Use one of the pre-existing
c subroutines in ies-io.f as a template and go on from there.  You are
c free to change or distribute my code as long as you 1) acknowlege me,
c Christine C. Molling (cmolling@facstaff.wisc.edu) as the author of
c the original code and 2) do not sell it.
c Of course if there is anything wrong with the code, or you somehow
c encounter damage to your reputation/wallet because of its use, you
c are forbidden to sue me or anyone else.
c Good Luck!
