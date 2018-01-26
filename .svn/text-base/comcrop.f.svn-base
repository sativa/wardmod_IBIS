c
c ---------------------------------------------------------------------
      module comcrop 
c     last edited by Bill Sacks 05.06.10
c ---------------------------------------------------------------------
c
      use comgrid
      use compar, only : npft

      implicit none
      save

c First define some parameters
c
      integer max_fert_years ! maximum number of years of data in the fertilizer maps
      integer max_cfert_years ! maximum number of years of data in the text file used to read the cfert* variables
      parameter (max_fert_years = 100)
      parameter (max_cfert_years = 2000)
c
c
      logical
     >  management_prescribed  ! true: use prescribed planting dates & cultivars, as read from maps; false: use prognostic planting dates & cultivars (but as of 5-16-10, days to harvest always comes from maps)

      integer
     >  irrigate_flag,   ! 0: no irrigation 1: irrigation strategy used everywhere 2: read map specifying where irrigation is used (ibis.infile)
     >  nstress,         ! 0: no nitrogen stress 1: apply nitrogen stress
     >  isoybean,        ! 0: soybeans not planted 1: soybeans planted (ibis.infile)
     >  imaize,          ! 0: maize not planted 1: maize not planted (ibis.infile)
     >  iwheat,          ! 0: wheat not planted 1: spring wheat planted 2: winter wheat (ibis.infile)
     >  iwheat_index,    ! index into wheat arrays: if iwheat > 0 then iwheat_index = iwheat;
     >                   !                          else iwheat_index = 1 (in this case its value doesn't matter)
     >  irotation,       ! 0: no crops rotated (monoculture) 2: two crops rotated 3: 3 crops rotated (see crops.f)
     >  isoilay,         ! soil layer for which drainage and leaching data is output
     >  overveg          ! 

      integer
     >  irrigate(npoi)          ! do we use irrigation in the given grid cell? (0 everywhere if irrigate_flag is 0; 1 everywhere if
                                ! irrigate_flag is 1; spatially-varying if irrigate_flag is 2)
c
      real      
     >  ztopmxsoy,     ! height maximum (m) for soybean
     >  ztopmxwht,     ! height maximum (m) for wheat
     >  ztopmxmze,     ! height maximum (m) for maize
     >  alphac,        ! small amount of N allowed to diffuse to roots under low uptake conditions
     >  gnmin,         ! minimum nitrogen fraction allowed in grain
     >  smax,          ! maximum stress factor allowed for nitrogen on vmax
     >  availn,        ! amount (fraction) of inorganic N pool available to plants
     >  cnmax,         ! maximum c:n allowed of aboveground crop residue
     >  consdays(npoi),  ! counter keeping track of number of days between freezing events
     >  maxcons(npoi),   ! maximum number of consecutive days between freezing events
     >  iniday(npoi),    ! last day of freeze in spring
     >  endday(npoi),    ! first day of freeze in fall
     >  gsdays(npoi),    ! number of days in the growing season  
     >  apar(npoi,npft)       ! total annual absorbed PAR
     
      real
     >  xirrig(npoi),       ! irrigated water application rate (mm/day) to crops
     >  xirriga(npoi),      ! irrigated application rate per timestep
     >  totirrig(npoi),     ! annual total irrigation applied (mm/yr)
     >  icropsum(npoi),     ! index - number of crop types planted in each grid cell
     >  df(npoi),           ! photoperiod factor for wheat crops
     >  vf(npoi),           ! vernalization factor for wheat crops 
     >  cumvd(npoi),        ! cumulative vernalization days
     >  hdidx(npoi),        ! cold-hardening index for wheat (varies from 0-2: stage 1 : 0-1, stage 2: 1-2)
     >  gddfzcorn(npoi),    ! seasonal gdd (base 8 C) in between last freeze and first freeze
     >  gddfzsoy(npoi),     ! seasonal gdd (base 10 C) in between last freeze and first freeze
     >  gddfzwht(npoi),     ! seasonal gdd (base 0 C) in between last freeze and first freeze
     >  conspdate(npoi),    ! constant planting date based on specified time period
     >  conshybrid(npoi),   ! constant hybrid (gdd) based on specified time period
     >  gddcorn(npoi,2100), ! hybrid planted for corn 
     >  gddsoy(npoi,2100),  ! hybrid planted for soy 
     >  gddwht(npoi,2100)   ! hybrid planted for wheat 
c
      integer 
     >  idop(npoi,npft),       ! day of year that crop was planted 
     >  iavepdate(npoi,npft),  ! average planting date over time
     >  harvdate(npoi,npft),   ! day of year that crop pft was harvested
     >  corndop(npoi,2100),    ! day that corn is planted
     >  soydop(npoi,2100),     ! day that soy  is planted
     >  whtdop(npoi,2100)      ! day that wheat is planted
c
      real
     >  cntops(npoi,npft),     ! cn ratio of plant residue
     >  cnroot(npoi,npft),     ! cn ratio of plant roots 
     >  croplive(npoi,npft),    ! 0 crops have been planted and living : 1 crops not living (WJS 05.16.10: I think the meaning here is reversed) (See also croppresent)
     >  grnfraccrop(npoi,npft),! green fraction of aboveground vegetation for crops 
     >  gddplant(npoi,npft),   ! accumulated growing degree days past planting date for crop = j
     >  gddtsoi(npoi,npft),    ! accumulated growing degree days past planting date based on top layer soil temperature  
     >  gddmaturity(npoi,npft),! accumulated growing degrees needed for plant to reach both vegetative and physiological maturity
     >  thrlai(npoi,npft),     ! lai threshold for crops when senescence begins 
     >  peaklai(npoi,npft),    ! 0: lai not at maximum allowable, 1: lai at peak value allowed
     >  lai_decl_per_day(npoi,npft), ! daily LAI decline from maturity to harvest
     >  hui(npoi,npft),        ! heat unit index
     >  phuf(npoi,npft),       ! previous day's value of the heat unit factor
     >  tlai(npoi,npft), 
     >  templai(npoi,npft), 
     >  harvidx(npoi,npft),    ! end of year harvest index for crop
     >  fnleaf(npoi,npft),     ! current fraction of nitrogen in leaf dry matter
     >  fnstem(npoi,npft),     ! current fraction of nitrogen in stem dry matter
     >  fnroot(npoi,npft),     ! current fraction of nitrogen in root dry matter
     >  fngrain(npoi,npft),    ! current fraction of nitrogen in grain dry matter
     >  fnplant(npoi,npft),    ! current fraction of nitrogen in entire plant 
     >  tnplant(npoi,npft),    ! total nitrogen in plant dry matter
     >  grainn(npoi,npft),     ! total nitrogen offtake by crop in grain 
     >  cumlvs(npoi,npft),     ! total number of leaves emerged
     >  idpp(npoi,npft),       ! number of days past planting
     >  dpgf(npoi,npft),       ! number of days past grain fill
     >  dpmature(npoi,npft),   ! number of days past maturity
     >  cropyld(npoi,npft),    ! crop yield in bu/ac
     >  dmleaf(npoi,npft),     ! leaf dry matter in Mg/ha
     >  dmstem(npoi,npft),     ! stem dry matter in Mg/ha
     >  dmresidue(npoi,npft),  ! aboveground leaf and stem residue dry matter in Mg/ha
     >  dmyield(npoi,npft),    ! yield dry matter in Mg/ha
     >  dmcrop(npoi,npft),     ! total crop dry matter production in Mg/ha
     >  cropn(npoi, npft),     ! nitrogen removed by crop in kg/ha/yr
     >  cropfixn(npoi,npft),   ! nitrogen fixation by crop in kg/ha/yr 
     >  nconcl(npoi,npft),     ! end of season N concentration in leaf (fraction)
     >  nconcs(npoi,npft),     ! end of season N concentration in stem (fraction)
     >  nconcr(npoi,npft),     ! end of season N concentration in root (fraction)
     >  nconcg(npoi,npft),     ! end of season N concentration in grain (fraction)
     >  leafout(npoi,npft),    ! gdd accumuation index for leaf emergence in crops
     >  cropplant(npoi,npft),  ! index keeping track of whether that crop has been planted during current year 
     >  croplaimx(npoi,npft),  ! maximum attained lai by crop during growing season 
     >  residuen(npoi, npft),  ! nitrogen contained in aboveground crop residue 
     >  dmroot(npoi,npft),     ! fine root dry matter (Mg/ha)
     >  matdate(npoi,npft),    ! day of year that plant reaches maturity
     >  hdate(npoi,npft),      ! harvest date (real value)
     >  pdate(npoi,npft),      ! planting date (real value) 
     >  crmclim(npoi,npft),    ! crop relative maturity rating (CRM) from Pioneer regression relationships
     >  crmact(npoi,npft),     ! best crop relative maturity rating (CRM) for that year based on total GDD accumulated 
     >  crmplant(npoi,npft),   ! crop relative maturity rating (CRM) based on GDD accumulated since planting of that year 
     >  ccdays(npoi,npft),     ! number of consecutive days with minimum temperatures below killing threshold
     >  fkill_gdd_missed(npoi,npft), ! fraction of the season (in a GDD sense) that the crop "missed out on" due to freeze kill (e.g., if gddmaturity=1500 but freeze kill happened at hui=1200, then this variable will be 0.2; this will be 0 if there was no freeze kill)
     >  dev_fraction_attained(npoi,npft), ! fraction of development attained at the time the crop reaches maturity (hui/gddmaturity); this will usually be 1 or slightly greater, but can be less than 1 if the crop is forced to reach "maturity" prematurely due to freeze kill, the max # days limit, etc.
     >  emerge_day(npoi,npft), ! day of year that leaf emergence occurs
     >  emerge_gdd(npoi,npft), ! GDD accumulated since planting at the time of leaf emergence (using air temp.)
     >  emerge_gddtsoi(npoi,npft), ! soil GDD accumulated since planting at the time of leaf emergence
     >  grainday(npoi,npft),   ! day of year that plant reaches grain fill stage
     >  graingdd(npoi,npft),   ! GDD accumulated since planting at the time the plant reached grain fill stage
     >  matgdd(npoi,npft),     ! GDD accumulated since planting at the time the plant reached maturity
     >  matgrnfraccrop(npoi,npft), ! green fraction of crops (i.e., value of grnfraccrop) at the time the plant reached maturity
     >  fertinput(npoi,npft),  ! annual fertilizer input for crop (kg-N per m-2)
     >  avehybrid(npoi,npft),  ! average hybrid planted (GDD) 
     >  arepr(npoi,npft),      ! fraction allocation to reproductive organs (grain/fruit) in crops
     >  astem(npoi,npft),      ! fraction allocation to stem in crops (non leaf/ non grain)
     >  astemi(npoi,npft),     ! fraction allocation to stem in crops at initial shift to grain 
     >  aleafi(npoi,npft)      ! fraction allocation to leaf in crops at initial shift to grain 
c

c croppresent: 1 if this crop is present here, 0 if not. The difference
c between croppresent and croplive is that, between maturity and
c harvest, croplive is 0 (crop is no longer photosynthesizing), but
c croppresent is 1 (there is still a crop in the field)
      integer croppresent(npoi, npft)

      integer
     >  fert_year1,             ! first year of data in fertilizer maps (read into fertmaize, fertsoy, fertwheat)
     >  fert_nyears,            ! number of years of data in fertilizer maps
     >  fert_lastyear,          ! last year of data in fertilizer maps (calculated from fert_year1 and fert_nyears)
     >  cfert_year1,            ! first year of data in the textfile used to read the cfert* variables 
     >  cfert_lastyear,         ! last year of data in the textfile used to read the cfert* variables 
     >  management_year1,       ! first year of data in planting date & cultivar-related maps
     >  management_lastyear     ! last year of data in planting date & cultivar-related maps


      real 
     >  fertmaize(npoi,max_fert_years), ! historical annual average N-fertilizer applied to maize
     >  cfertmaize(max_cfert_years),         ! estimated  annual average N-fertilizer applied to maize   crops across US Mississippi basin (1950-00) 
     >  fertsoy(npoi,max_fert_years),   ! historical annual average N-fertilizer applied to soybean
     >  cfertsoy(max_cfert_years),           ! estimated  annual average N-fertilizer applied to soybean crops across US Mississippi Basin (1950-00)
     >  fertwheat(npoi,max_fert_years), ! historical annual average N-fertilizer applied to wheat
     >  ndepfact(npoi,100),       ! historical annual average N-deposition factor applied to equations in biogeochem.f          (1950-00)  
     >  cfertwheat(max_cfert_years),         ! estimated  annual average N-fertilizer applied to wheat   crops across US Mississippi Basin (1950-00)  

c WJS 05.06.10: The following variables are used for prescribed-management runs (prescribing planting date & hybrid choice). Note
c that we only read in one file for each variable, regardless of the crop - i.e., unlike (say) the fertilization input, we do NOT
c have separate variables for each crop. Each variable stores the values for the current year.
     >  xinpdate(npoi),   ! planting date input 
     >  xinhybrid(npoi),  ! hybrid input (GDD required from planting to maturity)  
     >  xingrnfill(npoi), ! grnfill input (Fraction of total GDD requirement [given by xinhybrid] needed for grain-fill initiation)
     >  xindaystoharv(npoi), ! days between maturity and harvest input

     >  daygddc(npoi, 366),       ! gdd accumulation for each particular day of the year (corn) 
     >  daygdds(npoi, 366),       ! gdd accumulation for each particular day of the year (soy)  
     >  daygddw(npoi, 366)        ! gdd accumulation for each particular day of the year (wheat)  
c$$$     >  fert05maize(npoi),        ! 2005 average state level fertilizer for maize
c$$$     >  fert05soy(npoi),          ! 2005 average state level fertilizer for soybean
c$$$     >  fert05wheat(npoi),        ! 2005 average state level fertilizer for wheat 
c$$$     >  fert05sorg(npoi)          ! 2005 average state level fertilizer for sorghum

c
      real
     >  cfrac(npft),         ! fraction of dry matter production that is carbon
     >  baset(npft),         ! base temperature used to accumulate gdd for crops
     >  mxtmp(npft),         ! maximum gdd accumulation allowed for each crop pft (per day 0 C), or max temp. for gdd accumulation (meaning of this parameter depends on the GDD function used)
     >  tkill(npft),         ! temperature (K) at which crops are killed due to freeze
     >  hybgdd(npft),        ! maximum gdd for specified hybrids
     >  gddmin(npft),        ! minimum number of annual gdd needed to plant/grow a crop
     >  lfemerg(npft),       ! fraction of annual gdd (to reach physiological mat.) for leaf emergence
     >  grnfill(npft),       ! fraction of annual gdd (to reach physiological mat.) for grain fill initiation
     >  laicons(npft),       ! lai decline factor constant for crops
     >  allconsl(npft),      ! leaf allocation decline scaling factor after grain fill initiation
     >  allconss(npft),      ! stem allocation decline scaling factor after grain fill initiation
     >  laimx(npft),         ! maximum LAI of each crop allowed
     >  arooti(npft),        ! initial allocation of carbon to roots
     >  arootf(npft),        ! allocation of carbon to roots at end of growing season
     >  aleaff(npft),        ! allocation of carbon to leaves at end of growth cycle
     >  astemf(npft),        ! allocation of carbon to stems at end of growth cycle
     >  declfact(npft),      ! rate of LAI decline after grain fill inititation (dimensionless factor)
     >  mxgddgf(npft),       ! maximum gdd allowed past grain fill inititiation for phys. maturity
     >  mxdgfi(npft),        ! maximum number of days past planting allowed before auto-shift to grain fill
     >  mxmat(npft),         ! maximum number of days allowed past planting for physiological maturity to be reached
     >  fleafi(npft),        ! initial fraction of aboveground allocation going to leaf before grain fill begins
     >  fleaf(npft),         ! fraction of aboveground allocation going to leaf before grain fill begins
     >  cgrain(npft),        ! fraction of grain dry matter that is carbon
     >  convfact(npft),      ! factor converting kg/m2 of C in grain to bu/acre value
     >  maxhi(npft),         ! maximum harvest index
     >  fyield(npft),        ! fraction of C allocated to grain pool that is actually seed
     >  fnlfmx(npft),        ! maximum amount of N allowed in leaf at end of growing season
     >  fngrmx(npft),        ! maximum amount of N allowed in grain at end of growing season
     >  sratio(npft),        ! leaf:stem N allocation ratio
     >  rratio(npft),        ! leaf:root N allocation ratio
     >  fnopt(npft),         ! minimum leaf nitrogen content that doesn't experience N stress
     >  bfact(npft),         ! coefficient in LAI curve
     >  daystoharv(npft)     ! number of days between maturity and harvest
     
      real
     >  grnwht(2),            ! fraction of annual gdd (to reach phys. mat.) for wheat leaf emergence
     >  fleafiwht(2),         ! initial fraction of aboveground allocation going to leaf before grain fill begins
     >  fleafwht(2),          ! fraction of aboveground allocation going to leaf before grain fill begins
     >  mgddgf(2),            ! maximum gdd allowed past grain fill inititiation for phys. maturity
     >  mxmatwht(2),          ! maximum number of days allowed past planting for physiological maturity to be reached
     >  fnlfmxw(2),           ! leaf nitrogen maximum allowed for wheat
     >  fngrmxw(2),           ! grain nitrogen maximum allowed for wheat
     >  fnoptw(2)             ! optimum leaf nitrogen fraction for wheat
     
      end module comcrop
