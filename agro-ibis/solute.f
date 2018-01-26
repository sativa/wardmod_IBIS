c last modified C. Kucharik 10.02.02
c--------------------------------------------------------   
      subroutine leaching(irestart, irstyear, istep,iday,
     >                    imonth,iyear,iholdsoiln,iyear0) 
c--------------------------------------------------------
c CJK November 2000
c
c algorithms were derived from the ALEXI model (Anderson et al., 1999)
c it keeps track of inorganic nitrogen pools and its location
c and movement through the soil column profile.  Each layer has a total
c quantity of solute - where all is available to the plant for uptake -
c however, a buffering constant keeps a fixed % N in soil solution at any time while the rest
c is assumed bound to soil aggregates.  This allows for a buffering action to
c be built in where as solution N is decreased, N is allowed to move from
c the total soil pool to solution. 
c
c the inputs of nitrate to the model are : fertilizer, nitrogen deposition
c nitrogen fixation in soybean, nitrogen mineralization through biogeochemistry
c
c the outputs are the leaching through the bottom of the soil profile, and
c plant n uptake through carbon assimilation
c
c uses:
c
      use comgrid
      use compar
      use comsoi
      use comsum
      use comveg
      use comatm
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
     >         istep,
     >         imonth,
     >         iyear,
     >         iyear0,
     >         irestart,
     >         irstyear,
     >         iholdsoiln,
     >         i,
     >         j,
     >         m,
     >         k,
     >         l,
     >         n
c
      real     co,        ! initial total solute mass for the profile (kg-solute m-2) 
     >         cf,        ! buffering factor between inorganic N in soil vs. solution 
     >         tprod,     ! carbon assimilation for timestep by natural vegetation (mol co2 m-2)
     >         cassn,
     >         bal1,      ! solute balance calculation
     >         bal3,      ! solute balance calculation
     >         bal,       ! solute balance calculation after check 
     >         dif        ! difference calculation resulting from balance   
c
       real   cnwood,     ! c/n ratio of wood       - same value as in biogeochem.f
     >        cnfroot,    ! c/n ratio of fine roots - same value as in biogeochem.f
     >        cnleaf,     ! c/n ratio of leaf       - same value as in biogeochem.f 
     >        cndepth,
     >        frootavg,
     >        ypnuptake,
     >        grain,
     >        fyld,
     >        yld,
     >        hi,
     >        tnsoi,
     >        deficit,
     >        sum,
     >        rwork,  
     >        rwork2,
     >        sumcrop,
     >        soilwater,  ! CJK 6-8-16: sum of ice and water contents - and set a min value
     >        fnuptake,
     >        nloss,
     >        tupsoil,
     >        fnitrate    ! fraction of fertilizer addition that is nitrate, and will
                          ! be kept track of for leaching and concentration in solution
                          ! fraction by weight 
c
      real     depth(nsoilay),  ! cumulative depth to bottom of each soil layer
     >         snode(nsoilay)
c
      real     fin(npoi, nsoilay),            ! nitrogen input to layer 
     >         drn(npoi,nsoilay+1)            ! drainage out of layer (mm/timestep) 
c
c set constants
c
       rwork  = dtime * 12.e-03 
c
c these defined variables are also in biogeochem.f for decomposition
c to maintain proper relationship between nitrogen mineralization and
c nitrogen uptake here due to vegetation growth requirements, make sure
c these values are consistent for natural vegetation
c  
      cnwood  = 200.0
      cnfroot =  80.0
      cnleaf  =  60.0
c
      co = 0.0050      ! initial inorganic nitrogen in soil storage pools 
      cf = 0.07        ! buffering constant for inorganic N in soil solution
      fnitrate = 0.90  ! fraction of leached inorganic nitrogen that is nitrate - other portion is ammonium
                       ! values from kr brye
      nloss   = 0.50   ! fraction of immobile inorganic nitrogen that is lost annually to volatilization/denitrification 
      cndepth = 0.0
c 
c set day of year counter
c
      if (iday.eq.1 .and. imonth.eq.1 .and. istep.eq.1) idoy = 0 
      if (istep .eq. 1) idoy = idoy + 1
      
c------------------------------------------------------------------------------
c calculate cumulative depth to bottom of each layer
c this only needs to be called once per run
c depth could be local variable, where snode will be global
c------------------------------------------------------------------------------
c
         do 110 k = 1, nsoilay 
           if (k.eq.1) then
             depth(k) = hsoi(k)
           else
             depth(k) = depth(k-1) + hsoi(k)
           endif
 110     continue
c
         do 120 k = 1, nsoilay
           if (k.eq.1) then 
c            snode(k) = (depth(k+1) - depth(k)) / 2.0 
             snode(k) = (depth(k+1) - 0.) / 2.0 
           else if (k.eq.nsoilay) then
             snode(k) = (depth(k) - depth(k-1)) / 2.0 
           else
             snode(k) = (depth(k+1) - depth(k-1)) / 2.0             
           endif
           cndepth = cndepth + snode(k)
 120     continue
c
c------------------------------------------------------------------------------
c begin grid
c
       do 100 i = 1, npoi   
         if(cmask(i) == 0 ) cycle
c        
c set concentrations for first year of restart (beginning)
c
         if (irestart .eq. 1 .and. iyear .eq. irstyear .and.
     >       istep .eq. 1 .and. iday .eq. 1 .and. imonth .eq. 1) then
c
             ctoti(i) = 0.0
           do 105 k = 1, nsoilay
c
             if (iholdsoiln .eq. 0) then
               smsoil(i,k) = smsoil(i,k) * 0.0
               smsoln(i,k) = smsoln(i,k) * 0.0 
             endif
c
             csoln(i,k)  = smsoln(i,k) * (1.e+09)/ 
     >                     (1.e+04 * (wsoi(i,k) + wisoi(i,k)) * snode(k) * 10.)   
             ctoti(i)    = ctoti(i) + smsoil(i,k) + smsoln(i,k)
 105       continue
             ctot(i)     = ctoti(i)
         endif
c      
c set initial concentrations for first timestep of model run
c
         if (iyear.eq.iyear0 .and. istep.eq.1 .and. 
     >       iday.eq.1 .and. imonth.eq.1) then 
c
           ctoti(i)    = 0.0  ! total initial inorganic nitrogen in profile
           ctot(i)     = 0.0  ! total current inorganic nitrogen in soil profile (solution and soil) (kgm-2)
           drntot(i)   = 0.0  ! total drainage 
           ftot(i)     = 0.0  ! total leaching of inorganic  nitrogen
           yno3leach(i)= 0.0  ! total annual nitrate leaching
           assimn(i)   = 0.0   
           snbalance(i)= 0.0  ! soil nitrogen balance
c
           do 130 k = 1, nsoilay
c            smsoil(i,k) = co * (snode(k) / depth(nsoilay))
             smsoil(i,k) = co * (snode(k) / cndepth)
             smsoln(i,k) = cf * smsoil(i,k)
c
c adjust mass balance
c
             smsoil(i,k) = smsoil(i,k) - smsoln(i,k)
c
c check conversion factor from kg m-2 (ibis) to kg ha-1 (kris' equations) to mg liter-1
c concentration in solution
c have to account for ice fraction because wsoi can be 0.0 for much of the year 
c
             csoln(i,k) = smsoln(i,k) * (1.e+09)/ 
     >                     (1.e+04 * (wsoi(i,k) + wisoi(i,k)) * snode(k) * 10.)   
             ctoti(i)    = ctoti(i) + smsoil(i,k) + smsoln(i,k)
 130       continue
c      
         else if (iday.eq.1 .and. imonth.eq.1 .and. istep.eq.1) then
c
c initialize annual drainage and solute leached to zero at beginning of year
c
           drntot(i)   = 0.0     ! total drainage through profile (mm)
           ftot(i)     = 0.0     ! cumulative inorganic nitrogen leached (kg m-2)
           yno3leach(i)= 0.0     ! annual total nitrate leaching (kg no3 m-2)
           assimn(i)   = 0.0
           totnvegn(i) = 0.0     ! cumulative inorganic nitrogen uptake by natural vegetation (kg m-2 y-1)  
           taninp(i)   = 0.0     ! total annual inputs of inorganic nitrogen 
           snbalance(i)= 0.0     ! annual soil nitrogen balance
c
           ctoti(i)     = ctot(i) ! initial total inorganic N in profile at beginning of timestep 
           tsinp(i)     = 0.0     ! total daily inorganic N input
           tslay(i)     = 0.0     ! total daily inorganic N input to depth determined by soilay
           dtnleach(i)  = 0.0     ! daily total inorganic nitrogen solute (ammonium and nitrate) leached
           dnileach(i)  = 0.0     ! daily total nitrate-nitrogen leached
           tpnuptake(i) = 0.0     ! total daily plant nitrogen uptake
           ctot(i)      = 0.0     
           ddrn(i)      = 0.0     ! daily total drainage at specified depth in profile 
c
         else if (istep .eq. 1) then
c
c initialize daily sums used in mass balance calculation
c
           ctoti(i)     = ctot(i) ! initial total inorganic N in profile at beginning of timestep 
           tsinp(i)     = 0.0     ! total daily inorganic N input
           tslay(i)     = 0.0     ! total daily inorganic N input to depth determined by soilay
           dtnleach(i)  = 0.0     ! daily total inorganic nitrogen solute (ammonium and nitrate) leached
           dnileach(i)  = 0.0     ! daily total nitrate-nitrogen leached
           tpnuptake(i) = 0.0     ! total daily plant nitrogen uptake
           ctot(i)      = 0.0     
           ddrn(i)      = 0.0     ! daily total drainage at specified depth in profile 
c
c add denitrification/volatilization loss on daily timestep
c
c           do 125 k = 1, nsoilay
c
c              smsoil(i,k) = smsoil(i,k) * (1. - (nloss/ndaypy))
c
c 125       continue
c
         else  ! all other timesteps
c
           fin(i,1)  = 0.0     ! initialize input to top layer each timestep
           ctot(i)   = 0.0
c
         endif
c
c------------------------------------------------------------------------------
c calculating drainage (mm) through each node per timestep
c note that the drainage is from the above layer through the top of (k) - thus
c to get drainage out of the profile, have to assume that an additional hypothetical
c layer exists below, where bperm influences the gdrain calculated in soil.f  
c------------------------------------------------------------------------------
c 
           drn(i,1) = max(0.0, (fwtop(i) + fwpud(i))* dtime)  ! water infiltration (mm) into top layer
c
           do 140 k = 2, nsoilay+1      ! drainage out the bottom is an additional soil layer
c
             drn(i,k) = max(0.0, wflo(i,k) * dtime)  ! rate of drainage through each node -
c                                                   ! calculated in soil.f (kg m-2 s-1)
c                                                   ! drn(i,nsoilay+1) is equal to gdrain(i)  
 140       continue
c
c designate which layer in the profile you are tracking drainage 
c i.e. KRB measures drainage at ARL at 1.4 m with lysimeter
c
c          drntot(i) = drntot(i) + drn(i,nsoilay+1) 
           drntot(i) = drntot(i) + drn(i,isoilay) 
           ddrn(i)   = ddrn(i)   + drn(i,isoilay)
c
c determine nitrogen uptake by natural vegetation
c
             tprod = 0.0
c
             do 215 j = 1, 12     ! loop through natural vegetation types to check for carbon
                                  ! assimilation 
c
c total carbon assimilation for the timestep (tnpp has units of mol-co2 /m-2 / s)
c
               tprod = tprod + tnpp(i,j) * dtime  
 215         continue
c
c if natural vegetation has positive carbon assimilation for this timestep then calculate 
c how much nitrogen was required for growth based on carbon/nitrogen ratios dictated in biogeochem.f
c 
             cassn = 0.0
             if (tprod .gt. 0.0) then 
c
               do 220 j = 1, 12    ! loop through natural vegetation types
c
c calculate assimilated carbon for this timestep and nitrogen required
c to satisfy growth
c
c CJK 5/18/07               cassn = cassn + (tnpp(i,j) * rwork * awood(i,j) / cnwood) + 
c     >                         (tnpp(i,j) * rwork * aleaf(i,j) / cnleaf) +
c     >                         (tnpp(i,j) * rwork * aroot(i,j) / cnfroot) 
c
               cassn = cassn + (tnpp(i,j) * rwork * awood(i,j) / 180) + 
     >                         (tnpp(i,j) * rwork * aleaf(i,j) / 40) +
     >                         (tnpp(i,j) * rwork * aroot(i,j) / 60) 
c
 220           continue
c
               assimn(i) = assimn(i) + cassn
               tupsoil = 0.0 
c
c calculate total transpiration stream for both [u]pper and [l]ower
c canopies - calculated in canopy.f
c
               do 230 l = 1, nsoilay
                 tupsoil = tupsoil + (upsoiu(i,l) + upsoil(i,l)) * dtime 
 230           continue
c
c determine that the fraction of the total nitrogen required for uptake
c is proportional to the fraction of the total transpiration stream
c contribution of that layer 
c
c anuptake is same variable used in crops.f for actual nitrogen uptake for that
c timestep
c
               do 240 l = 1, nsoilay
c
                 if (tupsoil .gt. 0.0) then
                   fnuptake      = ((upsoiu(i,l) + upsoil(i,l)) * dtime) / tupsoil 
                 else
                   fnuptake = 0.0
                 endif
c
c add total nitrogen uptake for whole grid cell - could eventually
c be a combinatation of natural vegetation and crops
c anuptake is also being calculated in crops.f 
c 
                 anuptake(i,l) = min((smsoil(i,l)+smsoln(i,l)), fnuptake * cassn)
c
c update annual nitrogen uptake by natural vegetation (forests/grasses/shrubs) 
c
                 totnvegn(i)   = totnvegn(i) + anuptake(i,l) 
c
 240           continue
c
            endif  ! natural vegetation nitrogen uptake 
c
c------------------------------------------------------------------------------
           do 200 k = 1, nsoilay
c
c calculate the average root profile in the grid cell for upper and lower canopies
c combined
c
             frootavg    = (froot(k,1) + froot(k,2)) / 2.0
c
             if (k .eq. 1 .and. istep .eq. 1) then 
c
c first timestep/top layer - each day checked
c
               do 210 j = 1, npft
c               
c managed crop ecosystems: 
c human dependency - get at time of planting
c this could be changed in crops.f - to change timing/method of nitrogen fertilizer management
c value of fertnitro(j) to be 0 at other times of year
c
                 if (exist(i,j) .ne. 0.) then
                   if (fertnitro(i,j) .gt. 0.0 .and. croppresent(i,j) .eq. 1) then
                     fin(i,1) = fertnitro(i,j) 
                   else 
                     fin(i,1) = 0.
                   endif 
                 endif
 210           continue
c
c add nitrogen deposition to the amount coming into the top
c layer for that particular day - based on daily precipitation
c and calculated in biogeochem.f 
c only added on beginning timestep of each day
c
               fin(i,1) = fin(i,1) + deposn(i) 
c
             else
               fout(i,0) = 0.0
               nout(i,0) = 0.0
               fin(i,k)  = fout(i,k-1) 
             endif
c
c nitrogen movement - potential inputs are fertilizer, n-deposition, n-fixation,n-mineralization.  
c nitrogen mineralization is calculated in biogeochem.f as a daily rate and
c was converted to mole-N s-1
c
c have to assume where nitrogen mineralization is taking place in the profile -
c and how amount of atmopheric n-fixation (daily quantity) is being added by roots
c during this timestep (dtime)
c use a weighted average of fine root distribution profiles for lower and upper plant canopies
c fine root distribution is froot - initialized in initial.f  
c
c value for anuptake is calculated in nitrostress - which is applied to
c the vmax rate for each timestep 
c
c if crops are planted, no n-fixation used from biogeochem.f 
c that n-fixation is assumed from natural vegetation types that fix atmospheric N2
c
c n-fixation from soybeans is incorporated into soil inorganic N pools in crops.f
c 
c             if (icropsum(i) .gt. 0.) then
c                 fixsoin(i) = 0.0
c                 yfixsoin(i) = 0.0
c              endif
c
c assume natural nitrogen fixation = 0 
c cjk 11.18.01
c
             fixsoin(i) = 0.0
             yfixsoin(i) = 0.0
c
             smsoil(i,k) = smsoil(i,k) + fin(i,k) +
     >                     fixsoin(i) * frootavg * dtime / 86400. +       
     >                     (tnmin(i)  * frootavg * dtime * 0.014) -  
     >                     anuptake(i,k)  
c
c if total plant nitrogen uptake cannot be derived entirely from the immobile pool
c then remove it from the solution pool
c
             if (smsoil(i,k) .lt. 0.0) then
               deficit     = smsoil(i,k) 
               smsoln(i,k) = max(0.0, smsoln(i,k) + deficit) 
               smsoil(i,k) = 0.0
             endif
c
c keep sum of daily plant nitrogen uptake for mass balance calculation
c
             tpnuptake(i) = tpnuptake(i) + anuptake(i,k) 
             smsoil(i,k) = smsoil(i,k) - (cf*smsoil(i,k) - smsoln(i,k))
             smsoln(i,k) = smsoln(i,k) + (cf*smsoil(i,k) - smsoln(i,k))
c
c account for ice fraction when calculating concentration in solution
c
c CJK 6-8-2016 To prevent division by zero when wisoi and wsoi both
c happen to be equal to 0.0.  This happens in deeper layers on 
c occassion when driving the model with ZedX climate data.  Was not
c apparent with CRU/NCEP.
c
             soilwater = min(1.0,(max(0.001, wisoi(i,k)+wsoi(i,k))))
             csoln(i,k)  = smsoln(i,k) * (1.e+09)/ 
     >                     (1.e+04 * soilwater * snode(k) * 10.)   

C             if(wisoi(i,k) /= 0.0 ) then  !Y.Li
C               csoln(i,k)  = smsoln(i,k) * (1.e+09)/ 
C     >                       (1.e+04 * (wisoi(i,k) + wsoi(i,k)) * snode(k) * 10.)   
C             end if
c
c drainage is designated by subscript of layer it is going INTO 
c fout is for current layer 
c
c KRBs units are kg ha-1  - for fout
c since we are kg m-2     - need to reduce fout by larger constant
c
c we are interested in tracking nitrate only for leaching - add factor in these
c equations to account for the fact that csoln can contain both ammonium and nitrate 
c the fertnitro addition assumes that both nitrate and ammonium are added together
c
c            fout(i,k)   = min(smsoln(i,k), csoln(i,k) * drn(i,k+1) / 100.0)
             fout(i,k)   = min(smsoln(i,k), csoln(i,k) *
     >                     drn(i,k+1) / 1e+06)
             nout(i,k)   = min(smsoln(i,k), csoln(i,k) * fnitrate *
     >                     drn(i,k+1) / 1e+06)
             fout(i,k)   = max(0.0, fout(i,k))
             nout(i,k)   = max(0.0, nout(i,k))
c
c remove leached inorganic-N from layer reassign total nitrogen to each layer available to leach/total
c redistribute inorganic-N in layer between solution and total in soil
c
             smsoln(i,k) = smsoln(i,k) - fout(i,k)
             smsoil(i,k) = smsoil(i,k) - (cf * smsoil(i,k) - smsoln(i,k))
             smsoln(i,k) = smsoln(i,k) + (cf * smsoil(i,k) - smsoln(i,k))
c
             if (smsoil(i,k) .lt. 0.0) then
               smsoil(i,k) = 0.0
               smsoln(i,k) = 0.0
             endif
c
             ctot(i) = ctot(i) + smsoln(i,k) + smsoil(i,k)
c
c account for ice fraction in soil when calculating solute concentration
c 
c CJK 6-8-2016 To prevent division by zero when wisoi and wsoi both
c happen to be equal to 0.0.  This happens in deeper layers on 
c occassion when driving the model with ZedX climate data.  Was not
c apparent with CRU/NCEP.
c 
             csoln(i,k)  = smsoln(i,k) * (1.e+09)/ 
     >                     (1.e+04 * soilwater * snode(k) * 10.)   
c    >                     (1.e+04 * (wisoi(i,k)+ wsoi(i,k)) * snode(k) * 10.)   

c             if( wisoi(i,k) /= 0.0 ) then !Y.Li
c               csoln(i,k)  = smsoln(i,k) * (1.e+09)/ 
c     >                       (1.e+04 * (wisoi(i,k)+ wsoi(i,k)) * snode(k) * 10.)   
c             end if
c
c assign daily nitrate concentration for 1.4 m layer and
c total amount of daily drainage through that layer for arlington comparisons 
c
             if (istep .eq. 86400. / dtime) then
               daynconc(i,idoy) = fnitrate * adcsoln(i,isoilay)
               daydrn(i,idoy)   = ddrn(i)
             endif
c             
 200      continue  ! loop through all soil layers
c
c total inputs into the soil for this daily timestep
c
             tsinp(i) = tsinp(i) +  fin(i,1) + fixsoin(i) * dtime / 86400. +
     >                  tnmin(i) * dtime * 0.014        
c
c calculate the total daily inputs to top number of layers - determined by the
c soilay constant - to help in mass balance calculation to that depth
c
             do 280 k = 1, isoilay
               frootavg    = (froot(k,1) + froot(k,2)) / 2.0
               tslay(i)    = tslay(i) + tnmin(i) * frootavg * dtime * 0.014 +
     >                       fixsoin(i) * frootavg *  dtime / 86400.               
 280         continue
c
             tslay(i) = tslay(i) + fin(i,1)
c
c calculate total amount of total nitrogen leached out of entire profile
c for daily timestep
c
             dtnleach(i) = dtnleach(i) + fout(i,nsoilay) 
c
c convert to a rate (y-1) at end of day for daily output for Simon
c base on rate kg nitrate per hectare - for layer we are interested in designated
c as input to baseflow
c
             dnileach(i) = dnileach(i) + nout(i,isoilay) 
c             
c update annual total nitrate-nitrogen leaching - covert to kg/ha
c trying to compare to KRBs measurements at 1.4 m in profile
c or input to baseflow at soilay = 1.5 m for regional modeling 
c
             ftot(i)      = ftot(i)      + fout(i,isoilay)  * 1e+04
             yno3leach(i) = yno3leach(i) + nout(i,isoilay)  * 1e+04
c            ftot(i)      = ftot(i)      + nout(i,nsoilay) * 1e+04
c            ftot(i)      = ftot(i)      + nout(i,isoilay)  * 1e+04
c
c end of year calculation for flow-weighted mean nitrate concentration 
c
             if (istep .eq. 86400. / dtime .and.
     >           imonth .eq. 12 .and. iday .eq. 31) then 
c
c put in check for division by zero for annual drainage
c 
               if (drntot(i) .le. 0) then
                 concn(i) = 0.0
               else
                 sum = 0
                 do 250 l = 1, idoy
                   sum = sum + (daydrn(i,l)/drntot(i)) * daynconc(i,l)
 250             continue
                 concn(i) = sum
               endif
c
             endif
c
c------------------------------------------------------------------------------
c calculate mass balance approach at each day to make sure solute is being
c conserved 
c------------------------------------------------------------------------------
c
c calculate each timestep addition of nitrogen fixation from nitrostress routine
c in crops.f 
c
                  do 275 j = 1, npft
c
c only add fixed nitrogen to taninp for the top soil layers according to soilay
c
                    do 290 k = 1, isoilay
                      taninp(i) = taninp(i) + fixn(i,j) * froot(k,1)
 290                continue
c
c add fixation - which is a total for each timestep - to the total soil inputs
c 
                    tsinp(i)  = tsinp(i)  + fixn(i,j)
c
 275              continue
c
             if (istep .eq. 86400. / dtime) then
c
               taninp(i) = taninp(i) + tslay(i) 
c
               bal1  = ctoti(i)   + tsinp(i) -
     >                 dtnleach(i) - 
     >                 ctot(i)    - tpnuptake(i)
c
               if (bal1 .ne. 0.0) then      ! excess inputs vs. outputs 
                 dif = bal1
                 ctot(i) = 0.0
                 do 320 k = 1, nsoilay
c                  smsoil(i,k) = smsoil(i,k) + (dif*snode(k) / depth(nsoilay)) 
                   smsoil(i,k) = smsoil(i,k) + (dif*snode(k) / cndepth) 
                   smsoil(i,k) = smsoil(i,k) - (cf * smsoil(i,k) - smsoln(i,k))
                   smsoln(i,k) = smsoln(i,k) + (cf * smsoil(i,k) - smsoln(i,k))
                   ctot(i)     = ctot(i) + smsoil(i,k) + smsoln(i,k)
 320             continue
c
               endif
c
c recalculate the solute balance
c
c               bal = ctoti(i) + tsinp(i) - dtnleach(i) -
c     >               ctot(i)  - tpnuptake(i) 
c
c convert to a rate (y-1) at end of day for daily output for Simon
c base on rate kg per hectare
c
                dnileach(i) = dnileach(i) * 1.e+04 * 365. 
c
c cropn is in units of kg/ha 
c
                  ypnuptake = 0.0
                  do 335 j = 1, npft
                    ypnuptake = ypnuptake + cropn(i,j)
 335              continue
                    ypnuptake = ypnuptake + totnvegn(i) * 1.e+04
c
c calculate nitrogen balance for soil - crops - natural vegetation - inputs
c in kg ha-1
c
                  snbalance(i) = taninp(i)*1.e+04 - ypnuptake -
     >                           ftot(i)
c
           endif
 100    continue
c
        return 
        end
