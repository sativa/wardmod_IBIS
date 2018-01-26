c
c ---------------------------------------------------------------------
      module comhour
c     last update 10.02.02 C. Kucharik 
c ---------------------------------------------------------------------
c
c this module is to hold values for variables related to inputting
c a data file of meteorological quantities (hourly) to be used to drive
c IBIS using real weather data
c
      implicit none
      save

      integer
     > ndpts
c
      parameter (ndpts=9000) ! the number of data points in the file 
c
      integer  
     > iiyear(ndpts),
     > iimonth(ndpts),
     > iiday(ndpts),
     > iihour(ndpts)  
c
      real
     > var1(ndpts),
     > var2(ndpts), 
     > var3(ndpts), 
     > var4(ndpts), 
     > var5(ndpts), 
     > var6(ndpts) 
c
      end module comhour
