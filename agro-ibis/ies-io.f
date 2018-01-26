* ies-io.f last update 1/27/10 Bill Sacks, wsacks@gmail.com
* previous update 4/1/98 Christine C. Molling, cmolling@facstaff.wisc.edu
* University of Wisconsin - Madison; Institute for Environmental Studies;
* Center for Climatic Research; Climate, People, and Environment Program
*
* You are free to change or distribute my code as long as you
* 1) acknowlege me, Christine C. Molling (cmolling@facstaff.wisc.edu)
* as the author of the original code and 
* 2) do not sell it.
* Of course if there is anything wrong with the code, or you somehow
* encounter damage to your reputation/wallet because of its use, you
* are forbidden to sue me or anyone else.
*
* Read Optimization modifications by S. Slavin, 06 May 2016 
* Includes module with additional debugging functions, variables, etc.
* adds a dependency to: use rdopt
* The primary optimization is to not close files once they are opened,
* but to leave them open until close_logged_ncids() in rdopt.f is called.
* This behavior occurs when the variable flg_rdopt
* (declared in compar.f) is set.
* Slavin modifications tweeked by Cicada Dennis May 2016.
*
*     These subroutines can be used as a primary interface to the ies
* format netcdf files.  You can also use the slightly lower level functions
* found in ies.f or the low level netcdf commands (see the netcdf V3
* users guide).
*SUBROUTINE LIST:
* inifile - create a new file with dimensions and variables for
*  longitude, latitude, an optional 3rd real dimension, and time
*  Time must have units of days since a date!!!
* inifilec - create a new file with dimensions and variables for
*  longitude, latitude, a 3rd character-variable dimension, and time
*  Time must have units of days since a date!!!
* inivar - initialize a new variable
* endini - end initialization of file: sync file to disk, close file
* writevar - write a 2-, 3-, or 4-dimensional hyperslab of data
* readvar - read a 2-, 3-, or 4-dimensional hyperslab of data
* dimlen - get the length of a dimension
*
*.....................................................................
* inifile - create a new file with dimensions and variables for
*  longitude, latitude, an optional 3rd real dimension, and time.
*  Time must have units of days since a date!!!
*.....................................................................
      subroutine inifile(idies,filen,title,source,history,nlon,alons,
     . nlat,alats,name3rd,long3rd,units3rd,n3rd,vals3rd,pos3rd,tunits,
     . calendar,ierror)

*INPUT
*     filen - character*(*) - name for new file
*     title - character*(*) - title for file
*     source - character*(*) - source for file
*     history - character*(*) - date of creation, and possibly author
*     nlon - integer - number of point along longitude direction
*     alons - real(nlon) - values of longitudes
*     nlat - integer - number of point along latitude direction
*     alats - real(nlat) - values of latitudes
*     name3rd - character*(*) - name of 3rd dimension - use '' if you
*      don't want this dimension
*     long3rd - charqcter*(*) - long name for 3rd dimension variable
*      (ignored if nam3rd='')
*     n3rd - integer - size of 3rd dimension (ignored if name3rd='')
*     vals3rd - real(n3rd) - values along 3rd dimension (ignored if 
*      name3rd='')
*     pos3rd - character*(*) - if the 3rd dimension exists and refers to
*      a measured level or height, use 'up' if values increase with distance
*      above earth (as in height), or 'down' if values decrease with distance 
*      (as in pressure) (ignored if name3rd=''). If 3rd dimension does not
*      refer to a height or level, use ''.
*     tunits - character*(*) - units for time, must be in the form of days
*      since yyyy-mm-dd tt:tt:tt, where yyyy is a year or 0000, mm is a
*      month, dd is a day, and tt:tt:tt is (optional) time or 00:00:00
*     calendar - character*(*) - type of calendar.  Choose from 'noleap',
*      'gregorian','n ka BP', etc.  Use iescal if orbital parameters
*      differ from modern.
*OUTPUT
*     idies - integer - id number of new file for later use
*     ierror - integer - error code, 0 = no error, < 0 = an error occured

      integer idies, nlon, nlat, n3rd, ierror
      character*(*) filen, title, source, history, name3rd
      character*(*) long3rd, units3rd, pos3rd, tunits, calendar
      real alons(nlon), alats(nlat), vals3rd(n3rd)

      include 'netcdf.inc'

! S. Slavin added NF_WRITE_DISKLESS.
      integer :: NF_WRITE_DISKLESS = IOR( NF_WRITE, NF_DISKLESS )

      integer ierr, iddlon, idlon, iddlat, idlat, idd3rd, id3rd
      integer iddtime, idtime, iddlend, iddate, idtw, ll
      integer idims(2)

      ierror = 0

*     Open file
*     ---------
! S. Slavin changed NF_CREATE.
      ierr = NF_CREATE(filen,NF_WRITE_DISKLESS,idies)
! Cicada Dennis changed back to orig (below), then changed back to Slavin's
!     ierr = NF_CREATE(filen,NF_CLOBBER,idies)

      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Define global attributes
*     ------------------------
      ll = lenchr(title)
      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'title',ll,
     . title(1:ll))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ll = lenchr(source)
      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'source',ll,
     . source(1:ll))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ll = lenchr(history)
      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'history',ll,
     . history(1:ll))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ll = lenchr(calendar)
      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'calendar',ll,
     . calendar(1:ll))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'conventions',9,
     . 'NCAR-CSM'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Define dimensions
*     -----------------
      ierr = NF_DEF_DIM(idies,'longitude',nlon,iddlon)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_DEF_DIM(idies,'latitude',nlat,iddlat)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      idd3rd = -999
      if (name3rd .ne. '') then
         ierr = NF_DEF_DIM(idies,name3rd,n3rd,idd3rd)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in inifile'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
      end if
      ierr = NF_DEF_DIM(idies,'lengthd',10,iddlend)  
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_DEF_DIM(idies,'time',NF_UNLIMITED,iddtime)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Define dimension variables
*     --------------------------
      ierr = NF_DEF_VAR(idies,'longitude',NF_FLOAT,1,iddlon,idlon)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_DEF_VAR(idies,'latitude',NF_FLOAT,1,iddlat,idlat)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      if (idd3rd .ne. -999) then
         ierr = NF_DEF_VAR(idies,name3rd,NF_FLOAT,1,idd3rd,id3rd)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in inifile'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
      end if
      ierr = NF_DEF_VAR(idies,'time',NF_FLOAT,1,iddtime,idtime)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Define variables time_weights and date
*     --------------------------------------
      ierr = NF_DEF_VAR(idies,'time_weights',NF_FLOAT,1,iddtime,idtw)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      idims(1)=iddlend  
      idims(2)=iddtime
      ierr = NF_DEF_VAR(idies,'date',NF_CHAR,2,idims,iddate)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if


*     Attributes for dimension variables
*     ----------------------------------
      ierr = NF_PUT_ATT_TEXT(idies,idlon,'long_name',10,
     . 'longitude'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idlon,'units',13,
     . 'degrees_east'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idlat,'long_name',9,
     . 'latitude'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idlat,'units',14,
     . 'degrees_north'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      if (idd3rd .ne. -999) then
         ll = lenchr(long3rd)
         ierr = NF_PUT_ATT_TEXT(idies,id3rd,'long_name',
     .    ll,long3rd(1:ll))
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in inifile'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
         if (units3rd .eq. '') then
c handle empty string, which breaks the subsetting code in the 'else' clause (units3rd(1:ll))
c WJS 02.10.10: Ideally this check would be done for other input strings too, but for now I
c just do it for units3rd because that's the only one where the input string is sometimes ''.
           ierr = NF_PUT_ATT_TEXT(idies,id3rd,'units',
     .       0, '')
         else
           ll = lenchr(units3rd)
           ierr = NF_PUT_ATT_TEXT(idies,id3rd,'units',
     .       ll,units3rd(1:ll))
         end if
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in inifile'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idtime,'long_name',5,
     . 'time'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ll = lenchr(tunits)
      ierr = NF_PUT_ATT_TEXT(idies,idtime,'units',ll,
     . tunits(1:ll))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idtw,'long_name',29,
     . 'number of days per time step'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idtw,'units',5,
     . 'days'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,iddate,'long_name',25,   
     . 'label for each time step'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ll = lenchr(tunits)
      ierr = NF_PUT_ATT_TEXT(idies,iddate,'units',1,
     . char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     End define mode, enter data mode
*     --------------------------------
      ierr = NF_ENDDEF(idies)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Put dimension variables except for time
*     ---------------------------------------
      ierr = NF_PUT_VAR_REAL(idies,idlon,alons)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_VAR_REAL(idies,idlat,alats)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      if (idd3rd .ne. -999) then
         ierr = NF_PUT_VAR_REAL(idies,id3rd,vals3rd)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in inifile'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
      end if

*     Don't close file
*     ----------------

      return
      end
*.....................................................................
* inifilec - create a new file with dimensions and variables for
*  longitude, latitude, a 3rd character-variable dimension, and time
*  Time must have units of days since a date!!!
*.....................................................................
      subroutine inifilec(idies,filen,title,source,history,nlon,alons,
     . nlat,alats,name3rd,long3rd,units3rd,n3rd,len3rd,chars3rd,
     . tunits,calendar,ierror)

*INPUT
*     filen - character*(*) - name for new file
*     title - character*(*) - title for file
*     source - character*(*) - source for file
*     history - character*(*) - date of creation, and possibly author
*     nlon - integer - number of point along longitude direction
*     alons - real(nlon) - values of longitudes
*     nlat - integer - number of point along latitude direction
*     alats - real(nlat) - values of latitudes
*     name3rd - character*(*) - name of 3rd dimension
*     long3rd - charqcter*(*) - long name for 3rd dimension variable
*     n3rd - integer - size of 3rd dimension
*     len3rd - integer length of chracter strings in vals3rd
*     chars3rd - character*len3rd(n3rd) - values along 3rd dimension
*     tunits - character*(*) - units for time, must be in the form of days
*      since yyyy-mm-dd tt:tt:tt, where yyyy is a year or 0000, mm is a
*      month, dd is a day, and tt:tt:tt is (optional) time or 00:00:00
*     calendar - charater*(*) - type of calendar.  Choose from 'noleap',
*      'gregorian','n ka BP', etc.  Use iescal if orbital parameters
*      differ from modern.
*OUTPUT
*     idies - integer - id number of new file for later use
*     ierror - integer - error code, 0 = no error, < 0 = an error occured

      integer idies, nlon, nlat, n3rd, len3rd, ierror
      character*(*) filen, title, source, history, name3rd
      character*(*) long3rd, units3rd, tunits, calendar
      real alons(nlon), alats(nlat)
      character chars3rd(len3rd,n3rd)

      include 'netcdf.inc'

! S. Slavin added NF_WRITE_DISKLESS.
      integer :: NF_WRITE_DISKLESS = IOR( NF_WRITE, NF_DISKLESS )

      integer ierr, iddlon, idlon, iddlat, idlat, idd3rd, id3rd
      integer iddlen, iddtime, idtime, iddate, idtw, ll, idims(2)
      integer iddlend

      ierror = 0

*     Open file
*     ---------
! S. Slavin changed NF_CREATE.
      ierr = NF_CREATE(filen,NF_WRITE_DISKLESS,idies)
! Cicada Dennis changed back to orig (below), then changed back to Slavin's
!     ierr = NF_CREATE(filen,NF_CLOBBER,idies)

      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Define global attributes
*     ------------------------
      ll = lenchr(title)
      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'title',ll,
     . title(1:ll))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ll = lenchr(source)
      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'source',ll,
     . source(1:ll))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ll = lenchr(history)
      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'history',ll,
     . history(1:ll))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ll = lenchr(calendar)
      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'calendar',ll,
     . calendar(1:ll))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'conventions',9,
     . 'NCAR-CSM'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Define dimensions
*     -----------------
      ierr = NF_DEF_DIM(idies,'longitude',nlon,iddlon)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_DEF_DIM(idies,'latitude',nlat,iddlat)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_DEF_DIM(idies,name3rd,n3rd,idd3rd)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_DEF_DIM(idies,'length',len3rd,iddlen)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_DEF_DIM(idies,'lengthd',10,iddlend)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

      ierr = NF_DEF_DIM(idies,'time',NF_UNLIMITED,iddtime)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Define dimension variables
*     --------------------------
      ierr = NF_DEF_VAR(idies,'longitude',NF_FLOAT,1,iddlon,idlon)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_DEF_VAR(idies,'latitude',NF_FLOAT,1,iddlat,idlat)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      idims(1) = iddlen
      idims(2) = idd3rd
      ierr = NF_DEF_VAR(idies,name3rd,NF_CHAR,2,idims,id3rd)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_DEF_VAR(idies,'time',NF_FLOAT,1,iddtime,idtime)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Define variables time_weights and date
*     --------------------------------------
      ierr = NF_DEF_VAR(idies,'time_weights',NF_FLOAT,1,iddtime,idtw)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      idims(1)=iddlend
      idims(2)=iddtime
      ierr = NF_DEF_VAR(idies,'date',NF_CHAR,2,idims,iddate)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Attributes for dimension variables
*     ----------------------------------
      ierr = NF_PUT_ATT_TEXT(idies,idlon,'long_name',10,
     . 'longitude'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idlon,'units',13,
     . 'degrees_east'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idlat,'long_name',9,
     . 'latitude'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idlat,'units',14,
     . 'degrees_north'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      if (idd3rd .ne. -999) then
         ll = lenchr(long3rd)
         ierr = NF_PUT_ATT_TEXT(idies,id3rd,'long_name',
     .    ll,long3rd(1:ll))
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in inifilec'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
         if (units3rd .eq. '') then
c handle empty string, which breaks the subsetting code in the 'else' clause (units3rd(1:ll))
c WJS 02.10.10: Ideally this check would be done for other input strings too, but for now I
c just do it for units3rd because that's the only one where the input string is sometimes ''.
           ierr = NF_PUT_ATT_TEXT(idies,id3rd,'units',
     .       0, '')
         else
           ll = lenchr(units3rd)
           ierr = NF_PUT_ATT_TEXT(idies,id3rd,'units',
     .       ll,units3rd(1:ll))
         end if
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in inifilec'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idtime,'long_name',5,
     . 'time'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ll = lenchr(tunits)
      ierr = NF_PUT_ATT_TEXT(idies,idtime,'units',ll,
     . tunits(1:ll))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idtw,'long_name',29,
     . 'number of days per time step'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idtw,'units',5,
     . 'days'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,iddate,'long_name',25,
     . 'label for each time step'//char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ll = lenchr(tunits)
      ierr = NF_PUT_ATT_TEXT(idies,iddate,'units',1,
     . char(0))
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifile'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     End define mode, enter data mode
*     --------------------------------
      ierr = NF_ENDDEF(idies)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Put dimension variables except for time
*     ---------------------------------------
      ierr = NF_PUT_VAR_REAL(idies,idlon,alons)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_VAR_REAL(idies,idlat,alats)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_VAR_TEXT(idies,id3rd,chars3rd)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inifilec'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Don't close file, <inivar> needs access, <endini> will close it
*     ----------------

      return
      end
*.....................................................................
* inivar - initialize a new variable
*.....................................................................
      subroutine inivar(idies,varname,longname,units,ndims,dimnames,
     . valmissing,ierror)

*INPUT
*     idies - integer - id number of a new file from a previous call
*      to inifile, inifilec, or iesopen
*     varname - charcter*(*) - name for new variable
*     longname - character*(*) - descriptive name for variable
*     units - character*(*) - units for variable, compatible with udunits
*     ndims - integer - number of dimensions for this variable
*     dimnames - character*(*)(ndims) - name of each dimension, in order
*     valmissing - real - missing value, ignored if valmissing=0.
*OUTPUT
*     ierror - integer - error code, 0 = no error, < 0 = error occurred

      integer idies, ndims, ierror
      character*(*) varname, longname, units, dimnames(ndims)
      real valmissing

      include 'netcdf.inc'

      integer ierr, idvar, iddims(4)

      ierror = 0

*     Put into redefine mode, but don't fail if already in define mode
*     ----------------------------------------------------------------
      ierr = NF_REDEF(idies)
      if (ierr .ne. NF_NOERR) then
         print *, 'Warning in inivar'
         print *, NF_STRERROR(ierr)
      end if

*     Find id's of dimensions
*     -----------------------
      do i = 1, ndims
         ierr = NF_INQ_DIMID(idies,dimnames(i),iddims(i))
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in inivar'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
      end do

*     Define variable
*     ---------------
      ierr = NF_DEF_VAR(idies,varname,NF_FLOAT,ndims,iddims,idvar)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inivar'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Define attributes
*     -----------------
      ierr = NF_PUT_ATT_TEXT(idies,idvar,'long_name',
     . lenchr(longname),longname)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inivar'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      ierr = NF_PUT_ATT_TEXT(idies,idvar,'units',lenchr(units),
     . units)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inivar'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if
      if (valmissing .ne. 0.) then
         ierr = NF_PUT_ATT_REAL(idies,idvar,'missing_value',NF_FLOAT,
     .    1,valmissing)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in inivar'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if

         ierr = NF_PUT_ATT_REAL(idies,idvar,'_FillValue',NF_FLOAT,
     .    1,valmissing)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in inivar'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
      end if

*     Exit define mode
*     ----------------
      ierr = NF_ENDDEF(idies)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in inivar'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

      return
      end

*.....................................................................
* endini - end initialization of file: sync file to disk, close file
*.....................................................................
      subroutine endini(idies, ierror)

*INPUT
*     idies - integer - id of netcdf file (generated by inifile or inifilec)
*
*OUTPUT
*     ierror - integer - error code = 0 if no error, < 0 if error occurred

      integer idies, ierror

      include 'netcdf.inc'

      integer ierr

*     Close file, sync to disk
*     ------------------------
      ierr = NF_CLOSE(idies)

      if (ierr .ne. NF_NOERR) then
         print *, 'Error in endini'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

      return
      end

*.....................................................................
* writevar - write a 2-, 3-, or 4-dimensional hyperslab of data
*.....................................................................
      subroutine writevar(filen,varname,istart,icount,values,
     . times,tweights,dates,ierror)

*INPUT
*     filen - character*(*) - file name to which to write data
*     varname - character*(*) - name of variable to which to write
*     istart - integer(*) - starting points along each dimension of
*      the variable, for example the 1st point a 4-d variable is (1,1,1,1) 
*      and the 3rd level and 2nd time step of a 4d variable is (1,1,3,2).
*     icount - integer(*) - number of points along each dimension to write,
*      for example, to write a single lat/lon grid to a 4-d variable,
*      icount would be (nlon,nlat,1,1).
*     values - real(*) - real values of the designated hyperslab
*     times - real(*) - time vector vector for those points written (ignored
*      if variable is 2-d).
*     tweights - real(*) - number of days in each time step in times (ignored
*      if variable is 2-d).
*     dates - character*10(*) - character labels for each time step in times,
*      last character should be a null (dates ignored if variable 2-d).
*OUTPUT
*     ierror - integer - error code = 0 if no error, < 0 if error occurred

* Comment by Cicada Dennis: It seems like this subroutine could benefit from
* the same type of optimization as was implemented in readvar(), to not
* close the file at the end of the routine, but would probably need to 
* keep the list of files opened for writing separate from files opened
* for reading.
! S. Slavin, 06 May 2016 - Include rdopt module with additional debugging
!                          functions, variables, etc.
      use rdopt
      use compar,only: nc_file_bufrsize_KB,size_1KB,flg_detail

      character*(*) filen, varname
      character*10 dates(*)
      integer istart(*), icount(*), ierror
      real values(*), times(*), tweights(*)

      include 'netcdf.inc'

      integer idies, idvar, ndims, ierr, itype, idtime, idtw, iddate
      integer ist(2), ict(2)

      ierror = 0

*     Open file, but don't fail if already open
*     -----------------------------------------
      nc_file_bufrsize = nc_file_bufrsize_KB*size_1KB

! S. Slavin changed NF__OPEN() call to use a certain buffer size.
      write(*,*) 'From subroutine writevar():'
      write(*,*)
     &     'Call to NF__OPEN() with nc_file_bufrsize = ',
     &     nc_file_bufrsize
      write(*,'(a)') 'About to open file: ['//trim(filen)//']'
      ierr = NF__OPEN(filen,NF_WRITE,nc_file_bufrsize,idies)
! Cicada Dennis changed open call back to original to see if this fixes speed problem.
! but it did not seem to help. MMAP files seemed to be the problem.
!      ierr = NF_OPEN(filen,NF_WRITE,idies)
      write(*,'(a)') 'Opened file: ['//trim(filen)//
     &               '] with mode NF_WRITE'
      write(*,*) 'Call returned with nc_file_bufrsize = ',
     &           nc_file_bufrsize

      if (ierr .ne. NF_NOERR) then
         print *, 'Warning in writevar'
         print *, NF_STRERROR(ierr)
      end if

*     Get id of variable
*     ------------------
      ierr = NF_INQ_VARID(idies,varname,idvar)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in writevar'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Inquire number of dimensions
*     ----------------------------
      ierr = NF_INQ_VARNDIMS(idies,idvar,ndims)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in writevar'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Put variable
*     ------------
      ierr = NF_PUT_VARA_REAL(idies,idvar,istart,icount,values)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in writevar'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if

*     Put time value(s), time weights, dates if necesary
*     --------------------------------------------------
      if (ndims .gt. 2) then

         ierr = NF_INQ_VARID(idies,'time',idtime)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in writevar'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
         ierr = NF_PUT_VARA_REAL(idies,idtime,istart(ndims),
     .    icount(ndims),times)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in writevar'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
         ierr = NF_INQ_VARID(idies,'time_weights',idtw)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in writevar'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
         ierr = NF_PUT_VARA_REAL(idies,idtw,istart(ndims),
     .    icount(ndims),tweights)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in writevar'
            print *, NF_STRERROR(ierr)
            ierror = -1
            return
         end if
!         ierr = NF_INQ_VARID(idies,'date',iddate)
!         if (ierr .ne. NF_NOERR) then
!            print *, 'Error in writevar'
!            print *, NF_STRERROR(ierr)
!            ierror = -1
!            return
!         end if
!         ist(1)=1
!         ist(2)=istart(ndims)
!         ict(1)=10
!         ict(2)=icount(ndims)
!         ierr = NF_PUT_VARA_TEXT(idies,iddate,ist,
!     .    ict,dates)
!         if (ierr .ne. NF_NOERR) then
!            print *, 'Error in writevar'
!            print *, NF_STRERROR(ierr)
!            ierror = -1
!            return
!         end if

      end if

*     Close file
*     ----------
      ierr = NF_CLOSE(idies)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in writevar'
         print *, NF_STRERROR(ierr)
         ierror = -1
         return
      end if


      return
      end
*.....................................................................
* readvar - read a 2-, 3-, or 4-dimensional hyperslab of data
*.....................................................................

      subroutine readvar(filen,varname,name3d,istart,icount,values,
     . alons,alats,vals3d,times,ierror)

*INPUT
*     filen - character*(*) - file name from which to read data
*     varname - character*(*) - name of variable from which to read
*     name3d - character*(*) - name of 3rd, nontime dimension variable.
*      Ignored if varname variable is only 3-d.
*     istart - integer(*) - starting points along each dimension of
*      the variable, for example the 1st point a 4-d variable is (1,1,1,1) 
*      and the 3rd level and 2nd time step of a 4d variable is (1,1,3,2).
*     icount - integer(*) - number of points along each dimension to read,
*      for example, to read in a single lat/lon grid from a 4-d variable,
*      icount would be (nlon,nlat,1,1)
*OUTPUT
*     values - real(*) - returned real values of the designated hyperslab
*     alons - real(*) - longitude vector for those points read in
*     alats - real(*) - latitude vector for those points read in
*     vals3d - real(*) - vector of values of the 3rd dimension (unchanged if
*      variable is 2- or 3-d, or if 3rd dimension has character values).
*     times - real(*) - time vector vector for those points read in (unchanged
*      if variable is 2-d).
*     ierror - integer - error code = 0 if no error, < 0 if error occurred
  
! S. Slavin, 06 May 2016 - Include rdopt module with additional debugging
!                          functions, variables, etc.
      use rdopt
      use compar,only: nc_file_bufrsize_KB,size_1KB,flg_detail,flg_rdopt

      character*(*) filen, varname, name3d
      integer istart(*), icount(*), ierror
      real values(*), alons(*), alats(*), vals3d(*), times(*)
      character*20:: varname_alt,name3d_alt

      include 'netcdf.inc'

      integer idies, idvar, ndims, idlon, idlat, id3d, ierr, itype
! S. Slavin, 04 May 2016 - Set variable for file buffer size
      integer nc_file_bufrsize
      integer(kind=4) :: t_s,t_e                !clock counter, to check the elapsed time for section of interest
      integer(kind=4) :: time_s,time_e          !clock counter, to check the elapsed time for the subroutine/function
      integer(kind=4) :: clock_rate, clock_max  !see system_clock() for reference

      ierror = 0

      call system_clock(time_s,clock_rate,clock_max)

*     Open file, but don't fail if already open
*     -----------------------------------------
      call system_clock(t_s,clock_rate,clock_max)

!--------------------------------------------------------------------
! S. Slavin - Implement read optimization for NC files if flg_rdopt is set.

      if ( flg_rdopt == 0 ) then
!         ierr = NF_OPEN(filen,NF_NOWRITE,idies)
! Cicada Dennis changed NF_OPEN back to NF__OPEN call that was used in
! a prev. version. This uses a setable buffer size.
! When we set the size to 4MB the code seems to work faster.
         nc_file_bufrsize = nc_file_bufrsize_KB*size_1KB
         ierr = NF__OPEN(filen,NF_NOWRITE,nc_file_bufrsize,idies)
         if (ierr .ne. NF_NOERR) then
            print *, 'Warning in readvar ', trim(filen)
            print *, NF_STRERROR(ierr)
         end if

         write(*,'(a,i9.9)') 'Opened file: ['//trim(filen)//
     &        '] with mode NF_NOWRITE and default buffer size'
      else
         nc_file_bufrsize = nc_file_bufrsize_KB*size_1KB

! Call to see if idies is returned with a non-zero (presumably valid value)

         idies = 0
         call ncid_from_fname( filen, idies )

         if ( idies == 0 ) then
! d            write(*,'(a,i9.9)') 'Try to open file: ['//trim(filen)
! d     &           //'] with NF_MMAP and rdopt_mmap_bufrsize = ',
! d     &           rdopt_mmap_bufrsize
!            ierr = NF__OPEN(filen, NF_MMAP, rdopt_mmap_bufrsize, idies)
! Cicada Dennis changed NF__OPEN to not use memory mapping.
! NF_MMAP caused the node to run out of memory when many ibis processes were running.
           write(*,'(a,i9.9)') 'Try to open file: ['//trim(filen)
     &           //'] with NF_NOWRITE and nc_file_bufrsize = ',
     &           nc_file_bufrsize
            ierr = NF__OPEN(filen,NF_NOWRITE,nc_file_bufrsize,idies)

            if ( ierr /= NF_NOERR ) then

               if ( ierr == NF_ENOMEM ) then
                  write(*,'(a)') 'Error opening with NF_MMAP: '
     &                 //NF_STRERROR(ierr)
                  write(*,'(a)') 'Try to open file: ['//trim(filen)
     &              //'] with NF_NOWRITE!'

                  ierr = NF__OPEN(filen,NF_NOWRITE,nc_file_bufrsize,idies)

                  if ( ierr == NF_NOERR ) then
                     write(*,'(a)') 'Opened file: ['//trim(filen)//
     &                    '] with mode NF_NOWRITE'
                     write(*,'(a,i9.9)')
     &                    'Call returned with nc_file_bufrsize = ',
     &                    nc_file_bufrsize
                  else
                     write(*,*) 'Error in readvar: '//NF_STRERROR(ierr)
                     stop
                  end if
              else
                  write(*,'(a,i9.9)') 'Opened file: ['//trim(filen)//
     &                 '] with mode NF_MMAP and buffer size = ',
     &                 rdopt_mmap_bufrsize
               end if                        
            end if

! Log the open file in the ncid_map table of open netCDF units/ncids.
            call log_open_ncid( idies )

            if (ierr .ne. NF_NOERR) then
               print *, 'Warning in readvar ', trim(filen)
               print *, NF_STRERROR(ierr)
            end if

        else     ! Else do whatever happens if file is already open...
            write(*,'(a,i9.9)') 'File already open: ['//trim(filen)
     &           //'] on unit ', idies

         end if    ! End "if ( idies == 0 )"
      end if       ! End "if ( flg_rdopt == 0 )"

!--------------------------------------------------------------------

      call system_clock(t_e,clock_rate,clock_max)
      if(flg_detail == 1) then
        write(6,'(A,A,g16.8,A)') '(readvar)time to open file: ',
     &        trim(filen), real(t_e-t_s)/(real(clock_rate)),' seconds'
      end if

*     Get id of variable
*     ------------------
      ierr = NF_INQ_VARID(idies,varname,idvar)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in readvar ', trim(varname)
         print *, 'in file ', trim(filen)
         print *, NF_STRERROR(ierr)

         !some fix for issue of variable name
         if (trim(varname) == 'domtext') then !the name may be "soiltype"
           varname_alt = 'soiltype'
           write(6,*) 'Try variable name ',trim(varname_alt)
           ierr = NF_INQ_VARID(idies,varname_alt,idvar)

           write(6,*) NF_STRERROR(ierr)
           if(ierr .ne. NF_NOERR) then
             ierror = -1
             return
           end if
         end if !on if(trim(varname))
           
      end if

*     Inquire number of dimensions
*     ----------------------------
      ierr = NF_INQ_VARNDIMS(idies,idvar,ndims)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in readvar'
         print *, NF_STRERROR(ierr)
         ierror = -2
         return
      end if

*     Get variable
*     ------------
      call system_clock(t_s,clock_rate,clock_max)
      ierr = NF_GET_VARA_REAL(idies,idvar,istart,icount,values)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in readvar'
         print *, NF_STRERROR(ierr)
         ierror = -3
         return
      end if
      call system_clock(t_e,clock_rate,clock_max)
      if(flg_detail == 1) write(6,'(A,g16.8,A)') '(readvar)time to read var data: ',real(t_e-t_s)/(real(clock_rate)),' seconds'

*     Get values of dimension variables
*     ---------------------------------
      ierr = NF_INQ_VARID(idies,'longitude',idlon)
      if (ierr .ne. NF_NOERR) then  ! If we can't find a 'longitude' variable, look for 'lon'
        ierr = NF_INQ_VARID(idies,'lon',idlon)
      end if
      if (ierr .ne. NF_NOERR) then  ! try 'x'
        ierr = NF_INQ_VARID(idies,'x',idlon)
      end if
      if (ierr .ne. NF_NOERR) then
        print *, 'Error in readvar'
        print *, "Can't find longitude or lon variables in file"
        print *, NF_STRERROR(ierr)
        ierror = -4
      end if
      ierr = NF_GET_VARA_REAL(idies,idlon,istart(1),
     . icount(1),alons)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in readvar'
         print *, NF_STRERROR(ierr)
         ierror = -4
         return
      end if
      ierr = NF_INQ_VARID(idies,'latitude',idlat)
      if (ierr .ne. NF_NOERR) then  !If we can't find a 'latitude' variable, look for 'lat'
        ierr = NF_INQ_VARID(idies,'lat',idlat)
      end if
      if (ierr .ne. NF_NOERR) then  !try 'y'
        ierr = NF_INQ_VARID(idies,'y',idlat)
      end if
      if (ierr .ne. NF_NOERR) then
        print *, 'Error in readvar'
        print *, "Can't find latitude or lat variables in file"
        print *, NF_STRERROR(ierr)
        ierror = -5
      end if
      ierr = NF_GET_VARA_REAL(idies,idlat,istart(2),
     . icount(2),alats)
      if (ierr .ne. NF_NOERR) then
         print *, 'Error in readvar'
         print *, NF_STRERROR(ierr)
         ierror = -5
         return
      end if
      if (ndims .ge. 3) then
         ierr = NF_INQ_VARID(idies,'time',idtime)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in readvar'
            print *, NF_STRERROR(ierr)
            ierror = -6
            return
         end if
         ierr = NF_GET_VARA_REAL(idies,idtime,istart(ndims),
     .    icount(ndims),times)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in readvar'
            print *, NF_STRERROR(ierr)
            ierror = -7
            return
         end if
      end if

      if (ndims .eq. 4) then
         ierr = NF_INQ_VARID(idies,name3d,id3d)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in readvar(3rd variable): ', trim(name3d)
            write(6,*) 'in file: ', trim(filen)
            print *, NF_STRERROR(ierr)

            name3d_alt ='layer'
            write(6,*) "Try variable name ",trim(name3d_alt)
            ierr = NF_INQ_VARID(idies,name3d_alt,id3d)
            print *, NF_STRERROR(ierr)

            if (ierr .ne. NF_NOERR) then
               
              name3d_alt ='level'
              write(6,*) "Try variable name ",trim(name3d_alt)
              ierr = NF_INQ_VARID(idies,name3d_alt,id3d)
              print *, NF_STRERROR(ierr)

              if(ierr .ne. NF_NOERR) then
                ierror = -8
                return
              end if
            end if
         end if

         ierr = NF_INQ_VARTYPE(idies,id3d,itype)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in readvar'
            print *, NF_STRERROR(ierr)
            ierror = -9
            return
         end if
         if (itype .ne. NF_CHAR) ierr = NF_GET_VARA_REAL(idies,id3d,
     .    istart(3),icount(3),vals3d)
         if (ierr .ne. NF_NOERR) then
            print *, 'Error in readvar'
            print *, NF_STRERROR(ierr)
            ierror = -10
            return
         end if
      end if

*     Close file
*     ----------
! If rdopt routines are not being used, close netCDF file.
      if ( flg_rdopt == 0 ) then
         ierr = NF_CLOSE(idies)
      end if

      if (ierr .ne. NF_NOERR) then
         print *, 'Error in readvar'
         print *, NF_STRERROR(ierr)
         ierror = -11
         return
      end if

      call system_clock(time_e,clock_rate,clock_max)
      if(flg_detail == 1) then
         write(6,'(A,g16.8,A)')
     &     '(readvar)total time used: ',
     &     real(time_e-time_s)/(real(clock_rate)),' seconds'
      end if

      return
      end


*.....................................................................
* dimlen - get the length of a dimension
*.....................................................................

      subroutine dimlen(filen, dimname, length, ierror)

*INPUT
*     filen - character*(*) - file name from which to read data
*     dimname - character*(*) - name of dimension whose length we want
* OUTPUT
*     length - length of dimension
*     ierror - integer - error code = 0 if no error, < 0 if error occurred

* Comment by Cicada Dennis: It seems like this subroutine could benefit from
* the same read optimization as was implemented in readvar().
      use compar,only: nc_file_bufrsize_KB,size_1KB,flg_detail

      character*(*) filen, dimname
      integer length, ierror

      include 'netcdf.inc'

      integer idies, dimid, ierr

      ierror = 0

*     Open file, but don't fail if already open
*     -----------------------------------------

      write(*,*) 'From subroutine dimlen():'
      write(*,'(a)') 'About to open file: ['//trim(filen)//'] as NF_NOWRITE'
      ierr = NF_OPEN(filen,NF_NOWRITE,idies)

      if (ierr .ne. NF_NOERR) then
        print *, 'Warning in dimlen'
        print *, 'while opening file, with file:'
        print *, trim(filen)
        print *, NF_STRERROR(ierr)
      end if

*     Get dimension ID
*     ----------------
      ierr = NF_INQ_DIMID(idies,dimname,dimid)
      if (ierr .ne. NF_NOERR) then
        print *, 'Error in dimlen'
        print *, 'while trying to get dimid, with file:'
        print *, trim(filen)
        print *, NF_STRERROR(ierr)
        ierror = -1
        return
      end if

*     Get dimension length
*     --------------------
      ierr = NF_INQ_DIMLEN(idies,dimid,length)
      if (ierr .ne. NF_NOERR) then
        print *, 'Error in dimlen'
        print *, 'while trying to get dimension length, with file:'
        print *, trim(filen)
        print *, NF_STRERROR(ierr)
        ierror = -2
        return
      end if

*     Close file
*     ----------
      ierr = NF_CLOSE(idies)
      if (ierr .ne. NF_NOERR) then
        print *, 'Error in dimlen'
        print *, 'while trying to close file, with file:'
        print *, trim(filen)
        print *, NF_STRERROR(ierr)
        ierror = -3
        return
      end if

      return
      end      
