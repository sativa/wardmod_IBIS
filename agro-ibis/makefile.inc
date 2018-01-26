allocate.o: comatm.o

biogeochem.o: comatm.o combcs.o combgc.o comcrop.o comgrid.o comnitr.o compar.o
biogeochem.o: compft.o comsoi.o comsum.o comtex.o comveg.o

canopy.o: comsat.h
canopy.o: com1d.o comatm.o combcs.o comcrop.o comgrid.o comhyd.o compar.o
canopy.o: comsno.o comsoi.o comveg.o


climate.o: combcs.o comcrop.o comgrid.o comnitr.o compar.o compft.o comsum.o
climate.o: comveg.o comwork.o

com1d.o: comgrid.o

comatm.o: comgrid.o compar.o comwork.o

combcs.o: comgrid.o

combgc.o: comgrid.o

comcrop.o: comgrid.o compar.o




comhyd.o: comgrid.o compar.o



comnitr.o: comgrid.o compar.o

compar.o: comgrid.o

compft.o: compar.o

comsno.o: comgrid.o compar.o

comsoi.o: comgrid.o compar.o

comsum.o: comgrid.o compar.o


comveg.o: comgrid.o compar.o

comwork.o: comgrid.o compar.o

crops.o: com1d.o comatm.o combcs.o comcrop.o comgrid.o comnitr.o compar.o
crops.o: compft.o comsno.o comsoi.o comsum.o comveg.o comwork.o io-dyn.o

diagnostics.o: com1d.o comatm.o comdiag.o comgrid.o comhyd.o compar.o comsno.o
diagnostics.o: comsoi.o comsum.o comveg.o comwork.o


hourlymet.o: comsat.h
hourlymet.o: comatm.o combcs.o comcrop.o comgrid.o comhour.o comhyd.o comnitr.o
hourlymet.o: compar.o comsoi.o comsum.o comveg.o comwork.o

humid.o: comsat.h
humid.o: compar.o

ies-io.o: /opt/cray/netcdf/4.3.3.1/INTEL/140/include/netcdf.inc
ies-io.o: compar.o rdopt.o

initial.o: comsat.h
initial.o: comatm.o combcs.o comcrop.o comgrid.o comhour.o comhyd.o comnitr.o
initial.o: compar.o compft.o comsno.o comsoi.o comsum.o comtex.o comveg.o
initial.o: comwork.o

io-dyn.o: /opt/cray/netcdf/4.3.3.1/INTEL/140/include/netcdf.inc
io-dyn.o: clim_file_utils.o combcs.o comcrop.o comgrid.o cominput.o compar.o
io-dyn.o: comwork.o netcdf_wrappers.o

io.o: /opt/cray/netcdf/4.3.3.1/INTEL/140/include/netcdf.inc
io.o: clim_file_utils.o comatm.o combcs.o comcrop.o comgrid.o cominput.o
io.o: comnitr.o compar.o comsno.o comsoi.o comsum.o comveg.o comwork.o
io.o: netcdf_wrappers.o

io_hourly.o: com1d.o comatm.o combcs.o comcrop.o comgrid.o comhyd.o comnitr.o
io_hourly.o: compar.o comsoi.o comsum.o comveg.o comwork.o

main.o: com1d.o comatm.o combcs.o comcrop.o comdiag.o comgrid.o cominput.o
main.o: comnitr.o compar.o comsoi.o comsum.o comveg.o io-dyn.o rdopt.o

math.o: comgrid.o compar.o

netcdf_wrappers.o: comgrid.o compar.o comwork.o

physiology.o: comsat.h
physiology.o: com1d.o comatm.o combcs.o comcrop.o comgrid.o comnitr.o compar.o
physiology.o: compft.o comsno.o comsoi.o comsum.o comveg.o

radiation.o: com1d.o comatm.o combcs.o comgrid.o compar.o comsno.o comsoi.o
radiation.o: comveg.o

rdopt.o: /opt/cray/netcdf/4.3.3.1/INTEL/140/include/netcdf.inc

read_cfert.o: comcrop.o comgrid.o compar.o

read_input_description.o: clim_file_utils.o cominput.o compar.o

readpars.o: combgc.o comcrop.o comgrid.o comnitr.o compar.o compft.o comsoi.o
readpars.o: comtex.o comveg.o

snow.o: com1d.o combcs.o comgrid.o compar.o comsno.o comsoi.o comveg.o

soil.o: comsat.h
soil.o: com1d.o comatm.o combcs.o comgrid.o comhyd.o compar.o comsno.o comsoi.o
soil.o: comveg.o

solute.o: comatm.o combcs.o comcrop.o comgrid.o comnitr.o compar.o comsoi.o
solute.o: comsum.o comveg.o

stats.o: com1d.o comatm.o combcs.o comcrop.o comgrid.o comhyd.o comnitr.o
stats.o: compar.o comsno.o comsoi.o comsum.o comveg.o io-dyn.o

toa_radiation.o: comatm.o comgrid.o compar.o comwork.o

utilities.o: combcs.o comcrop.o comgrid.o comnitr.o compar.o comsoi.o comwork.o

vegetation.o: comatm.o combcs.o comcrop.o comgrid.o compar.o compft.o comsoi.o
vegetation.o: comsum.o comveg.o

weather.o: comsat.h
weather.o: comatm.o combcs.o comcrop.o comgrid.o comhour.o cominput.o comnitr.o
weather.o: compar.o comsoi.o comsum.o comveg.o comwork.o

daily_variables.mod: io-dyn.o
management_para.mod: io-dyn.o
