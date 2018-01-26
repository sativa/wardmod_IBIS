#!/bin/bash
#These will be archive.tar.gz
INPUTS='ibis ibis.infile yearsrun.dat  hist.fert.1945.1996 params.* restart'

#Run commands...
gnutar -zcf ibis_archive_sim17.tar.gz $INPUTS
