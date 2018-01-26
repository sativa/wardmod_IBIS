#!/bin/bash
# Written by Cicada Dennis
# 2016/05/25
#
# $1 is intended to be a listing of all subdirectories of the local process's output/ directory.
# (eg: myprocDir/output/*  - /tmp/mpi/proc#/output/*)
# $2 is intended to be the top level directory being used in shared memory (eg: /tmp/mpi)
# $3 is intended to be the local process's work directory (eg: myprocDir - /tmp/mpi/proc#)
echo ""
echo "====================================="
echo "*****    File System Report     *****"
echo ""
echo "ls -l $1:"
ls -l $1
echo ""
echo "df:"
df
echo ""
echo "du -d 0 $2:"
du -d 0 $2
echo ""
echo "du $3:"
du $3
echo ""
echo "====================================="
echo "***** End of File System Report *****"
echo ""
