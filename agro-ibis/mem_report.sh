#!/bin/bash
echo
echo "********* node:$1 process:$2 Start Memory Report *********"
echo 
echo "Free command:"
echo "------------"
free
echo 
echo "Top command:"
echo "------------"
top -b -u hbrokaw -n 1 | head -n 30
echo 
echo "ps command:"
echo "------------"
ps -eo "%p %P %C %z %t %x %c %U %c" --sort "-vsz" | head -n 30
echo
echo "********* node:$1 process:$2 End of Memory Report *********"
echo
