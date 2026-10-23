#!/bin/bash
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :
#
set -o nounset
set -o errexit
#command || { echo "command failed"; exit 1; }

if [ $# -lt 1 ]
then
	echo "Usage: $0 <PID> [time.interval=10s]"
	exit 2
fi > /dev/stderr

PID=$1
waittime=${2:-10s}

while true
do
	MEM=`pmap -x $PID | tail -n 1`
	if [ "$MEM" == "" ]
	then
		exit	## process does not exist, or error occurs
	else
		echo $MEM
		sleep $waittime
	fi
done

