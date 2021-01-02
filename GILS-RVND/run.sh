#!/bin/bash

larger=0
smaler=0

usage()
{
    echo "Usage:"
    echo "    ./run.sh -l	    Run the benchmark test just on the LARGERS instances (size >= 1500)."
    echo "    ./run.sh -s           Run the benchmark test just on the SMALLERS instances (size < 1500)."
}

larger=0
smaller=0

while getopts ":ls" opt; do
    case ${opt} in
	l)
	    $larger=1
	    ;;
	s)
	    $smaller=1
	    ;;
	*)
	    usage
	    exit 1
	    ;;
    esac
done

if [[ $large -eq 0  && $smaler -eq 0 ]]
then
    usage
    exit 1
fi

shift $((OPTIND-1))
exit 0

echo "--TSP Benchmark--"
FILENAME="./benchmark/results/bm-$(date +"%FT%T").txt"

make

k=1
for instance in instances/*; do
	echo $instance >> ${FILENAME}

	echo "Running $instance"
	echo "Instance $k of 67" 

	for i in {1..1}; do
		./tsp ${instance} | grep 'COST\|TIME' | awk "{print $1}" >> ${FILENAME}
	done

	k=$(($k + 1))
	#if (("$k" >"1")); then
		#break
	#fi
done

echo "-" >> ${FILENAME}

if [ `stat -c %A benchmark/summarycount.py | sed 's/...\(.\).\+/\1/'` != "x" ]; then
  chmod u+x benchmark/summarycount.py
fi

if [ `stat -c %A benchmark/bm.py | sed 's/...\(.\).\+/\1/'` != "x" ]; then
  chmod u+x benchmark/bm.py
fi

echo "Running bm.py to compute averages..."
./benchmark/bm.py

echo "Finishing up summary..."
./benchmark/summarycount.py

echo "Benchmark completed."
