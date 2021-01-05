#!/bin/bash

usage()
{
    echo "Usage: Range from ARG1 to ARG2"
    echo " "
    echo "      ./run.sh -s [ARG1] -l [ARG2] "
    echo " "
    echo "      -l      The larger instance size. "
    echo "      -s      The smaller instance size. "
    echo " "
}

larger=0
smaller=0

while getopts "l: s:" opt; do
    case ${opt} in
	l)
	    larger=$OPTARG
	    ;;
	s)
	    smaller=$OPTARG
	    ;;
	*)
	    usage
	    exit 1
	    ;;
    esac
done

if [[ $larger -eq 0  && $smaller -eq 0 ]]
then
    usage
    exit 1
fi

shift $((OPTIND-1))

echo "--TSP Benchmark--"
FILENAME="./benchmark/results/bm-$(date +"%FT%T").txt"

make

k=1
files=$(ls instances/ | wc -l)

for instance in instances/*; do

    instance_sz=${instance//[!0-9]/}

    if [[ $instance_sz -gt $larger && $larger -ne 0 ]]
    then
        continue
    elif [[ $instance_sz -lt $smaller && $smaller -ne 0 ]]
    then
        continue
    fi

	echo $instance >> ${FILENAME}
	echo "Running $instance"
	echo "Instance $k of $files" 

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

cd benchmark/

echo "Running bm.py to compute averages..."
./bm.py

echo "Finishing up summary..."
./summarycount.py

echo "Benchmark completed."
