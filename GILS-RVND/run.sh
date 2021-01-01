#!/bin/bash

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
	if (("$k" >"1")); then
		break
	fi
done

echo "-" >> ${FILENAME}

if [ `stat -c %A benchmark/summarycount.py | sed 's/...\(.\).\+/\1/'` != "x" ]; then
  chmod u+x benchmark/summarycount.py
fi

echo "Running bm.py to compute averages..."

cd benchmark
python3 bm.py

echo "Finishing up summary..."
python3 summarycount.py
cd ..

echo "Benchmark completed."
