#!/bin/bash
for i in {1..10000}
do
    ./build/main -v 1 -m assets/random-32-32-10.map -t 2 -N 30 -S tree.trace.yaml -s 0   >> result.txt
done
