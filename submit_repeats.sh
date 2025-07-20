

#!/bin/bash

# Get the number of repeats from input argument
if [ -z "$1" ]; then
    echo "Usage: $0 <number_of_repeats>"
    exit 1
fi
NUM_REPEATS=$1

# Loop to submit jobs for each repeat
for i in $(seq 0 $((NUM_REPEATS - 1)))
do
    sbatch one_job.sh $i
    echo "Submitted repeat $i"
done
