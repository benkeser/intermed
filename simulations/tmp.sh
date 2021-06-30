#! /bin/bash

total_jobs=9
max_jobs=3
analysis=tmp

# sub job that does merging
last_batch=$(
for number in $(seq $(($total_jobs - $max_jobs)) $total_jobs); do
    printf '%s\n' "${analysis}${number}"
done | paste -sd ',' -
)

echo $last_batch