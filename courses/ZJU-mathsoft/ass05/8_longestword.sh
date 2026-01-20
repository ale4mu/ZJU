#!/bin/bash

words=$(tr '[:upper:]' '[:lower:]' < "war_and_peace.txt" | \
grep -oE '\b[[:alpha:]]+\b')

longest_word=$(echo "$words" | \
awk '{ if (length($0) > max_length) { max_length = length($0); max_word = $0 } } END { print max_word }')

count=$(echo "$words" | grep -c "^$longest_word$")
echo "$longest_word:$count"