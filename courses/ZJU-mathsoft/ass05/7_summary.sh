#!/bin/bash

keyword="$1"
grep -i "\b$keyword\b" "war_and_peace.txt" > temp.tmp
counts=$(wc -l < temp.tmp)
echo "$keyword:$counts"
head -n 5 temp.tmp
rm -f temp.tmp