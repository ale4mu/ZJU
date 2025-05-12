#!/bin/bash
  
tr '[:upper:]' '[:lower:]' < war_and_peace.txt | \
tr -cs 'a-z' '\n' | \
grep -v '^$' | \
sort | \
uniq -c | \
sort -nr | \
awk '{print $2 " " $1}' | \
head -n 10