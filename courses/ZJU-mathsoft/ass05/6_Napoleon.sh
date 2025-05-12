#!/bin/bash

grep "Napoleon" "war_and_peace.txt" > "napoleon_lines.txt"
count=$(wc -l < "napoleon_lines.txt")
echo "包含Napoleon的行数: $count" 