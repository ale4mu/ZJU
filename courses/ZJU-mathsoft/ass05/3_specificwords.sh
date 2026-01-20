#!/bin/bash

wars=$(grep -o -i '\bwar\b' war_and_peace.txt | wc -l)
peaces=$(grep -o -i '\bpeace\b' war_and_peace.txt | wc -l)
echo "wars: $wars"
echo "peaces: $peaces"