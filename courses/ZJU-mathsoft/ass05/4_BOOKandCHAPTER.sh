#!/bin/bash

books=$(grep -E 'BOOK [A-Z]+:' war_and_peace.txt | wc -l)
chapters=$(grep -E 'CHAPTER[[:space:]]+[IVXLCDM][IVXLCDM]*' war_and_peace.txt | wc -l)
echo "Total Books: $books"
echo "Total Chapters: $chapters"