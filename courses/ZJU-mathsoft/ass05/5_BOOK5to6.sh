#!/bin/bash

mkdir -p "samples"

awk '
/BOOK FIVE/ && !in_block {
    book_count++
}

book_count == 2 && /CHAPTER I/ {
    in_block = 1
}

book_count == 2 && /CHAPTER II/ {
    in_block = 0
    exit
}

in_block {
    print
}
' "war_and_peace.txt" > "samples/book5.txt"