#!/bin/bash

read lines words characters <<< $(wc -lwm < war_and_peace.txt)
echo "Lines: $lines"
echo "Words: $words"
echo "Characters: $characters"