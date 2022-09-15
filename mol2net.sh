#!/bin/bash
python 1_2_mol2net.py
for i in ./*; do mv "$i" "$(echo "$i" | tr -d "'")"; done
python 3_netAna.py

