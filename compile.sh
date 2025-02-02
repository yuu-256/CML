#!/bin/zsh
gfortran -fcheck=all code/main.f90
./a.out
python3 utils/plot_2D.py
