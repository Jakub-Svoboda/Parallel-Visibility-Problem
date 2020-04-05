#!/bin/bash

param=$1;


#preklad cpp zdrojaku
mpic++ --prefix /usr/local/share/OpenMPI -o vid.exe vid.cpp


#spusteni
mpirun --prefix /usr/local/share/OpenMPI -np 10 vid.exe $param

