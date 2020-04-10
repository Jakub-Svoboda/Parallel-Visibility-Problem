#!/bin/bash


param=$1;



function logTwo {
    local x=0
    for (( y=$1-1 ; $y > 0; y >>= 1 )) ; do
        let x=$x+1
    done
    echo $x
}

#preklad cpp zdrojaku
#mpic++ --prefix /usr/local/share/OpenMPI -o vid.exe vid.cpp &&

number=$( grep -o ',' <<< "$param" | grep -c , )  &&
number=$((number+1))     &&
result=$(logTwo $number)
result=$((2**result))

procNum=$((result/2))

#spusteni
mpirun --prefix /usr/local/share/OpenMPI -np $procNum vid.exe $param

