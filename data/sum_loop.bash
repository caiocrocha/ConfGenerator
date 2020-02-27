#!/bin/bash
while getopts t:i: option
do
case "${option}"
in
t) TAM=${OPTARG};;
i) NUM_TIMES=${OPTARG};;
esac
done
T0=0; T1=0; T2=0
for (( I=0; I<$NUM_TIMES; I++ ))
do
    STR=$(python testes.py --tam $TAM);
    read -r AUX0 AUX1 AUX2 <<< $(echo $STR);
    T0=$(awk '{print $1+$2}' <<<"${T0} ${AUX0}");
    T1=$(awk '{print $1+$2}' <<<"${T1} ${AUX1}");
    T2=$(awk '{print $1+$2}' <<<"${T2} ${AUX2}");
done
echo "$T0"
echo "$T1"
echo "$T2"

