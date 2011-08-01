#!/bin/bash

CELLS=2000
GENS=1

RKLIM=50
RKSTEP=0.01
INTERVAL=1

MINRATE=0
MAXRATE=10
INITCONC=0


for HILL in {20,30,40}
do
	echo ""
	echo "Hill = $HILL ($CELLS cells)"
	echo "Begin time:`date`"
	echo "./EvoDevo --maxbasic 2 --maxptm 2 --maxprom 2 --cells $CELLS --gens $GENS --interval $INTERVAL --rklim $RKLIM --rkstep $RKSTEP --hill $HILL --minrate $MINRATE --maxrate $MAXRATE --csvData --csvCell --outputall >> ../output/scorehill$HILL"
	./EvoDevo --maxbasic 2 --maxptm 2 --maxprom 2 --cells $CELLS --gens $GENS --interval $INTERVAL --rklim $RKLIM --rkstep $RKSTEP --hill $HILL --minrate $MINRATE --maxrate $MAXRATE --csvData --csvCell --outputall >> ../output/scorehill$HILL &
	echo "End time: `date`"
done
