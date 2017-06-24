#!/bin/bash

version=2.1
RUNPATH=/home/lmingari/lidar_v$version
OUTPATH=$RUNPATH/Output

$RUNPATH/main.py Bariloche > $OUTPATH/Bariloche/last.log 2>&1
$RUNPATH/main.py Comodoro > $OUTPATH/Comodoro/last.log 2>&1
$RUNPATH/main.py Cordoba > $OUTPATH/Cordoba/last.log 2>&1
$RUNPATH/main.py Neuquen > $OUTPATH/Neuquen/last.log 2>&1
$RUNPATH/main.py Tucuman > $OUTPATH/Tucuman/last.log 2>&1
