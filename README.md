# AlignmentRepaC

The AlignmentRepaC repository is a fast C++ implementation of some of the *practicable inducers* described in the paper *The Theory and Practice of Induction by Alignment* at https://greenlake.co.uk/. The AlignmentRepaC repository depends on the [AlignmentC repository](https://github.com/caiks/AlignmentC) for the underlying *model* framework. 

## Installation

The `AlignmentRepaC` module requires [modern C++](https://en.cppreference.com/w/) version 17 or later to be installed.

For example in Ubuntu,
```
sudo apt-get update -y
sudo apt install -y g++
```
Then download the zip file or use git to get the underlying rapidjson and AlignmentC repositories, and the AlignmentRepaC repository -
```
git clone https://github.com/Tencent/rapidjson.git
git clone https://github.com/caiks/AlignmentC.git
git clone https://github.com/caiks/AlignmentRepaC.git
```

## Build

Ubuntu debug -
```sh
cd AlignmentC

g++ -I../rapidjson/include -std=gnu++17 -g -c AlignmentUtil.cpp Alignment.cpp AlignmentApprox.cpp AlignmentAeson.cpp 

cd ../AlignmentRepaC

g++ -I../rapidjson/include -I../AlignmentC -std=gnu++17 -g -o main main.cpp AlignmentRepa.cpp AlignmentAesonRepa.cpp AlignmentRandomRepa.cpp AlignmentPracticableRepa.cpp AlignmentPracticableIORepa.cpp ../AlignmentC/AlignmentUtil.o ../AlignmentC/Alignment.o ../AlignmentC/AlignmentApprox.o ../AlignmentC/AlignmentAeson.o

./main

```
Ubuntu release -
```sh
cd AlignmentC

g++ -I../rapidjson/include -std=gnu++17 -O3 -c AlignmentUtil.cpp Alignment.cpp AlignmentApprox.cpp AlignmentAeson.cpp 

cd ../AlignmentRepaC

g++ -I../rapidjson/include -I../AlignmentC -std=gnu++17 -O3 -o main main.cpp AlignmentRepa.cpp AlignmentAesonRepa.cpp AlignmentRandomRepa.cpp AlignmentPracticableRepa.cpp AlignmentPracticableIORepa.cpp ../AlignmentC/AlignmentUtil.o ../AlignmentC/Alignment.o ../AlignmentC/AlignmentApprox.o ../AlignmentC/AlignmentAeson.o

./main

```
Windows debug -
```sh
cd AlignmentC-master

cl -I../rapidjson-master/include /EHsc /DEBUG /Zi /c AlignmentUtil.cpp Alignment.cpp AlignmentApprox.cpp AlignmentAeson.cpp 

cd ..\AlignmentRepaC-master

cl -IC:../rapidjson-master/include -I../AlignmentC-master /EHsc /DEBUG /Zi main.cpp AlignmentRepa.cpp AlignmentAesonRepa.cpp AlignmentRandomRepa.cpp AlignmentPracticableRepa.cpp AlignmentPracticableIORepa.cpp ../AlignmentC-master/AlignmentUtil.obj ../AlignmentC-master/Alignment.obj ../AlignmentC-master/AlignmentApprox.obj ../AlignmentC-master/AlignmentAeson.obj 

main
```
Windows release -
```sh
cd AlignmentC-master

cl -I../rapidjson-master/include /EHsc /O2 /c AlignmentUtil.cpp Alignment.cpp AlignmentApprox.cpp AlignmentAeson.cpp 

cd ..\AlignmentRepaC-master

cl -IC:../rapidjson-master/include -I../AlignmentC-master /EHsc /O2 main.cpp AlignmentRepa.cpp AlignmentAesonRepa.cpp AlignmentRandomRepa.cpp AlignmentPracticableRepa.cpp AlignmentPracticableIORepa.cpp ../AlignmentC-master/AlignmentUtil.obj ../AlignmentC-master/Alignment.obj ../AlignmentC-master/AlignmentApprox.obj ../AlignmentC-master/AlignmentAeson.obj 

main
```
