# AlignmentRepaC

C++ implementation of practicable Aligned Induction 

Windows debug -
```sh
cd /d C:\zzz\caiks\AlignmentC-master

cl -I../rapidjson-master/include /EHsc /DEBUG /Zi /c AlignmentUtil.cpp Alignment.cpp AlignmentApprox.cpp AlignmentAeson.cpp 

cd /d C:\zzz\caiks\AlignmentRepaC-master

cl -IC:../rapidjson-master/include -I../AlignmentC-master /EHsc /DEBUG /Zi main.cpp AlignmentRepa.cpp AlignmentAesonRepa.cpp AlignmentRandomRepa.cpp AlignmentPracticableRepa.cpp AlignmentPracticableIORepa.cpp ../AlignmentC-master/AlignmentUtil.obj ../AlignmentC-master/Alignment.obj ../AlignmentC-master/AlignmentApprox.obj ../AlignmentC-master/AlignmentAeson.obj 

main
```
Windows release -
```sh
cd /d C:\zzz\caiks\AlignmentC-master

cl -I../rapidjson-master/include /EHsc /O2 /c AlignmentUtil.cpp Alignment.cpp AlignmentApprox.cpp AlignmentAeson.cpp 

cd /d C:\zzz\caiks\AlignmentRepaC-master

cl -IC:../rapidjson-master/include -I../AlignmentC-master /EHsc /O2 main.cpp AlignmentRepa.cpp AlignmentAesonRepa.cpp AlignmentRandomRepa.cpp AlignmentPracticableRepa.cpp AlignmentPracticableIORepa.cpp ../AlignmentC-master/AlignmentUtil.obj ../AlignmentC-master/Alignment.obj ../AlignmentC-master/AlignmentApprox.obj ../AlignmentC-master/AlignmentAeson.obj 

main
```
```
Ubuntu debug -
```sh
git clone https://github.com/Tencent/rapidjson.git

cd AlignmentC

g++ -I../rapidjson/include -std=gnu++17 -g -c AlignmentUtil.cpp Alignment.cpp AlignmentApprox.cpp AlignmentAeson.cpp 

cd ../AlignmentRepaC

g++ -I../rapidjson/include -I../AlignmentC -std=gnu++17 -g -o main main.cpp AlignmentRepa.cpp AlignmentAesonRepa.cpp AlignmentRandomRepa.cpp AlignmentPracticableRepa.cpp AlignmentPracticableIORepa.cpp ../AlignmentC/AlignmentUtil.o ../AlignmentC/Alignment.o ../AlignmentC/AlignmentApprox.o ../AlignmentC/AlignmentAeson.o

./main

```
Ubuntu release -
```sh
cd /home/cliff/Documents/projects/CAIKS4
git clone https://github.com/Tencent/rapidjson.git

cd AlignmentC

g++ -I../rapidjson/include -std=gnu++17 -O3 -c AlignmentUtil.cpp Alignment.cpp AlignmentApprox.cpp AlignmentAeson.cpp 

cd ../AlignmentRepaC

g++ -I../rapidjson/include -I../AlignmentC -std=gnu++17 -O3 -o main main.cpp AlignmentRepa.cpp AlignmentAesonRepa.cpp AlignmentRandomRepa.cpp AlignmentPracticableRepa.cpp AlignmentPracticableIORepa.cpp ../AlignmentC/AlignmentUtil.o ../AlignmentC/Alignment.o ../AlignmentC/AlignmentApprox.o ../AlignmentC/AlignmentAeson.o

./main


```
