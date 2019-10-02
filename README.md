# AlignmentRepaC
C++ implementation of practicable Aligned Induction 

Windows debug -
```sh
cd /d C:\zzz\caiks\AlignmentC-master

cl -I../rapidjson-master/include /EHsc /DEBUG /Zi /c AlignmentUtil.cpp Alignment.cpp AlignmentApprox.cpp AlignmentAeson.cpp 

cd /d C:\zzz\caiks\AlignmentRepaC-master

cl -IC:../rapidjson-master/include -I../AlignmentC-master /EHsc /DEBUG /Zi main.cpp AlignmentRepa.cpp ../AlignmentC-master/AlignmentUtil.obj ../AlignmentC-master/Alignment.obj ../AlignmentC-master/AlignmentApprox.obj ../AlignmentC-master/AlignmentAeson.obj 

main
```
Windows release -
```sh
cd /d C:\zzz\caiks\AlignmentC-master

cl -I../rapidjson-master/include /EHsc /O2 /c AlignmentUtil.cpp Alignment.cpp AlignmentApprox.cpp AlignmentAeson.cpp 

cd /d C:\zzz\caiks\AlignmentRepaC-master

cl -IC:../rapidjson-master/include -I../AlignmentC-master /EHsc /O2 main.cpp AlignmentRepa.cpp ../AlignmentC-master/AlignmentUtil.obj ../AlignmentC-master/Alignment.obj ../AlignmentC-master/AlignmentApprox.obj ../AlignmentC-master/AlignmentAeson.obj 

main
```
