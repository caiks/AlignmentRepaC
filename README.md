# AlignmentRepaC

The AlignmentRepaC repository is a fast C++ implementation of some of the *practicable inducers* described in the paper *The Theory and Practice of Induction by Alignment* at https://greenlake.co.uk/. The AlignmentRepaC repository depends on the [AlignmentC repository](https://github.com/caiks/AlignmentC) for the underlying *model* framework. 

## Installation

The `AlignmentRepaC` module requires [modern C++](https://en.cppreference.com/w/) version 17 or later to be installed.

For example, in Ubuntu bionic (18.04),
```
sudo apt-get update -y && sudo apt install -y git g++ cmake

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
mkdir AlignmentC_build AlignmentRepaC_build
cd AlignmentRepaC_build
cmake -DCMAKE_BUILD_TYPE=DEBUG ../AlignmentRepaC
make

```
Ubuntu release -
```sh
mkdir AlignmentC_build AlignmentRepaC_build
cd AlignmentRepaC_build
cmake -DCMAKE_BUILD_TYPE=RELEASE ../AlignmentRepaC
make

```
Windows debug -
```sh
mkdir AlignmentC_build AlignmentRepaC_build
cd /d AlignmentRepaC_build
"C:\Program Files\CMake\bin\cmake" -G "Visual Studio 14 2015" -A x64 ../AlignmentRepaC
"C:\Program Files\CMake\bin\cmake" --build . --config Debug --target AlignmentRepaC_test

```
Windows release -
```sh
mkdir AlignmentC_build AlignmentRepaC_build
cd /d AlignmentRepaC_build
"C:\Program Files\CMake\bin\cmake" -G "Visual Studio 14 2015" -A x64 ../AlignmentRepaC
"C:\Program Files\CMake\bin\cmake" --build . --config Release --target AlignmentRepaC_test

```

## Usage

Ubuntu -
```sh
cd ..
mkdir AlignmentRepaC_run
cd AlignmentRepaC_run
../AlignmentRepaC_build/AlignmentRepaC_test 

```
Windows debug -
```sh
cd ..
mkdir AlignmentRepaC_run
cd AlignmentRepaC_run
..\AlignmentRepaC_build\Debug\AlignmentRepaC_test.exe

```
Windows release -
```sh
cd ..
mkdir AlignmentRepaC_run
cd AlignmentRepaC_run
..\AlignmentRepaC_build\Release\AlignmentRepaC_test.exe 

```
