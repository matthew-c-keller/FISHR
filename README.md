# FISHR
#This is the Old version of FISHR. The new version is called FISHR2 located at 

- https://github.com/matthew-c-keller/FISHR2.git



FISHR2 includes significant changes from this verion of FISHR.



This is the intial README file for the FISHR program, which detects IBD segments between individuals. 

Note 
- FISHR needs the bmid, bsid, bmatch file from the germline to be able to run. 
- The wrapper program calls either one of FISHR or FISHR_Low_Ram. So you need to compile  FISHR and FISHR_Low_Ram seperately(below).

To Compile FISHR: Navigate to ErrorFinder23.3 folder. Let the binary file be in the ErrorFinder23.3 folder (same folder).
g++ -O2 ./ErrorFinderMain.cpp -o FISHR

TO Compile FISHR_Low_ram: Navigate to ErrorFinder23.3_Low_Ram folder. Let the binary file be in the ErrorFinder23.3_Low_Ram folder (same folder).
g++ -O2 -o FISHR_Low_Ram ./ErrorFinderMain.cpp

To Compile Wrapper: Navigate to ErrorFinder23.3 folder. Move the binary file to  ../Wrapper (i.e. FISHR folder)
g++ -O2 -o FISHR wrapper_fishr.cpp 

TO RUN:
1. Normal versioon of FISHR
./FISHR -ped-file ./Beagle.Phased.Group2.1k.ped  -bmatch GL_OUT.bmatch -bsid GL_OUT.bsid -bmid GL_OUT.bmid  -reduced 64  3 -window 50 -gap 100 -output-type finalOutput -count.gap.errors TRUE  -emp-pie-threshold 0.015  -ma-threshold 0.2  -log-file logs |gzip > FISHR_OUT.gz

2. Low ram version of FISHR (Notice the -low_ram flag)
./FISHR -ped-file ./Beagle.Phased.Group2.1k.ped  -bmatch GL_OUT.bmatch -bsid GL_OUT.bsid -bmid GL_OUT.bmid  -reduced 64  3 -window 50 -gap 100 -output-type finalOutput -count.gap.errors TRUE  -emp-pie-threshold 0.015  -ma-threshold 0.2  -log-file logs -low_ram |gzip > FISHR_OUT.gz
