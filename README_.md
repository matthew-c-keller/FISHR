FISHR-Project

TO COMPILE: Navigate to ErrorFinder23.3
g++ -O2 -o FISHR ./ErrorFinderMain.cpp

TO RUN:
./FISHR -ped-file ./Beagle.Phased.Group2.1k.ped  -bmatch GL_OUT.bmatch -bsid GL_OUT.bsid -bmid GL_OUT.bmid  -reduced 64  3 -window 50 -gap 100 -output-type finalOutput -count.gap.errors TRUE   -emp-pie-threshold 0.015  -ma-threshold 0.2  -log-file logs |gzip > FISHR_OUT.gz



Updated Code for moving window averages.
FISHR is an open-source project being developed at the University of Colorado, Boulder that is used in genetic analysis. It is written in C++

