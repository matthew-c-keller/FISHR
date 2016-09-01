# FISHR

This is the intial README file for the FISHR program, which detects IBD segments between individuals. 

Note - FISHR needs the bmid, bsid, bmatch file from the germline to be able to run. 

To compile:
g++ -O2 ./ErrorFinder23.3/ErrorFinderMain.cpp -o FISHR

TO RUN:
./FISHR -ped-file ./Beagle.Phased.Group2.1k.ped  -bmatch GL_OUT.bmatch -bsid GL_OUT.bsid -bmid GL_OUT.bmid  -reduced 64  3 -window 50 -gap 100 -output-type finalOutput -count.gap.errors TRUE -extendSNP 40  -emp-pie-threshold 0.015  -ma-threshold 0.2  -log-file logs |gzip > FISHR_OUT.gz
