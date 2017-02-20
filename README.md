# FISHR
# This is the old version of FISHR. The new version is called FISHR2 located at 

- https://github.com/matthew-c-keller/FISHR2.git

FISHR2 includes significant changes and improvements from this verion of FISHR, including: 

1) Can be used directly on SHAPEIT formatted (HAPS/SAMPLE) data in addition to the "phased PED" file format used by GERMLINE and FISHR v.1. 

2) Running FISHR2 is now a single-step process. The user runs FISHR2 directly on the phased data rather than running GERMLINE2 first on the phased data and then FISHR on the GERMLINE output.

3) IBD2 and IBD4  shared segments (where 2 or 4 IBD segments exist at the same location between two individuals) can now be detected by FISHR2.

4) use "./FISHR2 -help" command to view the list of flags and their usuage. 



# Below is the intial README file for the FISHR program

FISHR detects IBD segments between individuals using genomewide SNP data.  
- FISHR needs the bmid, bsid, bmatch file from the GERMLINE2 to be able to run. 
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
