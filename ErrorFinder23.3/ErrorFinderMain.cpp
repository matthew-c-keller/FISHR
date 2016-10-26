//#include "ErrorFinderManager.cpp"
#include <iostream>
#include <string>
#include<cstring>
#include <vector>
#include "ErrorFinderManager.hpp"

//using namespace std;

int main(int argc,char *argv[])
{

		
         ErrorFinderManager manager;
         manager.performConsolidation(argc,argv);
         return 0;
}


/*
 -bmatch GL.HE.Omincm3.bit20.err0.bmatch -bsid GL.HE.Omincm3.bit20.err0.bsid -bmid GL.HE.Omincm3.bit20.err0.bmid -reduced 100 2.0 -ped-file Beagle.Phased.Group2.1k.ped -window 50 -ma-threshold 0.2 -gap 2 -pct-err-threshold 0.8 -trueCM 6 -trueSNP 500 -count.gap.errors TRUE -PIE.dist.length 3 -output.type finalErrorsOutput -log.file s8v2t1 -extendSNP 60
 *
 * */
