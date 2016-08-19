#include "ErrorFinderManager.hpp"
#include <cstring>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#define GetCurrentDir getcwd




ErrorFinderManager::ErrorFinderManager():WINDOW(50),
                    MIN_SNP(30),GAP(1), MA_SNP_END(50), TRUESNP( 600 ),
                    MIN_CM(0.4), MA_ERR_THRESHOLD_START(0.08),
                    MA_ERR_THRESHOLD_END(0.08),PCT_ERR_THRESHOLD( 0.90 ),
                    HO_THRESHOLD( 0.98 ), TRUECM( 6 ),PIELENGTH( 3 ),
                    ISMOL( false), COUNTGAPERR( false ),MA_THRESHOLD(0.8),EMPIRICAL_MA_RESULT(-1.0), EMPIRICAL_PIE_RESULT(-1.0),EXTENDSNP(0)
{
}

void ErrorFinderManager::performConsolidation(int argc,char *argv[])
{
       bool goodParam=true;
       bool thresholdError = false; //if the user supplies both -empricial-ma-threshold and -ma-threshold, that is an error. This will be used to detect such an error.
       bool pieThresholdError = false;
       char currentPath[FILENAME_MAX];
       if(!GetCurrentDir(currentPath,sizeof(currentPath))){
	        cerr << "Error reading current directory" << endl;
		return;
       }
       for(int i=1;i<argc;i++)
        {

    	   /*Code for expanding window*/
    	   	    if (strcmp(argv[i],"-extendSNP")==0 && i<argc-1)
    	   	    {
    	   	    	EXTENDSNP=atoi(argv[++i]);
    	   	    	//cout<<"entered snp = "<<EXTENDSNP<<endl;
    	   	    }

    	   	    else if(strcmp(argv[i],"-bmatch")==0&&i<argc-1)
                {
                     BMATCHFILE=string(argv[++i]);
                }
                else if(strcmp(argv[i],"-bmid")==0&&i<argc-1)
                {
                     BMIDFILE=string(argv[++i]);
                }
                else  if(strcmp(argv[i],"-bsid")==0&&i<argc-1)
                {
                     BSIDFILE=string(argv[++i]);
                }
                 else if(strcmp(argv[i],"-reduced")==0&&i<argc-2)
                {
                     MIN_SNP=atoi(argv[++i]);
                     MIN_CM=atof(argv[++i]);
                }
                else  if(strcmp(argv[i],"-ped-file")==0&&i<argc-1)
                {
                     PEDFILE=string(argv[++i]);
                }
		            else  if(strcmp(argv[i],"-holdout-ped")==0&&i<argc-1)
                {
                     HPEDFILE=string(argv[++i]);
                }
		            else  if(strcmp(argv[i],"-holdout-map")==0&&i<argc-1)
                {
                     HMAPFILE=string(argv[++i]);
                }
                else  if(strcmp(argv[i],"-window")==0&&i<argc-1)
                {
                     WINDOW=atoi(argv[++i]);
                }
                else  if(strcmp(argv[i],"-holdout-threshold")==0&&i<argc-1)
                {
                     HO_THRESHOLD=atof(argv[++i]);
                }
                else  if(strcmp(argv[i],"-trueCM")==0&&i<argc-1)
                {
                     TRUECM=atof(argv[++i]);
                }
                else  if( strcmp( argv[i],"-trueSNP" )==0 && i < argc-1 )
                {
                     TRUESNP=atoi(argv[++i]);
                }
                else  if(strcmp(argv[i],"-holdout-missing")==0&&i<argc-1)
                {
                     HO_MISSING= string(argv[++i]);
                }
                else  if(strcmp(argv[i],"-gap")==0&&i<argc-1)
                {
                     GAP=atoi(argv[++i]);
                }
                else  if(strcmp(argv[i],"-ma-snp")==0&&i<argc-1)
                {
                     MA_SNP_END=atoi(argv[++i]);
                }
                else  if(strcmp(argv[i],"-pct-err-threshold")==0&&i<argc-1)
                {
                     if(pieThresholdError == true){
                      cerr << "ERROR: You have supplied both -emp-pie-threshold and -pct-err-threshold parameters, but only one is allowed. Please try again." << endl;
                      exit(1);
                     }
                     PCT_ERR_THRESHOLD=atof(argv[++i]);
                     pieThresholdError = true;
                }
                else if(strcmp(argv[i],"-emp-pie-threshold")==0&&i<argc-1)
                {
                     if(pieThresholdError == true){
                      cerr << "ERROR: You have supplied both -emp-pie-threshold and -pct-err-threshold parameters, but only one is allowed. Please try again." << endl;
                      exit(1);
                     }
                     EMPIRICAL_PIE_RESULT=atof(argv[++i]);
                     pieThresholdError = true;
                }
		            else  if(strcmp(argv[i],"-output-type")==0&&i<argc-1)
                {
                     OPTION=string(argv[++i]);
                     //cout<<OPTION<<endl;
                }
                else if(strcmp(argv[i], "-snpfile") ==0&&i<argc-1){
                     SNPWEIGHTFILE = string(argv[++i]);
                }
                else  if(strcmp(argv[i],"-log-file")==0&&i<argc-1)
                {
                     LOGFILE=string(argv[++i]);
                }
		else if(strcmp(argv[i],"-ma-threshold")==0&&i<argc-1)//adding new -ma-threshold argument
		{
		     if(thresholdError == true){ //user has already supplied an empirical-ma-threshold, so exit the program with an error message
			cerr << "ERROR: You have supplied both -empirical-ma-threshold and -ma-threshold parameters, but only one is allowed. Exiting program."<< endl;
                        exit(1);
		     }
		     MA_THRESHOLD=atof(argv[++i]);
		     thresholdError = true;
		}
		else if(strcmp(argv[i],"-empirical-ma-threshold")==0 && i<argc-1){
		     if(thresholdError == true){
			cerr << "ERROR: You have supplied both -empirical-ma-threshold and -ma-threshold parameters, but only one is allowed. Exiting program."<< endl;
			exit(1);
		     }
		     EMPIRICAL_MA_RESULT = atof(argv[++i]); //use the user supplied empirical ma threshold, instead of calculating it via true ibd segments
		     thresholdError = true;
		}
                else  if(strcmp(argv[i],"-PIE.dist.length")==0&&i<argc-1)
                {
                     string MOL=string(argv[++i]);
                     if( MOL.compare( "MOL" ) ==0 )
                     {
                        ISMOL = true;
                     }
                     else
                     {
                        PIELENGTH = atof( MOL.c_str() );//redundant string oeration
                     }
                }
                else  if(strcmp(argv[i],"-count.gap.errors")==0&&i<argc-1)
                {
                     string option=string(argv[++i]);
                     if( option.compare( "TRUE" ) ==0 )
                     {
                        COUNTGAPERR = true;
                     }
                }

                else
                {
                        wrongParam += " " + string(argv[i]);

                        goodParam=false;
                }
        }
        if((!goodParam)||BMATCHFILE.compare("")==0||BSIDFILE.compare("")==0||BMIDFILE.compare("")==0||PEDFILE.compare("")==0)
        {

                displayError( argv[0] );
                return;

        }
        if( OPTION.compare( "" ) == 0 )
        {
            cerr<< " please provide a valid output-type option " <<endl;
            exit( -1 );
        }
        if( LOGFILE.compare( "" ) == 0 )
        {
            cerr<< " default log file name is FISH " <<endl;
            LOGFILE = "FISH";
        }
        eCalculator.createLogFile( LOGFILE  );
        eCalculator.countGapErrors( COUNTGAPERR );
        time_t startTime;
       // cout<<"helloe5"<<endl;
        time (&startTime);
        string str_head = "****************************************************\n";
        str_head += "****************************************************\n";
        str_head += "FISHR LOG FILE INFORMATION\n\n";
        string str1 = " The program started at: " + string( ctime ( &startTime ) );
        string str = str_head + "Program working directory was: " + string(currentPath) +
		     " \nProgram version was: " + string(argv[0]) +
		     " \nProgram options:\n-bmatch file: " + BMATCHFILE +
                     " \n-bmid file: " + BMIDFILE +
                     " \n-bsid file: " + BSIDFILE +
                     " \n-ped file: " + PEDFILE;
		     if(HPEDFILE.compare( "" ) !=0 && HMAPFILE.compare( "" ) !=0){
		     str +=
                     " \n-holdout ped file: " + HPEDFILE +
                     " \n-holdoutmap file: " + HMAPFILE;
		     }
                     str += " \n-output type: " + OPTION +
                     " \n- missing SNP representation in pedfile: " + HO_MISSING +
                     " \n-log file: " + LOGFILE;
        str = str + " \nmin snp length : " + NumberToString( MIN_SNP )  +
                    " \nmin cm length : " + NumberToString( MIN_CM  ) +
                    " \ngap to consolidate : " + NumberToString( GAP ) +
                    " \nmoving averages window size : " + NumberToString(  WINDOW ) +
                    " \ndiscard ends to calculate pct err : " + NumberToString( MA_SNP_END );
                    if(EMPIRICAL_PIE_RESULT < 0.0){
                      str += " \nuser supplied percentage error threshold: " + NumberToString( PCT_ERR_THRESHOLD ) +
                      "\nuser did not supply an empirical-pie-threshold: NA";
                    } else {
                      str += "\nuser did not supply a percentage error threshold: NA";
                      str += "\nuser supplied empirical-pie-threshold: " + NumberToString( EMPIRICAL_PIE_RESULT );
                    }
                    if(EMPIRICAL_MA_RESULT < 0.0){
                      str += "\nuser supplied ma threshold: " + NumberToString(MA_THRESHOLD) +
                      "\nuser did not supply an empirical-ma-threshold: NA";
                    }else{
                      str += "\nuser did not supply an ma-threshold: NA";
                      str += "\nuser supplied empirical-ma-threshold: " + NumberToString(EMPIRICAL_MA_RESULT);
                    }
                    "\n****************************************************\n";

                  eCalculator.log( str );
        initiateErrorFinder(); //but bypass true calculations if empirical threshold is supplied
      if( HPEDFILE.compare( "" ) !=0 && HMAPFILE.compare( "" ) != 0 )
      {
           cerr<<"entering into HoldOut Mode Matching"<<endl;
/*           consolidator.performTrim(eCalculator,
                              WINDOW,MA_SNP_END,
                             MA_ERR_THRESHOLD_START,
                             MA_ERR_THRESHOLD_END,MIN_SNP,
                             MIN_CM,PCT_ERR_THRESHOLD, OPTION, HO_THRESHOLD, true  );*/
	   consolidator.performTrim(eCalculator, WINDOW, MA_SNP_END,MA_THRESHOLD,MIN_SNP,MIN_CM,PCT_ERR_THRESHOLD,OPTION,HO_THRESHOLD,true,EMPIRICAL_MA_RESULT,EMPIRICAL_PIE_RESULT,EXTENDSNP);//<piyush> added the param EXTENDSNP for calculating moving window
           cerr<<" Main Trim operation has completed "<< endl;
           cerr<< " Hold out trim has completed" <<endl;
           if( (OPTION.compare( "finalOutput" ) == 0) || (OPTION.compare( "Full" ) == 0 ) )
           {
                consolidator.finalOutPut( eCalculator, MIN_CM, MIN_SNP );
           }

      }
       else
       {
            if( OPTION.compare( "Error3") == 0 )
            {
               cerr<< " Error: You have provided option:Error3 " << endl
                   <<" you can use Error3 only if you provided"
                   <<" hold out ped and map file, program with not output anything" << endl;
              exit( -1 );
            }

           consolidator.performTrim(eCalculator, WINDOW, MA_SNP_END, MA_THRESHOLD, MIN_SNP, MIN_CM, PCT_ERR_THRESHOLD, OPTION, HO_THRESHOLD, false,EMPIRICAL_MA_RESULT,EMPIRICAL_PIE_RESULT,EXTENDSNP);//<piyush> added the param EXTENDSNP for calculating moving window
	   //Removing this option for now. Can add/remove later as need be
      /*     consolidator.performTrim(eCalculator,WINDOW,MA_SNP_END,
                                    MA_ERR_THRESHOLD_START,
                                    MA_ERR_THRESHOLD_END,MIN_SNP,MIN_CM,
                                    PCT_ERR_THRESHOLD, OPTION, HO_THRESHOLD, false );*/
          if( (OPTION.compare( "finalOutput" ) == 0) || (OPTION.compare( "Full" ) == 0 ) )
           {
		//cerr <<"DEBUG: ENTERING CONSOLIDATOR FINAL OUTPUT" << endl;
                consolidator.finalOutPut( eCalculator, MIN_CM, MIN_SNP );
		//cerr <<"DEBUG: EXITING CONSOLIDATOR FINAL OUTPUT" << endl;
           }
        }
        time_t endTime;
        time (&endTime);
        /*Provide the log file with information about start and end times*/
         string str2 = " The program ended at: " +
                         string( ctime ( &endTime ) );
         //remove trailing newline, for readability.
         str1.erase(std::remove(str1.begin(),str1.end(),'\n'),str1.end());
         str2.erase(std::remove(str2.begin(),str2.end(),'\n'),str2.end());

         str2 = str2 +  "           Total time ( in seconds): " +
                        NumberToString( ( endTime - startTime ) );

         //add start time
         str2 = "\n\n" + str1 + "        " + str2;
         //eCalculator.log( str1 );
          eCalculator.log( str2 );

}
void ErrorFinderManager::displayError(std::string argv)
{
                cerr<<"these parameters are not allowed "<<wrongParam<<endl;
                 cerr << "Usage: " << argv << " -bmatch [BMATCH FILE]  -bsid [BSID FILE] -bmid [BMID FILE] -reduced [min_snp] [min_cm] "
             <<" -ped-file [ped file] -window [window width to calculate moving averages] "
             <<" -gap [max gap to consolidate two matches]"
               <<" -pct-err-threshold [max percentage of errors in a match after the trim] OR -emp-pie-threshold" 
               <<" -ma-threshold [specifies percentile to be drawn from trulyIBD data for MA calculations] OR -empirical-ma-threshold"
               <<" Note that if both -emp-pie-threshold and empirical-ma-threshold are supplied, then -trueSNP and -trueCM will be ignored"
              <<"-output-type [ must provide any of these. it can be "
               << "MovingAverages  or Error1 or Error2 or Error3 or ErrorRandom1 " 
                << "or ErrorRandom2 or Error3 or ErrorRandom3 or Full "
                <<  "look at the description about how these works in wiki ]"
              << "(optional) -holdout-ped [new ped file path] -holdout-map [new map file] "
              << "-holdout-threshold [threshold to drop a match with new ped file ]"
              << " -holdout-missing [missing value representation in new ped file] "
              << " -log-file [log file name]"
              << " -trueCM [ true match maximum cm length] " 
              << " - trueSNP [ true match SNP length]"
              << " -PIE.dist.length [ can be MOL or any cm distance length "
              << "please refer wiki for more details on how to use this option"
              << "-count.gap.errors [ TRUE or FALSE to include gap errors in errors count ]"
            << endl;


}
void ErrorFinderManager::initiateErrorFinder()
{
        eCalculator.readBmidFile(BMIDFILE);
        cerr<<"Reading bmid file completed"<<endl;
        eCalculator.readBsidFile(BSIDFILE);
        cerr<<"Reading bsid file completed"<<endl;
        eCalculator.readPedFile(PEDFILE, HO_MISSING);
        cerr<<"Reading ped file completed"<<endl;
        int pers_count=eCalculator.getNoOfPersons();
        consolidator.readMatches(BMATCHFILE, pers_count, eCalculator, TRUESNP, TRUECM, EXTENDSNP,PEDFILE  );
        cerr<<"Reading bmatch file completed"<<endl;
        if( !(SNPWEIGHTFILE.empty())){ //verify that this check actually works
          consolidator.readUserSuppliedSnpWeights(SNPWEIGHTFILE);
        }
//testing to see what happens when consolidation is removed
        consolidator.performConsolidation(eCalculator,GAP,MIN_SNP,MIN_CM,EXTENDSNP);
        cerr<<"Consolidation completed"<<endl;
        if( HPEDFILE.compare( "" ) !=0 && HMAPFILE.compare( "" ) != 0 )
       {
        eCalculator.changeMapFile( HMAPFILE );
           eCalculator.readHPedFile( HPEDFILE, HO_MISSING );
           cerr<< " new map and ped File has read" <<endl;
        cerr<< " calculating true percentage errors" <<endl;
        if( ISMOL )
        {
          consolidator.findTruePctErrors( eCalculator, MA_SNP_END, true, WINDOW,MA_THRESHOLD, EMPIRICAL_MA_RESULT);

        }
        else
        {
          consolidator.findTrueSimplePctErrors( eCalculator, PIELENGTH, true, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT );
        }
         cerr<< " true percentage errors calculated "<<endl;
       }
       else
       {
          cerr<< " calculating true percentage errors" <<endl;
          if( ISMOL )
          {
            consolidator.findTruePctErrors( eCalculator, MA_SNP_END, false, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT );
          }
          else
          {
            consolidator.findTrueSimplePctErrors( eCalculator, PIELENGTH, false, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT );//<piyush>
          }
            cerr<< " true hold out percentage errors calculated "<<endl;
	
      }   
  
}
