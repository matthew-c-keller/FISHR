#include "Consolidator.hpp"
//#include "ErrorCalculator.cpp"
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <iostream>
#include <exception>
#include <algorithm>
#include <assert.h>
#include <sstream>
#include <iomanip>
using namespace std;
template <typename T>
  string NumberToString ( T Number )
  {
     ostringstream ss;
     ss << Number;
     return ss.str();
  }

bool Consolidator::compareFunction(SNP s1, SNP s2)
{

        return (s1.start<s2.start);

}

void Consolidator::sortMatches()
{
         for(int i=0;i<person_count;i++)
        {
                for(int j=i;j<person_count;j++)
                {


                        std::sort(m_matches[i][j].begin(),m_matches[i][j].end(),compareFunction);

                }
        }
}

void Consolidator::readMatches(string path,int pers_count, ErrorCalculator& eCalculator, int trueSNP, float trueCM )
{
     person_count=pers_count;
     if(pers_count<=0)
     {
       std::cerr<<"wrong BSID file, check it, reading ped file failed"<<std::endl;
       return;
     }
     try
     {
        person_count=pers_count;
        m_matches.resize(pers_count+1);
        m_trueMatches.resize( pers_count + 1 );
        for(int i=0;i<pers_count;++i)
        {
             m_matches[i].resize(pers_count+1);
             m_trueMatches[i].resize( pers_count + 1 );
        }
        unsigned int pid[2];
        unsigned int sid[2];
        unsigned int dif,hom[2];
        ifstream file_bmatch(path.c_str(),ios::binary);
        if( !file_bmatch )
         {
             cerr<<"unable to open the bmatch file, exiting the program" << endl;
             exit( -1 );

         }
   
        while ( !file_bmatch.eof())
        {
    pid[0] = -1;
                file_bmatch.read( (char*) &pid[0] , sizeof( unsigned int ) );
                if ( pid[0] == -1 ) continue;
                file_bmatch.read( (char*) &pid[1] , sizeof( unsigned int ) );
                file_bmatch.read( (char*) &sid[0] , sizeof( unsigned int ) );
                file_bmatch.read( (char*) &sid[1] , sizeof( unsigned int ) );
                file_bmatch.read( (char*) &dif , sizeof( int ) );
                file_bmatch.read( (char*) &hom[0] , sizeof( bool ) );
                file_bmatch.read( (char*) &hom[1] , sizeof( bool ) );
                 if(pid[0]>=pers_count||pid[1]>=pers_count)
                 {
                      cerr<<"problem with bsid file, check it please"<<endl;
                      return;

                 }
                 SNP snp;
                 snp.start=sid[0];
                 snp.end=sid[1];

                 if(pid[0]<=pid[1])
                       m_matches[(pid[0])][(pid[1])].push_back(snp);
                 else  
                       m_matches[(pid[1])][(pid[0])].push_back(snp);
                 if( ( eCalculator.getCMDistance( sid[ 1 ] ) - 
                                eCalculator.getCMDistance( sid[ 0 ] ) ) >= trueCM && 
                                 ( sid[ 1 ] - sid[ 0 ] ) >= trueSNP &&  pid[0] != pid[1] )
                 {
                     if(pid[0]<=pid[1])
                       m_trueMatches[(pid[0])][(pid[1])].push_back(snp);
                     else
                       m_trueMatches[(pid[1])][(pid[0])].push_back(snp);
   
                 }  

        }
        file_bmatch.close();
    }
    catch(exception &e)
    {
       cerr<<"Error:"<<e.what()<<endl;
       exit( -1 );
    }
        
}

void Consolidator::performConsolidation(ErrorCalculator& eCalculator, int gap,int min_snp,float min_cm)
{
         int consolidations = 0, removed = 0;
   for(int i=0;i<person_count;++i)//for each person
        {
                for(int j=i;j<person_count;++j)//compare with each other person
                {
                        int temp1=-1,temp2=-1;
                         for(int l=0;l<m_matches[i][j].size();++l)//for each match
                         {

                               temp1= m_matches[i][j][l].start;

                                temp2= m_matches[i][j][l].end;

        
                                if(temp2==-1||temp1==-1){continue;}

                                for(int k=l+1;k<m_matches[i][j].size();++k) //for each other match
                                {
                                    if((m_matches[i][j][k].start-temp2-1)<=gap)
                                    {
                                        ++consolidations;
                                        temp2=m_matches[i][j][k].end;
                                        m_matches[i][j][l].end=temp2;
                                        m_matches[i][j][k].end=-1;
                                    }
                                    else break;

                               }
            //this may be what is causing our initial drops to never show up...
            
                               if( ( (temp2-temp1)<min_snp) || ( (eCalculator.getCMDistance(temp2)-eCalculator.getCMDistance(temp1))<min_cm) )
                              {
               
                                       ++removed;
                          //             m_matches[i][j][l].end=m_matches[i][j][l].start=-1;
                              }
      }
    }
  }
        std::string str = " \n Number of Consolidations: " +
                            NumberToString( consolidations );
        str =  str + " \n Number of matches removed due to length: "
                     +  NumberToString( removed );
        eCalculator.log( str );
       
}
//overloaded version
void Consolidator::performTrim(ErrorCalculator& e_obj,int window,
                               int ma_snp_ends, float ma_threshold,
                               int min_snp,float min_cm,
                               float per_err_threshold, string option,
                               float hThreshold, bool holdOut,float empirical_threshold, string path, float empirical_pie_threshold)
{
  int removed1 =0, removed2 = 0, removed3 = 0, removed4 = 0;
  int not_removed = 0;
  int total_count = global_initial; //for keeping track of total #SH
  bool wrongOption = false;
  float per_err_threshold1;
  if(empirical_pie_threshold >= 0.0){
    per_err_threshold1 = empirical_pie_threshold;
  } else {
    per_err_threshold1 = getPctErrThreshold( per_err_threshold );
  }

  std::string str = " \n corresponding empirical error threshold for " + NumberToString( per_err_threshold ) + " is: " + NumberToString( per_err_threshold1 )  + " \n";
  float hThreshold1 = 0;
  if( holdOut )
  {
    hThreshold1 = getHoldOutThreshold( hThreshold );
    str = str + " \n corresponding holdout emperical error threshold for " + NumberToString( hThreshold ) + " is: " + NumberToString( hThreshold1 )  + " \n";
  }
  e_obj.log( str );
  per_err_threshold = per_err_threshold1;
  hThreshold = hThreshold1;

  //begin streaming in bmatch data
  unsigned int pid[2];
  unsigned int sid[2];
  unsigned int dif,hom[2];
  ifstream file_bmatch(path.c_str(),ios::binary);
  int pers_count = person_count;
  SNP snp;
  if( !file_bmatch )
  {
    cerr<<"unable to open the bmatch file, exiting the program" << endl;
    exit( -1 );
  }
  //will need to replace all instances of m_matches in this function with pid[0], pid[1], etc.
  //will probably encapsulate all of the below code into a while !eof loop
  while ( !file_bmatch.eof()){
    //read in bmatch data
    pid[0] = -1;
    file_bmatch.read( (char*) &pid[0] , sizeof( unsigned int ) );
    if ( pid[0] == -1 ) continue;
    file_bmatch.read( (char*) &pid[1] , sizeof( unsigned int ) );
    file_bmatch.read( (char*) &sid[0] , sizeof( unsigned int ) );
    file_bmatch.read( (char*) &sid[1] , sizeof( unsigned int ) );
    file_bmatch.read( (char*) &dif , sizeof( int ) );
    file_bmatch.read( (char*) &hom[0] , sizeof( bool ) );
    file_bmatch.read( (char*) &hom[1] , sizeof( bool ) );
    if(pid[0]>=pers_count||pid[1]>=pers_count)
    {
      cerr<<"problem with bsid file, check it please"<<endl;
      return;
    }

    snp.start=sid[0];
    snp.end=sid[1];
    if(pid[0] > pid[1]){
      int temp = pid[0];
      pid[0] = pid[1];
      pid[1] = temp;
    }

    int i = pid[0];
    int j = pid[1];

    //this next segment deals with truly ibd stuff...which I kind of want to ignore for now. It may be the case that IBD stuff needs to be accounted 
    //for back in the initiateErrorFinder() stage of things. Or we may simply need to store them in an actual vector...
    /*so for now...
    if( ( eCalculator.getCMDistance( sid[ 1 ] ) - 
    eCalculator.getCMDistance( sid[ 0 ] ) ) >= trueCM && 
    ( sid[ 1 ] - sid[ 0 ] ) >= trueSNP &&  pid[0] != pid[1] )
    {
    if(pid[0]<=pid[1])
    m_trueMatches[(pid[0])][(pid[1])].push_back(snp);
    else
    m_trueMatches[(pid[1])][(pid[0])].push_back(snp);

    }...*/
    //a bmatch line should have been read in by this point

    if(snp.end==-1)
    {
      continue; 
    }   
    int temp1=snp.start;
    int temp2=snp.end;
    int pers1 = pid[0], pers2 = pid[1];
    if( option.compare( "ErrorRandom1" ) == 0 || option.compare( "ErrorRandom2" ) == 0 || option.compare( "ErrorRandom3" ) == 0 )
    {
      pers1 = std::rand() % e_obj.getNoOfPersons();
      pers2 = std::rand() % e_obj.getNoOfPersons();
      if( pers1 > pers2 )
      {
        pers1 = pers1 + pers2;
        pers2 = pers1 - pers2;
        pers1 = pers1 - pers2;
      }
    }

    vector<vector<int> > errors=e_obj.checkErrors(pers1, pers2, temp1, temp2);
    vector<int>finalErrors=e_obj.getFinalErrors(errors);
    /*Inject implied error at start/end of SH here*/
    vector<int>::iterator it;
    it = finalErrors.begin(); //go to the start of the vector
    if(finalErrors[0] != 1){
      finalErrors.insert(it,1); //inject an error at position 1, if not already there
    }   
    /*End inject implied error section*/
    vector<int>trimPositions;
    vector<float>movingAverages;
    float threshold;


    //length 32 should survive this, since the reduced arg is 32 (although the reduced cM may be the culprit)

    if( (e_obj.isInitialCmDrop(temp1,temp2,min_cm)) || ((temp2-temp1) < min_snp) ){ //initial drop. Don't calculate MA
      trimPositions.push_back(temp1);
      trimPositions.push_back(temp2);
      trimPositions.push_back(1);
    }
    else
    {
        movingAverages = e_obj.getMovingAverages(finalErrors,temp1,temp2,window);
        if(empirical_threshold < 0.0){
          threshold = e_obj.getCutoff(); //this is the empirical average threshold for moving averages
        } else {
          threshold = empirical_threshold;
        }
        trimPositions = e_obj.getTrimPositions(movingAverages,temp1,temp2,threshold,min_cm); 
    }
    //-----------------
    int beforeTrimStart = temp1;        
    int beforeTrimEnd = temp2;
    snp.end = temp2 = temp1+trimPositions[1];
    snp.start = temp1 = temp1+trimPositions[0];
    int del0 = trimPositions[0];
    int del1 = trimPositions[1];
    //in order to update PIE calculation, may have to alter the finalErrors vector before passing it into getThreshold(). But the alteration
    //would need to occur at the start and end of the trimmed SH, not just at the start and end points of the initial SH.
    float per_err = e_obj.getThreshold(finalErrors,del0,del1,ma_snp_ends);
    //add new weighted option
    /*
     For this new option, we only output SH that are not dropped. So, the output is finalOutput + weighted column.
    */
    if( (option.compare("weightedOutput") == 0) ){
      int snp1 = 0, snp2 = 0, hlength = 0;
      float noOfOppHom = 0;
      /*if( holdOut )
      {
        snp1 = e_obj.getNewSnp( temp1 );
        snp2 = e_obj.getNewSnp( temp2 );
        hlength = snp2 - snp1;
        if( hlength <= 0 )
        {
          hlength = 1;
        }
        noOfOppHom = e_obj.getOppHomThreshold( pers1, pers2, m_matches[i][j][l].start, m_matches[i][j][l].end );
      }*/
      if( ( (beforeTrimEnd - beforeTrimStart) < min_snp) || ( (trimPositions.size() == 3) && (trimPositions[2] == 1) ) ){ 
        removed4++;
        continue;
      }
      if( (( temp2-temp1 ) < min_snp) || (trimPositions.size() == 3) ){ //removed2 a tpos.size of 3 indicates trimming due ot cM 
        removed2++;
        continue;
      }
      if( per_err > per_err_threshold){
        removed1++;
        continue;
      }
      if( holdOut && hThreshold < ( noOfOppHom ) / hlength ){
        removed3++;
        continue; 
      } //removed3

      not_removed++;
      m_weighted_sh.push_back(Weighted_SH(temp1,temp2,i,j)); //build the vector of SH that passed
      continue;
    }//end weghtedOutput
    if( (option.compare("FullPlusDropped") == 0) ){
      int snp1 = 0, snp2 = 0, hlength = 0;
      float noOfOppHom = 0;                
      if( holdOut )
      {
        snp1 = e_obj.getNewSnp( temp1 );
        snp2 = e_obj.getNewSnp( temp2 );
        hlength = snp2 - snp1;
        if( hlength <= 0 )
        {
          hlength = 1;
        }
        noOfOppHom = e_obj.getOppHomThreshold( pers1, pers2, snp.start, snp.end );
      }


      if( ( (beforeTrimEnd - beforeTrimStart) < min_snp) || ( (trimPositions.size() == 3) && (trimPositions[2] == 1) ) ){ 
        std::vector<float>movingAverages;
        temp1 = beforeTrimStart;
        temp2 = beforeTrimEnd;
        e_obj.fullPlusDroppedOutput(i,j,temp1,temp2,min_snp,min_cm,finalErrors,per_err,1);
        snp.start= snp.end=-1;
        removed4++;
        continue;
      }

      if( (( temp2-temp1 ) < min_snp) || (trimPositions.size() == 3) ){
        e_obj.fullPlusDroppedOutput(i,j,temp1,temp2,min_snp,min_cm,finalErrors,per_err,2); 
        removed2++;
        snp.start= snp.end=-1;
        continue;
      }
      if( per_err > per_err_threshold){
        e_obj.fullPlusDroppedOutput(i,j,temp1,temp2,min_snp,min_cm,finalErrors,per_err,3);
        removed1++;
        continue;
      }
      if( holdOut && hThreshold < ( noOfOppHom ) / hlength ){
        e_obj.fullPlusDroppedOutput(i,j,temp1,temp2,min_snp,min_cm,finalErrors,per_err,4);
        removed3++;
        snp.start= snp.end=-1;
        continue; 
      } 
      e_obj.finalOutPut(i,j,temp1,temp2,min_cm);
      continue;
    } //end FullPlusDropped

    //add finalOutput here
    if(option.compare("finalOutput") == 0){
    //only print if it passes all of the tests.
      int snp1 = 0, snp2 = 0, hlength = 0;
      float noOfOppHom = 0;
      if( holdOut )
      {
        snp1 = e_obj.getNewSnp( temp1 );
        snp2 = e_obj.getNewSnp( temp2 );
        hlength = snp2 - snp1;
        if( hlength <= 0 )
        {
          hlength = 1;
        }
        noOfOppHom = e_obj.getOppHomThreshold( pers1, pers2, snp.start, snp.end );
      }
      if( ( (beforeTrimEnd - beforeTrimStart) < min_snp) || ( (trimPositions.size() == 3) && (trimPositions[2] == 1) ) ){
        std::vector<float>movingAverages;
        temp1 = beforeTrimStart;
        temp2 = beforeTrimEnd;
        snp.start= snp.end=-1;
        removed4++;
        continue;
      }
      if( (( temp2-temp1 ) < min_snp) || (trimPositions.size() == 3) ){ 
        removed2++;
        snp.start= snp.end=-1;
        continue;
      }
      if( per_err > per_err_threshold){
        removed1++;
        continue;
      }
      if( holdOut && hThreshold < ( noOfOppHom ) / hlength ){
        removed3++;
        snp.start= snp.end=-1;
        continue;
      } 
      e_obj.finalOutPut(i,j,temp1,temp2,min_cm);
      continue;
    }//end finaloutput

    //Calculate Error1
    if( (option.compare("Error1") == 0 ) || (option.compare("ErrorRandom1") == 0) || (option.compare("Error") == 0) ){

      if( ( (beforeTrimEnd - beforeTrimStart) < min_snp) || ( (trimPositions.size() == 3) && (trimPositions[2] == 1) ) ){ //dropped before trimming
        //don't bother printing out ma for this one. But go back and change it so that it doesn't actually calc it
        std::vector<float>movingAverages;//null             
        //trying something special in this case. This can be removed once idrops aren't being trimmed
        //test code
        temp1 = beforeTrimStart;
        temp2 = beforeTrimEnd;
        //
        e_obj.errorOutput(i,j,temp1,temp2,min_snp,min_cm,movingAverages,finalErrors,per_err,temp1,temp2,beforeTrimStart,beforeTrimEnd,1);
        snp.start= snp.end=-1;
        removed4++; //seems ok
        continue;
      } 
      if( (( temp2-temp1 ) < min_snp) || ((trimPositions.size() == 3) && (trimPositions[2] == 2) ) ) //dropped after trimming
      {
        e_obj.errorOutput(i,j,temp1,temp2,min_snp,min_cm,movingAverages,finalErrors,per_err,temp1,temp2,beforeTrimStart,beforeTrimEnd,2);
        ++removed2;
        snp.start= snp.end=-1;
        continue;
      }
      if( per_err > per_err_threshold ) //dropped due to pie
      {
        e_obj.errorOutput(i,j,temp1,temp2,min_snp,min_cm,movingAverages,finalErrors,per_err,temp1,temp2,beforeTrimStart,beforeTrimEnd,3);
        ++removed1;
        snp.start= snp.end=-1;
        continue;
      }
      e_obj.errorOutput(i,j,temp1,temp2,min_snp,min_cm,movingAverages,finalErrors,per_err,temp1,temp2,beforeTrimStart,beforeTrimEnd,0);//no drop
      continue;
    }//end error1

    int snp1 = 0, snp2 = 0, hlength = 0;
    float noOfOppHom = 0;
    if( holdOut )
    {
      snp1 = e_obj.getNewSnp( temp1 );
      snp2 = e_obj.getNewSnp( temp2 );
      hlength = snp2 - snp1;
      if( hlength <= 0 )
      {
        hlength = 1;
      }
      noOfOppHom = e_obj.getOppHomThreshold( pers1, pers2, snp.start, snp.end );
    }
    //update drop order 2/26/14
    if( (( temp2-temp1 ) < min_snp) || (trimPositions.size() == 3) )
    {
      ++removed2;
      snp.start= snp.end=-1;
      continue;
    }
    if( per_err > per_err_threshold )
    {
      ++removed1;
      snp.start= snp.end=-1;
      continue;
    }
    if( option.compare("MovingAverages")==0 ) //make this ma2
    {
      if( holdOut)
      {
        e_obj.middleHoldOutPut(i,j,temp1,temp2, min_snp,min_cm,movingAverages,trimPositions,per_err, noOfOppHom, hlength );
      }  
      else
      {
        e_obj.middleOutPut(i,j,temp1,temp2, min_snp, min_cm,movingAverages, trimPositions,per_err );
      }
      continue;
    }
    if(option.compare("Error2")==0 || option.compare( "ErrorRandom2" ) == 0)
    {
      if( holdOut)
      { 
        e_obj.middleHoldOutPut(i,j,temp1,temp2, min_snp, min_cm, finalErrors, trimPositions, per_err, noOfOppHom, hlength );
      }
      else
      {
        e_obj.middleOutPut(i,j,temp1,temp2, min_snp, min_cm, finalErrors, trimPositions, per_err);
      }
      continue;
    }
    if ( holdOut && hThreshold < ( noOfOppHom ) / hlength )
    {
      ++removed3;
      snp.start= snp.end=-1;
      continue;
    }
    if( option.compare("Error3")==0 || option.compare( "ErrorRandom3" ) == 0  )
    {
      e_obj.middleHoldOutPut(i,j,temp1,temp2, min_snp, min_cm, finalErrors, trimPositions, per_err, noOfOppHom, hlength );
    }

  }//end while(!eof)
  /*Now, let's handle weighted output if need be*/
  if( option.compare("weightedOutput") == 0 ){
    float snp_average_count = 0.0;
    int start_position = find_genome_min();
    int end_position = find_genome_max();
    int genome_length = (end_position - start_position);
    genome_vector.resize(genome_length,0);
    for(int i = 0; i < m_weighted_sh.size(); i++){
      update_genome(m_weighted_sh[i].snp1, m_weighted_sh[i].snp2);
    }
    snp_average_count = average_snp_count();
    for(int i = 0; i < m_weighted_sh.size(); i++){
      m_weighted_sh[i].snp_weight = update_snp_weight(m_weighted_sh[i].snp1, m_weighted_sh[i].snp2);
    }
    for(int i = 0; i < m_weighted_sh.size(); i++){
      m_weighted_sh[i].final_weight = ( snp_average_count / (m_weighted_sh[i].snp_weight));
      e_obj.weightedOutput(m_weighted_sh[i].per1, m_weighted_sh[i].per2, m_weighted_sh[i].snp1, m_weighted_sh[i].snp2, m_weighted_sh[i].final_weight);
    }
  }
  /*End weighted output*/
  str = " \n No of matches removed due to percentage error: " + NumberToString( removed1 );
  str = str+ " \n No of matches removed due to length of trimming by moving averages: " + NumberToString( removed2 );
  str = str+ " \n No of matches removed due hold out ped file checking: "+ NumberToString( removed3 );
  str = str+ " \n No of matches removed due to initial length: " + NumberToString( removed4 );

  e_obj.log( str );
}

void Consolidator::performHoldOutTrim( ErrorCalculator& ecal, 
     float threshold, std::string hMiss,
     std::string option )
{
   threshold = getHoldOutThreshold( threshold );
   for(int i=0;i<person_count;i++)
   {
      for(int j=i;j<person_count;j++)
      {
          for(int l=0;l<m_matches[i][j].size();l++)
          {
              if(m_matches[i][j][l].end==-1)
              {
                  continue;
              }
              int temp1=m_matches[i][j][l].start;
              int temp2=m_matches[i][j][l].end;
              int  snp1 = ecal.getNewSnp( temp1 );
              int snp2 = ecal.getNewSnp( temp2 );
              int length = snp2 - snp1;
              float noOfOppHom = ecal.getOppHomThreshold( i, j, m_matches[i][j][l].start, m_matches[i][j][l].end );
              if( length != 0 && threshold < (noOfOppHom) / length )
              {
                 m_matches[i][j][l].start = m_matches[i][j][l].end = -1;
              }
         }
     }
  }
}
void Consolidator::finalOutPut(ErrorCalculator &e,float min_cm, int min_snp )const
{

     for(int i=0;i<person_count;i++)
        {
                for(int j=i;j<person_count;j++)
                {

                         for(int l=0;l<m_matches[i][j].size();l++)
                         { 
        if(m_matches[i][j][l].start==-1||m_matches[i][j][l].end==-1 || ( m_matches[i][j][l].end - m_matches[i][j][l].start ) < min_snp ) continue;
                                e.finalOutPut(i,j,m_matches[i][j][l].start,m_matches[i][j][l].end ,min_cm);                
                         }

    }

  }

}
/* Testing purposes June 9th 2014 at 2:30pm
void Consolidator::findTrueSimplePctErrors( ErrorCalculator &e_obj, float PIElength, bool holdOut,int window, float ma_threshold, float empirical_ma_threshold )
{
   for(int i=0;i<person_count;i++)
   {
      for(int j=i;j<person_count;j++)
      {
          for(int l=0;l<m_trueMatches[i][j].size();l++)
          {
              if(m_trueMatches[i][j][l].end==-1)
              {
                  continue;
              }
       //-------------------------------------------------------------------------------------------------
       int t1 = m_trueMatches[ i ][ j ][ l ].start +
                          ( m_trueMatches[ i ][ j ][ l ].end -
                            m_trueMatches[ i ][ j ][ l ].start ) * 0.25;
             int t2 = m_trueMatches[ i ][ j ][ l ].end -
                          ( m_trueMatches[ i ][ j ][ l ].end -
                            m_trueMatches[ i ][ j ][ l ].start ) * 0.25;

       vector<vector<int> > trueErrors=e_obj.checkErrors( i, j, t1, t2);
             vector<int>finalTrueErrors=e_obj.getFinalErrors( trueErrors );
       //This section handles finding the maximum moving averages amongst trulyIBD segments
       std::vector<float> av;
       float current_max;
       if(empirical_ma_threshold < 0.0){
             av = e_obj.getTrueMovingAverages(finalTrueErrors,t1,t2,window);
             current_max = av[0];
             for(int q = 1; q < av.size(); q++){
                   if(av[q] > current_max){
                           current_max = av[q];
                   }
              }
             e_obj.addMaxAverage(current_max);
       }
       //------------------------------------------------------------------------------
             int temp1 = m_trueMatches[i][j][l].start;
             int temp2 = m_trueMatches[i][j][l].end;
             float startCM = e_obj.getCMDistance( temp1 );
             float endCM = e_obj.getCMDistance( temp2 );
             float mid1CM = startCM + ( endCM - startCM ) / 2 - PIElength / 2;
             float mid2CM = startCM + ( endCM - startCM ) / 2 + PIElength / 2;
             while( e_obj.getCMDistance( temp1 ) <= mid1CM || e_obj.getCMDistance( temp2 ) >=mid2CM )
             {
                if( e_obj.getCMDistance( temp1 ) <= mid1CM )
                {
                  ++temp1;
                }
                if( e_obj.getCMDistance( temp2 ) >=mid2CM )
                {
                  --temp2;
                }
             }
              vector<vector<int> > errors=e_obj.checkErrors( i, j, temp1, temp2);

                  vector<int>finalErrors=e_obj.getFinalErrors( errors );
//                  float per_err = e_obj.getThreshold(finalErrors,temp1, temp2, 0 );
      float per_err = e_obj.getThreshold(finalErrors,temp1,temp2); //overload!
                  m_errors.push_back( per_err );
                 if( holdOut  )
                 {
                        float oppHom =
                                     ( e_obj.getOppHomThreshold( i, j,
                                               temp1,
                                             temp2 ) ) / ( temp2 -temp1 );
                        m_holdOutErrors.push_back( oppHom );
                 }

          }
     }
  }
   
   //this section actually handles the sorting of the max averages, and the setting of the user supplied percentile.
   vector<float>maxes;
   float cutoff = empirical_ma_threshold; //assume the user wanted to supply a value. This value will be overwritten shortly if they did not.
   if(empirical_ma_threshold < 0.0){
   maxes = e_obj.getMaxAverages();
   std::sort(maxes.begin(),maxes.end());
   e_obj.setMaxAverage(maxes);
   cutoff = e_obj.getXthPercentile(ma_threshold); //<-make that an actual user input value
   }
   
   e_obj.setCutoff(cutoff);//set the actual threshold to be used when calculating MA in all other SH

   if(empirical_ma_threshold < 0.0){
   std::string outt = "\n User supplied ma-threshold is: " + NumberToString(ma_threshold);
   outt = outt + "\n Moving Averages will be tested usign the empirical threshold: " + NumberToString(cutoff);
   e_obj.log(outt);
   } else {
   std::string outt = "\n User supplied empirical-ma-threshold is: " + NumberToString(cutoff);
   outt = outt + "\n Moving Averages will be tested usign the empirical threshold: " + NumberToString(cutoff);
   e_obj.log(outt);
   }
   //----------------------------------------


   std::sort( m_errors.begin(), m_errors.end() );
   std::sort( m_holdOutErrors.begin(), m_holdOutErrors.end() );
   std::string str =  " \n No of segments deemed to be IBD for finding empirical error threshold "
                      + NumberToString( m_errors.size() );
   str  = str + " \n No of segments deemed to be IBD for finding empirical error threshold in hold out are "
                      + NumberToString( m_holdOutErrors.size() );
        e_obj.log( str );

}

void Consolidator::findTruePctErrors( ErrorCalculator &e_obj,int ma_snp_ends, bool holdOut,int window,float ma_threshold, float empirical_ma_threshold )
{
   for(int i=0;i<person_count;i++)
   {
      for(int j=i;j<person_count;j++)
      {
          for(int l=0;l<m_trueMatches[i][j].size();l++)
          {
              if(m_trueMatches[i][j][l].end==-1)
              {
                  continue;
              }
        //handle moving averages calculation
        //
        int t1 = m_trueMatches[ i ][ j ][ l ].start +
                          ( m_trueMatches[ i ][ j ][ l ].end -
                            m_trueMatches[ i ][ j ][ l ].start ) * 0.25;
              int t2 = m_trueMatches[ i ][ j ][ l ].end -
                           ( m_trueMatches[ i ][ j ][ l ].end -
                            m_trueMatches[ i ][ j ][ l ].start ) * 0.25;    
        //now we have the positions of the first and last 25% of the truly ibd SH
        //all that's left to do is to pass them into the moving averages function, and obtain the max ma
        //then store that in a vector, sort them, and find the xth percentile of that vector. That will be
        //the ma that we use later
        //for that "finalErrors" parameters, need to get the number of errors along the truly IBD SH first...
        vector<vector<int> > trueErrors=e_obj.checkErrors( i, j, t1, t2);
              vector<int>finalTrueErrors=e_obj.getFinalErrors( trueErrors );

        //handles MA calculations
        std::vector<float> av;
        float current_max;
        if(empirical_ma_threshold < 0.0){
        av = e_obj.getTrueMovingAverages(finalTrueErrors,t1,t2,window);
              current_max = av[0];
              for(int q = 1; q < av.size(); q++){
                   if(av[q] > current_max){
                           current_max = av[q];
                   }
              }
              e_obj.addMaxAverage(current_max);
        }
        
        //
              int temp1 = m_trueMatches[ i ][ j ][ l ].start +
                          ( m_trueMatches[ i ][ j ][ l ].end -
                            m_trueMatches[ i ][ j ][ l ].start ) * 0.15; //Should probably stop doing this
              int temp2 = m_trueMatches[ i ][ j ][ l ].end - 
                           ( m_trueMatches[ i ][ j ][ l ].end - 
                            m_trueMatches[ i ][ j ][ l ].start ) * 0.15;
              int start =0, end =0, fend = ( temp2 -temp1 )  ;
               
      //since we are using MOL at this point, this will pick out a random SH from the set of non-truly IBD SH 
      //and use that length to define the region over which we find PIE. Unless you are changing something with MOL,
      //don't ever read this next block
                  int randPers1, randPers2, pos;
                  randPers1 = std::rand() % person_count;
                  randPers2 = std::rand() % person_count;
                  if( randPers1 > randPers2 )
                  {
                     randPers1 = randPers1 + randPers2;
                     randPers2 = randPers1 - randPers2;
                     randPers1 = randPers1 - randPers2;
                  }
                  while( m_matches[ randPers1 ][ randPers2 ].size() <= 0 )
                  {
                    randPers1 = std::rand() % person_count;
                    randPers2 = std::rand() % person_count;
                   if( randPers1 > randPers2 )
                   {
                      randPers1 = randPers1 + randPers2;
                      randPers2 = randPers1 - randPers2;
                      randPers1 = randPers1 - randPers2;
                   }

                  }
                  pos = std::rand() % m_matches[ randPers1 ][ randPers2 ].size();
                  int len = m_matches[ randPers1 ][ randPers2 ][ pos ].end 
                            - m_matches[ randPers1 ][ randPers2 ][ pos ].start;
                  if( len >= fend || len <= 0)
                  {
                      continue;
                  } 
                  temp1 = temp1;
                  temp2 = temp1 + len;
      //end crazy MOL stuff
                  vector<vector<int> > errors=e_obj.checkErrors( i, j, temp1, temp2);

                  vector<int>finalErrors=e_obj.getFinalErrors( errors );
//                  float per_err = e_obj.getThreshold(finalErrors,temp1,temp2,ma_snp_ends );
      float per_err = e_obj.getThreshold(finalErrors,temp1,temp2);//overload
                  m_errors.push_back( per_err );
                 if( holdOut  )
                 {
                        float oppHom = ( e_obj.getOppHomThreshold( i, j, temp1, temp2 ) ) / ( temp2 -temp1 );
                        m_holdOutErrors.push_back( oppHom );
                 }
          }
       }  
   }
   vector<float>maxes;
   float cutoff = empirical_ma_threshold;
   if(empirical_ma_threshold < 0.0){
   maxes = e_obj.getMaxAverages();
   std::sort(maxes.begin(),maxes.end());
   e_obj.setMaxAverage(maxes);
   cutoff = e_obj.getXthPercentile(ma_threshold); 
   }
   e_obj.setCutoff(cutoff);//set the actual threshold to be used when calculating MA in all other SH
   //
   std::sort( m_errors.begin(), m_errors.end() );
   std::sort( m_holdOutErrors.begin(), m_holdOutErrors.end() );
   std::string str =  " \n No of elements in error check are: "
                      + NumberToString( m_errors.size() );
   str  = str + " \n No of elements in hold  error check are: "
                      + NumberToString( m_holdOutErrors.size() );

        e_obj.log( str );

}
Remove this comment to restore*/
void Consolidator::findTrueSimplePctErrors( ErrorCalculator &e_obj, float PIElength, bool holdOut,int window, float ma_threshold, float empirical_ma_threshold, string path, int person_count, int trueSNP, float trueCM  )
{
  //begin streaming in bmatch data
  unsigned int pid[2];
  unsigned int sid[2];
  unsigned int dif,hom[2];
  ifstream file_bmatch(path.c_str(),ios::binary);
  int pers_count = person_count;
  SNP snp;
  if( !file_bmatch )
  {
    cerr<<"unable to open the bmatch file, exiting the program" << endl;
    exit( -1 );
  }
  //will need to replace all instances of m_matches in this function with pid[0], pid[1], etc.
  //will probably encapsulate all of the below code into a while !eof loop
  while ( !file_bmatch.eof()){
    //read in bmatch data
    pid[0] = -1;
    file_bmatch.read( (char*) &pid[0] , sizeof( unsigned int ) );
    if ( pid[0] == -1 ) continue;
    file_bmatch.read( (char*) &pid[1] , sizeof( unsigned int ) );
    file_bmatch.read( (char*) &sid[0] , sizeof( unsigned int ) );
    file_bmatch.read( (char*) &sid[1] , sizeof( unsigned int ) );
    file_bmatch.read( (char*) &dif , sizeof( int ) );
    file_bmatch.read( (char*) &hom[0] , sizeof( bool ) );
    file_bmatch.read( (char*) &hom[1] , sizeof( bool ) );
    if(pid[0]>=pers_count||pid[1]>=pers_count)
    {
      cerr<<"problem with bsid file in findTrueSimplePctErrors, check it please"<<endl;
      return;
    }

    snp.start=sid[0];
    snp.end=sid[1];
    if(pid[0] > pid[1]){
      int temp = pid[0];
      pid[0] = pid[1];
      pid[1] = temp;
    }
    int i = pid[0];
    int j = pid[1];
    if(snp.end==-1)
    {
      continue;
    }
    //check to see if this SH is trulyIBD or not. do not currently need this, as germline input will be trueIBD
    if( ( e_obj.getCMDistance( sid[ 1 ] ) - 
    e_obj.getCMDistance( sid[ 0 ] ) ) >= trueCM && 
    ( sid[ 1 ] - sid[ 0 ] ) >= trueSNP &&  pid[0] != pid[1] )
    {
      //this indicates a trulyIBD segment, so we can go to work on it with the following code. Otherwise, just continue.
      /*If you need to unroll these loops as well, you would add the stream here, and pluck out items that are 
      trulyIBD. This would result in streaming the data twice, as I don't think that there is a way to run non-IBD segments
      without first having run all of the trulyIBD segments.*/
      //-------------------------------------------------------------------------------------------------
      int t1 = snp.start +
      ( snp.end -
      snp.start ) * 0.25;
      int t2 = snp.end -
      ( snp.end -
      snp.start ) * 0.25;

      vector<vector<int> > trueErrors=e_obj.checkErrors( i, j, t1, t2);
      vector<int>finalTrueErrors=e_obj.getFinalErrors( trueErrors );
      //This section handles finding the maximum moving averages amongst trulyIBD segments
      std::vector<float> av;
      float current_max;
      if(empirical_ma_threshold < 0.0){ //this would indicate that the user did NOT supply an empirical-ma-threshold argument. So we need to calculate one.
        av = e_obj.getTrueMovingAverages(finalTrueErrors,t1,t2,window);
        current_max = av[0];
        for(int q = 1; q < av.size(); q++){
          if(av[q] > current_max){
            current_max = av[q];
          }
        }
        e_obj.addMaxAverage(current_max);
      }
      //------------------------------------------------------------------------------
      int temp1 = snp.start;
      int temp2 = snp.end;
      float startCM = e_obj.getCMDistance( temp1 );
      float endCM = e_obj.getCMDistance( temp2 );
      float mid1CM = startCM + ( endCM - startCM ) / 2 - PIElength / 2;
      float mid2CM = startCM + ( endCM - startCM ) / 2 + PIElength / 2;
      while( e_obj.getCMDistance( temp1 ) <= mid1CM || e_obj.getCMDistance( temp2 ) >=mid2CM )
      {
        if( e_obj.getCMDistance( temp1 ) <= mid1CM )
        {
          ++temp1;
        }
        if( e_obj.getCMDistance( temp2 ) >=mid2CM )
        {
          --temp2;
        }
      }
      vector<vector<int> > errors=e_obj.checkErrors( i, j, temp1, temp2);

      vector<int>finalErrors=e_obj.getFinalErrors( errors );
      //                  float per_err = e_obj.getThreshold(finalErrors,temp1, temp2, 0 );
      float per_err = e_obj.getThreshold(finalErrors,temp1,temp2); //overload!
      m_errors.push_back( per_err );
      if( holdOut  )
      {
        float oppHom = ( e_obj.getOppHomThreshold( i, j, temp1, temp2 ) ) / ( temp2 -temp1 );
        m_holdOutErrors.push_back( oppHom );
      }
    }else {//end check for trueIBD segment 
      continue; //this segment is not trueIBD, so go onto next SH
    }
  }//end while(!eof)

  //this section actually handles the sorting of the max averages, and the setting of the user supplied percentile.
  vector<float>maxes;
  float cutoff = empirical_ma_threshold; //assume the user wanted to supply a value. This value will be overwritten shortly if they did not.
  if(empirical_ma_threshold < 0.0){
    maxes = e_obj.getMaxAverages();
    std::sort(maxes.begin(),maxes.end());
    e_obj.setMaxAverage(maxes);
    cutoff = e_obj.getXthPercentile(ma_threshold); //<-make that an actual user input value
  }

  e_obj.setCutoff(cutoff);//set the actual threshold to be used when calculating MA in all other SH

  if(empirical_ma_threshold < 0.0){
    std::string outt = "\n User supplied ma-threshold is: " + NumberToString(ma_threshold);
    outt = outt + "\n Moving Averages will be tested usign the empirical threshold: " + NumberToString(cutoff);
    e_obj.log(outt);
  } else {
    std::string outt = "\n User supplied empirical-ma-threshold is: " + NumberToString(cutoff);
    outt = outt + "\n Moving Averages will be tested usign the empirical threshold: " + NumberToString(cutoff);
    e_obj.log(outt);
  }
  std::sort( m_errors.begin(), m_errors.end() );
  std::sort( m_holdOutErrors.begin(), m_holdOutErrors.end() );
  std::string str =  " \n No of segments deemed to be IBD for finding empirical error threshold " + NumberToString( m_errors.size() );
  str  = str + " \n No of segments deemed to be IBD for finding empirical error threshold in hold out are " + NumberToString( m_holdOutErrors.size() );
  e_obj.log( str );

}//end this function

void Consolidator::findTruePctErrors( ErrorCalculator &e_obj,int ma_snp_ends, bool holdOut,int window,float ma_threshold, float empirical_ma_threshold, string path, int person_count, int trueSNP, float trueCM )
{
   //begin streaming in bmatch data
        unsigned int pid[2];
        unsigned int sid[2];
        unsigned int dif,hom[2];
        ifstream file_bmatch(path.c_str(),ios::binary);
        int pers_count = person_count;
        SNP snp;
        if( !file_bmatch )
         {
             cerr<<"unable to open the bmatch file, exiting the program" << endl;
             exit( -1 );

         }
         //will need to replace all instances of m_matches in this function with pid[0], pid[1], etc.
         //will probably encapsulate all of the below code into a while !eof loop
         
         while ( !file_bmatch.eof()){
         //read in bmatch data
                pid[0] = -1;
                file_bmatch.read( (char*) &pid[0] , sizeof( unsigned int ) );
                if ( pid[0] == -1 ) continue;
                file_bmatch.read( (char*) &pid[1] , sizeof( unsigned int ) );
                file_bmatch.read( (char*) &sid[0] , sizeof( unsigned int ) );
                file_bmatch.read( (char*) &sid[1] , sizeof( unsigned int ) );
                file_bmatch.read( (char*) &dif , sizeof( int ) );
                file_bmatch.read( (char*) &hom[0] , sizeof( bool ) );
                file_bmatch.read( (char*) &hom[1] , sizeof( bool ) );
                if(pid[0]>=pers_count||pid[1]>=pers_count)
                {
                      cerr<<"problem with bsid file in findTruePctErrors, check it please"<<endl;
                      return;

                }
                 
                snp.start=sid[0];
                snp.end=sid[1];
                if(pid[0] > pid[1]){
                  int temp = pid[0];
                  pid[0] = pid[1];
                  pid[1] = temp;
                }

                int i = pid[0];
                int j = pid[1];




              if(snp.end==-1)
              {
                  continue;
              }

               if( ( e_obj.getCMDistance( sid[ 1 ] ) - 
                                e_obj.getCMDistance( sid[ 0 ] ) ) >= trueCM && 
                                 ( sid[ 1 ] - sid[ 0 ] ) >= trueSNP &&  pid[0] != pid[1] )
                 {

        //handle moving averages calculation
        //
        int t1 = snp.start +
                          ( snp.end -
                            snp.start ) * 0.25;
              int t2 = snp.end -
                           ( snp.end -
                            snp.start ) * 0.25;   
        //now we have the positions of the first and last 25% of the truly ibd SH
        //all that's left to do is to pass them into the moving averages function, and obtain the max ma
        //then store that in a vector, sort them, and find the xth percentile of that vector. That will be
        //the ma that we use later
        //for that "finalErrors" parameters, need to get the number of errors along the truly IBD SH first...
        vector<vector<int> > trueErrors=e_obj.checkErrors( i, j, t1, t2);
              vector<int>finalTrueErrors=e_obj.getFinalErrors( trueErrors );

        //handles MA calculations
        std::vector<float> av;
        float current_max;
        if(empirical_ma_threshold < 0.0){
        av = e_obj.getTrueMovingAverages(finalTrueErrors,t1,t2,window);
              current_max = av[0];
              for(int q = 1; q < av.size(); q++){
                   if(av[q] > current_max){
                           current_max = av[q];
                   }
              }
              e_obj.addMaxAverage(current_max);
        }
        
        //
              int temp1 = snp.start +
                          ( snp.end -
                            snp.start ) * 0.15; //Should probably stop doing this
              int temp2 = snp.end - 
                           ( snp.end - 
                             snp.start ) * 0.15;
              int start =0, end =0, fend = ( temp2 -temp1 )  ;
               
      //since we are using MOL at this point, this will pick out a random SH from the set of non-truly IBD SH 
      //and use that length to define the region over which we find PIE. Unless you are changing something with MOL,
      //don't ever read this next block
              /*
                  int randPers1, randPers2, pos;
                  randPers1 = std::rand() % person_count;
                  randPers2 = std::rand() % person_count;
                  if( randPers1 > randPers2 )
                  {
                     randPers1 = randPers1 + randPers2;
                     randPers2 = randPers1 - randPers2;
                     randPers1 = randPers1 - randPers2;
                  }
                  while( m_matches[ randPers1 ][ randPers2 ].size() <= 0 )
                  {
                    randPers1 = std::rand() % person_count;
                    randPers2 = std::rand() % person_count;
                   if( randPers1 > randPers2 )
                   {
                      randPers1 = randPers1 + randPers2;
                      randPers2 = randPers1 - randPers2;
                      randPers1 = randPers1 - randPers2;
                   }

                  }
                  pos = std::rand() % m_matches[ randPers1 ][ randPers2 ].size();
                  int len = m_matches[ randPers1 ][ randPers2 ][ pos ].end 
                            - m_matches[ randPers1 ][ randPers2 ][ pos ].start;
                  if( len >= fend || len <= 0)
                  {
                      continue;
                  } 
                
                  temp1 = temp1;
                  temp2 = temp1 + len;
                  */
      //end crazy MOL stuff
                 
                  vector<vector<int> > errors=e_obj.checkErrors( i, j, temp1, temp2);
                  vector<int>finalErrors=e_obj.getFinalErrors( errors );
//                  float per_err = e_obj.getThreshold(finalErrors,temp1,temp2,ma_snp_ends );
                  float per_err = e_obj.getThreshold(finalErrors,temp1,temp2);//overload
                  m_errors.push_back( per_err );
                 if( holdOut  )
                 {
                        float oppHom = ( e_obj.getOppHomThreshold( i, j, temp1, temp2 ) ) / ( temp2 -temp1 );
                        m_holdOutErrors.push_back( oppHom );
                 }

          }else{

            continue;
          }
    }//end while()
      
   vector<float>maxes;
   float cutoff = empirical_ma_threshold;
   if(empirical_ma_threshold < 0.0){
   maxes = e_obj.getMaxAverages();
   std::sort(maxes.begin(),maxes.end());
   e_obj.setMaxAverage(maxes);
   cutoff = e_obj.getXthPercentile(ma_threshold); 
   }
   e_obj.setCutoff(cutoff);//set the actual threshold to be used when calculating MA in all other SH
   //
   std::sort( m_errors.begin(), m_errors.end() );
   std::sort( m_holdOutErrors.begin(), m_holdOutErrors.end() );
   std::string str =  " \n No of elements in error check are: "
                      + NumberToString( m_errors.size() );
   str  = str + " \n No of elements in hold  error check are: "
                      + NumberToString( m_holdOutErrors.size() );

        e_obj.log( str );

}




float Consolidator::getHoldOutThreshold( float threshold )
{
 int pos = ( m_holdOutErrors.size() ) * threshold;
 if( pos >= m_holdOutErrors.size() || pos < 0 )
 {
    cerr<< "wrong threshold value. it should be between 0-99.99";
    exit( -1 );
 }
 return m_holdOutErrors[ pos ];  
}
float Consolidator::getPctErrThreshold( float threshold)
{
   int pos = ( m_errors.size()  ) * threshold;
 if( pos >= m_errors.size() || pos < 0 )
 {
    cerr<< "wrong threshold value. it should be between 0-99.99"
        << "errors list size is " << m_errors.size()<<endl ;
    exit( -1 );
 }
 return m_errors[ pos ];
}

/*Weighted output functions*/
void Consolidator::update_genome(int snp1,int snp2){
  int s1,s2;
  /*if(snp1 != 0){
    s1 = snp1-1;
  } else {
    cerr << "WE HAVE A SNP THAT IS 0!" << endl;
  }
  if(snp2 != 0){
    s2 = snp2-1;
  } else{
    cerr <<  "WE HAVE A SNP THAT IS 0" << endl;
  }*/

  s1 = snp1;
  s2 = snp2;
  for(int i = s1; i <= s2; i++){
    genome_vector[i] += 1;
  }
}
void Consolidator::print_genome(){
  cout << "Printing Geneome..." << endl;
  for(int i = 0; i < genome_vector.size(); i++){
    cout << genome_vector[i] << ", ";
  }
}
float Consolidator::update_snp_weight(int snp1,int snp2){
/*  if( (snp1 == 0) || (snp2 == 0)){
    cerr << "WE HAVE A SNP THAT IS 0" << endl;
  }
  snp1 = snp1 - 1;
  snp2 = snp2 - 1;*/
  int length = ((snp2 - snp1) + 1); //counts are inclusive, so [10,20] has 11 elements, not 10.
  float sum = 0.0;
  float weight = 0.0;
  for(int i = snp1; i <= snp2; i++){
    sum += genome_vector[i];
  }
  weight = sum / length;
  return weight;
}
int Consolidator::find_genome_min(){
  int min = m_weighted_sh[0].snp1;
  for(int i = 1; i < m_weighted_sh.size(); i++){
    if(m_weighted_sh[i].snp1 < min){
      min = m_weighted_sh[i].snp1;
    }
  }
  return min;
}
int Consolidator::find_genome_max(){
  int max = m_weighted_sh[0].snp2;
  for(int i = 1; i < m_weighted_sh.size(); i++){
    if(m_weighted_sh[i].snp2 > max){
      max = m_weighted_sh[i].snp2;
    }
  }
  return max;
}
float Consolidator::average_snp_count(){
  float sum = 0.0;
  float avg = 0.0;
  int zeros = 0;
  for(int i = 0; i < genome_vector.size(); i++){
    if(genome_vector[i] == 0){
      zeros += 1;
    } else {
      sum += genome_vector[i];
    }
  }
  avg = sum / ( genome_vector.size() - zeros);
  return avg;
}
/**/
