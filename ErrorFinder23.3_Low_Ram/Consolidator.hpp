#include "ErrorCalculator.cpp"

class SNP
{
 public:
 int start,end;

};
class Weighted_SH{
public:
    int snp1, snp2;
    int per1,per2;
    float final_weight;
    float snp_weight;
    Weighted_SH(int s1,int s2,int p1,int p2){
        snp1 = s1;
        snp2 = s2;
        per1 = p1;
        per2 = p2;
    }
};
class Consolidator
{
        public:
	//overloaded version to support new -ma-threshold argument and to depricate -ma-err-threshold-start and end
	void performTrim(ErrorCalculator& e_obj,int window,int ma_snp_ends,float ma_threshold,int min_snp,float min_cm,float per_err_threshold,string option,float hThreshold,bool holdout,float empirical_threshold, std::string path, float emp_pie_thresh);
 
       
	void performConsolidation(ErrorCalculator& ecal, 
              int gap, int min_snp, float min_cm);

        void performHoldOutTrim( ErrorCalculator& ecal, 
             float threshold, std::string hMiss, std::string );      
 
        void readMatches( std::string path, int pers_count, 
                     ErrorCalculator& eCalculator, int trueSNP, float trueCM );

        void finalOutPut(ErrorCalculator &e,float min_cm, int min_snp )const;
        void findTruePctErrors( ErrorCalculator &e, int ma_snp_ends, bool holdOut,int window,float ma_threshold, float empirical_ma_threshold, std::string path, int pr_count, int tsnp, float tcm);//overlaoding to add window argument and threshold as well
        void findTrueSimplePctErrors( ErrorCalculator &e, float PIElength, bool holdOut, int window, float ma_threshold, float empirical_ma_threshold, std::string path, int pr_count, int tsnp, float tcm); //removed int pers_count and path
        //void findTruePctErrors( ErrorCalculator &e, int ma_snp_ends, bool holdOut,int window,float ma_threshold, float empirical_ma_threshold);//overlaoding to add window argument and threshold as well
        //void findTrueSimplePctErrors( ErrorCalculator &e, float PIElength, bool holdOut, int window, float ma_threshold, float empirical_ma_threshold); //removed int pers_count and path
        float getPctErrThreshold( float threshold);
        float getHoldOutThreshold( float threshold );

        /*Weighted output functions*/
        void update_genome(int snp1,int snp2);
        void print_genome();
        float update_snp_weight(int snp1,int snp2);
        int find_genome_min();
        int find_genome_max();
        float average_snp_count();
        /**/

        void setPersonCount(int p_count){
            person_count = p_count;
        }

	//new methods to handle calculations of moving averages -nate 2/11/2014
	std::vector < std::vector < std::vector < SNP > > > getTrueMatches(){
		return m_trueMatches;
	}

	void setTrueMatches(std::vector < std::vector < std::vector < SNP > > > x){
		m_trueMatches = x;
	}
        private:
       int person_count;

        std::vector < std::vector < std::vector < SNP > > > m_matches;
        std::vector < std::vector < std::vector < SNP > > > m_trueMatches;
        std::vector< float > m_errors;
        std::vector< float > m_holdOutErrors;
        //new, for weighting algo
        std::vector < Weighted_SH > m_weighted_sh;
        //this should actually hold the ints for the genome array in the weighted algorithm
        std::vector <int> genome_vector;
        void sortMatches();

        static bool compareFunction(SNP s1, SNP s2);   
        std::string consolidated_str;
        std::string initial_drop_str;
        std::string emp_pie_thresh_str;
        std::string ma_drop_str;
        std::string pie_drop_str;
        std::string emp_ma_thresh_str;
        std::string ibg_str;
        std::string final_sh_str;
        std::string ma_thresh_str;
        int global_initial;
};
