#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <iomanip>
#include <string.h>
#include <cstdlib>

bool REDUCE=false, CM_FLAG=false, SNP_FLAG=false,INDEX=false,HAPLOID=true;
int MIN_SNP=0,pers_count=0;
float MIN_CM=0.0;
using namespace std;

struct Marker
{
	short	chr;
	string	rsid;
	float	cm_distance;
	long	bp_distance;
};

//This function simply calculates the distance in cM between two markers.
float get_distance( Marker& m_left , Marker& m_right , bool& genetic )
{
	if ( m_left.cm_distance == -1 || m_right.cm_distance == -1 )
	{
		genetic = false;
		return m_right.bp_distance - m_left.bp_distance;
	} else {
		genetic = true;
		return m_right.cm_distance - m_left.cm_distance;
	}
}

int main (int argc, char* argv[])
{
	if(argc < 4|| argc>8){
			cerr << "Usage: " << argv[0] << " [BMATCH FILE] [BSID FILE] [BMID FILE] or [BMATCH FILE] [BSID FILE] [BMID FILE] -reduced [max_snp] [max_cm] -make-index(optional)" << endl;
			return 0;
	}

	if(strcmp(argv[4],"-reduced")==0)
	{
		REDUCE=true;

		if(argc>=7){
			CM_FLAG=true;
			SNP_FLAG=true;
			MIN_SNP=atoi(argv[5]);
			MIN_CM=atof(argv[6]);
		}
		if(argc>=8&&strcmp(argv[7],"-no-haploid")==0)
			HAPLOID=false;
		else if(argc>=8&&strcmp(argv[7],"-make-index")==0)
			INDEX=true;

		if(argc>=9&&strcmp(argv[8],"-make-index")==0)
			INDEX=true;
	} else {
		 cerr << "Usage: " << argv[0] << " [BMATCH FILE] [BSID FILE] [BMID FILE] or [BMATCH FILE] [BSID FILE] [BMID FILE] -reduced -min_snp -min_cm -make-index(optional)" << endl;
	}


	string line, discard;
	ifstream file_bmatch( argv[1] , ios::binary );
	ifstream file_bsid( argv[2] );
	ifstream file_bmid( argv[3] );
	if(!file_bmatch || !file_bsid || !file_bmid ) { cerr << "file could not be opened" << endl; return 0; }

	stringstream ss;
	string item;
	// load samples
	vector< string > sample_id;
	while( getline(file_bsid , line) )
	{
		//cout<<"the line: "<<line<<endl;
		ss.clear();ss.str(line);
		while(getline(ss,item,' '))
		{
		//	cout<<"the item is: "<<item<<endl;
			sample_id.push_back( item );
			pers_count++;
		}
	}

	file_bsid.close();
	// cout<<"number of persons are "<<pers_count<<endl;
	// load markers
	vector< Marker > marker_id;
	Marker cur_marker;
	while ( getline(file_bmid , line) )
	{
		ss.clear(); ss.str( line );
        ss >> cur_marker.chr >> cur_marker.rsid >> cur_marker.cm_distance >> cur_marker.bp_distance;
		//cout<<"the marker_id"<<cur_marker.bp_distance<<endl;
		marker_id.push_back( cur_marker );
	}
	file_bmid.close();

	// load matches
	unsigned int pid[2];
	unsigned int sid[2],min1;
	int dif;
	float min2;
	bool hom[2] , genetic;

	while ( !file_bmatch.eof() )
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
		min1=sid[1]-sid[0]+1;
		min2=get_distance( marker_id[sid[0]] , marker_id[sid[1]] , genetic );
		if(SNP_FLAG&&CM_FLAG)if(min1<MIN_SNP||min2<MIN_CM) continue;
		if(!HAPLOID&&(pid[0]*2>=sample_id.size()||pid[1]*2>=sample_id.size()))
		{
			cerr<<"the person ids are not correct, please check the bsid file pid0="<<pid[0]<<" pid1="<<pid[1]<<endl;
			return -1;
		}
		if(sid[0]>=marker_id.size()||sid[1]>=marker_id.size())
		{
			cerr<<"the marker file bmid is not matching to bmatch file, please check it"<<endl;
			return -1;
		}

		string sub;
		if(!HAPLOID)
			sub=sample_id[ pid[0]*2];
		else
			sub=sample_id[pid[0]];
		cout<< sub << '\t';
		if(!HAPLOID)
			sub=sample_id[pid[1]*2];
		else sub=sample_id[pid[1]];
			cout << sub << '\t';
		if(!REDUCE)
			cout << marker_id[ sid[0] ].chr << '\t';

		cout << marker_id[ sid[0] ].bp_distance << "\t" << marker_id[ sid[1] ].bp_distance << "\t";

		if(!REDUCE){
			cout << marker_id[ sid[0] ].rsid << '\t' << marker_id[ sid[1] ].rsid << '\t';
		}
		cout <<min1 << '\t' << setiosflags(ios::fixed) << setprecision(2) <<min2 << '\t';
		if(!REDUCE){
			if ( genetic ) cout << "cM\t"; else cout << "MB\t";
			cout << dif << '\t';
			if ( hom[0] ) cout << "1\t"; else cout << "0\t";
			if ( hom[1] ) cout << "1\t"; else cout << "0\t";
		}
		cout << endl;
		if(INDEX)
			cout <<pid[0]<<"\t"<<pid[1]<<"\t"<<sid[0]<<"\t"<<sid[1]<< endl;
	}
	file_bmatch.close();
}
