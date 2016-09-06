#include <iostream>
#include <sstream>
#include <unistd.h>
#include <sstream>
#include <string>
#include <cstring>
#include <cstdlib>
using namespace std;

int main(int argc, char * argv[])
{	
	std::string  line ="";
	std::string low_ram = "-low_ram";
	std::stringstream ss("");
	bool flag=false;
	//char** tokens;
	
	//*tokens = strtok(*tokens, " ");
	//execvp(tokens[0], tokens);

//std::cout<<argc <<std::endl;
ss.clear();
//ss<<argv[0]<<" ";
for(int i=1;i<argc;i++)
{
	if (!strcmp(low_ram.c_str(),argv[i]))
	{
		//-low_ram flag is passed
		flag = true;
		continue;
	}
	else
	{
		//std::cout<<argv[i]<<std::endl;
		ss<< argv[i]<<" ";
		//std::cout<<ss.str()<<std::endl;
	}
}



if (flag == true)//lowram
{
	flag = system(("./ErrorFinder23.3_Low_Ram/FISHR_Low_Ram " +  ss.str()).c_str());
}
else if(flag == false)
{
	flag = system(("./ErrorFinder23.3/FISHR " + ss.str()).c_str());
}

else
{

}



/*
tokens = new char*[ss.str().size() + 1];	//http://www.cplusplus.com/forum/unices/51539/
strcpy(*tokens, ss.str().c_str());

if (flag == true)
{
	system("/home/piyush/Documents/IBG-FISHR-GITHUB/Project-FISHR-Low-Ram/FISHR_Low_Ram",tokens);
}
else if(flag == false)
{
	system("/home/piyush/Documents/IBG-FISHR-GITHUB/Project-FISHR/FISHR",tokens);
}

else
{

}
*/
	return 0;
}
