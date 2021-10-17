#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
using namespace std;


//======================= I/O related ==========================//
//-------- utility ------//
void getBaseName(string &in,string &out,char slash,char dot)
{
	int i,j;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	i++;
	for(j=len-1;j>=0;j--)
	{
		if(in[j]==dot)break;
	}
	if(j==-1)j=len;
	out=in.substr(i,j-i);
}
void getRootName(string &in,string &out,char slash)
{
	int i;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	if(i<=0)out=".";
	else out=in.substr(0,i);
}

//================= decompose '#' commented file into sub files =============//

void Parse_Species(string &in,string &out)
{
	//init
	out="";
	//process
	int i;
	int length=(int)in.length();
	int pos1=in.find("OS=");
	if(pos1>=0 && pos1<length)
	{
		int found=0;
		int pos2=-1;
		pos1+=3;
		for(i=pos1;i<length;i++)
		{
			if(in[i]=='=')
			{
				found=1;
				pos2=i;
				break;
			}
		}
		if(found==1)
		{
			int found2=0;
			int pos3=-1;
			for(i=pos2;i>=pos1;i--)
			{
				if(in[i]==' ')
				{
					found2=1;
					pos3=i-1;
					break;
				}
			}
			if(found2==1)
			{
				out=in.substr(pos1,pos3-pos1+1);
			}
		}
	}
}


void Decompose_Commented_File(string &input_file)
{
	//start
	ifstream fin;
	string buf,temp;
	//read
	fin.open(input_file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"input_file %s not found!\n",input_file.c_str());
		exit(-1);
	}
	//load
	int count=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		if(buf=="")continue;
		string out;
		Parse_Species(buf,out);
		if(out!="")printf("%s\n",out.c_str());
//		if(buf[0]=='>')printf("%s\n",buf.c_str());
//		count++;
	}
//	printf("count=%d\n",count);
}

//----------- main -------------//
int main(int argc,char **argv)
{
	//---- Decompose_Commented_File ----//
	{
		if(argc<2)
		{
			fprintf(stderr,"OutFile_Species <input_file> \n");
			exit(-1);
		}
		//argument
		string input_file=argv[1];
		//process
		Decompose_Commented_File(input_file);
		exit(0);		
	}
}
