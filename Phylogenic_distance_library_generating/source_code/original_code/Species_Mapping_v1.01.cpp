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

//=================== upper and lower case ====================//
//----------upper_case-----------//
void toUpperCase(char *buffer)
{
	for(int i=0;i<(int)strlen(buffer);i++)
	if(buffer[i]>=97 && buffer[i]<=122) buffer[i]-=32;
}
void toUpperCase(string &buffer)
{
	for(int i=0;i<(int)buffer.length();i++)
	if(buffer[i]>=97 && buffer[i]<=122) buffer[i]-=32;
}
//----------lower_case-----------//
void toLowerCase(char *buffer)
{
	for(int i=0;i<(int)strlen(buffer);i++)
	if(buffer[i]>=65 && buffer[i]<=90) buffer[i]+=32;
}
void toLowerCase(string &buffer)
{
	for(int i=0;i<(int)buffer.length();i++)
	if(buffer[i]>=65 && buffer[i]<=90) buffer[i]+=32;
}

//--------- Parse_Str_Str --------//
int Parse_Str_Str(string &in,vector <string> &out, char separator)
{
	istringstream www(in);
	out.clear();
	int count=0;
	for(;;)
	{
		string buf;
		if(!getline(www,buf,separator))break;
		out.push_back(buf);
		count++;
	}
	return count;
}


//=================== load species mapping ==============//
//-> main file: names.dmp_species_taxid
/*
all     1
root    1
Bacteria        2
Monera  2
Procaryotae     2
Prokaryota      2
Prokaryotae     2
...
*/

int Load_Species_TaxID(string &in_file, 
	map<string, int > &nam_mapping,
	map<int, int > &taxid_mapping)
{
	map<string, int>::iterator iter;
	map<int, int>::iterator iter_;
	//--- list for mapping ---//
	taxid_mapping.clear();
	nam_mapping.clear();
	ifstream fin;
	string buf,temp,name;
	//read
	fin.open(in_file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"in_file %s not found!!\n",in_file.c_str());
		exit(-1);
	}
	//proc
	int count=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		istringstream www(buf);
		//-> read
		string species;
		int taxid;
		getline(www,species,'\t');
		getline(www,temp,'\t');
		taxid=atoi(temp.c_str());
		//-> map
		iter = nam_mapping.find(species);
		if(iter != nam_mapping.end())
		{
			fprintf(stderr,"duplicated species string\t%s\r",species.c_str());
//			continue;
		}
		else
		{
			nam_mapping.insert(map < string, int >::value_type(species, taxid));
			//-> lowercase map
			string lowerspec=species;
			toLowerCase(lowerspec);
			nam_mapping.insert(map < string, int >::value_type(lowerspec, taxid));
		}
		//-> non-redundant taxid
		iter_ = taxid_mapping.find(taxid);
		if(iter_ != taxid_mapping.end())continue;
		count++;
		taxid_mapping.insert(map < int, int >::value_type(taxid, count));
		//printf
		fprintf(stderr,"cur-> %d\r",count);
	}
	fprintf(stderr,"load Species_TaxID file %s\n",in_file.c_str());
	//return
	return count;
}

//============= load species node =============//
//-> traditional taxonomy category
/*
Kingdom -> Jie
Phylum  -> Men
Class   -> Gang
Order   -> Mu
Family  -> Ke
Genus   -> Shu
Species -> Zhong
*/

//-> current available taxonomy category in 'nodes.dmp'
/*
//              rank_id    number
no              -> 1   -> (201582)  
superkingdom    -> 2   -> (5)
kingdom         -> 3   -> (3)
subkingdom      -> 4   -> (1)
superphylum     -> 5   -> (2)
phylum          -> 6   -> (214)
subphylum       -> 7   -> (25)
superclass      -> 8   -> (5)
class           -> 9   -> (304)
subclass        -> 10  -> (130)
infraclass      -> 11  -> (16)
superorder      -> 12  -> (49)
order           -> 13  -> (1390)       => level 3
suborder        -> 14  -> (313)
infraorder      -> 15  -> (96)
parvorder       -> 16  -> (10)
superfamily     -> 17  -> (796)
family          -> 18  -> (8388)       => level 2
subfamily       -> 19  -> (2619)
tribe           -> 20  -> (1886)
subtribe        -> 21  -> (465)
genus           -> 22  -> (79470)      => level 1
subgenus        -> 23  -> (1238)
species         -> 24  -> (1194362)    => level 0
subspecies      -> 25  -> (19778)
varietas        -> 26  -> (6531)
forma           -> 27  -> (455)
*/

//---- Int_To_Rank ----//
void Int_To_Rank(int rank,string &rank_str)
{
	rank_str="";
	if(rank==1) rank_str="no"          ;
	if(rank==2) rank_str="superkingdom";
	if(rank==3) rank_str="kingdom"     ;
	if(rank==4) rank_str="subkingdom"  ;
	if(rank==5) rank_str="superphylum" ;
	if(rank==6) rank_str="phylum"      ;
	if(rank==7) rank_str="subphylum"   ;
	if(rank==8) rank_str="superclass"  ;
	if(rank==9) rank_str="class"       ;
	if(rank==10)rank_str="subclass"    ;
	if(rank==11)rank_str="infraclass"  ;
	if(rank==12)rank_str="superorder"  ;
	if(rank==13)rank_str="order"       ;
	if(rank==14)rank_str="suborder"    ;
	if(rank==15)rank_str="infraorder"  ;
	if(rank==16)rank_str="parvorder"   ;
	if(rank==17)rank_str="superfamily" ;
	if(rank==18)rank_str="family"      ;
	if(rank==19)rank_str="subfamily"   ;
	if(rank==20)rank_str="tribe"       ;
	if(rank==21)rank_str="subtribe"    ;
	if(rank==22)rank_str="genus"       ;
	if(rank==23)rank_str="subgenus"    ;
	if(rank==24)rank_str="species"     ;
	if(rank==25)rank_str="subspecies"  ;
	if(rank==26)rank_str="varietas"    ;
	if(rank==27)rank_str="forma"       ;
	return;
}

int Species_NonRed_Number;
int *Species_Taxonomy_ID;
int *Species_Taxonomy_Rank;
int *Species_Taxonomy_Parent;

//------ load TaxID node ---------//
//-> main file: nodes.dmp_taxid_prev_rank_better
/*
1       1       1       no
2       131567  2       superkingdom
6       335928  22      genus
7       6       24      species
9       32199   24      species
10      1706371 22      genus
11      1707    24      species
13      203488  22      genus
14      13      24      species
16      32011   22      genus
....
*/
int Load_Species_NodeRank(string &in_file,
	map<int, int > &taxid_mapping,
	int *Species_Taxonomy_ID,
	int *Species_Taxonomy_Rank,
	int *Species_Taxonomy_Parent)
{
	map<int, int>::iterator iter;
	//--- list for mapping ---//
	ifstream fin;
	string buf,temp,name;
	//read
	fin.open(in_file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"in_file %s not found!!\n",in_file.c_str());
		exit(-1);
	}
	//proc
	int count=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		istringstream www(buf);
		//-> read
		int taxid_node;
		int taxid_parent;
		int taxid_rank;
		getline(www,temp,'\t');
		taxid_node=atoi(temp.c_str());
		getline(www,temp,'\t');
		taxid_parent=atoi(temp.c_str());
		getline(www,temp,'\t');
		taxid_rank=atoi(temp.c_str());
		//-> map1
		iter = taxid_mapping.find(taxid_node);
		if(iter == taxid_mapping.end())
		{
			fprintf(stderr,"taxid_node %d not found in taxid_mapping \n",taxid_node);
			exit(-1);
		}
		int taxid_node_=taxid_mapping[taxid_node]-1;
		//-> map2
		iter = taxid_mapping.find(taxid_parent);
		if(iter == taxid_mapping.end())
		{
			fprintf(stderr,"taxid_parent %d not found in taxid_mapping \n",taxid_parent);
			exit(-1);
		}
		int taxid_parent_=taxid_mapping[taxid_parent]-1;
		//-> record
		Species_Taxonomy_ID[taxid_node_]=taxid_node;        //-> record the Real-World TaxID
		Species_Taxonomy_Rank[taxid_node_]=taxid_rank;      //-> record the Taxonomy Rank
		Species_Taxonomy_Parent[taxid_node_]=taxid_parent_; //-> record the parent node
		//-> count
		count++;
		fprintf(stderr,"node->%d\r",count);
	}
	fprintf(stderr,"load Species_NodeRank file %s\n",in_file.c_str());
	//return
	return count;
}

//========== bad species check ==========//
//-> example of Bad Species:
/*
Chlamydia trachomatis serovar L2 (strain 434/Bu / ATCC VR-902B)
Human herpesvirus 1 (strain KOS)
Human adenovirus E serotype 4
Human adenovirus B serotype 7
Propionibacterium freudenreichii subsp. shermanii (strain ATCC 9614 / CIP 103027 / CIRM-BIA1)
Influenza A virus (strain A/Swine/Hong Kong/81/1978 H3N2)
Deinococcus gobiensis (strain DSM 21396 / JCM 16679 / CGMCC 1.7299 / I-0)
Ficus variegata Blume, 1825
Ziziphus jujuba witches'-broom phytoplasma
Fusarium solani subsp. pisi
...
*/
//-> our strategy is starting from the first string block that hit thhe mapping to the last, 
//   and select the lowest species rank
void Bad_Species_Check(string &in_file,
	map<string, int > &nam_mapping,
	map<int, int > &taxid_mapping,
	int *Species_Taxonomy_Rank)
{
	map<string, int>::iterator iter;
	map<int, int>::iterator iter_;
	//--- list for mapping ---//
	ifstream fin;
	string buf,temp,name;
	//read
	fin.open(in_file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"in_file %s not found!!\n",in_file.c_str());
		exit(-1);
	}
	//proc
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		//-> load species strings
		vector <string> out_str;
		int str_num=Parse_Str_Str(buf,out_str,' ');
		//-> record the least species_rank
		int least_species_rank=0;
		int least_species_taxid=-1;
		string species_str="";
		for(int i=0;i<str_num;i++)
		{
			species_str=species_str+out_str[i];

//printf("species_str=%s\n",species_str.c_str());

			//-> map1
			iter = nam_mapping.find(species_str);
			if(iter == nam_mapping.end())
			{
				species_str=species_str+" ";
				continue;
			}
			int ori_num=nam_mapping[species_str];

//printf("ori_num=%d\n",ori_num);

			//-> map2
			int cur_num=taxid_mapping[ori_num]-1;
			int cur_rank=Species_Taxonomy_Rank[cur_num];

//printf("cur_rank=%d\n",cur_rank);

			if(cur_rank>least_species_rank)
			{
				least_species_rank=cur_rank;
				least_species_taxid=ori_num;
			}

			//final add
			species_str=species_str+" ";
		}
		//printf
		if(least_species_taxid==-1)
		{
			fprintf(stderr,"failed for '%s' \n",buf.c_str());
			continue;
		}
		printf("%s\t%d\t%d\n",buf.c_str(),least_species_taxid,least_species_rank);
	}
}

//=========== main process for BAD SPECIES ============//
void Main_Process_For_Bad_Species(
	string &species_taxid_file, 
	string &taxid_rank_file, 
	string &bad_speices_file)
{
	//load species_taxid_file
	map<string, int > nam_mapping;
	map<int, int > taxid_mapping;
	Species_NonRed_Number=Load_Species_TaxID(species_taxid_file, 
		nam_mapping,taxid_mapping);
	//create main data structure
	Species_Taxonomy_ID=new int[Species_NonRed_Number];
	Species_Taxonomy_Rank=new int[Species_NonRed_Number];
	Species_Taxonomy_Parent=new int[Species_NonRed_Number];
	//load taxid_rank_file
	Load_Species_NodeRank(taxid_rank_file,taxid_mapping,
		Species_Taxonomy_ID,Species_Taxonomy_Rank,Species_Taxonomy_Parent);
	//process bad_speices_file
	Bad_Species_Check(bad_speices_file,
		nam_mapping,taxid_mapping,Species_Taxonomy_Rank);
}



//------------ main -------------//
int main(int argc, char** argv)
{
	//------- Proc_Bad_Species -----//
	{
		if(argc<4)
		{
			fprintf(stderr,"Proc_Bad_Species <species_taxid_file> <taxid_rank_file> <bad_speices_file>\n");
			exit(-1);
		}
		string species_taxid_file=argv[1];
		string taxid_rank_file=argv[2];
		string bad_speices_file=argv[3];
		//process
		Main_Process_For_Bad_Species(species_taxid_file,taxid_rank_file,bad_speices_file);
		//exit
		exit(0);
	}
}

