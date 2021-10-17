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

//--------- remove specific char -----//
void Remove_Char(string &in, char c, string &out)
{
	int len=(int)in.length();
	out="";
	for(int i=0;i<len;i++)
	{
		if(in[i]==c)continue;
		out.push_back(in[i]);
	}
}

//--------- remove specific char -----//
void Replace_Char(string &in, char c, char p, string &out)
{
	int len=(int)in.length();
	out="";
	for(int i=0;i<len;i++)
	{
		if(in[i]==c)out.push_back(p);
		else out.push_back(in[i]);
	}
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
kingdom         -> 3   -> (3)          => level 6
subkingdom      -> 4   -> (1)
superphylum     -> 5   -> (2)
phylum          -> 6   -> (214)        => level 5
subphylum       -> 7   -> (25)
superclass      -> 8   -> (5)
class           -> 9   -> (304)        => level 4
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



//========== specific bad species =======//
/*
// bad                              correct
Human astrovirus-X               -> Human astrovirus X
'Asterotremella humicola'        -> Asterotremella humicola
[Clostridium] sordellii          -> Clostridium sordellii
Human T-cell leukemia virus 1    -> Human T-cell leukemia virus type I
Human T-cell leukemia virus 2    -> Human T-cell leukemia virus type II
Maedi visna virus                -> Maedi visna virus MVV
Bovine herpesvirus 1.1           -> Bovine herpesvirus type 1.1
Mouse intracisternal a-particle  -> Mouse intracisternal A-particle
Ross river virus (strain 213970) -> Ross river virus (STRAIN 213970)
*/

void Specific_Bad_Species_Proc(string &in, string &out)
{
	out=in;
	vector <string> out_str;
	int str_num=Parse_Str_Str(in,out_str,' ');
	int in_len=(int)in.length();
	int inin_len=(int)out_str[0].length();
	//-> 1. Human astrovirus-X
	vector <string> human_astrovirus;
	Parse_Str_Str(out_str[1],human_astrovirus,'-');
	if(out_str[0]=="Human" && human_astrovirus[0]=="astrovirus")
	{
		Replace_Char(in,'-',' ',out);
		return;
	}
	//-> 2. 'Asterotremella humicola'
	if(in[0]=='\'' && in[in_len-1]=='\'')
	{
		Remove_Char(in,'\'',out);
		return;
	}
	//-> 3. [Clostridium] sordellii
	if(out_str[0][0]=='[' && out_str[0][inin_len-1]==']')
	{
		Remove_Char(in,'[',out);
		string temp=out;
		Remove_Char(temp,']',out);
		return;
	}
	//-> 4. Human T-cell leukemia virus X
	if(str_num>=5)
	{
		if(out_str[0]=="Human" && out_str[1]=="T-cell" && out_str[2]=="leukemia" && out_str[3]=="virus")
		{
			if(out_str[4]=="2")out="Human T-cell leukemia virus type II ";
			if(out_str[4]=="1")out="Human T-cell leukemia virus type I ";
			for(int i=5;i<str_num;i++)
			{
				if(i<str_num-1)out+=out_str[i]+" ";
				else out+=out_str[i];
			}
			return;
		}
	}
	//-> 5. Maedi visna virus
	if(str_num>=3)
	{
		if(out_str[0]=="Maedi" && out_str[1]=="visna" && out_str[2]=="virus")
		{
			out="Maedi visna virus MVV ";
			for(int i=3;i<str_num;i++)
			{
				if(i<str_num-1)out+=out_str[i]+" ";
				else out+=out_str[i];
			}
			return;
		}
	}
	//-> 6. Bovine herpesvirus 1.1
	if(str_num>=3)
	{
		if(out_str[0]=="Bovine" && out_str[1]=="herpesvirus")
		{
			out="Bovine herpesvirus type "+out_str[2]+" ";
			for(int i=3;i<str_num;i++)
			{
				if(i<str_num-1)out+=out_str[i]+" ";
				else out+=out_str[i];
			}
			return;
		}
	}
	//-> 7. Mouse intracisternal a-particle
	if(str_num>=3)
	{
		if(out_str[0]=="Mouse" && out_str[1]=="intracisternal" && out_str[2]=="a-particle")
		{
			out="Mouse intracisternal A-particle ";
			for(int i=3;i<str_num;i++)
			{
				if(i<str_num-1)out+=out_str[i]+" ";
				else out+=out_str[i];
			}
			return;
		}
	}
	//-> 8. Ross river virus (strain
	if(str_num>=4)
	{
		if(out_str[0]=="Ross" && out_str[1]=="river" && out_str[2]=="virus" && out_str[3]=="(strain")
		{
			out="Ross river virus (STRAIN ";
			for(int i=4;i<str_num;i++)
			{
				if(i<str_num-1)out+=out_str[i]+" ";
				else out+=out_str[i];
			}
			return;
		}
	}
}


//---------- bad species check single ----//
void Bad_Species_Check_Single(string &buf_,
	map<string, int > &nam_mapping,
	map<int, int > &taxid_mapping,
	int &least_species_rank,
	int &least_species_taxid)
{
	map<string, int>::iterator iter;
	map<int, int>::iterator iter_;
	//-> to lower string
	string buf=buf_;
	toLowerCase(buf);
	//-> load species strings
	vector <string> out_str;
	int str_num=Parse_Str_Str(buf,out_str,' ');
	//-> record the least species_rank
	least_species_rank=0;
	least_species_taxid=-1;
	string species_str=out_str[0]+" ";
	for(int i=1;i<str_num;i++)
	{
		species_str=species_str+out_str[i];
		//-> map1
		iter = nam_mapping.find(species_str);
		if(iter == nam_mapping.end())
		{
			species_str=species_str+" ";
			continue;
		}
		int ori_num=nam_mapping[species_str];
		//-> map2
		int cur_num=taxid_mapping[ori_num]-1;
		int cur_rank=Species_Taxonomy_Rank[cur_num];
		if(cur_rank>least_species_rank)
		{
			least_species_rank=cur_rank;
			least_species_taxid=ori_num;
		}
		//final add
		species_str=species_str+" ";
	}
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
	int least_species_rank;
	int least_species_taxid;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		//-> first round check
		Bad_Species_Check_Single(buf,nam_mapping,taxid_mapping,
			least_species_rank,least_species_taxid);
		//printf
		if(least_species_taxid==-1)
		{
			//--- try to save !!! ----//start
			string species_str;
			Specific_Bad_Species_Proc(buf,species_str);
			Bad_Species_Check_Single(species_str,nam_mapping,taxid_mapping,
				least_species_rank,least_species_taxid);
			//--- try to save !!! ----//end
			if(least_species_taxid==-1)
			{
				fprintf(stderr,"failed for '%s' \n",buf.c_str());
				continue;
			}
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

////////////////// ================= Main Process of Contacenate Two MSAs ============= ////////////
//-> data structure
/*

|  species_1  |  species_2  |  species_3  |  species_4  | ... |  species_N  |
<oriID,taxID>_1
<oriID,taxID>_2
...
<oriID,taxID>_k

where,
oriID is the original index from the input MSA
taxID is the canonical taxonomy ID

*/



//=================== given a species string, return Spec (0), Genus (1), and Family (2), etc ========//
void Return_Spec_Genu_Fami(string &species_str, 
	int &species_taxid, int &species_taxid_r,
	int &Spec, int &Genu, int &Fami, int &Orde, int &Clas, int &Phyl, int &King, 
	int &Spec_r, int &Genu_r, int &Fami_r, int &Orde_r, int &Clas_r, int &Phyl_r, int &King_r,
	map<string, int > &nam_mapping, map<int, int > &taxid_mapping,
	int *Species_Taxonomy_ID,
	int *Species_Taxonomy_Rank,
	int *Species_Taxonomy_Parent)
{
	//init
	Spec=-1;
	Genu=-1;
	Fami=-1;
	Orde=-1;
	Clas=-1;
	Phyl=-1;
	King=-1;
	Spec_r=-1;
	Genu_r=-1;
	Fami_r=-1;
	Orde_r=-1;
	Clas_r=-1;
	Phyl_r=-1;
	King_r=-1;
	species_taxid=-1;
	species_taxid_r=-1;
	//map
	map<string, int>::iterator iter;
	map<int, int>::iterator iter_;
	iter = nam_mapping.find(species_str);
	if(iter == nam_mapping.end())
	{
		fprintf(stderr,"WARNING -> failed to map the species_str '%s' to a taxid \n",
			species_str.c_str());
		return;
	}
	int ori_num=nam_mapping[species_str];
	int cur_num=taxid_mapping[ori_num]-1;
	int cur_rank=Species_Taxonomy_Rank[cur_num];
	species_taxid=Species_Taxonomy_ID[cur_num];
	species_taxid_r=cur_rank;
	//trace back on the Evolutionary Tree
	int curr_node=cur_num;
	int next_node=-1;
	int curr_rank=cur_rank;
	int next_rank=-1;
	for(int i=0;i<30;i++)
	{
		//check break
		if(curr_node==0)break;
		//record
		if(curr_rank<=27 &&curr_rank>=24)
		{
			int Spec_=Species_Taxonomy_ID[curr_node];
			int Spec_r_=curr_rank;
			if(Spec_r!=24)
			{
				Spec_r=Spec_r_;
				Spec=Spec_;
			}
		}
		if(curr_rank<=23 && curr_rank>=20)
		{
			int Genu_=Species_Taxonomy_ID[curr_node];
			int Genu_r_=curr_rank;
			if(Genu_r!=22)
			{
				Genu_r=Genu_r_;
				Genu=Genu_;
			}
		}
		if(curr_rank<=19 && curr_rank>=17)
		{
			int Fami_=Species_Taxonomy_ID[curr_node];
			int Fami_r_=curr_rank;
			if(Fami_r!=18)
			{
				Fami_r=Fami_r_;
				Fami=Fami_;
			}
		}
		if(curr_rank<=16 && curr_rank>=12)
		{
			int Orde_=Species_Taxonomy_ID[curr_node];
			int Orde_r_=curr_rank;
			if(Orde_r!=13)
			{
				Orde_r=Orde_r_;
				Orde=Orde_;
			}
		}
		if(curr_rank<=11 && curr_rank>=8)
		{
			int Clas_=Species_Taxonomy_ID[curr_node];
			int Clas_r_=curr_rank;
			if(Clas_r!=9)
			{
				Clas_r=Clas_r_;
				Clas=Clas_;
			}
		}
		if(curr_rank<=7  && curr_rank>=5)
		{
			int Phyl_=Species_Taxonomy_ID[curr_node];
			int Phyl_r_=curr_rank;
			if(Phyl_r!=6)
			{
				Phyl_r=Phyl_r_;
				Phyl=Phyl_;
			}
		}
		if(curr_rank<=4  && curr_rank>=2)
		{
			int King_=Species_Taxonomy_ID[curr_node];
			int King_r_=curr_rank;
			if(King_r!=3)
			{
				King_r=King_r_;
				King=King_;
			}
		}
		//get next node
		next_node=Species_Taxonomy_Parent[curr_node];
		next_rank=Species_Taxonomy_Rank[next_node];
		//check rank
		if(next_rank!=1 && curr_rank!=1)
		{
			if(next_rank>curr_rank)
			{
				fprintf(stderr,"WARNING %s -> next_rank %d larger or equal to curr_rank %d \n",
					species_str.c_str(),next_rank,curr_rank);
			}
		}
		//update
		curr_node=next_node;
		curr_rank=next_rank;
	}
}

//=========== Main Process (trace species) ===========//
void Main_Process_For_Trace_Species(
	string &species_taxid_file, 
	string &taxid_rank_file, 
	string &trace_speices_file)
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
	//process trace_speices_file
	string in_file=trace_speices_file;
	ifstream fin;
	string buf;
	fin.open(in_file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"in_file %s not found!!\n",in_file.c_str());
		exit(-1);
	}
	int Spec;
	int Genu;
	int Fami;
	int Orde;
	int Clas;
	int Phyl;
	int King;
	int Spec_r;
	int Genu_r;
	int Fami_r;
	int Orde_r;
	int Clas_r;
	int Phyl_r;
	int King_r;
	int species_taxid;
	int species_taxid_r;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		Return_Spec_Genu_Fami(buf,species_taxid,species_taxid_r,
			Spec,Genu,Fami,Orde,Clas,Phyl,King,
			Spec_r,Genu_r,Fami_r,Orde_r,Clas_r,Phyl_r,King_r,
			nam_mapping,taxid_mapping,
			Species_Taxonomy_ID,Species_Taxonomy_Rank,Species_Taxonomy_Parent);
		printf("%s\t%d (%d)\t%d(%d) %d(%d) %d(%d) %d(%d) %d(%d) %d(%d) %d(%d)\n",
			buf.c_str(),species_taxid,species_taxid_r,
			Spec,Spec_r,Genu,Genu_r,Fami,Fami_r,Orde,Orde_r,Clas,Clas_r,Phyl,Phyl_r,King,King_r);
	}
}



//------------ main -------------//
int main(int argc, char** argv)
{
	//------- Trace_Species -----//
	{
		if(argc<4)
		{
			fprintf(stderr,"Trace_Species <species_taxid_file> <taxid_rank_file> <trace_speices_file>\n");
			exit(-1);
		}
		string species_taxid_file=argv[1];
		string taxid_rank_file=argv[2];
		string trace_speices_file=argv[3];
		//process
		Main_Process_For_Trace_Species(species_taxid_file,taxid_rank_file,trace_speices_file);
		//exit
		exit(0);
	}
}

