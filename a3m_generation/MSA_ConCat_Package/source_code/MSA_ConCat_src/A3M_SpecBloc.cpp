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

//================= decompose '#' commented file into sub files =============//
//-> example
/*
>sp|Q8WZ42-12|TITIN_HUMAN Isoform 12 of Titin OS=Homo sapiens GN=TTN
*/

//Todo: Need to modify
void Parse_Species(string &in,string &out)
{
	//init
	out="";
	//process
	int i;
	int length=(int)in.length();
	int pos1=in.find("Tax=");
	if(pos1>=0 && pos1<length)
	{
		int found=0;
		int pos2=-1;
		pos1+=4;
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
		/*
		else  //-> if not found '=', then we shall try '|' which is created by Sheng Wang
		{
			int foundws=0;
			int pos2ws=-1;
			for(i=pos1;i<length;i++)
			{
				if(in[i]=='|')
				{
					foundws=1;
					pos2ws=i;
					break;
				}
			}
			if(foundws==1)
			{
				int found2ws=0;
				int pos3ws=-1;
				for(i=pos2ws;i>=pos1;i--)
				{
					if(in[i]==' ')
					{
						found2ws=1;
						pos3ws=i-1;
						break;
					}
				}
				if(found2ws==1)
				{
					out=in.substr(pos1,pos3ws-pos1+1);
				}
			}
		}
		*/
	}
}

//============= Load Multiple Sequence Alignment =========//
//-> for UniProt only
// header example
/*
>sp|Q8WZ42-12|TITIN_HUMAN Isoform 12 of Titin OS=Homo sapiens GN=TTN
*/
int Multi_FASTA_Input_UniProt_Spec(string &multi_fasta,
	vector <string> &nam_list, vector <string> &fasta_list,
	vector <string> &nam_orig, vector <string> &species_rec,
	string &pivot_header, string &pivot_fasta)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(multi_fasta.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"file %s not found!\n",multi_fasta.c_str());
		exit(-1);
	}
	//load
	int header=1;
	int first=1;
	int count=0;
	int number=0;
	string name;
	string seq;
	nam_list.clear();
	nam_orig.clear();
	species_rec.clear();
	fasta_list.clear();
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		if(buf=="")continue;
		if(buf.length()>=1 && buf[0]=='>')
		{
			if(first!=1)
			{
				//-> the following foramt is based on uniprot_human_nr

				//Todo: Need to change here to read the new Uniref format
				int start=-1;
				for(int i=0;i<(int)buf.length();i++)
				{
					if(buf[i]=='_')
					{
						start=i+1;
						break;
					}
				}
				//check
				if(start<0 || start>=(int)buf.length())
				{
					fprintf(stderr,"MSA header format error here !! %s \n",buf.c_str());
					exit(-1);
				}
				//-> get end
				int end=-1;
				for(int i=start;i<(int)buf.length();i++)
				{
					if(buf[i]==' ')
					{
						end=i-1;
						break;
					}
				}
				//check
				if(end<0 || end>=(int)buf.length())
				{
					fprintf(stderr,"MSA header format error here !! %s \n",buf.c_str());
					exit(-1);
				}
				//-> get name string
				name=buf.substr(start,end-start+1);
				nam_list.push_back(name);
				//-> get original name
				temp=buf.substr(1,buf.length()-1);
				nam_orig.push_back(temp);
				//-> get species
				string spec_str;
				Parse_Species(buf,spec_str); //Need to modify
				species_rec.push_back(spec_str);
				//count
				count++;
			}
			else
			{
				name=buf.substr(1,buf.length()-1);
				pivot_header=name;
			}

			//----- record sequence ----//
			if(first!=1)
			{
				if(header==1)
				{
					header=0;
					pivot_fasta=seq;
				}
				else
				{
					fasta_list.push_back(seq);
					number++;
				}
			}
			first=0;
			seq="";
		}
		else
		{
			if(first!=1)seq+=buf;
		}
	}
	//final
	if(first!=1)
	{
		if(header==1)
		{
			header=0;
			pivot_fasta=seq;
		}
		else
		{
			fasta_list.push_back(seq);
			number++;
		}
	}
	//check
	if(number!=count)
	{
		fprintf(stderr,"%s -> num %d != count %d \n",multi_fasta.c_str(),number,count);
		exit(-1);
	}
	return count;
}

//=================== load Species_Taxonomy_Tree ========//
//-> example
/*
Schistosoma mansoni     6183 (24)       6183(24) 6181(22) 31245(18) 6180(13) 6178(9) 6157(6) 33208(3)
Sepia officinalis       6610 (24)       6610(24) 6609(22) 6608(18) 551287(13) 6605(9) 6447(6) 33208(3)
Sus scrofa      9823 (24)       9823(24) 9822(22) 9821(18) 314145(12) 40674(9) 7711(6) 33208(3)
Ovis aries      9940 (24)       9940(24) 9935(22) 9895(18) 314145(12) 40674(9) 7711(6) 33208(3)
...

[note]:
each row has three columns, separated by '\t',
the first column is the full name of the species,
the second column is the lowest rank in the taxonomy tree
the third column is the detailed rank order in 7 levels.
*/

int Load_SpecTaxTree(string &spec_tax_tree_file,
	vector <string> &spec_name, map<string, int > &ws_mapping,
	vector <int> &bot_digit, vector <vector <int> > &spec_tax_tree)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(spec_tax_tree_file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"file %s not found!\n",spec_tax_tree_file.c_str());
		exit(-1);
	}
	//load
	map<string, int >::iterator iter;
	ws_mapping.clear();
	spec_name.clear();
	bot_digit.clear();
	spec_tax_tree.clear();
	int retv;
	int count=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		//---------- get data process --------//
		//-> get whole row input
		vector <string> init_out;
		retv=Parse_Str_Str(buf,init_out, '\t');
		if(retv!=3)
		{
			fprintf(stderr,"ERROR -> each row shall contain only 3 columns separated with '\t' -> %s \n",buf.c_str());
			exit(-1);
		}
		//-> get bottom digit 
		vector <string> bot_out;
		retv=Parse_Str_Str(init_out[1],bot_out, ' ');
		if(retv!=2)
		{
			fprintf(stderr,"ERROR -> each bottom_digit shall contain only 2 columns separated with ' ' -> %s \n",init_out[1].c_str());
			exit(-1);
		}
		//-> get 7-level taxonomy tree
		vector <string> tree_out;
		retv=Parse_Str_Str(init_out[2],tree_out, ' ');
		if(retv!=7)
		{
			fprintf(stderr,"ERROR -> each bottom_digit shall contain only 7 columns separated with ' ' -> %s \n",init_out[2].c_str());
			exit(-1);
		}
		//---------- record data process --------//
		//-> record species name
		count++;
		iter = ws_mapping.find(init_out[0]);
		if(iter != ws_mapping.end())
		{
			fprintf(stderr,"duplicated mapping !! %s \n",init_out[0].c_str());
			exit(-1);
		}
		ws_mapping.insert(map < string, int >::value_type(init_out[0], count));
		spec_name.push_back(init_out[0]);
		//-> record bottom digit
		retv=atoi(bot_out[0].c_str());
		bot_digit.push_back(retv);
		//-> record 7-level taxonomy tree
		vector <int> seven_tmp;
		for(int i=0;i<7;i++)
		{
			vector <string> tree_out_tmp;
			retv=Parse_Str_Str(tree_out[i],tree_out_tmp, '(');
			if(retv!=2)
			{
				fprintf(stderr,"ERROR -> each tree_out[%d] shall contain only 2 columns separated with '(' -> %s \n",
					i,tree_out[i].c_str());
				exit(-1);
			}
			retv=atoi(tree_out_tmp[0].c_str());
			seven_tmp.push_back(retv);
		}
		spec_tax_tree.push_back(seven_tmp);
		//------ printf -------//
		if(count%5000==0)fprintf(stderr,"num->%d\r",count);
	}
	//return
	return count;
}


//---------- bad species check single ----//
void Bad_Species_Check_Single(string &buf,
	map<string, int > &nam_mapping, string &out_spec)
{
	out_spec="";
	map<string, int>::iterator iter;
	//-> load species strings
	vector <string> out_str;
	int str_num=Parse_Str_Str(buf,out_str,' ');
	//-> record the least species_rank
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
		out_spec=species_str;
		return;
	}
}

//================= output all species into Species_Assembly_Block ==============//
int Species_To_Block(vector <string> &species_rec, vector <pair<int,int> > &species_taxid,
	map<int, int > &block_mapping, vector <int> &block_taxid, vector <vector <int> > &block_occur,
	map<string, int > &ws_mapping, vector <int> &bot_digit, vector <vector <int> > &spec_tax_tree)
{
	//init
	int spec_count=0;
	map<int, int >::iterator iter_;
	block_mapping.clear();
	block_taxid.clear();
	block_occur.clear();
	//proc
	int i;
	int size=(int)species_rec.size();
	species_taxid.clear();
	map<string, int >::iterator iter;
	for(i=0;i<size;i++)
	{
		//check species
		string spec=species_rec[i];
		if(spec=="")
		{
			fprintf(stderr,"cur %d -> species name blank !! \n",i);
			species_taxid.push_back(pair<int,int>(-1,-1));
			continue;
		}
		iter = ws_mapping.find(spec);
		if(iter == ws_mapping.end())
		{
			fprintf(stderr,"cur %d -> species name %s not found !! \n",i,spec.c_str());
			//-- try to save !!! --//
			string out_spec;
			Bad_Species_Check_Single(spec,ws_mapping,out_spec);
			if(out_spec=="")
			{
				species_taxid.push_back(pair<int,int>(-1,-1));
				continue;
			}
			else
			{
				spec=out_spec;
				fprintf(stderr,"sav %d -> species name %s saved !! \n",i,spec.c_str());
			}
		}
		//push bot and spec
		int key=ws_mapping[spec]-1;
		int bot=bot_digit[key];
		int tax=spec_tax_tree[key][0];
		species_taxid.push_back(pair<int,int>(bot,tax));
		//record species
		iter_ = block_mapping.find(tax);
		if(iter_ == block_mapping.end()) //-> new species
		{
			spec_count++;
			block_mapping.insert(map < int, int >::value_type(tax, spec_count));
			block_taxid.push_back(tax);
			vector <int> tmp_rec;
			tmp_rec.push_back(i);
			block_occur.push_back(tmp_rec);
		}
		else                             //-> old species
		{
			int key_=block_mapping[tax]-1;
			block_occur[key_].push_back(i);
		}
	}
	//return
	return spec_count;
}

//=============== from A3M_NoGap to A3M_AssemBlock =========//
//-> [note]: the input A3M format must be processed by 'A3M_NoGap', say
//           1) the first input should be the query sequence
//           2) all other sequences should be sorted in descending order of BLOS62 score to the query sequence
void A3M_To_AssemBlock(string &a3m_input, string &spec_tax_tree_file, string &output)
{
	//input raw A3M
	vector <string> nam_list;
	vector <string> fasta_list;
	vector <string> nam_orig;
	vector <string> species_rec;
	string pivot_header;
	string pivot_fasta;
	int msa_num=Multi_FASTA_Input_UniProt_Spec(a3m_input,nam_list,fasta_list,nam_orig,species_rec,
		pivot_header,pivot_fasta);

fprintf(stderr,"msa_num=%d\n",msa_num);

	//input species_taxomony_tree
	vector <string> spec_name;
	map<string, int > spec_mapping;
	vector <int> bot_digit;
	vector <vector <int> > spec_tax_tree;
	int spec_num=Load_SpecTaxTree(spec_tax_tree_file,spec_name,spec_mapping,bot_digit,spec_tax_tree);

fprintf(stderr,"spec_num=%d\n",spec_num);

	//A3M to assembly block
	vector <pair<int,int> > species_taxid;  //-> <bottom_taxid,species_taxid>
	map<int, int > block_mapping;
	vector <int> block_taxid;
	vector <vector <int> > block_occur;
	int spec_count=Species_To_Block(species_rec,species_taxid,block_mapping,block_taxid,block_occur,
		spec_mapping,bot_digit,spec_tax_tree);

fprintf(stderr,"spec_count=%d\n",spec_count);

	//output assebmly block
	//-> output a new A3M file
	FILE *fp=fopen(output.c_str(),"wb");
	fprintf(fp,">%s\n",pivot_header.c_str());
	fprintf(fp,"%s\n",pivot_fasta.c_str());
	for(int i=0;i<spec_count;i++)
	{
		int spec_id=block_taxid[i];
		int spec_size=(int)block_occur[i].size();
		//bad species
		if(spec_id<0)
		{
			printf("null block %5d -> number %5d -> %d !! \n",i+1,spec_size,spec_id);
			continue;
		}
		//normal species
		for(int k=0;k<spec_size;k++)
		{
			int pos=block_occur[i][k];
			int pos_bot_id=species_taxid[pos].first;
			int pos_spec_id=species_taxid[pos].second;
			//-> output the same species
			if(pos_bot_id==spec_id)
			{
				fprintf(fp,">%s ! + %d %d %d \n",nam_orig[pos].c_str(),i+1,pos_bot_id,spec_id);
				fprintf(fp,"%s\n",fasta_list[pos].c_str());
			}
			else
			{
				fprintf(fp,">%s ! - %d %d %d \n",nam_orig[pos].c_str(),i+1,pos_bot_id,spec_id);
				fprintf(fp,"%s\n",fasta_list[pos].c_str());
			}
		}
		//---- printf output species assembly block ---//
		printf("block %5d -> number %5d -> species %d \n",i+1,spec_size,spec_id);
	}
}

//------------ main -------------//
int main(int argc, char** argv)
{
	//------- A3M_To_AssemBlock -----//
	{
		if(argc<4)
		{
			fprintf(stderr,"A3M_SpecBloc <a3m_nogap> <spec_tax_tree_file> <a3m_specbloc>\n");
			exit(-1);
		}
		string a3m_input=argv[1];
		string spec_tax_tree_file=argv[2];
		string a3m_assemblock=argv[3];
		//process
		A3M_To_AssemBlock(a3m_input,spec_tax_tree_file,a3m_assemblock);
		//exit
		exit(0);
	}
}

