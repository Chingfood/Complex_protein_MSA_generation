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

//=============== load MSA_SpecBloc =============//
//-> header example
/*
>tr|A0A0P0IJ73|A0A0P0IJ73_9CAUD Tail fiber protein OS=Acinetobacter phage IME-200 PE=4 SV=1 | 1751 / 719 -> 2.435327  ! + 1 1735582 1735582
*/
int Multi_FASTA_Input_SpecBloc(string &multi_fasta,
	vector <string> &nam_list, vector <string> &fasta_list,
	vector <string> &nam_orig, vector <vector <int> > &species_block,
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
	species_block.clear();
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
				//Todo:Need to change here if a3m file format is changed
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
				//-> get species block
				//--| first round
				vector <string> spec_str_out;
				int spec_bloc_num=Parse_Str_Str(buf,spec_str_out, '!');
				if(spec_bloc_num!=2)
				{
					fprintf(stderr,"MSA header spec_bloc %d format error !! %s \n",
						spec_bloc_num,buf.c_str());
					exit(-1);
				}
				//--| second round
				string spec_bloc_str=spec_str_out[1];
				vector <string> spec_str_out_fin;
				int spec_bloc_num_fin=Parse_Str_Str(spec_bloc_str,spec_str_out_fin, ' ');
				if(spec_bloc_num_fin!=5)
				{
					fprintf(stderr,"MSA header spec_bloc_fin %d format error !! %s \n",
						spec_bloc_num_fin,spec_bloc_str.c_str());
					exit(-1);
				}
				//--| third round
				vector <int> spec_str_out_fin_int (4,0);
				if(spec_str_out_fin[1]=="+")spec_str_out_fin_int[0]=1;
				else spec_str_out_fin_int[0]=0;
				for(int i=2;i<=4;i++)
				{
					int tmp_int=atoi(spec_str_out_fin[i].c_str());
					spec_str_out_fin_int[i-1]=tmp_int;
				}
				species_block.push_back(spec_str_out_fin_int);
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

//--------- parse species_block ---------//
//-> format
/*
+ 1 1735582 1735582

[note]:
1. +/- : indicates has sub-species or not
2.  1  : indicates species block id
3. 1735582 : sub-species tax_id
4. 1735582 : species tax_id
*/

int Parse_Species_Block(vector <vector <int> > &species_block,
	map<int, int > &spec_mapping,vector <int> &spec_taxid, vector <vector <int> > &spec_taxid_rec,
	map<int, int > &subs_mapping,vector <int> &subs_taxid, vector <vector <pair<int,int> > > &subs_taxid_rec)
{
	map<int, int >::iterator iter;
	//init spec
	int spec_count=0;
	spec_mapping.clear();
	spec_taxid.clear();
	spec_taxid_rec.clear();
	//init subs
	int subs_count=0;
	subs_mapping.clear();
	subs_taxid.clear();
	subs_taxid_rec.clear();
	//proc
	int i;
	int size=(int)species_block.size();
	int cur_pos=0;
	for(i=0;i<size;i++)
	{
		//-> record species
		int tax_id=species_block[i][3];
		iter = spec_mapping.find(tax_id);
		if(iter == spec_mapping.end()) //-> new species
		{
			spec_count++;
			spec_mapping.insert(map < int, int >::value_type(tax_id, spec_count));
			spec_taxid.push_back(tax_id);
			vector <int> tmp_rec;
			tmp_rec.push_back(i);
			spec_taxid_rec.push_back(tmp_rec);
			cur_pos=0;
		}
		else                           //-> old species
		{
			int key=spec_mapping[tax_id]-1;
			spec_taxid_rec[key].push_back(i);
			cur_pos++;
		}
		//-> record sub-species
		if(species_block[i][2]!=species_block[i][3])
		{
			int sub_id=species_block[i][2];
			iter = subs_mapping.find(tax_id);
			if(iter == subs_mapping.end()) //-> new species
			{
				subs_count++;
				subs_mapping.insert(map < int, int >::value_type(sub_id, subs_count));
				subs_taxid.push_back(sub_id);
				vector <pair<int,int> > tmp_rec;
				int first=spec_count-1;
				int second=cur_pos;
				tmp_rec.push_back(pair<int,int>(first,second));
				subs_taxid_rec.push_back(tmp_rec);
			}
			else                           //-> old species
			{
				int key=subs_mapping[tax_id]-1;
				int first=spec_count-1;
				int second=cur_pos;
				subs_taxid_rec[key].push_back(pair<int,int>(first,second));
			}
		}
	}
	//return
	return spec_count;
}

//========== concatenate two MSAs ==============//
int Concatenate_Two_MSAs(
	vector <vector <int> > &species_block1, vector <vector <int> > &spec_taxid_rec1,
	vector <vector <int> > &species_block2,
	map<int, int > &spec_mapping2, vector <vector <int> > &spec_taxid_rec2,
	map<int, int > &subs_mapping2, vector <vector <pair<int,int> > > &subs_taxid_rec2_,
	vector <pair<int,int> > &concatenate_result)
{
	//init
	int concatenate_num=0;
	map<int, int >::iterator iter;
	map<int, int >::iterator iter_;
	vector <vector <pair<int,int> > > subs_taxid_rec2=subs_taxid_rec2_;
	concatenate_result.clear();
	//set 1st MSA as pivot
	int spec_size1=(int)spec_taxid_rec1.size();
	for(int i=0;i<spec_size1;i++)
	{
		//-> init
		int size=(int)spec_taxid_rec1[i].size();
		if(size==0)continue;
		int init_pos=spec_taxid_rec1[i][0];
		int spec_id=species_block1[init_pos][3];
		iter = spec_mapping2.find(spec_id);
		if(iter == spec_mapping2.end())continue; //-> not found the same species in MSA2
		int init_key=spec_mapping2[spec_id]-1;
		int init_size=spec_taxid_rec2[init_key].size();
		vector <int> tmp_rec(init_size,1);
		vector <int> sel_rec(size,1);
		//-> search for sub-spec first
		for(int j=0;j<size;j++)
		{
			int pos=spec_taxid_rec1[i][j];
			if(species_block1[pos][2]!=species_block1[pos][3]) //contain sub-spec !!
			{
				int sub_spec=species_block1[pos][2];
				iter = subs_mapping2.find(sub_spec);
				if(iter != subs_mapping2.end()) //-> found the same sub-species in MSA2 !!!
				{
					//--> get position
					int key2=subs_mapping2[sub_spec]-1;
					int size2=subs_taxid_rec2[key2].size();
					//--> add to result
					if(size2>0)
					{
						int first=subs_taxid_rec2[key2][0].first;
						int second=subs_taxid_rec2[key2][0].second;
						int pos2=spec_taxid_rec2[first][second];
						concatenate_result.push_back(pair<int,int>(pos,pos2)); //pos from MSA1, and pos2 from MSA2
						concatenate_num++;
						//delete
						subs_taxid_rec2[key2].erase(subs_taxid_rec2[key2].begin());
						tmp_rec[second]=0;
						sel_rec[j]=0;
					}
				}
			}
		}
		//generate new species list
		vector <int> rel_rec;
		int rel_rec_size=0;
		for(int j=0;j<init_size;j++)
		{
			if(tmp_rec[j]==1)
			{
				rel_rec.push_back(spec_taxid_rec2[init_key][j]);
				rel_rec_size++;
			}
		}
		if(rel_rec_size==0)continue;
		//-> search for spec second
		for(int j=0;j<size;j++)
		{
			int pos=spec_taxid_rec1[i][j];
			if(sel_rec[j]==1)
			{
				int pos2=rel_rec[0];
				concatenate_result.push_back(pair<int,int>(pos,pos2));
				concatenate_num++;
				//delete
				tmp_rec.erase(tmp_rec.begin());
				rel_rec_size--;
				//check
				if(rel_rec_size==0)break;
			}
		}
	}
	//return 
	return concatenate_num;
}

//------------ main -------------//
int main(int argc, char** argv)
{
	//------- MSA_ConCat -----//
	{
		if(argc<4)
		{
			fprintf(stderr,"MSA_ConCat <a3m_specbloc_1> <a3m_specbloc_2> <a3m_concat>\n");
			fprintf(stderr,"[note]: pivot MSA is the first input a3m file \n");
			exit(-1);
		}
		string a3m_specbloc_1=argv[1];
		string a3m_specbloc_2=argv[2];
		string a3m_concat=argv[3];

		//------ step 1: load input two MSAs ------//
		//load MSA_1
		vector <string> nam_list_1;
		vector <string> fasta_list_1;
		vector <string> nam_orig_1;
		vector <vector <int> > species_block_1;
		string pivot_header_1;
		string pivot_fasta_1;
		int msa1_num=Multi_FASTA_Input_SpecBloc(a3m_specbloc_1,
			nam_list_1,fasta_list_1,nam_orig_1,species_block_1,
			pivot_header_1,pivot_fasta_1);
		printf("msa1_num=%d\n",msa1_num);
		//load MSA_2
		vector <string> nam_list_2;
		vector <string> fasta_list_2;
		vector <string> nam_orig_2;
		vector <vector <int> > species_block_2;
		string pivot_header_2;
		string pivot_fasta_2;
		int msa2_num=Multi_FASTA_Input_SpecBloc(a3m_specbloc_2,
			nam_list_2,fasta_list_2,nam_orig_2,species_block_2,
			pivot_header_2,pivot_fasta_2);
		printf("msa2_num=%d\n",msa2_num);

		//------ step 2: processs species block ------//
		//process species_block_1
		map<int, int > spec_mapping_1;
		vector <int> spec_taxid_1;
		vector <vector <int> > spec_taxid_rec_1;
		map<int, int > subs_mapping_1;
		vector <int> subs_taxid_1;
		vector <vector <pair<int,int> > > subs_taxid_rec_1;
		int spec1_num=Parse_Species_Block(species_block_1,
			spec_mapping_1,spec_taxid_1,spec_taxid_rec_1,
			subs_mapping_1,subs_taxid_1,subs_taxid_rec_1);
		printf("spec1_num=%d\n",spec1_num);
		//process species_block_2
		map<int, int > spec_mapping_2;
		vector <int> spec_taxid_2;
		vector <vector <int> > spec_taxid_rec_2;
		map<int, int > subs_mapping_2;
		vector <int> subs_taxid_2;
		vector <vector <pair<int,int> > > subs_taxid_rec_2;
		int spec2_num=Parse_Species_Block(species_block_2,
			spec_mapping_2,spec_taxid_2,spec_taxid_rec_2,
			subs_mapping_2,subs_taxid_2,subs_taxid_rec_2);
		printf("spec2_num=%d\n",spec2_num);

		//------ step 3: concatenate two MSAs ------//
		vector <pair<int,int> > concatenate_result;
		int concate_num=Concatenate_Two_MSAs(species_block_1,spec_taxid_rec_1,
			species_block_2,spec_mapping_2,spec_taxid_rec_2,subs_mapping_2,subs_taxid_rec_2,
			concatenate_result);
		printf("concate_num=%d\n",concate_num);
		FILE *fp=fopen(a3m_concat.c_str(),"wb");
		//-> output header
		fprintf(fp,"> %s & %s \n",pivot_header_1.c_str(),pivot_header_2.c_str());
		fprintf(fp,"%s%s\n",pivot_fasta_1.c_str(),pivot_fasta_2.c_str());
		//-> output others
		for(int i=0;i<concate_num;i++)
		{
			int first=concatenate_result[i].first;
			int second=concatenate_result[i].second;
			fprintf(fp,"> %s & %s \n",nam_orig_1[first].c_str(),nam_orig_2[second].c_str());
			fprintf(fp,"%s%s\n",fasta_list_1[first].c_str(),fasta_list_2[second].c_str());
		}
		fclose(fp);
		//exit
		exit(0);
	}
}

