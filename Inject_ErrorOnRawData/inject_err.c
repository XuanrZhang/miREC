#include <cstring>
#include <string>
#include <fstream>
#include <algorithm>
#include <vector>
#include <iostream>
#include <unistd.h>
//AIM: Generate Simulated Data 

using namespace std;

const int MAX_CHAR_NUM = 1<<28;

// induce errors in raw sequencing dataset

/*
Output File:
simuD1.fa (simulated data, reads with errors)
simuD1_true.fa (Groundtruth, correct reads, with errinfo in ID)
err.txt (only reads contain errors, with err info)
*/


// awk '{if((NR%2)==1)print $1;else print $0}' /home/xuanzhan/Data/miRNA/salmon/cleandata/cut_573.fastq > input.fq
// awk '{if(NR%4!=0)ORS=" ";else ORS="\n"}1' input.fq | awk '{print $1 " " $(NF-2) " " $NF}' > ./id_read.txt;
// awk '{print $2}' ./id_read.txt |sort |uniq -c| sort -r -nk1 > ./expreLevel_cor.txt

// awk '{if(NR%4!=0)ORS=" ";else ORS="\n"}1' witherr.fq | awk '{print $1 " " $(NF-2) " " $NF}' | awk '{print $2}' |sort |uniq -c| sort -r -nk1 > ./expreLevel_cor_orignal_err.txt
// awk '{if(NR%4!=0)ORS=" ";else ORS="\n"}1' witherr_573_3.fq | awk '{print $1 " " $(NF-2) " " $NF}' | awk '{print $2}' |sort |uniq -c| sort -r -nk1 > ./expreLevel_cor_orignal_err_573_3.txt


//./induce_err -s 2 -c expreLevel_cor.txt -f id_read.txt -o witherr_573_T.fq

// ./miREC.sh -f /home/xuanzhan/Data/miRNA/salmon/cleandata/cut_573.fastq  -s 8 -e 25 -t 26  -o original_cor.fastq
// ./miREC.sh -f /home/xuanzhan/Data/miRNA/code/witherr_605_1.fq  -s 8 -e 25 -t 26  -o Correct_0825_605_1.fastq 
//nohup ./miREC.sh -f /home/xuanzhan/Data/miRNA/code/witherr.fq  -s 8 -e 25 -t 26  -o Correct_0818.fastq &


const char* cs1p = "ATCG";

unsigned int seed;
std::string c_FN;
std::string f_FN;
std::string o_FN;


hash<string> str_hash;
int size_str_hash = MAX_CHAR_NUM; //2^30

struct KmerNode
{
    string kmer;
    int freq;
};

//store info from k-mer frequency file
vector<KmerNode> *kmerfreq; //2^30

//store info from fastq file,includes ID,raw read, phreh quality
vector<string> F_id;
vector<string> F_read;
vector<string> F_quality;

int totalread = 0;

std::ostream *v1logger;

//setting parameters, don't induce error in read_copy lower than "low_setting"
int low_setting = 5; 
// int err_rate_z = 1000000;
// int err_rate_m = 3;

// per read error rate is 3/100000;
inline void displayHelp(const char* prog) {
	printf("Induce a few man-made errors data v1.1, by XuanZhang, March 2021.\n");
	printf("Usage: %s -c <copynumber_FileName> -f <id_read_quality_FileName> -o <output_fileName> \n", prog);
	printf("\t\t -s is seed number\n");
	printf("\t\t -c is real read expression-level file, refer its copy number\n");
	printf("\t\t -f is id_read.txt\n");
	printf("\t\t -o is output file saved\n");
	printf("\t\t -z -m  is per read error rate (m/z) default(3/100000) -z 100000 -m 3\n");
	printf("Example:\n\t\t");
	printf("./induce_err -s 1 -c expreLevel_cor.txt -f id_read.txt -o witherr.fq\n\n");
	}

//show current params
inline void displayParams() {
	printf("seed is s = %d\n", seed);
	printf("copy count is c = %s\n",c_FN.c_str() );
	printf("miRNA read is from f: %s\n", f_FN.c_str());
	printf("Saved output file(withError) is: %s\n", o_FN.c_str());
}

//get params
inline void getPars(int argc, char* argv[]) {
	//v1logger = &std::cout;
	bool is1 = false, is2 = false, is3 = false;
	//bool iskmer = false; //four
	int oc;
	while ((oc = getopt(argc, argv, "s:c:f:o:hf")) >= 0) {
		switch (oc) {
			case 's':
				seed = atoi(optarg);
			// case 'z':
			// 	err_rate_z = atoi(optarg);
			// case 'm':
			// 	err_rate_m = atoi(optarg);
			case 'c':
				c_FN = optarg;
				is1 = true;
				break;
			case 'f':
				f_FN = optarg;
				is2 = true;
				break;
			case 'o':
				o_FN = optarg;
				is3 = true;
				break;
			case 'h':
				displayHelp(argv[0]);
				exit(0);
			case '?':
				std::cerr << "Error parameters.\n Please run 'miREC -h'\n";
				exit(1);
				break;
		}
	}

	if (!is1 || !is2 || !is3 ) {
		fprintf(stderr, "Required parameters are not provided!!\n\n");
		exit(1);
	}
	
	std::ifstream f;

	f.open(c_FN);
	if (f.fail()) {
		fprintf(stderr, "k-mer frequency file  '%s' does not exist.\n", c_FN.c_str());
		exit(1);
	}
	f.close();

	f.open(f_FN);
	if (f.fail()) {
		fprintf(stderr, "read expression-level file '%s' does not exist.\n", f_FN.c_str());
		exit(1);
	}
	f.close();

	//assert(kmer <= L && kmer%4==0);
}

inline void initial() {

	kmerfreq = new vector<KmerNode>[MAX_CHAR_NUM];

}

//string hash: transfer mer/read to integer
int str_hash_index(string s_seq)
{
    //upperCase
    transform(s_seq.begin(), s_seq.end(), s_seq.begin(), ::toupper);
    
    size_t b = str_hash(s_seq);

    int index = b % size_str_hash;

    return index;
}

//check freq_FIle
void readfreq_File(const char *refFile, vector<KmerNode> *kmerfreq) { // processing reference file

    FILE *fp = fopen(refFile, "r");
    if (NULL == fp) {
        printf("fail to open file %s\n", refFile);
        return;
    }

    int r_count = 0;
   
    int freq;
    char *kmer;
    kmer = new char[MAX_CHAR_NUM];
    int hash_index;

    KmerNode node;
    while (fscanf(fp, "%d %s",&freq,kmer) != EOF) {   
        
	        node.kmer = kmer;
	        node.freq = freq;

	        hash_index = str_hash_index(kmer);
	        kmerfreq[hash_index].push_back(node);
	        r_count = r_count + freq;

    }
    fclose(fp);
    cout<<"total frequency"<<r_count<<endl;
}

//check whether errorous kmer
int read_frecheck(string read){

	int hash_index = str_hash_index(read);

	if( kmerfreq[hash_index][0].freq > low_setting )
	{
		//might have protential errors
		return kmerfreq[hash_index][0].freq;
	}
	else
	{
		//no protential error in the read
		return 0;
	}
    		
}

//read fastq :ID_read_quality.txt(preprocess from row_fastq file)
void readrow_fastqFile(const char *refFile) { // processing reference file

	FILE *fp = fopen(refFile, "r");
	if (NULL == fp) {
		printf("fail to open file %s\n", refFile);
		return;
	}

	char *id, *read, *quality;
	id = new char[MAX_CHAR_NUM];
	read = new char[MAX_CHAR_NUM];
	quality = new char[MAX_CHAR_NUM];	

	while (fscanf(fp, "%s %s %s",id, read, quality) != EOF) {
			F_id.push_back(id);
			F_read.push_back(read);
			// F_flag.push_back(0);
			// Freq_b.push_back(0);
			// Freq_a.push_back(0);
			// Err_posi.push_back(0);
			// Err_posi_tmp.push_back(0);
			F_quality.push_back(quality);
			totalread++;		
	}
	cout<<"totalread"<<totalread<<endl;
	fclose(fp);
}


int main(int argc, char *argv[]) 
{

	if (argc < 3) {
		displayHelp(argv[0]);
		return 1;
	}
	
	getPars(argc, argv);

	displayParams();
	initial();
	
	readfreq_File(c_FN.c_str(),kmerfreq);
	readrow_fastqFile(f_FN.c_str());
	ofstream outfile(o_FN.c_str());
	ofstream outfileerr("err_list.txt");
	//ofstream outfile("correct_read.fastq");

	//write simulated reads with errors
	if(!outfile){
    cout << "Unable to open file to write simulated reads with errors";
        exit(1); // terminate with error
    }

    //write simulated reads with errors
	if(!outfileerr){
    cout << "Unable to open file to write simulated reads with errors";
        exit(1); // terminate with error
    }


    int ran,err_ran,ran_num,copy_num,indel_ran;//store random number
    int posi, lett;

    string read,qual;
    int sum_sub_error = 0;
    int sum_ins_error = 0;
    int sum_del_error = 0;

    srand(seed);

	int flag = 0;
	for(unsigned int i=0;i<F_id.size();i++) 
	{
		ran = rand();
		ran_num = ran%100000;
		
		if((ran_num<3) || (flag == 1))// induce 30 in total reads
		{

			copy_num = read_frecheck(F_read.at(i));
			if(  copy_num > 0 )
			{
			    
			    //induce error
			 
				read = F_read.at(i);
				//details of error bases
				err_ran = rand();
				posi = err_ran % read.size();
				lett = err_ran % 4;

				indel_ran = rand() % 10;
				if( indel_ran < 2 )
				{
					//error rate 0.02 indels
					if( indel_ran == 0 )
					{
						//induce insertion errors
						read.insert(posi,&cs1p[lett],1);
						// cout<<">"<<i<<"-"<<copy_num<<"-"<<copy_tmp<<"insert"<<endl;
						
						//revise quality value
						qual = F_quality.at(i);
						qual.insert(posi,&cs1p[2],1);
						F_quality.at(i) = qual;


						//Write errinfo
						outfileerr<<F_id.at(i)<<" "<<copy_num<<" "<<F_read.at(i)<<" "<<read<<" "<<posi<<" + "<<cs1p[lett]<<endl;
						
						//write witherror data.fastq
						outfile<<F_id.at(i)<<endl<<read<<endl<<"+"<<endl<<F_quality.at(i)<<endl;
						
						sum_ins_error++;
						flag = 0;
						
						continue;
					}
					else
					{
						//induce deletion errors
						read.erase(posi,1);
						// cout<<">"<<i<<"-"<<copy_num<<"-"<<copy_tmp<<"deletion"<<endl;

						//revise quality value
						qual = F_quality.at(i);
						qual.erase(posi,1);
						F_quality.at(i) = qual;
					

						//Write errinfo
						outfileerr<<F_id.at(i)<<" "<<copy_num<<" "<<F_read.at(i)<<" "<<read<<" "<<posi<<" - "<<cs1p[lett]<<endl;
						
						//write witherror data.fastq
						outfile<<F_id.at(i)<<endl<<read<<endl<<"+"<<endl<<F_quality.at(i)<<endl;
						
						sum_del_error++;
						flag = 0;

						continue;
					}
				}
				else
				{
					//subsitution error ~1%(0.08)
					if( read[posi] == cs1p[lett] ){lett = (err_ran+1) % 4;}

					read.replace(posi, 1, &cs1p[lett], 0, 1);;
					// cout<<">"<<i<<"-"<<copy_num<<"-"<<copy_tmp<<"subsitution"<<endl;
					
					//Write errinfo
					outfileerr<<F_id.at(i)<<" "<<copy_num<<" "<<F_read.at(i)<<" "<<read<<" "<<posi<<"  "<<cs1p[lett]<<endl;
					
					//write witherror data.fastq
					outfile<<F_id.at(i)<<endl<<read<<endl<<"+"<<endl<<F_quality.at(i)<<endl;
					
				
					sum_sub_error++;
					flag = 0;
					continue;
				}
			
			}
			else
			{
				flag = 1;
			
				//write witherror data.fastq
				outfile<<F_id.at(i)<<endl<<F_read.at(i)<<endl<<"+"<<endl<<F_quality.at(i)<<endl;
			}
		}
		else // don't induce error
		{	flag = 0;
			//write witherror data.fastq
			outfile<<F_id.at(i)<<endl<<F_read.at(i)<<endl<<"+"<<endl<<F_quality.at(i)<<endl;
		}	
		// cout << "finish" << i <<endl;
	}


	cout << "sub: "<<sum_sub_error<<"ins: "<<sum_ins_error<<"del: "<<sum_del_error<<endl;
	outfile.close();
	outfileerr.close();
}









