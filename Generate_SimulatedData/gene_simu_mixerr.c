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

/*
Output File:
simuD1.fa (simulated data, reads with errors)
simuD1_true.fa (Groundtruth, correct reads, with errinfo in ID)
err.txt (only reads contain errors, with err info)
*/

// ./gene_simu_sub -c /home/xuanzhan/Data/miRNA/salmon/cleandata/read_explevel_573.txt -f /home/xuanzhan/Data/miRNA/ref/mature.fa -o simuD1.fa -g simuD1_ture.fa
//./gene_simu_sub -s 2 -c /home/xuanzhan/Data/miRNA/salmon/cleandata/read_explevel_573.txt -f /home/xuanzhan/Data/miRNA/ref/mature.fa -o /home/xuanzhan/Data/miRNA/simu/data/s_d/s_d2/simusD2.fa -g /home/xuanzhan/Data/miRNA/simu/data/s_d/s_d2/simusD2_ture.fa

vector<int> K_freq;
vector<string> R_read;

const char* cs1p = "ATCG";

unsigned int seed;
std::string c_FN;
std::string f_FN;
std::string o_FN;
std::string g_FN;

inline void displayHelp(const char* prog) {
	printf("generate simulated data v1.1, by XuanZhang, March 2020.\n");
	printf("Usage: %s -c <copynumber_FileName> -f <matureRNA_FileName> -o <output_fileName> \n", prog);
	printf("\t\t -s is seed number\n");
	printf("\t\t -c is real read expression-level file, refer its copy number\n");
	printf("\t\t -f is mature.fa\n");
	printf("\t\t -o is output file saved\n");
	printf("Example:\n\t\t");
	printf("./gene_simu -s 2 -c read_explevel_573.txt -f /home/xuanzhan/Data/miRNA/ref/mature.fa -o simulated.fa\n\n");
	//./gene_simu -c /home/xuanzhan/Data/miRNA/salmon/cleandata/read_explevel_573.txt -f /home/xuanzhan/Data/miRNA/ref/mature.fa -o simulated.fa
}

//show current params
inline void displayParams() {
	printf("seed is s = %d\n", seed);
	printf("copy count is c = %s\n",c_FN.c_str() );
	printf("miRNA read is from f: %s\n", f_FN.c_str());
	printf("Saved output file(withError) is: %s\n", o_FN.c_str());
	printf("Saved output file(GroundTruth '_true.fa') is: %s\n", g_FN.c_str());
}

//get params
inline void getPars(int argc, char* argv[]) {
	//v1logger = &std::cout;
	bool is1 = false, is2 = false, is3 = false, is4 = false;
	//bool iskmer = false; //four
	int oc;
	while ((oc = getopt(argc, argv, "s:c:f:o:g:hf")) >= 0) {
		switch (oc) {
			case 's':
				seed = atoi(optarg);
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
			case 'g':
				g_FN = optarg;
				is4 = true;
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

	if (!is1 || !is2 || !is3 || !is4) {
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


//read copy_number file
void Copynum_File(const char *refFile) { // processing reference file

	FILE *fp = fopen(refFile, "r");
	if (NULL == fp) {
		printf("fail to open file %s\n", refFile);
		return;
	}

	int r_count = 0;

	int copy_num;
	char *kmer;
	kmer = new char[MAX_CHAR_NUM];

	while (fscanf(fp, "%s %d",kmer, &copy_num) != EOF ){
		
		if (copy_num < 500 && copy_num > 5) 
		{
			K_freq.push_back(copy_num);
			r_count++;
		}

		//printf("%s %d\n", K_mer[r_count].c_str(), K_freq[r_count]);
		
	}
	fclose(fp);
}

//read mature.fa (/home/xuanzhan/Data/miRNA/ref/mature.fa)
void MiRNA_File(const char *refFile) { // processing reference file

	FILE *fp = fopen(refFile, "r");
	if (NULL == fp) {
		printf("fail to open file %s\n", refFile);
		return;
	}

	int freq = 1;
	char *kmer;
	kmer = new char[MAX_CHAR_NUM];

	while (fscanf(fp, "%s",kmer) != EOF) {
		
		if (freq%6 == 0) {
			R_read.push_back(kmer);		
		}
		else{

		}
		freq ++;	
	}
	fclose(fp);
}

// string induce_suberr(string read, int posi, int lett )
// {



// 	return read;
// }

int main(int argc, char *argv[]) {

	if (argc < 3) {
		displayHelp(argv[0]);
		return 1;
	}
	
	getPars(argc, argv);

	displayParams();
	//initial();
	
	Copynum_File(c_FN.c_str());
	MiRNA_File(f_FN.c_str());
	ofstream outfile(o_FN.c_str());
	ofstream outfiletrue(g_FN.c_str());
	ofstream outfileerr("err.txt");
	//ofstream outfile("correct_read.fastq");

	//write simulated reads with errors
	if(!outfile){
    cout << "Unable to open file to write simulated reads with errors";
        exit(1); // terminate with error
    }

    //write Ground Truth reads
	if(!outfiletrue){
    cout << "Unable to open file to write Ground Truth reads";
        exit(1); // terminate with error
    }

    unsigned int copy_tmp;
    int ran,err_ran,ran_num,copy_num,indel_ran;//store random number
    int posi, lett;

    string read,err_read;
    int sum = 0;
    int sum_sub_error = 0;
    int sum_ins_error = 0;
    int sum_del_error = 0;

    srand(seed);

	for(unsigned int i=1;i<K_freq.size();i++) 
	{
		copy_num = K_freq.at(i);
		copy_tmp = copy_num;
		sum = copy_tmp + sum;

		//Simulated Data with errors
		while(copy_tmp > 0){

			read = R_read.at(i);
			ran = rand();
			ran_num = ran%100;

			//induce errors	
			if(ran_num <= 1){

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

						//Write errinfo
						outfileerr<<">"<<i<<"-"<<copy_num<<"-"<<copy_tmp<<" "<<R_read.at(i)<<" "<<read<<" "<<posi<<" + "<<cs1p[lett]<<endl;
						
						//write Ground Truth reads
						outfiletrue<<">"<<i<<"-"<<copy_num<<"-"<<copy_tmp<<" "<<posi<<" + "<<cs1p[lett]<<" "<<read<<endl;
						outfiletrue<<R_read.at(i)<<endl;

						//write simulated data
						outfile<<">"<<i<<"-"<<copy_num<<"-"<<copy_tmp<<endl;
						outfile<<read<<endl;

						sum_ins_error++;
						copy_tmp--;
						continue;
					}
					else
					{
						//induce deletion errors
						read.erase(posi,1);
						// cout<<">"<<i<<"-"<<copy_num<<"-"<<copy_tmp<<"deletion"<<endl;
						//Write errinfo
						outfileerr<<">"<<i<<"-"<<copy_num<<"-"<<copy_tmp<<" "<<R_read.at(i)<<" "<<read<<" "<<posi<<" - "<<R_read.at(i)[posi]<<endl;
						//write Ground Truth reads
						outfiletrue<<">"<<i<<"-"<<copy_num<<"-"<<copy_tmp<<" "<<posi<<" - "<<R_read.at(i)[posi]<<" "<<read<<endl;
						outfiletrue<<R_read.at(i)<<endl;
						//write simulated data
						outfile<<">"<<i<<"-"<<copy_num<<"-"<<copy_tmp<<endl;
						outfile<<read<<endl;

						sum_del_error++;
						copy_tmp--;
						continue;
					}
				}
				else
				{
					//subsitution error ~1%(0.08)
					if( read[posi] == cs1p[lett] )
					{	
						lett = (err_ran+1) % 4;	
					}

					read.replace(posi, 1, &cs1p[lett], 0, 1);;
					// cout<<">"<<i<<"-"<<copy_num<<"-"<<copy_tmp<<"subsitution"<<endl;
					//Write errinfo
					outfileerr<<">"<<i<<"-"<<copy_num<<"-"<<copy_tmp<<" "<<R_read.at(i)<<" "<<read<<" "<<posi<<" "<<R_read.at(i)[posi]<<" "<<cs1p[lett]<<endl;
					
					//write Ground Truth reads
					outfiletrue<<">"<<i<<"-"<<copy_num<<"-"<<copy_tmp<<" "<<posi<<" "<<R_read.at(i)[posi]<<" "<<cs1p[lett]<<" "<<read<<endl;
					outfiletrue<<R_read.at(i)<<endl;

					//write simulated data
					outfile<<">"<<i<<"-"<<copy_num<<"-"<<copy_tmp<<endl;
					outfile<<read<<endl;
					
					sum_sub_error++;
					copy_tmp--;
					continue;
				}
			}
			else
			{	
				//no error induced
				//write simulated data
				outfile<<">"<<i<<"-"<<copy_num<<"-"<<copy_tmp<<endl;
				outfile<<R_read.at(i)<<endl;
				
				//write Ground Truth reads
				outfiletrue<<">"<<i<<"-"<<copy_num<<"-"<<copy_tmp<<endl;
				outfiletrue<<R_read.at(i)<<endl;

				copy_tmp --;
			}	
		}
	}
	cout<<"sub: "<<sum_sub_error<<"ins: "<<sum_ins_error<<"del: "<<sum_del_error<<endl;

	outfile.close();
	outfiletrue.close();
	outfileerr.close();
}







