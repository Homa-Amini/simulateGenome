/***read Simulation ****/
/***Copy to work*******/
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <popt.h>
#include <stdexcept>
#include "faidx.h"
#include <time.h>
#include <random>
#include <math.h>
#include <string.h>
#include <limits>
#include <functional>   // std::plus
// Boost library
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/detail/seed_impl.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_real.hpp>
// ART library for seq-error
#include "empdist.h"

using namespace std;

boost::random::mt19937 gen;
double empdist::prob_err[HIGHEST_QUAL];

enum substitution_Type{
	Deaminate = 0,
	Sequencing_Error = 1 ,
	Divergence = 2,
	Mismatch = 3,
};


int nextBlock(faidx_t *REFIndex,string ChrName,  int &Begin,  int Length,string &Ref)
{
        int len[1]={0};
	char ChromosomeName[ChrName.size()+1];//as 1 char space for null is also required
	strcpy(ChromosomeName, ChrName.c_str());
        Ref=string(faidx_fetch_seq (REFIndex,ChromosomeName,Begin,Begin+Length,len));
	return 1;
}

int Uniform_randomPosition(signed int start,signed  int end) { //closed interval [start, end]
	boost::random::uniform_int_distribution<> dist( start , end );
	return dist(gen);
}

double Uniform_randomNumber(double start,double end){
	boost::uniform_real<> uni_dist(start, end);
	return uni_dist(gen);
}

 
bool makeMatrix(unsigned long int (&Matrix)[5][5], char first,char second)
{
	string Bases="ACGTN";
	int row=Bases.find(first);
	int col=Bases.find(second);
	if((row !=std::string::npos) && (col != std::string::npos))
	{
		Matrix[row][col]++;
		return true;
	}
	return false;
}


string Generate_BacteriaRead( long int length){
        string Bases="ACGT";
	string bacteriaRead="";
        for(int index=0; index < length; index++){
		int randomNumber = Uniform_randomPosition(0,3);
		char BacteriaBase=Bases.at(randomNumber);
		bacteriaRead.append(1,BacteriaBase);
           }
	return bacteriaRead;
}



int Substitute(string& Read,substitution_Type sType,double Substitution_Ratio5Prime,double Substitution_Ratio3Prime,double Substitution_RatioMiddle,double success_fraction){
  try{
	string Bases="ACGT";
	char seenBase = ' ';
	char SubstituteBase = Bases[Uniform_randomPosition(0,Bases.length()-1)];
	int SubstitutePosition = Uniform_randomPosition(0,Read.length()-1);
	int overhang_length= 20000000; // a large number bigger that maximum read length 
	double randNb = 0;
	int numberOfSubstitution = 0 ;
//	int nbOfMissMatches = int(Substitution_Ratio5Prime);
	vector <int> randomPositions;
  	vector<int>::iterator it;
	switch(sType){
		case Deaminate :
			seenBase = 'C';
			SubstituteBase = 'T';
			randNb = Uniform_randomNumber(0,1);
			if(randNb<0.5){
				double random_number= Uniform_randomNumber(0,1);
				while(overhang_length > Read.length()){ //almost always true caze we don't have reads longer than 6milion bp
					overhang_length=(log(1-random_number)/log(1-success_fraction));
					random_number= Uniform_randomNumber(0,1);
				}
				for(int deaminateIndex=0; deaminateIndex < overhang_length ; deaminateIndex++){
					if( (Read[deaminateIndex] == seenBase) && (Uniform_randomNumber(0,1) < Substitution_Ratio5Prime) ){
				        	Read.replace(deaminateIndex,1,1,SubstituteBase);
					}
		                        if( (Read[(Read.length()-1)-deaminateIndex] == seenBase) && ( Uniform_randomNumber(0,1) < Substitution_Ratio3Prime) ){
						Read.replace((Read.length()-1)-deaminateIndex,1,1,SubstituteBase);
					}
				}
			} else{
				for(int deaminateIndex=0; deaminateIndex < Read.length() ; deaminateIndex++){
					if( (Read[deaminateIndex] == seenBase) && ( Uniform_randomNumber(0,1) < Substitution_RatioMiddle ) ){
						Read.replace(deaminateIndex,1,1,SubstituteBase);
					}
				}
			}
			break;
		case Sequencing_Error :
			seenBase = Read[0] ;
			while (SubstituteBase == seenBase){  //if substitute base is the same base as seen base go for another random base
				SubstituteBase = Bases[Uniform_randomPosition(0,Bases.length()-1)];
			}
			if( Uniform_randomNumber(0,1) < Substitution_Ratio5Prime) {
				Read.replace(0,1,1,SubstituteBase);
			}
			break;
		case Divergence :
			for(int readIndex=0; readIndex < Read.length() ; readIndex++){
				seenBase = Read[readIndex];
				while (SubstituteBase == Read[readIndex]){  //if substitute base is the same base as seen base go for another random base
					SubstituteBase = Bases[Uniform_randomPosition(0,Bases.length()-1)];
				}
				if( Uniform_randomNumber(0,1) < Substitution_Ratio5Prime ){
					Read.replace(readIndex,1,1,SubstituteBase);
					numberOfSubstitution ++;
				}
			}
			break;
		case Mismatch :
			seenBase = Read[SubstitutePosition];
			while (SubstituteBase == seenBase){  //if substitute base is the same base as seenBase go for another random base
				SubstituteBase = Bases[Uniform_randomPosition(0,Bases.length()-1)];
			}
			Read.replace(SubstitutePosition,1,1,SubstituteBase);
			return SubstitutePosition;
			break;	
			
		}
	}
// }
 catch (const out_of_range& oor) {
	cerr << "Out of range: "<<oor.what() << " happend on read " << Read << endl;
 }
}


string Reverse_Complement(string sequence){

 try{ 
	replace(sequence.begin(), sequence.end(),'A','a');
	replace(sequence.begin(), sequence.end(),'C','c');
	replace(sequence.begin(), sequence.end(),'G','g');
	replace(sequence.begin(), sequence.end(),'T','A');
	replace(sequence.begin(), sequence.end(),'a','T');
	replace(sequence.begin(), sequence.end(),'c','G');
	replace(sequence.begin(), sequence.end(),'g','C');

	return string(sequence.rbegin(),sequence.rend());
 }
 catch (const out_of_range& oor) {
	cerr << "Out of range: "<<oor.what() << " happemd on read " << sequence <<endl;
 }
} 


int Diverge_Genome(const string& REFaddr,const string& diverged_genome,double divRate){
	cerr << "reading modelled Refernce .... " << endl;
	cerr << "Diverge it in to "<< diverged_genome << "......" << endl;
	string IndexREFFile=REFaddr+".fai";
	ifstream REF(REFaddr.c_str());
	
	faidx_t *REFIndex = fai_load(IndexREFFile.c_str());
	if (!REF.is_open())
	{
		cerr << "Unable to open Reference file" << endl;
                REF.close();
                return 0; 
	}
	if(!REFIndex)
        {
                cerr << "Unable to open Reference index file"<< endl;
		cerr  <<"Creating it now ..." <<endl;
                REFIndex = fai_load(REFaddr.c_str());
                if(!REFIndex){
                        cerr << "Try to creat it but it seems it's not its day!" <<endl;
                        REF.close();
                        fai_destroy(REFIndex);
                        cerr << "Unable to open Reference Index file" << endl <<"Try to creat it but it seems it's not its day!" << endl;
                        return 0;
                }
        }
	int nbOFchromosomes= faidx_nseq(REFIndex);
	cerr << "found " << nbOFchromosomes << " chromosoms" << endl;
	ofstream readFile(diverged_genome.c_str());
	streambuf *placeHolder = cout.rdbuf();
	cout.rdbuf(readFile.rdbuf());
	string ChrName = "";
	for(int chrNB=0; chrNB < nbOFchromosomes ; chrNB ++ ) {
		ChrName = faidx_iseq(REFIndex, chrNB);
		cout<<">"<<ChrName<<endl;
		int Chromosome_Length = faidx_seq_len(REFIndex,faidx_iseq(REFIndex,chrNB));	
		for( int numberOfBlocks=0; numberOfBlocks < ((Chromosome_Length/100)+1); numberOfBlocks++){
			int StartBlock =  100 * numberOfBlocks ; //Uniform_randomPosition(0, GenomeLength-100);
			string tmpSequence="";
			if(nextBlock(REFIndex,ChrName, StartBlock , 99 ,tmpSequence)){ // 99 = 100 -1 means block size
				Substitute(tmpSequence,Divergence,divRate,0,0,0);//This substitution is divergence
			}
			cout << tmpSequence << endl;
		}
	}
	readFile.close();
	cout.rdbuf(placeHolder);
	return 1;	

}

bool SeqError(empdist qdist,short readLength, const string& sequence, const string& quality,bool highQuality){
	
	if(!highQuality){
		vector<short> qual ;
		qdist.get_read_qual(qual, readLength, true);
     		////sequencing error substitution
		for(int len_index=0;len_index<readLength;len_index++){
			quality.append(1,char(qual[len_index]+33));  //calculating Base Quality score
			double P = pow(10,(-1)*(qual[len_index])/10); //calculating error probability base on Quality score
			string tmpChar(1,sequence[len_index]);
			Substitute(tmpChar,Sequencing_Error,P,0,P,0);//This substitution is sequencing error
			sequence[len_index] = tmpChar[0];
		}
	} 
	else{
		quality=string(readLength,'I');
	}
	return true;
}

bool loadSeqError(const string& inputError, short q_shift_up,short readLength,empdist& qdist){
	// ART library
	short q_shift_up_2 = q_shift_up ;
	short max_qual_s=64;//muist be 64
	empdist::set_err_prob();
	//empdist qdist;
	if(ErrorProfileAddress.empty()){
		qdist.setdist("HS10", false, readLength);//HS10
	} else {
		qdist.setdist(ErrorProfileAddress,ErrorProfileAddress,false, readLength);
	}
	// Some sanity check
	int profile_size=qdist.qual_dist_first.size();
	if(readLength > profile_size){
		if(profile_size == 0){
			cerr << "Fatal Error: " <<  ErrorProfileAddress << ", is not a valid profile." << endl << endl;
		} else {
			cerr<<"Fatal Error: The read length, "<<readLength<<", exceeds the maximum first read profile length, "<<profile_size<< "." <<endl<<endl;
		}
		return 0;
	}
	// increase overal base quality value
	for(size_t i=0; i<qdist.qual_dist_first.size(); i++){
		for(map<unsigned int, unsigned short>::iterator it=qdist.qual_dist_first[i].begin(); it!=qdist.qual_dist_first[i].end(); it++){
			if(q_shift_up<0 && (-q_shift_up>it->second)){ it->second=0; }
			else{ it->second+=q_shift_up; if(it->second>max_qual_s) it->second=max_qual_s; }
		}
	}
	for(size_t i=0; i<qdist.qual_dist_second.size(); i++){
		for(map<unsigned int, unsigned short>::iterator it=qdist.qual_dist_second[i].begin(); it!=qdist.qual_dist_second[i].end(); it++){
			if(q_shift_up_2<0 && (-q_shift_up_2>it->second)){ it->second=0; }
			else{ it->second+=q_shift_up_2; if(it->second>max_qual_s) it->second=max_qual_s; }
		}
	}
	return true;	
		// End seq error
}

string GenerateRead(const string& REFaddr,unsigned long int NumberOfGenomicReads){

	cerr << "reading model refernce ......" << endl;
	string IndexREFFile=REFaddr+".fai";
	ifstream REF(REFaddr.c_str());
	
	faidx_t *REFIndex = fai_load(IndexREFFile.c_str());
	if (!REF.is_open())
	{
		cerr << "Unable to open Reference file" << endl;
                REF.close();
                return 0; 
	}
	if(!REFIndex)
        {
                cerr << "Unable to open Reference index file"<< endl <<"Create it now ..." <<endl;
                REFIndex = fai_load(REFaddr.c_str());
                if(!REFIndex){
                        cerr << "Try to creat it but it seems it's not its day!" <<endl;
                        REF.close();
                        fai_destroy(REFIndex);
                        cerr << "Unable to open Reference Index file" << endl <<"Try to creat it but it seems it's not its day!" << endl;
                        return 0;
                }
        }

}

int Generate_reads_Post_Process(const string& REFaddr,unsigned long int NumberOfGenomicReads, int BeginReadLength, int EndReadLength, int Step, unsigned long int GenomeLength, unsigned long int Chromosome_Length_Threshold , double DeaminationRatio5Prime,double DeaminationRatio3Prime, double LowDeaRatioMiddle,double success_fraction,unsigned long int NumberOfBacteriaReads, bool highQuality, const string& ErrorProfileAddress, short q_shift_up, int nbOfMismatches){       
	short q_shift_up_2 = q_shift_up ;
	short max_qual_s=64;//muist be 64
	vector<int> V_MismatchPosition(nbOfMismatches+1);
	std::fill (V_MismatchPosition.begin(),V_MismatchPosition.end(),-1);
	
	cerr << "reading model refernce ......" << endl;
	cerr << "First quality shift:  " << q_shift_up << endl;
	cerr << "Second quality shift:  " << q_shift_up_2 << endl;
	cerr << "Error Profile address is " << ErrorProfileAddress << endl;
	string IndexREFFile=REFaddr+".fai";
	ifstream REF(REFaddr.c_str());
	
	faidx_t *REFIndex = fai_load(IndexREFFile.c_str());
	if (!REF.is_open())
	{
		cerr << "Unable to open Reference file" << endl;
                REF.close();
                return 0; 
	}
	if(!REFIndex)
        {
                cerr << "Unable to open Reference index file"<< endl <<"Create it now ..." <<endl;
                REFIndex = fai_load(REFaddr.c_str());
                if(!REFIndex){
                        cerr << "Try to creat it but it seems it's not its day!" <<endl;
                        REF.close();
                        fai_destroy(REFIndex);
                        cerr << "Unable to open Reference Index file" << endl <<"Try to creat it but it seems it's not its day!" << endl;
                        return 0;
                }
        }
	// Seq error
	// ART library
	bool nRef = false;
    	int StartBlock=0;
	for( int readLength=BeginReadLength; readLength<=EndReadLength;readLength+=Step){
		// Seq error just loading or reading it ,No substitution.
		// ART library
		empdist::set_err_prob();
		empdist qdist;
		if(ErrorProfileAddress.empty()){
			qdist.setdist("HS10", false, readLength);//HS10
		} else {
			qdist.setdist(ErrorProfileAddress,ErrorProfileAddress,false, readLength);
		}
		// Some sanity check
		int profile_size=qdist.qual_dist_first.size();
		if(readLength > profile_size){
			if(profile_size == 0){
				cerr << "Fatal Error: " <<  ErrorProfileAddress << ", is not a valid profile." << endl << endl;
			} else {
				cerr<<"Fatal Error: The read length, "<<readLength<<", exceeds the maximum first read profile length, "<<profile_size<< "." <<endl<<endl;
			}
			return 0;
		}
		// increase overal base quality value
		for(size_t i=0; i<qdist.qual_dist_first.size(); i++){
			for(map<unsigned int, unsigned short>::iterator it=qdist.qual_dist_first[i].begin(); it!=qdist.qual_dist_first[i].end(); it++){
				if(q_shift_up<0 && (-q_shift_up>it->second)){ it->second=0; }
				else{ it->second+=q_shift_up; if(it->second>max_qual_s) it->second=max_qual_s; }
			}
		}
		for(size_t i=0; i<qdist.qual_dist_second.size(); i++){
			for(map<unsigned int, unsigned short>::iterator it=qdist.qual_dist_second[i].begin(); it!=qdist.qual_dist_second[i].end(); it++){
				if(q_shift_up_2<0 && (-q_shift_up_2>it->second)){ it->second=0; }
				else{ it->second+=q_shift_up_2; if(it->second>max_qual_s) it->second=max_qual_s; }
			}
		}
		
		// End seq error
		if(NumberOfGenomicReads>0){ //because we don't always want Genomic reads, sometimes just bacterial one.
			cerr << "Choping genomic reads in length of " << readLength << endl;
			string lengthString = static_cast<ostringstream*>( &(std::ostringstream() << (readLength)) )->str();
	//		string nbOfMissMatchesString = static_cast<ostringstream*>( &(std::ostringstream() << (nbOfMisMatches)) )->str();
		//	ofstream readFile(("/mnt/scratch/Homa/ChoppedReads/GenomicRead_"+lengthString+"_"+nbOfMisMatchesString+".sam").c_str());
			ofstream readFile(("/mnt/scratch/Homa/ChoppedReads/GenomicRead_"+lengthString+".sam").c_str());
			streambuf *placeHolder = cout.rdbuf();
			cout.rdbuf(readFile.rdbuf());
			cout << "@HD\tVN:1.5\tSO:coordinate" << endl;
			cout << "@SQ\tSN:fakeGenome\tLN:"<<GenomeLength << endl;
			int nbOFchromosomes= faidx_nseq(REFIndex);
			string tmpSequence="";
			//int Coverage = 1;
			//for( long int numberOfReads=0; numberOfReads<(((GenomeLength/readLength)*Coverage)+1); numberOfReads++){
			for(int numberOfReads=0; numberOfReads < NumberOfGenomicReads; ++numberOfReads)
			{
				int RandomChromosomeNB = Uniform_randomPosition(0, nbOFchromosomes-1);
				string ChrName = faidx_iseq(REFIndex,RandomChromosomeNB);
			 	int Chromosome_Length = faidx_seq_len(REFIndex,faidx_iseq(REFIndex,RandomChromosomeNB));
				if(readLength > Chromosome_Length){
					if(Chromosome_Length == 0){
						cerr << "Fatal Error: 0 , is not a valid chromosome length." << endl << endl;
					} else {
						cerr<<"Fatal Error: The read length, "<<readLength<<", is bigger than Chromosome length, "<<  Chromosome_Length << "." <<endl<<endl;
					}
					return 0;
				}
					
				StartBlock =Uniform_randomPosition(0,Chromosome_Length-readLength);
				nextBlock(REFIndex,ChrName,StartBlock,readLength-1,tmpSequence);
				std::size_t foundN = tmpSequence.find('N');
				if(foundN == std::string::npos)
				{
					//string keepOriginalTMPsequence = tmpSequence ;
					string oneTolatesttmpSequence = tmpSequence ;
					for(int MismatchIndex=0; MismatchIndex != (nbOfMismatches + 1); ++MismatchIndex){
						int newSubstitutedPosition = -1 ;
						vector<int>::iterator it = V_MismatchPosition.end() ;

						if(MismatchIndex != 0){
							newSubstitutedPosition = Substitute(tmpSequence,Mismatch,0,0,0,0);//This substitution adds miss matches to reads
							it = find (V_MismatchPosition.begin(), V_MismatchPosition.end(), newSubstitutedPosition);
						}

						if (it != V_MismatchPosition.end()){
							//Element has found in the vector, go for a another position
							MismatchIndex--;
							//tmpSequence = oneTolatesttmpSequence;
						}
						else{
							oneTolatesttmpSequence = tmpSequence ;

							tmpSequence = (numberOfReads % 2)? tmpSequence : Reverse_Complement(tmpSequence) ;
							Substitute(tmpSequence,Deaminate,DeaminationRatio5Prime,DeaminationRatio3Prime,LowDeaRatioMiddle,success_fraction);//This substitution is deamination
							tmpSequence = (numberOfReads % 2)? tmpSequence : Reverse_Complement(tmpSequence) ;
							string ReadName="Genomic";
							int flag=(numberOfReads % 2)? 0 : 16;
							string Quality = "";
							if(!highQuality){
								vector<short> qual ;
								qdist.get_read_qual(qual, readLength, true);
     								////sequencing error substitution
								for(int len_index=0;len_index<readLength;len_index++){
									Quality.append(1,char(qual[len_index]+33));  //calculating Base Quality score
									double P = pow(10,(-1)*(qual[len_index])/10); //calculating error probability base on Quality score
									string tmpChar(1,tmpSequence[len_index]);
									Substitute(tmpChar,Sequencing_Error,P,0,P,0);//This substitution is sequencing error
									tmpSequence[len_index] = tmpChar[0];
								}
							} 
							else{
									Quality=string(readLength,'I');
							}
							//cout << "Quality is:" << Quality <<endl;
								
							V_MismatchPosition[MismatchIndex]=newSubstitutedPosition;	
							cout << ReadName <<"_"<<numberOfReads<<"_"<<MismatchIndex<< "\t"<< flag <<"\t*\t"<<StartBlock+1<<"\t0\t"<<readLength<<"M\t*\t0\t0\t"<< tmpSequence << "\t" << Quality<<"\tZP:i:"<<StartBlock+1<<"\tZF:i:"<<flag<<"\tZC:Z:"<<ChrName<< "\tZM:i:"<<MismatchIndex<<endl;
							//tmpSequence = oneTolatesttmpSequence;
						}
						tmpSequence = oneTolatesttmpSequence;
					}
				} else{
					--numberOfReads;
					//cerr << "Throw away read as it contain N(s)" << endl;
				}
			}
			readFile.close();
			cout.rdbuf(placeHolder);
		}
		//Genrate Bacteria
		if(NumberOfBacteriaReads>0){
			cerr << "Generating Bacteria Reads with length " << readLength << endl;
			string lengthString = static_cast<ostringstream*>( &(std::ostringstream() << (readLength)) )->str();
			ofstream readFile(("/mnt/scratch/Homa/ChoppedReads/Bacteria_"+lengthString+".sam").c_str());
			streambuf *placeHolder = cout.rdbuf();
			cout.rdbuf(readFile.rdbuf());
			cout << "@HD\tVN:1.5" << endl;
			cout << "@SQ\tSN:BacteriaReads\tLN:"<<NumberOfBacteriaReads << endl;
			string tmpSequence="";
			for(int numberOfReads=0;numberOfReads != NumberOfBacteriaReads; ++numberOfReads){
				tmpSequence = (numberOfReads % 2)? Generate_BacteriaRead(readLength) : Reverse_Complement(Generate_BacteriaRead(readLength)) ;
				Substitute(tmpSequence,Deaminate,DeaminationRatio5Prime,DeaminationRatio3Prime,LowDeaRatioMiddle,success_fraction);//This substitution is deamination
				tmpSequence = (numberOfReads % 2)? Generate_BacteriaRead(readLength) : Reverse_Complement(Generate_BacteriaRead(readLength)) ;
				string ReadName="Bacteria";
				int flag=(numberOfReads % 2)? 0 : 16;
				string Quality = "";
				vector<short> qual ;
				qdist.get_read_qual(qual, readLength, true);
				for(int len_index=0;len_index<readLength;len_index++){
					Quality.append(1,char(qual[len_index]+33));
					double P = pow(10,(-1)*(qual[len_index])/10);
					string tmpChar(1,tmpSequence[len_index]);
					Substitute(tmpChar,Sequencing_Error,P,0,P,0);//This substitution is sequencing error
					tmpSequence[len_index] = tmpChar[0];
				}
				cout << ReadName <<"_"<<numberOfReads<< "\t"<<flag<<"\tBacterialRead\t"<<-1<<"\t0\t"<<readLength<<"M\t*\t0\t0\t"<< tmpSequence << "\t" << Quality <<"\tZP:i:-1\tZF:i:"<<flag <<"\tZC:Z:Bacteria"<<endl;
			}
			readFile.close();
			cout.rdbuf(placeHolder);
		}
		
	}
	return 1;
}




int reading(unsigned long int (&Matrix)[5][5],const string& REFaddr,int BlockSize)
{
	cerr << "reading Refernce ......" << endl;
	string IndexREFFile=REFaddr+".fai";
	ifstream REF(REFaddr.c_str());
	faidx_t *REFIndex = fai_load(IndexREFFile.c_str());
	if (!REF.is_open())
	{
		cerr << "Unable to open Reference file" << endl;
                REF.close();
                return 0; 
	}
	if(!REFIndex)
        {
                cerr << "Unable to open Reference index file"<< endl;
		cerr << "Will create it now ..." <<endl;
                REFIndex = fai_load(REFaddr.c_str());
                if(!REFIndex){
                        cerr << "Try to creat it but it seems it's not its day!" <<endl;
                        REF.close();
                        fai_destroy(REFIndex);
                        cerr << "Unable to open Reference Index file" << endl <<"Try to creat it but it seems it's not its day!" << endl;
                        return 0;
                }
        }
	
	cerr << "Reference index has been generated correctly" << endl;
	int nbOFchromosomes= faidx_nseq(REFIndex);
	cerr << "found "<< nbOFchromosomes << " chromosomes" << endl;


	for(int chrNB=0; chrNB < nbOFchromosomes ; chrNB ++ )
	{
		string ChrName = faidx_iseq(REFIndex, chrNB);
		bool nRef = false;
		string Ref;
		int StartBlock=0;
		if(!nextBlock(REFIndex,ChrName,StartBlock,BlockSize,Ref)){
			cerr << "Error in reading Reference" <<endl;
			return 0;
		}
		char first='X',second='X';
		for(int Index=0;Index<=Ref.length();)
		{
			if(nRef)
			{
				StartBlock+=BlockSize;
				Ref="";
				if(!nextBlock(REFIndex,ChrName,StartBlock,BlockSize,Ref)){
					cerr << "Error in reading Reference" <<endl;
					return 1;
				} else
				{
					nRef = false;
					Index=0;
					second=(Ref.c_str())[Index];
					makeMatrix(Matrix,first,second);
					Index++;
				}
			} else
			{
				first=(Ref.c_str())[Index];
				if((++Index)<Ref.length())
				{
					second=(Ref.c_str())[Index];
					makeMatrix(Matrix,first,second);
				} else
				{
					second='X';
					nRef=true;
				}
			}
		
		}
	}
	return 1;

}

void PrintMatrix(unsigned long int (&Matrix)[5][5],const string& matrixFilename)
{
	cerr << "writting matrix to file "<< matrixFilename <<endl;
	ofstream matrixFile(matrixFilename.c_str());
	streambuf *placeHolder = cout.rdbuf();
	cout.rdbuf(matrixFile.rdbuf());
	char Base[5]={'A','C','G','T','N'};
	cout << "\tA\tC\tG\tT\tN"<<endl;
	
	for (int index=0 ; index<5 ;index++)
	{
		cout << Base[index] <<"\t"<< Matrix[index][0]<<"\t"<< Matrix[index][1]<<"\t"<< Matrix[index][2]<<"\t"<< Matrix[index][3]<<"\t"<< Matrix[index][4] << endl;	
	}
	matrixFile.close();
	cout.rdbuf(placeHolder);
}

void PrintFraction(unsigned long int (&Matrix)[5][5],long double (&TransitoinMatrix)[5][5],const string& transitionFilename)
{
	cerr << "writting transition matrix to file "<< transitionFilename <<endl;
	ofstream transitionFile(transitionFilename.c_str());
	streambuf *placeHolder = cout.rdbuf();
	cout.rdbuf(transitionFile.rdbuf());
	char Base[5]={'A','C','G','T','N'};
	cout << "\tA\tC\tG\tT\tN"<<endl;
	
	for (int index=0 ; index<5 ;index++)
	{
		unsigned long int sum = Matrix[index][0]+Matrix[index][1]+ Matrix[index][2]+Matrix[index][3]+Matrix[index][4];
		TransitoinMatrix[index][0]=(Matrix[index][0]*1.0)/sum;
		TransitoinMatrix[index][1]=(Matrix[index][1]*1.0)/sum;
		TransitoinMatrix[index][2]=(Matrix[index][2]*1.0)/sum;
		TransitoinMatrix[index][3]=(Matrix[index][3]*1.0)/sum;
		TransitoinMatrix[index][4]=(Matrix[index][4]*1.0)/sum;
		cout << Base[index] << "\t" << TransitoinMatrix[index][0] << "\t" << TransitoinMatrix[index][1] << "\t" << TransitoinMatrix[index][2] << "\t" << TransitoinMatrix[index][3] << "\t" << TransitoinMatrix[index][4] << endl;
	}
	transitionFile.close();
	cout.rdbuf(placeHolder);
}

void generate_Genome(long double (&TransitionMatrix)[5][5], unsigned long int GenomeLength, unsigned long int Chromosome_Length_Threshold,string ChrName,const string& readFilename)
{
 double startP[]={0.25,0.25,0.25,0.25};
 int pos=0;
 char Base[]={'A','C','G','T','N'};
 cerr << "generating Genome with length "<< GenomeLength <<endl;
 ofstream readFile(readFilename.c_str());
 streambuf *placeHolder = cout.rdbuf();
 cout.rdbuf(readFile.rdbuf());
 cout<<">"<<ChrName<<"1"<<endl;
 for(unsigned long int Base_Produced=1;Base_Produced <= GenomeLength;Base_Produced++)
 {
	double resultP[]={0,0,0,0};
	double random=Uniform_randomNumber(0,1);//(rand()%100)/100.0 ;
	bool ex=false;
	double tmpSum=0;
	pos=0;
	for( pos=0;!ex && pos <5;)	
	{
		double difference = tmpSum - random ;
		if( (difference) >= (-1)*0.0001 ){// || (pos == 4)) {
			pos+=(pos==0)?1:0;
			resultP[--pos]=1;
			ex=true;
		} else
		{
			tmpSum+=startP[pos];
			pos++;
		}				
	}
	if( pos >=0 &&  pos<= 3)
	{
		cout<<Base[pos];
	} else {
		cerr << "Pos is: "<< pos << " tmpsum : "<< tmpSum << " random: " << random << " Base_Production: "<< Base_Produced<< endl;
		cerr << startP[0] << "," << startP[1] << "," << startP[2] << "," << startP[3] << ","<< startP[4] <<  endl;
		return ;
	}
	if( (Base_Produced != GenomeLength) && ((Base_Produced%100 == 0) || (Base_Produced%Chromosome_Length_Threshold == 0)) )
	{
		cout<<endl;
		if (Base_Produced%Chromosome_Length_Threshold == 0)
		{
			cout<<">"<<ChrName<<(Base_Produced/Chromosome_Length_Threshold)+1<<endl;
			for (int index=0;index<4;index++)
			{
				startP[index]=0.25;
			}
			continue;
		}
	}
	std::copy(resultP, resultP+4,startP);
	for (int i =0 ; i<4 ;i++)
	{
		double sum=0;
		for(int j=0;j<4;j++)
		{
			sum+=(startP[j]*TransitionMatrix[j][i]);
		} 
		resultP[i]=sum;
	}
	std::copy(resultP, resultP+4,startP);
 }
 readFile.close();
 cout.rdbuf(placeHolder);
}


int main (int argc, const char* argv[])
{

	const char *inputFile , *outputFile, *ChrName, *inputError;	
	int lenB=20,lenE=50,Step=1,cov=1,ran=0;
	unsigned long int NumberOfBacteriaReads=1000000, NumberOfGenomicReads=1000000;
	unsigned long int GenomeLength=1000000000;
 	unsigned long int Chromosome_Length_Threshold=250000000; // unwise sam default chromosome length format
	unsigned long int Matrix[5][5]={0};
	long double TransitionMatrix[5][5]={0};
	int q_shift_up = 0;
	int nbOfMismatches = 0;
	bool highQuality = false;
	bool newGenome = false;
	double Sigma5Prime=0.9, Sigma3Prime=0.93, delta=0.02, p=0.3, divRate=0.00;
	inputError="";
	outputFile="Genome.Simulated";
	ChrName="";	
	int chromosomeNumber=0;
	poptContext optCon; /* context for parsing command-line options */
	enum option_tags { opt_none,opt_highQuality , opt_version , opt_newGenome} ;

        struct poptOption allSwitches[] =
	 {
			{ "inputFile",                          'i', POPT_ARG_STRING, &inputFile, opt_none, "Sequence FILE", "FILE" },
                        { "outputFile",                          'o', POPT_ARG_STRING, &outputFile, opt_none, " Output FILE", "FILE" },
			{ "NbBacterialReads",                         'b', POPT_ARG_LONG | POPT_ARGFLAG_SHOW_DEFAULT, &NumberOfBacteriaReads, 0,"Number  of bacterial reads", "n" },
			{ "NbGenomicReads",                         'g', POPT_ARG_LONG | POPT_ARGFLAG_SHOW_DEFAULT, &NumberOfGenomicReads, 0,"Number of wanted Genomic reads", "n" },
			{ "lenBegin",                          'f', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT, &lenB, 0,"Generate reads from this length ", "n" },
			{ "lenEnd",                         't', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT, &lenE, 0,"Generate reads up to this length", "n" },
			{ "Step",                         's', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT, &Step, 0,"Increament Steps between lenBegin and lenEnd of reads", "n" },
		  	{ "numberOfMismatches",                          'z', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT, &nbOfMismatches, 0, "introduced mismatches per read", "n" },
                        { "GenomeLength",                         'l', POPT_ARG_LONG | POPT_ARGFLAG_SHOW_DEFAULT, &GenomeLength, 0,"simulated GenomeLength", "n" },
		  	{ "newGenome",                          'n', POPT_ARG_NONE , 0,opt_newGenome, "Create a new Genome", "n" },
                        { "QualityShifUp",                         'q', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT, &q_shift_up, 0,"the amount to shift every read quality score by", "n" },
			{ "ChromosomName",                         'k', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &ChrName, opt_none,"read name", "FILE" },
			{ "Divergence",                         'd', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &divRate, 0,"Divergence rate in %", "n" },
			{ "DeaminationRatio5Prime",                         '5', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &Sigma5Prime, 0,"Deamination ratio at 5Prime", "n" },
			{ "DeaminationRaio3Primet",                         '3', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &Sigma3Prime, 0,"Deamination ratio at 3Prime", "n" },
		        { "Deamination_Middle",                         'm', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &delta, 0,"Deamination ratio at middle", "n" },
			{ "ProbabilityOverhang",                         'p', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &p, 0,"Success Probability at Overhang", "n" },
		        { "ErrorProgile",                         'e', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &inputError, opt_none,"illumina Error Profile", "FILE" },
			{ "HighestQuality",                  'h',POPT_ARG_NONE,0,opt_highQuality, "reads with highest quality score, without any sequencing error"},

			{ "Chromosom_Length",                         'c', POPT_ARG_LONG | POPT_ARGFLAG_SHOW_DEFAULT, &Chromosome_Length_Threshold, 0,"Chromosome Length", "n" },
			 { "Version", 'v', POPT_ARG_NONE,0,opt_version, "Print version control"},
                        POPT_AUTOHELP   { NULL, 0, 0,0,0,"Help Options", NULL } 
       };

 	poptContext opt = poptGetContext("readSim", argc, argv, allSwitches, 0);
	poptSetOtherOptionHelp(opt, "[OPTIONS]*");
	

	if( argc <= 1 ) {
		poptPrintHelp(opt, stderr, 0 ); 
		return 1; 
	}
        int rc = -1 ;

	for( rc = poptGetNextOpt( opt ) ; rc > 0 ; rc = poptGetNextOpt(opt)){
                 switch( rc )
                 {
                         case opt_version:
                                 std::cerr << poptGetInvocationName(opt) << " revision 0.0.1" << std::endl ;
                                 return 0;

                         case opt_highQuality:
                                 std::cerr << "Producing read with highest qualityScore " << std::endl ;
                                 highQuality = true ;
                                 break;
			case opt_newGenome:
				newGenome = false;
				break;
                 }
         }

        
	if ( rc != -1) {
		std::cerr << "need argument" <<std::endl; 
		poptPrintUsage(opt, stderr, 0); 
		return 0;
	}

	
	ios_base::sync_with_stdio(false);
	if(p < 0 || p >1)
	{
		cerr<<"P must be between 0 and 1"<<endl;
		return 1;
	}

	if (Step == 0)
	{	
		cerr<<"Step must be bigger than or equal to one (Step >=1)" <<endl;
		return 1;
	}	
	string transitionFilename=string(outputFile)+".Transition";
	string matrixFilename=string(outputFile)+".Matrix";
	string readFilename="";
	if (newGenome)
	{
		readFilename=string(outputFile)+".tmp.fa";
		if (!(reading(Matrix,string(inputFile),4000))) 
		{
			cerr << "Error in reading reference file" <<endl;
			return 1;
		} else 
		{
			PrintMatrix(Matrix,matrixFilename);
			PrintFraction(Matrix,TransitionMatrix,transitionFilename);
			generate_Genome(TransitionMatrix,GenomeLength,Chromosome_Length_Threshold,ChrName,readFilename);
			cerr <<"Done by generate temporary genome " << endl;
		}
	} else 
	{
		readFilename=string(inputFile);
	}

	//string diverged_genome="Diverged_"+readFilename;
	
/*	if(!(Diverge_Genome(readFilename,diverged_genome,divRate))){
		cerr << "Error in Divergence function" << endl;
		return 1;
	}*/

	//Divide my genome to different Chromosoms of length 250 milion bp
	if (!(Generate_reads_Post_Process(readFilename,NumberOfGenomicReads,lenB,lenE,Step,GenomeLength,Chromosome_Length_Threshold,Sigma5Prime,Sigma3Prime,delta,p,NumberOfBacteriaReads,highQuality,inputError,q_shift_up,nbOfMismatches)) ){
		cerr << "Error in Generate_reads function" <<endl;
		return 1;
	}

	for( int readLength=BeginReadLength; readLength<=EndReadLength;readLength+=Step){
	//for lenB to lenE Step
	//Producing read
	//GenerateRead(ReadFileName,NumberOfGenomicRead);
	//Mismatch
	//MisMatch(Read,nbOfMismatches);
	//Damage
	//Damage(Read,Sigma5Prime,Sigma3Prime,delta,p)
	//SeqError
	//SeqError(inputError,q_shift_up,highQuality);
	}


	cerr << "I'm done! move on" << endl;
	return 0;

}
