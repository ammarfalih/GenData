// General Data classes.
// GenData version 1.08
// programmed by Ammar Falih A. Hasan
// 2006-2008
//for more information see GenData_1_08.pdf

#include<stdio.h>

#ifndef GENDATA_H
#define GENDATA_H

const bool ADD_CORRECTION=false;

//helper functions
bool add(bool a, bool b);//add logical bools

//Interleaver types
enum Interleaver {SAME,REVERSE,RANDOM,ROWCOLUMN,CIRCLEDIST};

//Convert : converts between binary and integer numbers.
class Convert{
	int memSize;//measured by No. of bits
public:
	//public members
	bool* array;
	int state;
	//constructors
	Convert(int memorySize);
	//access functions
	Convert& operator=(const Convert &copy);
	int getSize() const;
	//helper functions
	bool* inttobool();//returns array from state
	bool* inttobool(int newState);//assigns state then calculates and returns array
	int booltoint();//returns state from array
	int booltoint(bool *copy);// assigns array and returns state
	void shiftRight();//shifts the array to the right
	void shiftLeft();//shifts the array to the left
	//destructors
	~Convert();
};

//SqMatrix : assists in building matrix related interleavers
class SqMatrix{
	int **mat;
	int rows,columns;
	void setValues();
public:
	//public members
	int *inter;
	//constructor
	SqMatrix(int size);
	//access functions
	int getRows();
	int getColumns();
	int getSize();
	void print();
	//helper functions
	void setRowColumn();
	//destructor
	~SqMatrix();
};

//InterMat : formats interleaver matrix array.
class InterMat{
	int *matrix;
	int matSize;
public:
	//constructor
	InterMat(Interleaver l, int blockSize, int firstVal=1);//fills the matrix map with the given interleaver technique
	InterMat(Interleaver l, int blockSize, int p, bool similed);//performs a simile interleaver
											// p is the delay unit given by 2^(delay elements)-1
											// make sure that blockSize is divisible by p
	//access functions
	int* getMat()const;//returns a handle to the matrix from this object.
	//helper functions
	bool isConsistent();//retrurns true if the interleaver is consitent
							// used for diagnostic operations
	//file handling
	bool saveMat(const char* fileName);
	bool loadMat(const char* fileName);
	//destructor
	~InterMat();
};

//StateTran : 
//State transitions for a 4-bit memory encoder with the optimum transfer function by default
// this transfer function can be modified by giving a suitable Convert class equation
class StateTran{
	int state;
	bool *tailMat;//zero state leading path matrix
	// each state have 4 direction transitions, two previous ones(for 0 or 1 input) and two next ones.
	// [state][direction]
	int **direction;
	// each state have a parity according to input
	//[state][input]
	bool **parity;
	void fillTailTable();
	void createTable(bool g1[], bool g2[]);
public:
	//constructors
	StateTran(int s=0);// default 23/33 encoder
	StateTran(bool g1[], bool g2[]);//make sure the sizes of the given arrays
					// are 5 elements, also make sure to supply a feedback element 1 for g1[0] 
	//access functions
	void reset(int newState=0);
	int getState();
	//Helper functions
	bool getParity(int input);
	int getNextState(int input);
	int getPreviousState(int input);
	double tailer(int state);//returns a value that leads to the nearest path to
							// the zero state according to the given state
	//destructor
	~StateTran();
};

//DataStr : Data Stream class.
class DataStr{
	double *bits;
	int datSize;
public:
	//constructors
	DataStr(int s);//creates unintialized DataStr of size s
	DataStr(const DataStr &dataStream);// default assignment constructor
	DataStr(DataStr *dataStream);//copies objects
	//access functions
	DataStr& operator=(const DataStr &dataStream);//assignment function
	void setBits();//assigns random bits values
	void setBits(double choice);//assigns the given double to all the bits
	void setBits(const double* copy);//assigns a copy of the given values
	double* getBits();//returns a pointer to the bits.
	double* getCBits() const;
	int getSize() const;
	//helper functions
	DataStr permute(const InterMat *mat,bool option);//returns permuted version 
										// option = true->interleave ,or  false->deinterleave
	DataStr encode(StateTran *encoder);	//returns parity bits from this DataStr 
										//object according to the given encoder
	DataStr hard();//returns hard limited version of this data(+1.0,-1.0)
	void tail(StateTran *encoder);//tails last bits of this object to lead to the zero state
	int compare(DataStr *original);//compares this object with original, returns No. of errors
	void halfPuncture(bool position);// sets data with 0 between one value and the other
									// position=false-> start from first value,position=true->otherwise
	//file handling
	bool saveData(const char* fileName);
	bool loadData(const char* fileName);
	//destructor
	~DataStr();
};

// BCJR : Decoder component algorithm.
class BCJR{
	DataStr *msg,*prty,*LeIn,*LeOut;
	StateTran *states;
	double Lc;//Lc=4*SNR_channel
	int N;// size of the Stream of Data(DataStr) Handled in the operations
	bool LeInSet;
	double** A;//alpha k of state s : Ak(s) : A[k][s]
	double** B;//beta k of s : Bk(s) : B[k][s]
	double*** G;// gamma k of s' and s : Gk(s',s) : G[k][s][s']. Note that s' have only two values
	double*** Ge;// Extrinsic Gamma function, same as above.
	void createTables();// creates alpha, beta,gamma and gammaE according to inputs.
public:
	//constructor
	BCJR(DataStr *message,DataStr *parity,double Lchannel,StateTran *encoder);
	//access functions
	void setLeIn(DataStr *extrinsic);// sets LeIn, to be used in successive iterations
	DataStr getLeOut();// calculate LeOut bits and return a copy from the result
	DataStr getFinal();// adds (input message+input extrinsic+output extrinsic)
	//destructor
	~BCJR();
};

//Channel : simulates a channel with a given noise.
class Channel{
	DataStr *message,*parity1,*parity2;
	double SNR_c;
	double gaussian(double variance);
	bool punctured;
	bool singleData;
public:
	//constructor
	Channel(DataStr *msg,DataStr *par1,DataStr *par2,
			double SNR_channel,bool punctured=false);
	Channel(DataStr *msg, double SNR_channel);
	//access functions
	DataStr getDat(int choice);
	bool isPunctured();
	//Helper functions
	void applyNoise();//applies noise to all streams
	DataStr decode(int iterations,InterMat *interleaver,StateTran *encoder);
								//decodes the three DataStr objects with iterative BCJR algorithm
	//destructor
	~Channel();
};

//Property Options enumeration
enum PropOp {DESCRIPTION,INT_SIZE,PACKETS,SAMPLES,START_DB,END_DB,
				INTERLEAVER};

//PropertyUse: Holds data necessary for implementation by other classes
class PropertyUse{
	char fileName[16];
	char description[512];
	char size[8];
	char packets[8];
	char samples[8];
	char startDB[8];
	char endDB[8];
	char interleaver[16];
public:
	//constructors
	PropertyUse();// use default: properties.txt
	PropertyUse(char* fName);

	//access functions
	void setOption(PropOp op,char* opDat);
	//void setSize(int s);
	//void setPackets(int p);
	//void setSamples(int smp);
	//void setStartEndDB(int start,int end);
	//void setInterleaver(Interleaver i);
	char* getOption(PropOp op);
	//int getSize();
	//int getPackets();
	//int getSamples();
	//int getStartDB();
	//int getEndDB()
	//Interleaver getInterleaver();

	//file handling
	bool saveData();
	bool loadData();

	//destructors
	~PropertyUse();
};

//Use: an implementation of GenData library
class Use{
	FILE *file;
	DataStr *msg,*par1,*par2,*inter;
	InterMat *mat;
	StateTran *encoder;
	int SIZE;//size of block
	bool fileOpened;
	bool TYPE;//type of Use, true=uncoded, false=coded
	bool punctured;
	bool warningPrinted;
	bool encoderInitialized;
	int processCoded(double SNR, int noOfPackets);
	int iIter;
public:
	//constructors
	Use(int packetSize);//uncoded version
	Use(int interleaverSize, bool puncturedState, Interleaver intType=RANDOM,
				bool SIMILED=false, int p=15, bool defaultEncoder=true, int iIterations=3);//coded version
	//helper functions
	double processAndSave(double SNR_dB,int noOfPackets);//Processes single Pe event,
										//saves it to the opened stream and returns the result.
	bool initializeEncoder(bool g1[],bool g2[]);//initialize encoder
	//file handling
	bool openStream(const char* fileName);
	bool closeStream();
	//destructors
	~Use();
};

#endif
