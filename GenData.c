// implementations of the GenData.h V1.08 classes
// 2006-2008

#include "GenData.h"
#include "random.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stddef.h>

//create a static Random class object
//rndRandom created static becuase we don't want it to be initialized each time we go in
Random rndNumber;

bool add(bool a, bool b){
   //return a==b ? false : true; //<-- old one
   return a^b;// a XOR b <-- faster
}
double correct(double MaxVal,double vals[],int n){
	double Es = 0;//exponentials sum
	if(MaxVal==HUGE_VAL||MaxVal==-HUGE_VAL)
		return 0;
	int i;
	for(i=0;i<n;i++)
		if(vals[i]!=-HUGE_VAL)
			Es += exp(vals[i]-MaxVal);
	return log(Es);
}
double max(double vals[],int n,bool AddCorrect){
	double MaxVal = -HUGE_VAL;
	int i;
	for(i=0;i<n;i++)
		if(MaxVal<vals[i])
			MaxVal=vals[i];
	if(MaxVal==-HUGE_VAL||MaxVal==HUGE_VAL)
		return MaxVal;
	if(AddCorrect)
		return MaxVal+correct(MaxVal,vals,n);
	return MaxVal;
}

//Convert class implementation
Convert::Convert(int m){
	if(m>0&&m<33)
		memSize = m;
	else{
		printf("err in Convert: value for bits out of bounds\n");
		exit(0);
	}
	//bool array[memSize]; <--  this doesn't work, as the compiler refuses to assign variable array size
	array = new bool[memSize];
	int x;
	for(x=0;x<memSize;x++){
		array[x]=0;
	}
	state = 0;
}
Convert& Convert::operator=(const Convert &copy){
	if(memSize==copy.getSize()){
		state = copy.state;
		int i;
		for(i=0;i<memSize;i++)
			array[i]=copy.array[i];
	}
	else
		printf("err in Convert: assignment with different sizes\n");
	return *this;
}
int Convert::getSize() const{
	return memSize;
}
bool* Convert::inttobool(){
	int st = state;int next;
	int i;
	for (i=0;i<memSize;i++){
		next = st >> 1;// shift the state to the right
		if ((next << 1) == st)
			array[i] = false;
		else
			array[i] = true;
		st = next;
	}
	return array;
}
bool* Convert::inttobool(int s){
	state = s;
	return inttobool();
}
int Convert::booltoint(){
	state = 0;// state = 0000
	int i;
	for (i=0;i<memSize;i++)
		if (array[i] == true)
			state |= (1 << i);// state = state | (0001<<x)
	return state;
}
int Convert::booltoint(bool *copy){
	int i;
	for (i=0;i<memSize;i++){
		array[i] = copy[i];
	}
	return booltoint();
}
void Convert::shiftRight(){
	bool temp = array[0];
	int i;
	for(i=0;i<(memSize-1);i++)
		array[i]=array[i+1];
	array[memSize-1]=temp;
	booltoint();
}
void Convert::shiftLeft(){
	bool temp = array[memSize-1];
	int i;
	for(i=(memSize-1);i>0;i--)
		array[i]=array[i-1];
	array[0]=temp;
	booltoint();
}

Convert::~Convert(){
	delete[] array;
}

//SqMatrix class implementation
SqMatrix::SqMatrix(int size){
	double dSize=(double)size;
	int iRoot=(int)(sqrt(dSize));
	int r=iRoot+1,iDev=1;
	double dR,dDev;
	while((r--)>1){
		dR=(double)r;
		dDev=dSize/dR;
		iDev=(int)dDev;
		if(dDev==((double)iDev))
			break;
	}
	rows=r;
	columns=iDev;
	mat = new int*[rows];
	int i;
	for(i=0;i<rows;i++)
		mat[i]=new int[columns];
	inter = new int[size+1];
	setValues();
}
void SqMatrix::setValues(){
	int i,j,count=1;
	for(i=0;i<rows;i++){
		for(j=0;j<columns;j++){		
			mat[i][j]=count;
			inter[count]=count;
			count++;
		}
	}
}
int SqMatrix::getRows(){
	return rows;
}
int SqMatrix::getColumns(){
	return columns;
}
int SqMatrix::getSize(){
	return rows*columns;
}
void SqMatrix::print(){
	int i,j;
	for(i=0;i<rows;i++){
		for(j=0;j<columns;j++)
			printf("%d ",mat[i][j]);
		printf("\n");
	}
	//setRowColumn();
	for(i=1;i<=getSize();i++)
		printf("%d ",inter[i]);
	printf("\n");
}
void SqMatrix::setRowColumn(){
	int i,j,count=1;
	for(j=0;j<columns;j++)
		for(i=0;i<rows;i++)
			inter[count++]=mat[i][j];
}	
SqMatrix::~SqMatrix(){
	delete[] inter;
	int i;
	for(i=0;i<rows;i++)
		delete[] mat[i];
	delete[] mat;
}

//InterMat class  implementation
int compareDist(int val1,int val2,int SIZE){
	if(val1==val2)
		return 0;
	if(val1>SIZE||val2>SIZE){
		printf("err in compareDist: vals>SIZE\n");
		return 0;
	}
	int inner,outer;
	inner = abs(val1-val2);
	if(val2>val1)
		outer = SIZE-val2+val1;
	else
		outer = SIZE-val1+val2;
	if(inner<outer)
		return inner;
	else
		return outer;
}
InterMat::InterMat(Interleaver l, int blockSize, int x){
	matSize = blockSize;
	matrix = new int[blockSize+1];
	// Do Not Use Switch-Case, the compiler refuses to pass Local variable initialization
	
	int i,j;
	//put all interleaver cases here :
	//case SAME: default operation
	for(i=0;i<=blockSize;i++)
		matrix[i]=i;
	//case REVERSE:Reverse the block 
	if(l==REVERSE){
		for(int i=1;i<=blockSize;i++)
			matrix[i]=(blockSize+1)-i;
	}
	// case RANDOM:generate a random interleaver
	if(l==RANDOM){
		bool *locations=new bool[blockSize+1];
		for(i=1;i<=blockSize;i++)
			locations[i]=false;
		for(i=1;i<=blockSize;i++){
			int count=0;
			// r is a number initially set between 0 -> blockSize-1, and continues dcreasing
			double r = ((double)(blockSize+1-i))*rndNumber.doublerandom();
			//intR is initially between 1->blockSize and decreases
			int intR = ((int)r)+1;
			while(intR!=0){
				count++;
				if(!locations[count])
					intR--;
			}
			matrix[i]=count;
			locations[count]=true;
		}
		delete[] locations;
	}
	//case ROWCOLUMN
	if(l==ROWCOLUMN){
		SqMatrix m(blockSize);
		m.setRowColumn();
		for(i=1;i<=blockSize;i++)
			matrix[i]=m.inter[i];
		for(i=1;i<=x;i++){
			for(j=1;j<=blockSize;j++){
				int v=matrix[j];
				v++;if(v>blockSize) v=1;
				matrix[j]=v;
			}
		}
	}
	//case CIRCLEDIST
	if(l==CIRCLEDIST){
		const int N = blockSize;
		if(N>2&&x<=N){
			for(i=0;i<=N;i++){
				matrix[i]=0;
			}
			int position; int L,R,highest,hPosition;
			int* distances = new int[N+1];
			if(x>N) x=1;//just in case the user provided an x larger than N
			matrix[1]=x;x++;if(x>N) x=1;
			matrix[2]=x;x++;if(x>N) x=1;
			for(i=3;i<=N;i++){
				//printf("\n%d\n",x);
				position = 1;
				for(j=0;j<=N;j++)
					distances[j]=0;
				while(position!=i){
					L=compareDist(x,matrix[position],N);
					if(matrix[position+1]!=0)
						R=compareDist(x,matrix[position+1],N);
					else
						R=compareDist(x,matrix[1],N);
					distances[position]=L+R;
					position++;				
				}
				highest=0;
				for(j=1;j<=N;j++){
					if(highest<distances[j])
						highest=distances[j];
				}
				hPosition=1;
				for(j=N;j>=1;j--){
					if(highest==distances[j]){
						hPosition=j;
						break;
					}
				}
				j=N-1;
				while(j!=hPosition){
					matrix[j+1]=matrix[j];
					j--;
				}
				matrix[hPosition+1]=x;
				x++;
				if(x>N)
					x=1;
			}
			delete[] distances;
		} else{
			printf("err in InterMat: CIRCLEDIST minimum interleaver size is 3\n");
			printf("                 and/or x>N\n");
			for(i=0;i<=1;i++){
				matrix[i]=i;
			}
		}
	}
}
InterMat::InterMat(Interleaver l,int N, int p, bool similed){
	matSize = N;
	matrix = new int[N+1];	
	int i,j;
	if(N%p==0){
		int** seq;//seq[sequence][position]
		seq = new int*[p];
		for(i=0;i<p;i++){
			seq[i] = new int[N/p];
		}
		for(i=0;i<N;i++){
			seq[i%p][i/p]=i+1;
		}
		InterMat partials(l,N/p);
		int* partialPositionArray=partials.getMat();
		int partialPosition =1;
		for(i=0;i<p;i++){
			double* seqVals;
			seqVals = new double[N/p+1];
			for(j=1;j<=N/p;j++){
				seqVals[j]=(double)seq[i][j-1];
			}
			DataStr seqDat(N/p);
			seqDat.setBits(seqVals);
			InterMat seqIntType(l,N/p,partialPositionArray[partialPosition]);
			partialPosition++; if(partialPosition>(N/p)) partialPosition=1;
			DataStr seqInterleaved = seqDat.permute(&seqIntType,true);
			double* vals;
			vals = seqInterleaved.getBits();
			for(j=0;j<N/p;j++){
				seq[i][j]=(int)vals[j+1];
			}
			delete[] seqVals;
		}
		int n = 1;
		for(i=0;i<N/p;i++){
			for(j=0;j<p;j++){
				matrix[n]=seq[j][i];
				n++;
			}
		}
		//cleanup
		for(i=0;i<p;i++){
			delete[] seq[i];
		}
		delete[] seq;
	} else{
		printf("err in StateTran: SIMILE works only on block sizes\n");
		printf("                  divisible by constraint lengths K\n");
		for(i=0;i<=N;i++)
			matrix[i]=i;
	}
}
int* InterMat::getMat()const{
	return matrix;
}
//Helper functions
bool InterMat::isConsistent(){
	int i,j;
	bool found;
	for(i=1;i<=matSize;i++){
		found = false;
		for(j=1;j<=matSize;j++){
			if(matrix[j]=i){
				found = true;
				break;
			}
		}
		if(!found)
			return false;
	}
	return true;
}
//file handling
bool InterMat::saveMat(const char* fileName){
	FILE* fileMat;
	if((fileMat=fopen(fileName,"wb"))==NULL){
		printf("err in InterMat: cannot create file\n");
		return false;
	}
	fwrite(&matSize,sizeof(int),1,fileMat);
	fwrite(matrix,sizeof(int),matSize+1,fileMat);
	fclose(fileMat);
	return true;
}
bool InterMat::loadMat(const char* fileName){
	FILE* fileMat;
	if((fileMat=fopen(fileName,"rb"))==NULL){
		printf("err in InterMat: cannot read file %s\n",fileName);
		return false;
	}
	int size=0;
	fread(&size,sizeof(int),1,fileMat);
	if(size!=matSize){
		printf("err in InterMat: file size doesn't match\n");
		printf("                 sizes given are %d and %d \n",size,matSize);
		fclose(fileMat);
		return false;
	}
	fread(matrix,sizeof(int),matSize+1,fileMat);
	fclose(fileMat);
	return true;
}
InterMat::~InterMat(){
	delete[] matrix;
}

//StateTran class implementation
StateTran::StateTran(int st){
	// create the table
	//Encoder type : 4-bit memory RSC with trans. func. in octal (g1,g2)=(23,33)  or in binary (g1,g2)=(10011,11011)
	//this is the optimum encoder which minimizes the bit error probability by maximizing the effective free distance
	bool g1[5],g2[5];//feedback = g1, forward = g2
	g1[0]=true ;g1[1]=false;g1[2]=false;g1[3]=true ;g1[4]=true;//feedback equation
	g2[0]=true ;g2[1]=true ;g2[2]=false;g2[3]=true ;g2[4]=true;//feedforward equation
	
	createTable(g1,g2);
	
	tailMat = new bool[16];
	fillTailTable();
	state = st;
}
StateTran::StateTran(bool g1[], bool g2[]){
	
	createTable(g1,g2);
	
	tailMat = new bool[16];
	fillTailTable();
	state = 0;
}
void StateTran::fillTailTable(){
	bool valueSet[16];// status of the values of tailMat matrix
	int** tailTable;//tailTable[level][value]
	tailTable = new int*[6];
	int i,j;
	for(i=0;i<6;i++){// assumed all tailMat values will be filled by less than 6 levels
		tailTable[i]=new int[32];
		for(j=0;j<32;j++)
			tailTable[i][j]=0;
	}
	for(i=0;i<16;i++)
		valueSet[i]=false;
	int level,N=1;
	tailTable[0][0]=0;//<==this is the final required state
	for(level=0;level<=4;level++){
		for(i=0;i<N;i++){
			state = tailTable[level][i];
			int prevZero = getPreviousState(0);
			tailTable[level+1][i*2]=prevZero;
			if(!valueSet[prevZero]){
				tailMat[prevZero] = false;
				valueSet[prevZero] = true;
			}
			int prevOne = getPreviousState(1);
			tailTable[level+1][i*2+1] = prevOne;
			if(!valueSet[prevOne]){
				tailMat[prevOne] = true;
				valueSet[prevOne] = true;
			}
		}
		N*=2;//N=N*2
	}
	for(i=0;i<16;i++){
		if(!valueSet[i])
			printf("err in StateTran: couldn't create tail table:%d\n",i);
	}
	//cleanup
	for(i=0;i<6;i++)
		delete[] tailTable[i];
	delete[] tailTable;
}
void StateTran::createTable(bool g1[], bool g2[]){
	//initialize
	direction = new int*[16];// [state]
	int i;
	for (i=0;i<16;i++){
		direction[i]= new int[4];//[direction]
	}
	parity = new bool*[16];// [state]
	for (i=0;i<16;i++){
		parity[i] = new bool[2];// [parity]
	}
	Convert c(4),t(4);//c=current state, t=temporary state
	for(int st=0;st<16;st++){
		c.inttobool(st);// note that array bit state is handled in reverse as in
						// x3,x2,x1,xo not to be confused with the encoder's equations
		bool y1[5],y2[5];
		for(int i=1;i<5;i++){
			y1[i]=(c.array[4-i]&g1[i]);
			y2[i]=(c.array[4-i]&g2[i]);
		}
		for(int intIn=0;intIn<2;intIn++){
			bool boolIn = (intIn==1)? true:false;
			
			bool first = false;//the first result after the input adder
			y1[0]=(boolIn&g1[0]);// note that g1[0] must always be true, so as to include the input
			int i;
			for(i=0;i<5;i++){
				first = add(first,y1[i]);
			}
			bool parityIn = false;
			y2[0]=(first&g2[0]);
			for(i=0;i<5;i++){
				parityIn = add(parityIn,y2[i]);
			}
			parity[st][intIn]=parityIn;
			t=c;t.shiftRight();t.array[3]=first;
			int nextState = t.booltoint();
			//printf("%d -%d- > %d\n",st,intIn,nextState);
			direction[st][intIn+2]=nextState;
			//if all states are to be assured to have previous state
			//then the following line shall be changed
			//instead of saving the estimations for the next state
			//the calculations must be saved for the previous states
			//the idea is to shiftLeft!!, and fill t.array[0] with
			//????????
			//direction[nextState][intIn]=st;
			// this is a try!!
			Convert p(4);
			p=c;p.shiftLeft();
			bool yp[5];//,yp2[5];
			yp[0]=p.array[0];
			for(i=1;i<=3;i++){
				yp[i]=(p.array[4-i]&g1[i]);
			//	y2[i]=(c.array[4-i]&g2[i]);
			}
			//bools addition
			yp[4]=boolIn^yp[0]^yp[1]^yp[2]^yp[3];
			p.array[0]=g1[4]&yp[4];
			int previousState=p.booltoint();
			direction[st][intIn]=previousState;
		}
	}
}
void StateTran::reset(int newState){
	if(newState>=0&&newState<16){
		state = newState;
	}
	else
		printf("err in StateTran: newState out of bounds\n");
}
int StateTran::getState(){//this is used for special purposes
	return state;
}
bool StateTran::getParity(int in){
	return parity[state][in];
}
int StateTran::getNextState(int in){
	return direction[state][in+2];
}
int StateTran::getPreviousState(int in){
	return direction[state][in];
}
double StateTran::tailer(int state){
	return (tailMat[state]==true)?1.0:-1.0;
}
StateTran::~StateTran(){
	int i;
	for (i=0;i<16;i++){
		delete[] direction[i];
	}
	delete[] direction;
	for (i=0;i<16;i++){
		delete[] parity[i];
	}
	delete[] parity;
	delete[] tailMat;
}

//DataStr class implementation
DataStr::DataStr(int s){// doesn't assign the bits values
	datSize = s;
	bits = new double[datSize+1];
	//setBits(); //<-- old
}
DataStr::DataStr(const DataStr& dataStream){// default constructor
	datSize = dataStream.getSize();
	bits = new double[datSize+1];// this is very important, bits must be initialized
	setBits(dataStream.getCBits() );
}
DataStr::DataStr(DataStr *dataStream){
	datSize = dataStream->getSize();
	bits = new double[datSize+1];
	setBits(dataStream->getBits() );
}
DataStr& DataStr::operator=(const DataStr& dataStream){// default assignment
	if(datSize == dataStream.getSize()){
		//bits = new double[datSize];// this is not a constructor, this is an access function
		setBits(dataStream.getCBits() );
	}
	else
		printf("err in DataStr: assignment with different sizes\n");
	return *this;
}
int DataStr::getSize() const{
	return datSize;
}
void DataStr::setBits(){// sets random values to the bits
	int i;
	for (i=1;i<=datSize;i++)
		bits[i] = (rndNumber.boolrandom())? 1.0 : -1.0;
}
void DataStr::setBits(double choice){
	int i;
	for (i=1;i<=datSize;i++)
		bits[i]=choice;
}
void DataStr::setBits(const double* copy){
	int i;
	for(i=1;i<=datSize;i++)
		bits[i]=copy[i];
}
double* DataStr::getBits(){// return a pointer to this object data
	return bits;
}
double* DataStr::getCBits() const{// the same but constant for compatibility with defaults
	return bits;
}
DataStr DataStr::permute(const InterMat *mat,bool option){
	double* manipulated = new double[datSize+1];
	const int* pMat = mat->getMat();
	int matLocation;
	int i;
	for(i=1;i<=datSize;i++){
		matLocation = pMat[i];
		if(option)
			manipulated[matLocation]=bits[i];
		else
			manipulated[i]=bits[matLocation];
	}
	DataStr permuted(datSize);
	permuted.setBits(manipulated);
	delete[] manipulated;
	return permuted;
}
DataStr DataStr::encode(StateTran *encoder){
	double *parBit = new double[datSize+1];
	int decision;
	encoder->reset(0);
	int i;
	for(i=1;i<=datSize;i++){
		decision = (bits[i]>0)? 1 : 0;
		parBit[i]=(encoder->getParity(decision))? 1.0 : -1.0;
		encoder->reset(encoder->getNextState(decision) );
	}
	//printf("encoder last state is %d\n",encoder->getState());
	DataStr parity(datSize);
	parity.setBits(parBit);
	delete[] parBit;
	return parity;
}
void DataStr::tail(StateTran *s){
	//StateTran s;<--not recomended 
	s->reset(0);
	//double *dats = getBits();<--no need as bits is fully accessed from within the class
	int tailState;int decision;
	int i;
	for(i=1;i<=datSize-4;i++){
		decision = (bits[i]>0)?1:0;
		tailState = s->getNextState(decision);
		s->reset(tailState);
	}
	for(i=datSize-3;i<=datSize;i++){
		bits[i]=s->tailer(tailState);
		decision = (bits[i]>0)?1:0;
		tailState = s->getNextState(decision);
		s->reset(tailState);
	}
}
int DataStr::compare(DataStr *original){
	int errors=0;
	const double* origin=original->getBits();
	int i;
	for(i=1;i<=datSize;i++){
		errors += ((bits[i]*origin[i])<0.0)? 1 : 0;
	}
	return errors;
}
DataStr DataStr::hard(){
	double* limit = new double[datSize+1];
	int i;
	for (i=1;i<=datSize;i++){
		limit[i] = (bits[i]>0)? 1.0 : -1.0;
	}
	DataStr hardened(datSize);
	hardened.setBits(limit);
	delete[] limit;
	return hardened;
}
void DataStr::halfPuncture(bool position){
	int i,p;
	if(position)
		p=0;
	else
		p=1;
	for(i=p;i<=datSize;i+=2){
		bits[i]=0;
	}
}
//file handling
bool DataStr::saveData(const char* fileName){
	FILE* fileData;
	if((fileData=fopen(fileName,"wb"))==NULL){
		printf("error in DataStr: cannot create file\n");
		return false;
	}
	fwrite(&datSize,sizeof(int),1,fileData);
	fwrite(bits,sizeof(double),datSize+1,fileData);
	fclose(fileData);
	return true;
}
bool DataStr::loadData(const char* fileName){
	FILE* fileData;
	if((fileData=fopen(fileName,"rb"))==NULL){
		printf("err in DataStr: cannot open file %s\n",fileName);
		return false;
	}
	int x=0;
	fread(&x,sizeof(int),1,fileData);
	if(x!=datSize){
		printf("err in DataStr: size of file doesn't match\n");
		fclose(fileData);
		return false;
	}
	fread(bits,sizeof(double),datSize+1,fileData);
	fclose(fileData);
	return true;
}
DataStr::~DataStr(){
	delete[] bits;
}

//BCJR class implementation
BCJR::BCJR(DataStr *message,DataStr *parity ,
			double Lchannel,StateTran *encoder){
	LeInSet=false;
	states = encoder;
	Lc = Lchannel;//Lc = 4*p*SNR=4*SNR_channel
	if(message->getSize()==parity->getSize())
		N = message->getSize();
	else{
		printf("err in BCJR: invalid message and parity size match\n");
		exit(0);
	}
	msg = new DataStr(message);
	prty = new DataStr(parity);
	LeIn = new DataStr(N);
	LeOut = new DataStr(N);
	//variables initialization
	A = new double*[N+1];
	B = new double*[N+1];
	G = new double**[N+1];
	Ge = new double**[N+1];
	for(int k=0;k<=N;k++){
		A[k]=new double[16];
		B[k]=new double[16];
		G[k]=new double*[16];
		Ge[k]=new double*[16];
		for(int s=0;s<16;s++){
			A[k][s]=0;
			B[k][s]=0;
			G[k][s]=new double[2];
			Ge[k][s]=new double[2];
			for(int in=0;in<=1;in++){
				G[k][s][in]=0;//G[k][s][s']
				Ge[k][s][in]=0;//Ge[k][s][s']
			}
		}
	}
}
void BCJR::createTables(){
	double Uk,Cp,Ys,Yp;//Uk = Cs
	double *msgP=msg->getBits();
	double *prtyP=prty->getBits();
	double *LeInP=LeIn->getBits();
	//Gamma Tabel
	int k;
	for(k=1;k<=N;k++){
		Ys=msgP[k];
		Yp=prtyP[k];
		for(int s=0;s<16;s++){
			for(int in=0;in<=1;in++){
				states->reset(s);
				int p = states->getPreviousState(in);
				states->reset(p);
				Cp =(states->getParity(in))? 1.0:-1.0;
				Uk = (in==1)? 1.0:-1.0;
				Ge[k][s][in]=0.5*Lc*Yp*Cp;
				G[k][s][in]=0.5*Uk*(LeInP[k]+Lc*Ys)+Ge[k][s][in];
			}
		}
	}
	int s,so;
	//Alpha Table
	for(s=1;s<16;s++)
		A[0][s]=-HUGE_VAL;
	A[0][0]=0;//Ao(0)=0, the rest are = -INF
	for(k=1;k<=N;k++){
		for(s=0;s<16;s++){
			states->reset(s);
			double first[2];
			for(int in=0;in<2;in++){
				so = states->getPreviousState(in);
				first[in]=A[k-1][so]+G[k][s][in];
			}
			A[k][s]=max(first,2,ADD_CORRECTION);
		}
	}
	//Beta Table
	for(s=1;s<16;s++)
		B[N][s]=-HUGE_VAL;
	B[N][0]=0;//Bn(0)=0,the rest =-INF
	for(k=N;k>=2;k--){
		for(so=0;so<16;so++){
			states->reset(so);
			double first[2];
			for(int in=0;in<2;in++){
				s = states->getNextState(in);
				first[in]=B[k][s]+G[k][s][in];
			}
			B[k-1][so]=max(first,2,ADD_CORRECTION);
		}
	}
}
void BCJR::setLeIn(DataStr *extrinsic){
	if(extrinsic->getSize()==N){
		LeIn->setBits(extrinsic->getBits() );
		LeInSet=true;
	} else
		printf("err in BCJR: invalid LeIn size\n");
}
DataStr BCJR::getLeOut(){
	if(!LeInSet){
		printf("err in BCJR: LeIn not yet set\n");
		exit(0);
	}
	createTables();
	int s,s0,s1;
	double numerator,denomenator;
	double *LeOutData = LeOut->getCBits();
	double *LeOutP = new double[N+1];
	int k;
	for(k=1;k<=N;k++)
		LeOutP[k]=LeOutData[k];
	for(k=1;k<=N;k++){
		double positive[16],negative[16];int count=0;
		for(s=0;s<16;s++){
			states->reset(s);
			s0=states->getPreviousState(0);
			s1=states->getPreviousState(1);
			positive[count]=A[k-1][s1]+B[k][s]+Ge[k][s][1];
			negative[count]=A[k-1][s0]+B[k][s]+Ge[k][s][0];
			count++;
		}
		numerator = max(positive,16,ADD_CORRECTION);
		denomenator = max(negative,16,ADD_CORRECTION);
		LeOutP[k]=numerator-denomenator;
	}
	DataStr Le(N);
	Le.setBits(LeOutP);
	delete[] LeOutP;
	return Le;	
}
DataStr BCJR::getFinal(){
	double *final=new double[N+1];
	double *msgP=msg->getBits();
	double *LeInP=LeIn->getBits();
	double *LeOutP=LeOut->getBits();
	int k;
	for(k=1;k<=N;k++)
		final[k]=Lc*msgP[k]+LeInP[k]+LeOutP[k];
	DataStr finalD(N);
	finalD.setBits(final);
	delete[] final;
	return finalD;
}
BCJR::~BCJR(){
	delete msg;
	delete prty;
	delete LeIn;
	delete LeOut;
	int k;
	for(k=0;k<=N;k++){
		for(int s=0;s<16;s++){
			delete[] G[k][s];
			delete[] Ge[k][s];
		}
		delete[] A[k];
		delete[] B[k];
		delete[] G[k];
		delete[] Ge[k];
	}
	delete[] A;
	delete[] B;
	delete[] G;
	delete[] Ge;
}	

//channel class implementation
Channel::Channel(DataStr *msg,DataStr *par1,DataStr *par2,
										double SNR_channel,bool punctured){
	SNR_c = SNR_channel;
	message = new DataStr(msg);
	parity1 = new DataStr(par1);
	parity2 = new DataStr(par2);
	singleData = false;
	this->punctured=punctured;
}
Channel::Channel(DataStr *msg, double SNR_channel){
	message = new DataStr(msg);
	parity1 = 0;
	parity2 = 0;
	SNR_c = SNR_channel;
	punctured = false;
	singleData = true;
}
double Channel::gaussian(double variance){
	double R=0;//24 Uniform Random variables
	double k;
	k = sqrt(variance/2.0);
	// add 24 uniform RV to obtain a simulation of normality
	for (int i=0;i<24;i++)
		R+=rndNumber.doublerandom();
	return (k*(R-0.5*24));
}
bool Channel::isPunctured(){
	return punctured;
}
DataStr Channel::getDat(int choice){
	if(!singleData){
		switch(choice){
			case 0:
				return DataStr(message);
			case 1:
				return DataStr(parity1);
			case 2:
				return DataStr(parity2);
			default:
				return DataStr(message);
		}
	} else
		return DataStr(message);
}
void Channel::applyNoise(){
	if(!singleData){
		int messageSize=message->getSize();
		int parity1Size=parity1->getSize();
		int parity2Size=parity2->getSize();

		int count;
		if(punctured)
			count=2;
		else
			count=1;

		int i;
		double variance=1/(2*SNR_c);
		double *messageBits=message->getBits();
		for(i=0;i<=messageSize;i++)
			messageBits[i] += gaussian(variance);
		double *parity1Bits = parity1->getBits();
		for(i=count;i<=parity1Size;i+=count)
			parity1Bits[i] += gaussian(variance);
		double *parity2Bits = parity2->getBits();
		for(i=count-1;i<=parity2Size;i+=count)
			parity2Bits[i] += gaussian(variance);
	} else{
		int i;
		int messageSize=message->getSize();
		double variance=1/(2*SNR_c);
		double *messageBits=message->getBits();
		for(i=0;i<=messageSize;i++)
			messageBits[i] += gaussian(variance);		
	}
}
DataStr Channel::decode(int iterations,InterMat *interleaver,StateTran *encoder){
	if(!singleData){
		double Lc=4*SNR_c;//Lc = 4*p/No, p is the coding rate
		int N = message->getSize();
		DataStr m = getDat(0);
		DataStr p1 = getDat(1);
		DataStr p2 = getDat(2);
		DataStr msgInter = m.permute(interleaver,true);

		if(punctured){
			p1.halfPuncture(false);
			p2.halfPuncture(true);
		}

		BCJR bcjr1(&m,&p1,Lc,encoder);
		BCJR bcjr2(&msgInter,&p2,Lc,encoder);
		DataStr LeIn1(N),LeIn2(N);
		LeIn1.setBits(0.0);
		if(iterations<=0)
			iterations = 1;
		while((iterations--)>0){
			bcjr1.setLeIn(&LeIn1);
			LeIn2 = (bcjr1.getLeOut()).permute(interleaver,true);
			bcjr2.setLeIn(&LeIn2);
			LeIn1 = (bcjr2.getLeOut()).permute(interleaver,false);
		}
		DataStr result(N);
		result = ((bcjr2.getFinal()).permute(interleaver,false));
		return result;
	} else{
		return DataStr(message);
	}
}
Channel::~Channel(){
	delete message;
	delete parity1;
	delete parity2;
}

//ProperyUse class implementation
//constructors
PropertyUse::PropertyUse(){
	char* temp="properties.gdc";
	strcpy(fileName,temp);
	saveData();
}
PropertyUse::PropertyUse(char *fName){
	strcpy(fileName,fName);
	saveData();
}	
void PropertyUse::setOption(PropOp op,char* opDat){	
	switch(op){
	case DESCRIPTION:
		strcpy(description,opDat);
		break;
	case INT_SIZE:
		strcpy(size,opDat);
		break;
	case PACKETS:
		strcpy(packets,opDat);
		break;
	case SAMPLES:
		strcpy(samples,opDat);
		break;
	case START_DB:
		strcpy(startDB,opDat);
		break;
	case END_DB:
		strcpy(endDB,opDat);
		break;
	case INTERLEAVER:
		strcpy(interleaver,opDat);
		break;
	}
	saveData();
}
char* PropertyUse::getOption(PropOp op){
	switch(op){
	case DESCRIPTION:
		return description;
		break;
	case INT_SIZE:
		return size;
	case PACKETS:
		return packets;
	case SAMPLES:
		return samples;
	case START_DB:
		return startDB;
	case END_DB:
		return endDB;
	case INTERLEAVER:
		return interleaver;
	}
	return description;
}
bool PropertyUse::saveData(){
	FILE* f;
	if(f=fopen(fileName,"wb")){
		fwrite(this,sizeof(PropertyUse),1,f);
		return true;
	} else{
		return false;
	}
}
bool PropertyUse::loadData(){
	FILE* f;
	if(f=fopen(fileName,"rb")){
		fread(this,sizeof(PropertyUse),1,f);
		return true;
	} else{
		return false;
	}
}
PropertyUse::~PropertyUse(){
}


//Use class implementation
//uncoded constructor
Use::Use(int packetSize){
	TYPE=true;SIZE=packetSize;
	msg = new DataStr(SIZE);
	encoder=0;mat=0;par1=0;inter=0;par2=0;
	punctured = false;warningPrinted=false;fileOpened=false;
}
//coded constructor
Use::Use(int interleaverSize, bool puncturedState, Interleaver intType, 
			bool SIMILED, int p, bool defaultEncoder, int iIterations){
	iIter=iIterations;
	TYPE=false;SIZE=interleaverSize;
	//Interleaver intType = iType;
	msg = new DataStr(SIZE);
	// Encoder initialization
	if(defaultEncoder){
		encoder = new StateTran();
		encoderInitialized=true;
	} else{
		encoderInitialized=false;
	}
	if(SIMILED)
		mat = new InterMat(intType,SIZE,p,true);
	else
		mat = new InterMat(intType,SIZE);
	par1 = new DataStr(SIZE);
	inter = new DataStr(SIZE);
	par2 = new DataStr(SIZE);
	punctured=puncturedState;warningPrinted=false;fileOpened=false;
}
double Use::processAndSave(double SNR_dB, int noOfPackets){
	int j,err=0;
	double SNR = pow(10,SNR_dB/10);
	if(TYPE){
		for(j=1;j<=noOfPackets;j++){
			msg->setBits();
			Channel chanUncode(msg,SNR);
			chanUncode.applyNoise();
			DataStr msgNo = chanUncode.getDat(0);
			err+=msgNo.compare(msg);
		}
	} else{
		if(encoderInitialized)
			err = processCoded(SNR,noOfPackets);
		else
			printf("err in Use: Encoder not initialized\n");
	}
	double Pe =((double)err)/(((double)SIZE)*((double)noOfPackets));
	if(fileOpened){
		fwrite(&Pe,sizeof(double),1,file);
		fwrite(&SNR_dB,sizeof(double),1,file);
		fflush(file);
	} else{
		if(!warningPrinted){
			printf("Warning in Use: cannot save processed data\n");
			warningPrinted = true;
		}
	}
	return Pe;
}
bool Use::initializeEncoder(bool g1[], bool g2[]){
	if(!encoderInitialized){
		encoder = new StateTran(g1,g2);
		encoderInitialized = true;
		return true;
	} else{
		printf("err in Use: encoder already initialized\n");
		return false;
	}
}
int Use::processCoded(double SNR, int noOfPackets){
	int j,err=0;
	for(j=1;j<=noOfPackets;j++){
		double p;
		if(punctured)
			p=0.5;//p=1/2 punctured
		else
			p=0.33333;//p=1/3 Not punctured
		double SNR_coded=SNR*p;
		msg->setBits();
		msg->tail(encoder);
		*par1 = msg->encode(encoder);
		*inter = msg->permute(mat,true);
		*par2 = inter->encode(encoder);
		Channel chanCoded(msg,par1,par2,SNR_coded,punctured);
		chanCoded.applyNoise();
		DataStr decoded = chanCoded.decode(iIter,mat,encoder);
		err+=decoded.compare(msg);
	}
	return err;
}
//file handling
bool Use::openStream(const char* fileName){
	if(!fileOpened){
		if((file=fopen(fileName,"wb"))==NULL){
			printf("err in Use: couldn't open file for writing\n");
			return false;
		}
		fileOpened=true;
		return true;
	}
	printf("err in Use: file already opened\n");
	return false;
}
bool Use::closeStream(){
	if(fileOpened){
		fclose(file);
		fileOpened = false;
		return true;
	}
	printf("err in Use: file not opened\n");
	return false;
}
Use::~Use(){
	if(fileOpened){
		closeStream();
	}
	delete msg;
	delete par1;
	delete par2;
	delete inter;
	delete mat;
	delete encoder;
}
//END
