#include<bits/stdc++.h>
#include <iostream>
#include<fstream>
#include <string>
#include <limits>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
using namespace std;

#define _USE_MATH_DEFINES
#include <math.h>
#include "cgal.h"

#define MAX_REC_LEN 1024
  ///  edited  ///
#define readSize 101

typedef struct qEntry{
    long long int pos;
    double prob;
}qE;

typedef struct information
{
    double prob;
    long int pos;
    char *contigField;
    long int readLen;
    char* conName;
} info;

typedef struct buildString
{
    long int pos;
    string cigar;
    string read;
} bs;

typedef struct readCigar
{
    string cigar;
    string read;
    int contigNo;
} CR;

 /// edited  ///

int noContigs=0;
int noReads=0;
long int contigLength=0;
int maxReadLength=0;

long int totalContigLength=0;
vector<char*> contigs;
vector<char*> contigNames;
vector<long int> contigLengths;


FILE *contigFile;
FILE *mapFile;
FILE *summaryFile;
FILE *outFile;
ofstream infoFile; /// edited

char *contigFileName;
const char *mapFileName;
static int iteration = 0;

long int *insertCounts;
int maxInsertSize=0;
int MAX_INSERT_SIZE=100000;
double insertSizeMean;
double insertSizeVar;


long int errorTypes[5][5];
long int baseCounts[5];
long int *errorPos;
long int *inPos;
long int *inLengths;
long int *delPos;
long int *delLengths;
long int *readLengths;

long int *effectiveLengths;

double errorTypeProbs[5][5];
double baseErrorRates[5];
double *errorPosDist;
double *inPosDist;
double *inLengthDist;
double *delPosDist;
double *delLengthDist;
double *insertLengthDist;

double *noErrorProbs;
int tmpCount=0;
int toAlign;

long int erroredReads=0;
long int uniqueMappedReads=0;
long int discardedReads=0;
long int totalCount,unCount;

char tempCigar[500], tempMD[500];
char noErrorCigar[500], noErrorMD[500];
///edited
void writeInfoToFile(vector<info> multiMap, map<long int, string> &, map<long int, string> &);
void writeBsToFile(vector<bs> &bsCollection, int count);
void createBs(long int pos, char* readString, char* cigar, vector<bs> &B);
void createInfo(double errorProb1, long int pos1, char* contigField1, char* readString1, vector<info> &M);
void unixSort();
void calculateProbability();


void initInsertCounts(int max)
{
    maxInsertSize=max;
    insertCounts=new long int[maxInsertSize];
    for(int i=0; i<maxInsertSize; i++)
    {
        insertCounts[i]=1;
    }
}

void updateInsertCounts(int index)
{

    if(index<=0)
        return;
    if(index<maxInsertSize)
    {
        insertCounts[index]++;
    }
    else
    {

        if(index>MAX_INSERT_SIZE)
        {
            discardedReads++;
            return;
        }
        int tempInsertSize=max(maxInsertSize*2,index);
        long int *tempCounts=new long int[maxInsertSize];
        for(int i=0; i<maxInsertSize; i++)
        {
            tempCounts[i]=insertCounts[i];
        }
        insertCounts=new long int[tempInsertSize];
        for(int i=0; i<maxInsertSize; i++)
        {
            insertCounts[i]=tempCounts[i];
        }
        for(int i=maxInsertSize; i<tempInsertSize; i++)
        {
            insertCounts[i]=1;
        }

        insertCounts[index]++;
        maxInsertSize=tempInsertSize;
        delete []tempCounts;

    }

}

void initErrorTypes(int readLength)
{
    for(int i=0; i<5; i++)
        for(int j=0; j<5; j++)
            errorTypes[i][j]=1;

    for(int i=0; i<5; i++)
        baseCounts[i]=1;

    errorPos=new long int[readLength];
    inPos=new long int[readLength];
    inPos=new long int[readLength];
    inPos=new long int[readLength];
    inPos=new long int[readLength];
    inLengths=new long int[readLength];
    delPos=new long int[readLength];
    delLengths=new long int[readLength];
    readLengths=new long int[readLength];

    for(int i=0; i<readLength; i++)
    {
        errorPos[i]=1;
        inPos[i]=1;
        inLengths[i]=1;
        delPos[i]=1;
        delLengths[i]=1;
        readLengths[i]=1;
    }
}


int getLength(char *read)
{

    int i=0;
    while(read[i])
    {
        if(read[i]=='A')
            baseCounts[0]++;
        else if(read[i]=='C')
            baseCounts[1]++;
        else if(read[i]=='G')
            baseCounts[2]++;
        else if(read[i]=='T')
            baseCounts[3]++;
        else
            baseCounts[4]++;

        i++;
    }

    return i;
}


void processErrorTypes(char *cigar, char *md, char *read, int strandNo)
{

    int readLength=getLength(read);
    readLengths[readLength-1]++;

    if(strcmp(md,noErrorCigar)!=0)
        erroredReads++;
    else
        return;


    int mdLength=strlen(md)-5;
    int tempLength=0;

    char *temp;
    int index=0,totalLength=0;


    int curIndex=0;
    int *inserts=new int[readLength];

    for(int i=0; i<readLength; i++)
    {
        inserts[i]=0;
    }

    int cigarLength=strlen(cigar);
//	char *tempCigar=new char[cigarLength];
    char cigarChar;

    strcpy(tempCigar,cigar);


    temp=strtok(tempCigar,"IDMS^\t\n ");

    while(temp!=NULL)
    {

        tempLength=atoi(temp);
        totalLength+=strlen(temp);
        cigarChar=cigar[totalLength];

        if(cigarChar=='M')
        {
            index+=tempLength;
            curIndex+=tempLength;
        }
        else if(cigarChar=='I' || cigarChar=='S')
        {
            if(strandNo==0)
            {
                inPos[index]++;
                inLengths[tempLength-1]++;

            }
            else
            {

                inPos[readLength-index-1]++;
                inLengths[tempLength-1]++;
            }

            inserts[curIndex]=tempLength;

            index+=tempLength;
        }
        else if(cigarChar=='D' )
        {
            if(strandNo==0)
            {
                delPos[index]++;
                delLengths[tempLength-1]++;

            }
            else
            {

                delPos[readLength-index-1]++;
                delLengths[tempLength-1]++;
            }
        }
        totalLength++;
        temp=strtok(NULL,"IDMS^\t\n ");
    }


    strcpy(tempMD,md);

    strtok(tempMD,":");
    strtok(NULL,":");


    index=0,totalLength=0,tempLength=0;

    int f, t;

    while((temp=strtok(NULL,"ACGTN^\t\n "))!=NULL)
    {
        tempLength=strlen(temp);


        totalLength+=tempLength;


        if(totalLength<mdLength)
        {
            char from=md[5+totalLength];


            if(from=='^')
            {
                totalLength++;
                index+=atoi(temp);
                for(int i=totalLength; i<mdLength; i++)
                {
                    from=md[5+totalLength];
                    if(from=='A' || from=='C' || from=='G' || from=='T'|| from=='N')
                        totalLength++;
                    else
                        break;

                }
            }
            else if(from=='A' || from=='C' || from=='G' || from=='T'|| from=='N')
            {
                totalLength++;
                index+=atoi(temp)+1;



                curIndex=0;
                for(int i=0; i<index; i++)
                {
                    curIndex+=inserts[i];
                }
                char to=read[index-1+curIndex];

                if(strandNo==0)
                    errorPos[index-1+curIndex]++;
                else
                    errorPos[readLength-index-curIndex]++;


                switch(from)
                {
                case 'A':
                    f=0;
                    break;
                case 'C':
                    f=1;
                    break;
                case 'G':
                    f=2;
                    break;
                case 'T':
                    f=3;
                    break;
                default:
                    f=4;
                }

                switch(to)
                {
                case 'A':
                    t=0;
                    break;
                case 'C':
                    t=1;
                    break;
                case 'G':
                    t=2;
                    break;
                case 'T':
                    t=3;
                    break;
                default:
                    t=4;
                }

                if(f==t)
                {

                }
                else

                    errorTypes[f][t]++;

            }
            else
                break;
        }

    }
    delete []inserts;

}

void computeProbabilites()
{

    int errorCount=0;

    for(int i=0; i<5; i++)
    {
        errorCount=0;
        for(int j=0; j<5; j++)
        {
            errorCount+=errorTypes[i][j];
        }
        for(int j=0; j<5; j++)
        {
            errorTypeProbs[i][j]=(double)errorTypes[i][j]/errorCount;
        }

        baseErrorRates[i]=errorCount/(double)baseCounts[i];
    }

    double sum=0;
    for(int i=0; i<4; i++)
        sum+=baseErrorRates[i];

    for(int i=0; i<4; i++)
    {
        baseErrorRates[i]=4*baseErrorRates[i]/sum;///???
    }

    baseErrorRates[4]=1;

    for(int i=maxReadLength-1; i>0; i--)
    {
        readLengths[i-1]=readLengths[i]+readLengths[i-1];
    }

    errorPosDist=new double[maxReadLength];

    for(int i=0; i<maxReadLength; i++)
    {
        errorPosDist[i]=(double)errorPos[i]/readLengths[i];
    }

    inPosDist=new double[maxReadLength];

    for(int i=0; i<maxReadLength; i++)
    {
        inPosDist[i]=(double)inPos[i]/readLengths[i];
    }

    inLengthDist=new double[maxReadLength];

    int inCount=0;

    for(int i=0; i<maxReadLength; i++)
    {
        inCount+=inLengths[i];
    }

    for(int i=0; i<maxReadLength; i++)
    {
        inLengthDist[i]=(double)inLengths[i]/inCount;
    }

    delPosDist=new double[maxReadLength];

    for(int i=0; i<maxReadLength; i++)
    {
        delPosDist[i]=(double)delPos[i]/readLengths[i];
    }

    delLengthDist=new double[maxReadLength];

    int delCount=0;

    for(int i=0; i<maxReadLength; i++)
    {
        delCount+=delLengths[i];
    }

    for(int i=0; i<maxReadLength; i++)
    {
        delLengthDist[i]=(double)delLengths[i]/delCount;
    }


    insertLengthDist=new double[maxInsertSize];

    long int insCount=discardedReads;

    sum=0;

    for(int i=0; i<maxInsertSize; i++)
    {
        insCount+=insertCounts[i];
        sum+=i*insertCounts[i];
    }
    insertSizeMean=sum/insCount;

    sum=0;

    for(int i=0; i<maxInsertSize; i++)
    {
        insertLengthDist[i]=(double)insertCounts[i]/insCount;

        sum+=insertCounts[i]*(insertSizeMean-i)*(insertSizeMean-i);
    }

    insertSizeVar=sum/insCount;


    noErrorProbs=new double[maxReadLength];

    double noErrorProb=1.0;

    for(int i=0; i<maxReadLength; i++)
    {
        noErrorProb*=(1-errorPosDist[i]-inPosDist[i]-delPosDist[i]);
        noErrorProbs[i]=noErrorProb;
    }

    effectiveLengths=new long int[maxInsertSize];

    for(int i=0; i<maxInsertSize; i++)
    {
        effectiveLengths[i]=-1;
    }

    long int totalContigLength=0;
    for(int i=0; i<contigLengths.size(); i++)
    {
        totalContigLength+=contigLengths[i];
    }
    effectiveLengths[0]=totalContigLength;

}

double dnorm(double x,double mean, double variance)
{
    double val=1/sqrt(M_PI*2*variance);
    val*=exp(-((x-mean)*(x-mean))/(2*variance));
    return val;
}

void processMapping(char *line)
{

    char * temp;
    char *qname, *rname, *mapq, *contigField;
    int	pos,flag;
    char * cigar, * readString; // * md, *nhstring;

    char md[500];
    char nhstring[500];

    int nh;

    int strandNo=0;
    char *temp1,*temp2,*temp3;

    qname=strtok(line,"\t");

    temp=strtok(NULL,"\t");
    flag=atoi(temp);


    strandNo=(flag&16)>>4;
    contigField =strtok(NULL, "\t");
    temp=strtok(NULL,"\t");
    pos=atoi(temp);



    cigar=strtok(NULL,"\t");




    temp=strtok(NULL,"\t");

    readString=strtok(NULL,"\t");

    int insertSize=atoi(temp);

    while((temp=strtok(NULL,"\t\n"))!=NULL)
    {
        if(temp[0]=='M' && temp[1]=='D')
        {
            strcpy(md,temp);
        }
        else if(temp[0]=='I' && temp[1]=='H')
        {
            strcpy(nhstring,(temp+5));
            nh=atoi(nhstring) ;
        }

    }


    if(nh==1 && md[5]!='^')
    {
        updateInsertCounts(insertSize);
        processErrorTypes(cigar,md,readString,strandNo);
        uniqueMappedReads++;
    }


}

long int getEffectiveLength(int insertSize)
{
    if(insertSize<0)
        return effectiveLengths[0];

    if(insertSize>=maxInsertSize)
    {
        long int effectiveLength=0;
        for(int i=0; i<contigLengths.size(); i++)
        {
            if(contigLengths[i]>=insertSize)
                effectiveLength+=(contigLengths[i]-insertSize+1);
        }
        return effectiveLength;

    }
    if(effectiveLengths[insertSize]==-1)
    {
        long int effectiveLength=0;
        for(int i=0; i<contigLengths.size(); i++)
        {
            if(contigLengths[i]>=insertSize)
                effectiveLength+=(contigLengths[i]-insertSize+1);
        }
        effectiveLengths[insertSize]=effectiveLength;
    }
    return effectiveLengths[insertSize];
}

double computeErrorProb(char *cigar, char *md, char *read, int strandNo)
{

    int readLength=strlen(read);



    double errorProb=noErrorProbs[readLength-1];

    if(md[5]=='^')
        return errorProb;


    char tempMD[1000], tempCigar[1000];

    int mdLength=strlen(md)-5;
    int tempLength=0;

    char *temp;
    int index=0,totalLength=0;


    int curIndex=0;
    int *inserts=new int[readLength];

    for(int i=0; i<readLength; i++)
    {
        inserts[i]=0;
    }

    int cigarLength=strlen(cigar);
    char cigarChar;

    strcpy(tempCigar,cigar);

    temp=strtok(tempCigar,"IDM^\t\n ");

    while(temp!=NULL)
    {

        tempLength=atoi(temp);
        totalLength+=strlen(temp);
        cigarChar=cigar[totalLength];

        if(cigarChar=='M')
        {
            index+=tempLength;
            curIndex+=tempLength;
        }
        else if(cigarChar=='I')
        {
            int i;
            if(strandNo==0)
            {
                //look up insert probs
                i=index;

            }
            else
            {
                i=readLength-index-1;
            }

            errorProb=errorProb*inPosDist[i]*inLengthDist[tempLength-1]/(1-errorPosDist[i]-inPosDist[i]-delPosDist[i]);

            inserts[curIndex]=tempLength;

            index+=tempLength;
        }
        else if(cigarChar=='D')
        {
            int i;
            if(strandNo==0)
            {
                i=index;
                //	look up delete probs
            }
            else
            {
                i=readLength-index-1;
            }

            errorProb=errorProb*delPosDist[i]*delLengthDist[tempLength-1]/(1-errorPosDist[i]-inPosDist[i]-delPosDist[i]);
        }
        totalLength++;
        temp=strtok(NULL,"IDM^\t\n ");
    }


    strcpy(tempMD,md);

    strtok(tempMD,":");
    strtok(NULL,":");

    index=0,totalLength=0,tempLength=0;

    int f, t;

    while((temp=strtok(NULL,"ACGTN^\t\n "))!=NULL)
    {
        tempLength=strlen(temp);

        totalLength+=tempLength;

        if(totalLength<mdLength)
        {
            char from=md[5+totalLength];

            if(from=='^')
            {
                totalLength++;
                index+=atoi(temp);
                for(int i=totalLength; i<mdLength; i++)
                {
                    from=md[5+totalLength];
                    if(from=='A' || from=='C' || from=='G' || from=='T'|| from=='N')
                        totalLength++;
                    else
                        break;
                }
            }
            else if(from=='A' || from=='C' || from=='G' || from=='T'|| from=='N')
            {
                totalLength++;
                index+=atoi(temp)+1;


                curIndex=0;
                for(int i=0; i<index; i++)
                {
                    curIndex+=inserts[i];
                }
                char to=read[index-1+curIndex];

                int i;
                if(strandNo==0)
                    i=index-1+curIndex;
                else
                    i=readLength-index-curIndex;


                errorProb=errorProb*errorPosDist[i]/(1-errorPosDist[i]-inPosDist[i]-delPosDist[i]);


                switch(from)
                {
                case 'A':
                    f=0;
                    break;
                case 'C':
                    f=1;
                    break;
                case 'G':
                    f=2;
                    break;
                case 'T':
                    f=3;
                    break;
                default:
                    f=4;
                }

                switch(to)
                {
                case 'A':
                    t=0;
                    break;
                case 'C':
                    t=1;
                    break;
                case 'G':
                    t=2;
                    break;
                case 'T':
                    t=3;
                    break;
                default:
                    t=4;
                }

                if(f==t)
                {


                }
                else
                {
                    //errorTypeProb
                    errorProb*=baseErrorRates[f]*errorTypeProbs[f][t];
                }
            }
            else
                break;

        }

    }

    delete []inserts;
    return errorProb;
}


void printHelp()
{

    cout<<"cgal v0.9.9-beta"<<endl;
    cout<<"----------------"<<endl;
    cout<<endl;
    cout<<"cgal - computes likelihood"<<endl;
    cout<<"Usage:"<<endl;
    cout<<"cgal [options] <contigfile>"<<endl;
    cout<<endl;
    cout<<"Required arguments:"<<endl;
    cout<<"<contigfile.sam>\t Assembly file in FASTA format"<<endl;
    cout<<endl;
    cout<<"Options:"<<endl;
    cout<<"-h [--help]\t\t Prints this message"<<endl;
    cout<<endl;
    cout<<"Output: "<<endl;
    cout<<"(In file out.txt) <numberContigs> <totalLikelihood> <mappedLikelihood> <unmappedLikelihood> <noReads> <noReadsUnmapped>"<<endl;
    cout<<"<numberContigs>\t\t Number of contigs"<<endl;
    cout<<"<totalLikelihood>\t Total log likelihood value"<<endl;
    cout<<"<mappedLikelihood>\t Likelihood value of reads mapped by the mapping tool"<<endl;
    cout<<"<unmappedLikelihood>\t Likelihood value corresponding to reads not mapped by alignment tool"<<endl;
    cout<<"<noReads>\t\t Total number of paired-end reads"<<endl;
    cout<<"<noReadsUnmapped>\t Number of reads not mapped by the alignment tool"<<endl;
    cout<<endl;
    exit(1);

}
