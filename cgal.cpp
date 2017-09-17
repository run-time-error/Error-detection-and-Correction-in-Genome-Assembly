
#include "cgalHeader.h"
static int countInfoFile=0;
void printGlobalMapData();
map<long int, vector<CR> > globalMap;
map <int , vector<double> > mapOfProbListToPosition;
int startCentralTendency = 232000;
int endCentralTendency = 233000;

void checkCentralTendency(deque <qE> Q2, int position) {


    cout << "Here in the function " << endl;
    vector <double> probList;
    for(int i=0;i<Q2.size();i++){
        probList.push_back(Q2[i].prob);
    }
    mapOfProbListToPosition[position] = probList;
}
void writeFileForCentralTendencyCheck(){
    cout << "Here" << endl;
    ofstream fileWriter ("centralTendencyCheck.txt");
    map<int , vector<double> >:: iterator it;
    for(it = mapOfProbListToPosition.begin(); it != mapOfProbListToPosition.end();it++ ){
        vector <double> probList = it->second;
        int position = it->first;
        for(int i=0;i<probList.size();i++){

            if(position >= startCentralTendency && position <= endCentralTendency){
                fileWriter << position << probList[i] << endl;
            }
        }
        fileWriter << "~~~~~~~~~~~~~~~" << endl;
    }
    fileWriter.close();

}


double computeLikelihood(const char *file)
{

    mapFile=fopen(file, "r");
    char *line1= new char[MAX_REC_LEN];
    char *line2= new char[MAX_REC_LEN];

    char *qname1,*qname2,preqname1[500],preqname2[500];

    int MAX_FILE_READ=MAX_REC_LEN/sizeof(line1[0]);

    double sum=0.0;
    double logsum=0.0;

    char *temp,*contigField1, *contigField2;
    char *rname1, *rname2;
    int	pos1,pos2,flag,strandNo1, strandNo2, insertSize1, insertSize2;
    char *cigar1, *cigar2, *readString1, *readString2, md1[1000], md2[1000];


    double insertSizeProb;

    double errorProb1, errorProb2;

    int noUnmappedReads=0;
    int noUniqueMappedReads=0;
    int noMultMappedReads=0;

    preqname1[0]=preqname2[0]='*';
    preqname1[1]=preqname2[1]=0;


    vector<info> multiMapProb1, multiMapProb2;
    vector<bs> bsCollection1,bsCollection2;
    map<long int, string> cigarMap1,cigarMap2;
    map<long int, string> readMap1,readMap2;

    int fileCount=0;

    int it=0;

    while(fgets(line1, MAX_FILE_READ, mapFile)!=NULL)
    {

        if(line1[0]=='@')
            continue;
//????
        if(fgets(line2, MAX_FILE_READ, mapFile)==NULL)
            break;

        ///r001 99 ref 7 30 8M2I4M1D3M = 37 39 TTAGATAAAGGATACTG *
        ///r003 0 ref 9 30 5S6M * 0 0 GCCTAAGCTAA * SA:Z:ref,29,-,6H5M,17,0;

        qname1=strtok(line1,"\t");
        temp=strtok(NULL,"\t");
        flag=atoi(temp);


        strandNo1=(flag&16)>>4;///checking the 4th bit
        contigField1 = strtok(NULL, "\t");

        temp=strtok(NULL,"\t");
        pos1=atoi(temp);


        cigar1=strtok(NULL,"\t");///cigar string
        //cout<<cigar1<<endl;

        temp=strtok(NULL,"\t");



        insertSize1=atoi(temp);


        readString1=strtok(NULL,"\t");///sequence string


        while((temp=strtok(NULL,"\t\n"))!=NULL)
        {
            if(temp[0]=='M' && temp[1]=='D')
            {
                strcpy(md1,temp);
            }
        }

//		cout<<insertSize1<<" "<<cigar1<<" "<<md1<<endl;


///second of the pair

        qname2=strtok(line2,"\t");
        temp=strtok(NULL,"\t");
        flag=atoi(temp);


        strandNo2=(flag&16)>>4;
        contigField2 = strtok(NULL,"\t");

        temp=strtok(NULL,"\t");
        pos2=atoi(temp);



        cigar2=strtok(NULL,"\t");
        //cout<<cigar2<<endl;
        cout<<contigField1<<endl;
        temp=strtok(NULL,"\t");

        insertSize2=atoi(temp);

        readString2=strtok(NULL,"\t");


        while((temp=strtok(NULL,"\t\n"))!=NULL)
        {
            if(temp[0]=='M' && temp[1]=='D')
            {
                strcpy(md2,temp);
            }
        }

//		cout<<insertSize2<<" "<<cigar2<<" "<<md2<<endl;



        int insertSize=max(insertSize1, insertSize2);


        insertSizeProb=0;

        if(insertSize>=0 && insertSize<maxInsertSize)
        {
            insertSizeProb=insertLengthDist[insertSize];
        }

        if(insertSizeProb==0)
        {
            insertSizeProb=1/(double)uniqueMappedReads;
        }


        errorProb1=computeErrorProb(cigar1,md1,readString1,strandNo1);


        errorProb2=computeErrorProb(cigar2,md2,readString2,strandNo2);


        long int totalEffectiveLength=getEffectiveLength(insertSize);


        double prob=(1/(double)(totalEffectiveLength))*insertSizeProb*errorProb1*errorProb2;


//		cout<<errorProb1<<" "<<errorProb2<<" "<<insertSizeProb<<" "<<prob<<endl;




        if(strcmp(qname1,preqname1)==0 && strcmp(qname2,preqname2)==0)
        {
            sum+=prob;
            //cout<<prob<<endl;
            ///need to add code here

            createBs(pos1, readString1, cigar1, bsCollection1);
            createBs(pos2, readString2, cigar2, bsCollection2);

            cigarMap1[pos1] = cigar1;
            cigarMap2[pos2] = cigar2;
            readMap1[pos1] = readString1;
            readMap2[pos2] = readString2;

            createInfo(errorProb1, pos1, contigField1, readString1, multiMapProb1);
            createInfo(errorProb2, pos2, contigField2, readString2, multiMapProb2);

        }
        else if(strcmp("*",preqname1)!=0 && strcmp("*",preqname2)!=0)
        {
            if(sum<1e-320 || std::isnan(sum))
            {
                sum=1e-320;
            }
            logsum+=log(sum);
            ///edited
            writeInfoToFile(multiMapProb1, cigarMap1, readMap1);
            writeInfoToFile(multiMapProb2, cigarMap2, readMap2);
            if(bsCollection1.size()>10)
                writeBsToFile(bsCollection1,fileCount++);
            it++;
            multiMapProb1.clear();
            multiMapProb2.clear();

            cigarMap1.clear();
            cigarMap2.clear();
            readMap1.clear();
            readMap2.clear();

            bsCollection1.clear();
            bsCollection2.clear();

            cigarMap1[pos1] = cigar1;
            cigarMap2[pos2] = cigar2;
            readMap1[pos1] = readString1;
            readMap2[pos2] = readString2;

            createInfo(errorProb1, pos1, contigField1, readString1, multiMapProb1);
            createInfo(errorProb2, pos2, contigField2, readString2, multiMapProb2);

            createBs(pos1, readString1, cigar1, bsCollection1);
            createBs(pos2, readString2, cigar2, bsCollection2);

            ///edited
            sum=prob;
        }

        else
        {
            sum=prob;

            cigarMap1[pos1] = cigar1;
            cigarMap2[pos2] = cigar2;
            readMap1[pos1] = readString1;
            readMap2[pos2] = readString2;

            createInfo(errorProb1, pos1, contigField1, readString1, multiMapProb1);
            createInfo(errorProb2, pos2, contigField2, readString2, multiMapProb2);

            createBs(pos1, readString1, cigar1, bsCollection1);
            createBs(pos2, readString2, cigar2, bsCollection2);

        }

        strcpy(preqname1,qname1);
        strcpy(preqname2,qname2);
        //it++;


        if(std::isinf( logsum ))
        {
            cout<<it<<endl;
            exit(1);
        }


    }
    if(sum!=0)
        logsum+=log(sum);

    cout<<"$$$"<<it<<endl;
    fclose(mapFile);

    return logsum;
}
/// edited
void writeInfoToFile(vector<info> multiMap, map<long int, string> &cigar, map<long int, string> &read)
{
    vector<double> cdf;
    vector<info> data(multiMap);
    //cout<<data.size()<<"  "<<multiMap.size()<<endl;
    double s=0,s1=0;

    for(int i=0; i<multiMap.size(); i++)
    {
        info temp = multiMap[i];
        //fprintf(infoFile, "(%ld, %ld) --> %lf\n",temp.pos1,temp.pos2,temp.prob);
        s += temp.prob;

    }
    srand(time(NULL));
    for(int i=0; i<multiMap.size(); i++)
    {
        info temp = multiMap[i];
        //cout<<temp.prob<<endl;
        multiMap[i].prob = temp.prob/s;
    }
    double r = ((double) rand() / (RAND_MAX));
    for(int i=0; i<multiMap.size(); i++)
    {
        info temp = multiMap[i];
        s1 += temp.prob;
        info temp1 = data[i];
        if(r<=s1)
        {
            iteration++;
            infoFile << temp1.contigField<< " " << temp1.pos << " " << temp1.prob << endl;
            CR ob;
            int number = atoi(temp1.contigField);
            ob.contigNo = temp1.pos;
            ob.cigar = cigar[temp1.pos];
            ob.read = read[temp1.pos];
            globalMap[number].push_back(ob);
            //cout<<iteration<<endl;
            break;
        }
    }
    /*for(int i=0;i<multiMap.size();i++){
    info temp = multiMap[i];
    //cout<<temp.prob<<endl;
    temp.prob = temp.prob/s;

    s1 += temp.prob;
    cdf.push_back(s1);
    }
    srand(time(NULL));
    double r = ((double) rand() / (RAND_MAX));
    for(int i=cdf.size()-1;i>=0;i--){
    if(cdf[i]<r)
    {
        info temp = multiMap[i+1];
        //cout<<temp.prob<<endl;

        fprintf(infoFile, "%s %s %ld %ld %e\n",temp.contigField1,temp.contigField2,temp.pos1,temp.pos2,temp.prob);
        break;
    }
    }*/
    countInfoFile++;
    cout<<"Write info to file"<<endl;
}
void printGlobalMapData()
{
    map<long int, vector<CR> >:: iterator it;
    ofstream gOut("globalOutput.txt");
    for(it = globalMap.begin(); it!=globalMap.end(); it++)
    {
        vector<CR> tempCr = it->second;
        cout << "itpos" << it->first << " " << tempCr.size() << endl;
        for(int i=0;i<tempCr.size();i++)
        {
            gOut<<it->first<<" "<<tempCr[i].contigNo<<" "<<tempCr[i].cigar<<" "<<tempCr[i].read<<endl;
        }
    }
    gOut.close();
    system("sort -n -k1,1 -k2,2 globalOutput.txt > globalOutSort.txt");
}

void createBs(long int pos, char* readString, char* cigar, vector<bs> &B)
{
    bs bstemp1;
    bstemp1.pos = pos;
    bstemp1.read = readString;
    bstemp1.cigar = cigar;
    B.push_back(bstemp1);
}

void createInfo(double errorProb1, long int pos1, char* contigField1, char* readString1, vector<info> &M)
{
    info temp1;
    temp1.prob = errorProb1;
    temp1.pos = pos1;
    temp1.contigField = contigField1;
    temp1.readLen = strlen(readString1);
    M.push_back(temp1);
}

void writeBsToFile(vector<bs> &bsCollection, int count)
{
    if(count<10){
        string s = "bsFolder2/bsOutput_";
        string t = ".txt";
        stringstream oss;
        oss<< s<<count<<t;
        cout<<oss.str()<<endl;
        //FILE *fp;
        //fp = fopen(oss.str().c_str(),"w");
        ofstream bsFile(oss.str().c_str());
        for(int i=0;i<bsCollection.size();i++)
        {
            bs b = bsCollection[i];
            bsFile <<b.pos<<" "<< b.cigar<<" "<<b.read<<endl;
            //fprintf(fp,"%ld %s %s\n",b.pos,b.cigar,b.read);
        }
        //fclose(fp);
        bsFile.close();
    }
}


void unixSort()
{
    system("sort -n -k1,1 -k2,2 info.txt > infoOutput1.txt");
}

void calculateProbability()
{

    int lenIt = 0;
    ifstream inFile("infoOutput1.txt");

    qE q;
    string contigName, temp;
    int pos,contigLine=0;
    int contig,t1,t2,prevContig=0,tempPos,curContig;
    pair<bool,double> tempPair;
    int loopNo=0;
    while(!inFile.eof())
    {

        //cout<<prevContig<<"\t";
        //string fileName="outCtg"+std::to_string(prevContig);
        stringstream ss;
        ss << "ProbFolder2/outCtg" << prevContig << ".txt";
        ofstream outFile(ss.str().c_str());
        map<int,pair<bool,double> >M;
        if(prevContig>0)
             M[tempPos] = tempPair;
        while(inFile>>contig>>q.pos>>q.prob)
        {
            if(contig == prevContig)
                M[q.pos] = make_pair(true, q.prob);
            else
            {
                tempPair = make_pair(true,q.prob);
                tempPos = q.pos;
                curContig = prevContig;
                prevContig = contig;
                break;
            }
        }

        if(curContig < noContigs)
        {
            /// loop for each contig c
            double probability = 0.0;
            ofstream probf("probtest.txt");
            queue<qE> Q;
            deque <qE> deq;
            for(long int i=0;i<contigLengths[curContig];i++){
                map<int,pair<bool, double> > :: iterator it = M.find(i);
                if(it != M.end()){
                    // map e ache
                    double temp = fabs(log(it->second.second));
                    cout<<log(it->second.second)<<"\t"<<it->second.second<<endl;
                    probability = probability + temp;
                    probf<<probability<<"\t"<<log(it->second.second)<<endl;
                    qE qe;
                    qe.pos = i;
                    qe.prob = it->second.second;

                    Q.push(qe);
                    deq.push_back(qe);

                }
                double p;

                if(!Q.empty()){
                    int sz = (int) Q.size();
                    checkCentralTendency(deq, i);
                    p = probability/sz;
                    p = (-1)*p;
                }
                else{
                    p = (-1)*probability;
                }

                outFile<< i << " " << p <<" "<<Q.size()<< endl;
                loopNo++;
                if(!Q.empty()){
                    qE qeTemp = Q.front();
                    if( i - qeTemp.pos== (readSize-1)){
                        Q.pop();
                        deq.pop_front();
                        probability -= fabs(log(qeTemp.prob));

                        if(probability<-1E-320)
                            probability=0;
                    }
                }
            }
            if(curContig == 178)
            {
                writeFileForCentralTendencyCheck();
            }

        }
        outFile.close();

    }

    inFile.close();
    cout<<"LoopNo: "<<loopNo;
}

int main(int argc, char *argv[])
{
    /*	input contig file name, read file name
    	contig file - fasta format
    	read file - fastq format
    */


    if(argc<2)
        printHelp();

    if(strcmp(argv[1],"--help")==0 || strcmp(argv[1],"-h")==0)
        printHelp();



    contigFileName=argv[1];

    /// edited
    infoFile.open("info.txt");

    contigFile=fopen(contigFileName, "r");
    outFile=fopen("out.txt", "w");


    if (contigFile == NULL)
    {
        printf("Can't open contig file\n");
        exit(1);
    }
    char *line= new char[MAX_REC_LEN];

    char *line1= new char[MAX_REC_LEN];
    char *line2= new char[MAX_REC_LEN];

    int read;
    int MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);


    long int bufferLength=1024;

    char *contig=new char[bufferLength];
    contig[0]='\0';
    char *newcontig;
    char *contigName;
    contigLength=0;


    long int tempContigLength=0;


    while(fgets(line, MAX_FILE_READ, contigFile)!=NULL)///parsing genome.fastq file
    {
        if(line[0]==';')///if any comment found
        {
            continue;
        }
        else if(line[0]=='>')///if contigname found
        {
            contigName=new char[strlen(line)];
            strcpy(contigName,line+1);
            contigName[strlen(contigName)-1]='\0';
            contigNames.push_back(contigName);///contignames holds all contigname
            if(contigLength>0) /// init at 1078 with 0
            {
                noContigs++;
                contigs.push_back(contig);
                contigLengths.push_back(contigLength);
                totalContigLength+=contigLength;
                contigLength=0;
                bufferLength=1024;
                contig=new char[bufferLength];
                contig[0]='\0';
            }
        }
        else
        {
            read=strlen(line);
            tempContigLength=contigLength;///initially 0
            if(read<MAX_FILE_READ-1)///init with MAX_REC_LEN = 1024
            {
                contigLength+=(read-1);/// ???
            }
            else
            {
                contigLength+=MAX_FILE_READ-1;
                read++;/// ???

            }
            if(contigLength>bufferLength)
            {
                bufferLength=max(bufferLength*2,contigLength+1);
                newcontig=new char[bufferLength];
                strcpy(newcontig,contig);
                line[read-1]='\0';
                strcat(newcontig, line);
                delete []contig;
                contig=newcontig;
            }
            else
            {
                line[read-1]='\0';
                strcpy(contig+tempContigLength, line);
            }

        }

    }
    noContigs++; /// increment for last contig
    contigs.push_back(contig);
    contigLengths.push_back(contigLength);
    totalContigLength+=contigLength;



    fclose(contigFile);


//	cout<<"after reading contig file"<<endl;

    /*
    	use bfast or some other tool to map reads and save mapping
    */

    mapFileName="myout.sam";

    mapFile=fopen(mapFileName, "r");

    summaryFile=fopen("stat.txt","r");

    if (mapFile == NULL)
    {
        printf("Can't open map file\n");
        exit(1);
    }

    int count=0;

    fscanf(summaryFile,"%ld %ld %d %d %d",&totalCount, &unCount, &toAlign, &maxReadLength, &MAX_INSERT_SIZE);


    initInsertCounts(5000);


    initErrorTypes(maxReadLength);


    itoa(maxReadLength, noErrorCigar, 10);
    strcpy(noErrorMD,"MD:Z:");
    strcat(noErrorMD,noErrorCigar);
    strcat(noErrorCigar,"M");

//	cout<<"after allocation"<<endl;

    int noMatches=0;

    while(fgets(line1, MAX_FILE_READ, mapFile)!=NULL)
    {

        if(line1[0]=='@')
            continue;

        processMapping(line1);

        count+=1;

    }

    fclose(mapFile);
    fclose(summaryFile);

//	cout<<"after parameter calculation"<<endl;

    computeProbabilites();

//	cout<<"after prob calculation"<<endl;


    /*
    	for(int i=0;i<maxReadLength;i++)
    	{
    		cout<<i<<" "<<errorPos[i]<<" "<<errorPosDist[i]<<endl;
    	}
    	for(int i=0;i<5;i++)
    	{
    		for(int j=0;j<5;j++)
    		{
    			cout<<errorTypes[i][j]<<" ";
    		}
    		cout<<endl;
    	}
    	for(int i=0;i<5;i++)
    	{
    		for(int j=0;j<5;j++)
    		{
    			cout<<errorTypeProbs[i][j]<<" ";
    		}
    		cout<<endl;
    	}
    	for(int i=0;i<maxReadLength;i++)
    	{
    		cout<<i<<" "<<inPos[i]<<" "<<inPosDist[i]<<endl;
    	}

    	for(int i=0;i<maxReadLength;i++)
    	{
    		cout<<i<<" "<<inLengths[i]<<" "<<inLengthDist[i]<<endl;
    	}

    	for(int i=0;i<maxReadLength;i++)
    	{
    		cout<<i<<" "<<delPos[i]<<" "<<delPosDist[i]<<endl;
    	}

    	for(int i=0;i<maxReadLength;i++)
    	{
    		cout<<i<<" "<<delLengths[i]<<" "<<delLengthDist[i]<<endl;
    	}
    */

    /*
    	cout<<erroredReads<<endl;

    	for(int i=0;i<500;i++)
    	{
    		cout<<insertLengthDist[i]<<endl;
    	}


    	cout<<noErrorProb<<endl;
    	cout<<uniqueMappedReads<<endl;
    */



    double val1=computeLikelihood(mapFileName);
    /// edited ///
    cout << "balchal" << endl;
    infoFile.close();
    unixSort();
    calculateProbability();
    printGlobalMapData();
//	cout<<"after val 1"<<endl;

    double val2=0;//computeLikelihood("unmappedOut.sam");

//	cout<<"after val 2"<<endl;


    double value=val1;
    /*if(unCount==0)
    {
    	value=val1;
    }
    else if(unCount>toAlign)
    {
    	val2=val2/toAlign*unCount;
    	value=val1+val2;
    }
    else
    {
    	value=val1+val2;
    }
    */
    cout<<value<<endl;
    cout<<contigLength<<endl;

    fprintf(outFile,"%d\t%f\t%f\t%f\t%ld\t%ld\n",noContigs,value,val1,val2,totalCount,unCount);

    fclose(outFile);
    FILE *cFile = fopen("ContigName","w");
    for(int i=0; i<contigNames.size(); i++)
    {
        char* name = contigNames[i];
        fprintf(cFile,"%d %s\n",i,name);
    }
    fclose(cFile);
    FILE *cFile1 = fopen("ContigLength","w");
    for(int i=0; i<contigNames.size(); i++)
    {
        long int len = contigLengths[i];
        fprintf(cFile1,"%ld\n",len);
    }
    fclose(cFile1);
    ///edited


    return 0;
}
