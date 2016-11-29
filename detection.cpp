/// trying to draw a threshold for detecting errors ///

#include<bits/stdc++.h>
#define MAX_POS 1154134
#define contigSize 101
#define blockSize 500
//#define percentage 2

using namespace std;
double globalMin = 99999;
class data{
    public:
        int contigNo;
        int pos;
        double prob;
        bool operator<(const data &rhs) const { return prob < rhs.prob;}
        void print(ostream &s)
        {
            s<<contigNo<<"\t"<<pos<<"\t"<<prob<<endl;
        }
};

vector <data> dataBlock;

main()
{
    double percentage;
    ifstream inFile("outCtg311.txt");
    ifstream inFile2("outCtg311.txt");
    ofstream outFile("error311.txt");
    double t1,t2,p;
    while((inFile2>>t1>>p>>t2))
    {
        if(globalMin > p)
            globalMin =p;
        //cout<<p<<endl;
    }
    cout<<globalMin<<endl;

    int loopCount = (MAX_POS / blockSize)+1;
    for(int i=0;i<loopCount;i++)
    {
        for(int j=0;j<blockSize;j++)
        {
            data d;
            inFile >> d.contigNo >> d.pos >> d.prob;
            dataBlock.push_back(d);
            if(inFile.eof())
                break;
        }
        sort(dataBlock.begin(), dataBlock.end());
        //cout<<dataBlock[0].prob<<"\t"<<dataBlock[1].prob<<"\t"<<dataBlock[2].prob<<endl;
        //cout<<globalMin<< endl;
        double cmp = 100.0 * abs(dataBlock[0].prob - globalMin)/globalMin;
        if(cmp < 5)
            percentage = 10;
        else if(cmp < 25)
            percentage = 5;
        else
            percentage = 2;

        int dataCount = dataBlock.size()*percentage/100;
        for(int p=0;p<dataCount;p++)
        {
            //dataBlock[p].print(cout);
            dataBlock[p].print(outFile);
        }
        dataBlock.clear();
        outFile<<endl<<endl;
        //cout<<i<<endl;
    }
    inFile.close();
    outFile.close();

    return 0;

}
