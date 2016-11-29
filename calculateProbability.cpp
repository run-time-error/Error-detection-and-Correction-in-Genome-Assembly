#include<bits/stdc++.h>
#include<string>
#define MAX_POS 1154134
#define contigSize 101
using namespace std;

typedef struct qEntry{
    long long int pos;
    double prob;
}qE;


main()
{
    /*ifstream cLen("ContigLength.txt");
    //vector<int> contigLens;
    int x;
    while(!cLen.eof())
    {
        cLen>>x;
        cout<<x<<" ";
        //contigLens.push_back(x);
    }*/
    ifstream inFile("infoOutput.txt");
    //ifstream inFile("contig311data.txt");
    ifstream contigF("contigName.txt");
    qE q;
    string contigName, temp;
    int pos,contigLine=0;
    int contig,t1,t2,prevContig=0,tempPos;
    pair<bool,double> tempPair;
    int loopNo=0;
    while(!inFile.eof())
    {
        while(contigLine<=prevContig)
        {
            contigF>>pos>>contigName>>temp;
            contigLine++;
        }
        cout<<prevContig<<endl;
        string fileName="outCtg"+std::to_string(prevContig);
        ofstream outFile("ProbFolder/"+fileName+".txt");
        ofstream out2 ("posFile.txt");
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
                prevContig = contig;
                break;
            }
        }
        for(map<int,pair<bool, double> >::iterator it=M.begin(); it!= M.end();it++){
            //cout << it->first << ' ' <<  contig << ' ' << it->second.second << endl;
        }
        //cout<<contigLens[prevContig]<<endl;

        /// loop for each contig c
        double probability = 0.0;
        ofstream probf("probtest.txt");
        queue<qE> Q;
        for(int i=0;i<MAX_POS;i++){
//            cout << probability << endl;

            map<int,pair<bool, double> > :: iterator it = M.find(i);
            if(it != M.end()){
                // map e ache
                double temp = fabs(log(it->second.second));
                //cout<<log(it->second.second)<<"\t"<<it->second.second<<endl;
                probability = probability + temp;
                probf<<probability<<"\t"<<log(it->second.second)<<endl;
                qE qe;
                qe.pos = i;
                qe.prob = it->second.second;

                Q.push(qe);
//                cout<<"%%"<<Q.size()<<endl;
//                out2<<"push "<<qe.pos<<" "<<Q.size()<<endl;

            }
            //write into output file
            double p;
            //if(probability < 1e-320) probability = 1e-320;

            //int sz = (int) Q.size();

            if(!Q.empty()){
                int sz = (int) Q.size();
                p = probability/sz;
                p = -p;
                //cout << "q size not 0 inf"<<Q.size() << endl;
                //cout<<p<<endl;
            }
            else{
                p = -probability;
                //cout<<"q size 0"<<endl;
                //cout<<p<<endl;
            }
            //cout<<"$$"<<probability<<endl;

            outFile<< i << " " << p <<" "<<Q.size()<< endl;
            loopNo++;
            if(!Q.empty()){
                qE qeTemp = Q.front();
                if( i - qeTemp.pos== (contigSize-1)){
                    Q.pop();

                    probability -= fabs(log(qeTemp.prob));

                    if(probability<-1E-320)
                        probability=0;
    //                out2<<"pop "<<qeTemp.pos<<" "<<probability<<" "<<Q.size()<<endl;
                }
            }


        }
        outFile.close();
        out2.close();
        probf.close();

    }

    inFile.close();
    contigF.close();
    cout<<"LoopNo: "<<loopNo;
    ///end loop for contig c
    return 0;
}
