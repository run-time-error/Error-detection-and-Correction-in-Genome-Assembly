#include<bits/stdc++.h>
//#include "errorStruct.h"
#define readLength 101
using namespace std;
/// M = 0,  I = 1,  D = 2
/// A = 0,  C = 1,  G = 2,  T = 3,  N = 4
struct E{
    int type[3];
    int mChar[5];
    int iChar[5];
    vector<string> insertedS;
    E()
    {
        for(int i=0;i<3;i++)
            type[i] = 0;
        for(int i=0;i<5;i++)
        {
            mChar[i] = 0;
            iChar[i] = 0;
        }
    }
};
int start = 232251;
int finish = 232840;
map<char, int> my;

pair< vector<int>, vector<char> > processCigar(string s)
{
    vector<int> n;
    vector<char> c;
    int p = 0;
    for(int i=0;i<s.size();i++)
    {
        if(isalpha(s[i]))
        {
            c.push_back(s[i]);
            n.push_back(p);
            p=0;
        }
        else
        {
            p =10*p +(s[i] - '0');
        }
    }
    return make_pair(n,c);
}

main()
{
    ifstream in("bsFolder2/a.txt");
    ofstream out("bsFolder2/saureus3.txt");
    out<<"scaffold 4.2\t"<<start<<"\t"<<finish<<endl;
    string str[4];
    my['A'] = 0;
    my['C'] = 1;
    my['G'] = 2;
    my['T'] = 3;
    int eLength = finish-start+1;

    E arr[eLength];

    pair<vector <int>, vector<char> > vp;
    vector<int> n;
    vector<char> c;
    int iterator1 = 0;
    while(in>>str[0]>>str[1]>>str[2] >> str[3])
    {
//        cout<<str[0]<<endl;
//        cout<<str[1]<<endl;
//        cout<<str[2]<<endl;
        cout<<iterator1++<<endl;
        int epos = atoi(str[1].c_str()) - start;

//        cout << epos<< endl;
        string cigar = str[2];
        vp = processCigar(cigar);
        n = vp.first;
        c = vp.second;
        int iter = 0;
        for(int i=0;i<n.size() && epos<eLength;i++){
            int ni = n[i];
            char ci = c[i];
            //cout<<ni<<"\t"<<ci<<endl;
            if(ci == 'M')
            {
                //cout<<ni<<endl;
                for(int j = epos; j <= (epos + ni) && j<eLength ; j++)
                {
                    arr[j].type[0]++;
                    char nt = str[3][iter];

                    iter++;
                    //cout<< iter++<<"\t"<<nt<<endl;
                    arr[j].mChar[my[nt]]++;
                    //cout<< j << " " << arr[j].type[0]<<endl;

                }
                epos += ni;
            }
            else if( ci == 'D')
            {
                //int ni = n[i];
                for(int j = epos; j <= (epos + ni) && j<eLength; j++)
                {
                    arr[j].type[2]++;
                }
                epos += ni;
            }
            else if(ci == 'I')
            {
               //int ni = n[i];

               arr[epos].type[1]++;

               string temp = str[3].substr(iter,ni);
               //cout<<epos<<endl;
               arr[epos].insertedS.push_back(temp);

               iter += ni;
            }

        }



    }
    cout<<"\tM  I  D\t\tA  C  G  T"<<endl;
    out<<"\tM  I  D\t\tA  C  G  T"<<endl;
    for(int i=0;i<eLength;i++)
    {
        cout<<i<<"\t";
        out<<i<<"\t";
        for(int j=0;j<3;j++){
            cout<<arr[i].type[j]<<"  ";
            out<<arr[i].type[j]<<"  ";
        }
        out<<"\t";
        cout<<"\t";
        int m = 0,ind = 0;
        for(int j=0;j<4;j++)
        {
            cout<<arr[i].mChar[j]<<"  ";
            out<<arr[i].mChar[j]<<"  ";
            if(m < arr[i].mChar[j])
            {
                m = arr[i].mChar[j];
                ind = j;
            }
        }
        if(arr[i].insertedS.size()>0)
        {
            cout<<"\t";
            out<<"\t";
            for(int j=0;j<arr[i].insertedS.size();j++){
                cout<<arr[i].insertedS[j]<<"  ";
                out<<arr[i].insertedS[j]<<"  ";
            }
        }
        cout<<endl;
        out<<endl;
    }

    return 0;
}
