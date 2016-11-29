#include<iostream>
#include<fstream>
#include<vector>
#include<algorithm>
using namespace std;

typedef struct information{
    double prob;
    long int pos1;
    long int pos2;
    long int contigField1;
	long int contigField2;
	
	bool operator < (const struct information & i1)const{
		if(contigField1 == i1.contigField1){
			return pos1 < i1.pos1;
		}
		return contigField1 < i1.contigField1;
	}	
	
}info;


int main()
{
	ofstream outFile;
	ifstream inFile;
	
	inFile.open("info.txt");
	outFile.open("sortedInfo.txt");
	
	info temp;
	vector <info> infoData;
	
	while(inFile >>  temp.contigField1 >>
		temp.contigField2 >> temp.pos1 >>
		temp.pos2 >> temp.prob)
	
	{
		infoData.push_back(temp);		
	}
	cout << infoData.size() << endl;
	sort(infoData.begin(), infoData.end());	
	outFile << "CF1\tPos1\tProb\n"; 
	for(int i=0;i<infoData.size();i++){
		info temp = infoData[i];
		outFile << temp.contigField1 << "\t" << temp.pos1 << "\t\t"<< temp.prob << endl;
	}
	cout << "The End" << endl;
	outFile.close();
	inFile.close();
}
