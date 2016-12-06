#include<bits/stdc++.h>

using namespace std;

main()
{
    ifstream lenFile("ContigLength");
    int x;
    while(!lenFile.eof())
    {
        lenFile>>x;
        cout<<x<<endl;
    }

    return 0;
}
