
#include <iostream>
#include <fstream>
#include <string>


using namespace std;

int convert_to_good_ply(string filein, string fileout, int v_count_int, int num_faces)
{
    ifstream badmeshfile(filein);
    ofstream goodteapot(fileout);

    if(!badmeshfile || !goodteapot)
    {
        cout << "error opening files" << endl;
        return 1;
    }

 
    string line;
    getline(badmeshfile,line);
    goodteapot << "kcply" << endl;
    getline(badmeshfile,line);
    goodteapot << "element vertex "<< v_count_int*2+2 << endl;
    
    
    getline(badmeshfile,line);
    goodteapot << line << endl;
    int s1,s2,s3;
    string v1v2v3,v4v5v6,newline;
    for(int i = 0; i<v_count_int+1; i++)
    {
        getline(badmeshfile,line);
        s1 = line.find(" ");
        s2 = line.find(" ", s1+1);
        s3 = line.find(" ",s2+1);
        v1v2v3 = line.substr(0,s3); // maybe + 1 here.
        v4v5v6 = line.substr(s3+1,line.length()+1);
        newline = v1v2v3 + "\n" + v4v5v6;
        goodteapot << newline << endl;
    }
    for(int i=0; i<num_faces; i++)
    {
        getline(badmeshfile,line);
        newline = line;
        goodteapot << line << "\n";
    }
    return 0;
}



