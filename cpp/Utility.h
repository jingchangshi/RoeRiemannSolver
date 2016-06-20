#include <iostream>
#include <fstream>
#include <string>
using namespace std;

void write1D(string FileName, double* InArray1D_ptr, \
        int InArray1DN, int Precision){
    ofstream OutFile(FileName, ofstream::out);
    if(!OutFile) {cerr << "Error: Open output file!" << endl;}
    OutFile.precision(Precision);
    for(int i = 0; i < InArray1DN; i++){
        OutFile << scientific << *(InArray1D_ptr+i) << endl;
    }
}

void write2D(string FileName, double* InArray2D_ptr, \
    int InArray2DRol, int InArray2DCol, int Precision, string Separator){
    ofstream OutFile(FileName, ofstream::out);
    if(!OutFile) {cerr << "Error: Open output file!" << endl;}
    OutFile.precision(Precision);
    for(int i = 0; i < InArray2DRol; i++){
        for(int j = 0; j < InArray2DCol-1; j++){
            OutFile << scientific << *(InArray2D_ptr+i*InArray2DCol+j) << Separator;
        }
        OutFile << \
            *(InArray2D_ptr+i*InArray2DCol+InArray2DCol-1) << endl;
    }
}
