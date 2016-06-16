#include <iostream>
#include "Solution.h"
#include "Utility.h"
using namespace std;

int main(int argc, char* argv[]){
    int N = 501;
    int Flag_IC = 1;
    int Precision = 18;
    double CFL = 0.8;
    Solution Sol_Obj(N, Flag_IC, CFL);
    double* X_ptr = new double[N];
    double XL = 0.0, XH = 1.0;
    int ItN;
    double Time, TimeEnd;

    for(int i = 0; i < N; i++){
        *(X_ptr+i) = (XH-XL) / (N-1) * i;
    }
    Sol_Obj.setPts(X_ptr);
    write1D("PtsSol.txt", Sol_Obj.PtsSol_ptr, Sol_Obj.SolN, Precision);
    Sol_Obj.setIC();
    Sol_Obj.W2U_EntireField();

    ItN = 10000;
    Time = 0.0;
    TimeEnd = 0.2;
    for(int i = 0; i < ItN; i++){
        cout << "Iteration " << i << ": " << Time << "s" << endl;
        Sol_Obj.calcTimeStep();
        Time = Time + Sol_Obj.TimeStep;
        if(Time > TimeEnd){
            Sol_Obj.TimeStep = TimeEnd - Time;
            Time = TimeEnd;
            break;
        }
        Sol_Obj.marchTime();
    }
    // 
    write2D("SolW.txt", Sol_Obj.SolW_ptr, Sol_Obj.SolN, 3, Precision, " ");
    write2D("SolU.txt", Sol_Obj.SolU_ptr, Sol_Obj.SolN, 3, Precision, " ");
    return 0;
}
