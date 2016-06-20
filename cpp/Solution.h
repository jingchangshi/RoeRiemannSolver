#include <iostream>
#include <cmath>
using namespace std;

class Solution {
    public:
        Solution() = default;
        Solution(const int InN, const int InIC, double InCFL)
        :
            SolN(InN-1),
            FluN(InN),
            Flag_IC(InIC),
            CFL(InCFL)
        {}

        const int SolN;
        const int FluN;

        double* PtsSol_ptr = new double[SolN];
        double* PtsFlu_ptr = new double[FluN];
        double PtsStep;
        // Sol:   1   2   ...        SolN
        // Flu: 1   2   3 ... FluN-1      FluN
        double* SolW_ptr = new double[SolN*3];
        double* SolU_ptr = new double[SolN*3];
        double* Flu_ptr = new double[FluN*3];
        int Flag_IC;

        double TimeStep;
        double CFL;
        const double Gamma = 1.4;
        const double Gamma_Ext1 = (Gamma-1.0)/(Gamma+1.0);
//=========================================================
        void setPts(double* InPts_ptr){
            PtsFlu_ptr = InPts_ptr;
            PtsStep = ( *(PtsFlu_ptr+FluN-1) - *PtsFlu_ptr ) / SolN;
            for(int i = 0; i <= SolN-1; i++){
                *(PtsSol_ptr+i) = *(PtsFlu_ptr+i) + PtsStep/2.0;
            }
            cout << "Solution points and flux points are set!" << endl;
        }
        void setIC(){
            if (Flag_IC == 1){
                for(int i = 0; i < (int)ceil(SolN*0.3); i++){
                    *(SolW_ptr+i*3) = 1.0;
                    *(SolW_ptr+i*3+1) = 0.75;
                    *(SolW_ptr+i*3+2) = 1.0;
                }
                for(int i = (int)ceil(SolN*0.3); i < SolN; i++){
                    *(SolW_ptr+i*3) = 0.125;
                    *(SolW_ptr+i*3+1) = 0.0;
                    *(SolW_ptr+i*3+2) = 0.1;
                }
                cout << "IC is set to SOD test, Only SolW." << endl;
            }
        }
        void W2U(double* inW, double* outU){
            *(outU) = *(inW);
            *(outU+1) = *(inW) * *(inW+1);
            *(outU+2) = *(inW+2) / (Gamma-1.0) \
                + 0.5 * *(inW) * *(inW+1) * *(inW+1);
        }
        void U2W(double* inU, double* outW){
            *(outW) = *(inU);
            *(outW+1) = *(inU+1) / *(inU);
            *(outW+2) = (Gamma-1.0) \
                * ( *(inU+2) - 0.5 * *(outW) * *(outW+1) * *(outW+1) );
        }
        void W2U_EntireField(){
            for(int i = 0; i < SolN; i++){
                W2U(SolW_ptr+i*3, SolU_ptr+i*3);
            }
        }
        void U2W_EntireField(){
            for(int i = 0; i < SolN; i++){
                U2W(SolU_ptr+i*3, SolW_ptr+i*3);
            }
        }
        void calcFlux(double* inSolW, double* inSolU, double* outFlux){
            *(outFlux) = *(inSolU+1);
            *(outFlux+1) = *(inSolW+2) \
                + *(inSolW) * *(inSolW+1) * *(inSolW+1);
            *(outFlux+2) = *(inSolW+1) \
                * ( *(inSolU+2) + *(inSolW+2));
        }
        void calcRoeFlux(double* WL, double* WR, double* UL, double* UR, double* RoeFlux){
            double FluxL[3], FluxR[3];
            double EigVal[3], EigVec[9];
            double RoeVel, RoeEtha, HL, HR, RoeSS;
            double WaveA[3], UJump[3];
            // Calc Roe average
            HL = ( *(UL+2) + *(WL+2) ) / *(WL);
            HR = ( *(UR+2) + *(WR+2) ) / *(WR);
            RoeVel = ( sqrt(*(WL)) * *(WL+1) + sqrt(*(WR)) * *(WR+1) ) \
                / ( sqrt(*(WL)) + sqrt(*(WR)) );
            RoeEtha = ( sqrt(*(WL)) * HL + sqrt(*(WR)) * HR ) \
                / ( sqrt(*(WL)) + sqrt(*(WR)) );
            RoeSS = sqrt( (Gamma-1.0) * (RoeEtha-0.5*RoeVel*RoeVel) );
            //
            EigVal[0] = RoeVel - RoeSS;
            EigVal[1] = RoeVel;
            EigVal[2] = RoeVel + RoeSS;
            //
            EigVec[0*3+0] = 1.0;
            EigVec[1*3+0] = RoeVel - RoeSS;
            EigVec[2*3+0] = RoeEtha-RoeVel*RoeSS;
            EigVec[0*3+1] = 1.0;
            EigVec[1*3+1] = RoeVel;
            EigVec[2*3+1] = 0.5*RoeVel*RoeVel;
            EigVec[0*3+2] = 1.0;
            EigVec[1*3+2] = RoeVel + RoeSS;
            EigVec[2*3+2] = RoeEtha+RoeVel*RoeSS;
            //
            for(int i = 0; i < 3; i++){
                UJump[i] = *(UR+i) - *(UL+i);
            }
            WaveA[1] = (Gamma-1.0) / (RoeSS*RoeSS) \
                * ( UJump[0]*(RoeEtha-RoeVel*RoeVel) + RoeVel*UJump[1] - UJump[2] );
            WaveA[0] = (UJump[0]*(RoeVel + RoeSS)-UJump[1]-RoeSS*WaveA[1]) / (2.0*RoeSS);
            WaveA[2] = UJump[0] - WaveA[0] - WaveA[1];
            // Harten-Hyman Entropy fix
            double SS_L = sqrt(Gamma*WL[2]/WL[0]), \
                   SS_R = sqrt(Gamma*WR[2]/WR[0]);
            // Direct approach
            double Rho_StarL = WL[0] + WaveA[0], \
                   Rho_StarR = WR[0] - WaveA[2];
            double U_Star = \
                   ( WL[0] * WL[1] + WaveA[0]*(RoeVel-RoeSS) ) \
                    / ( WL[0] + WaveA[0] );
            double P_Star = (Gamma-1.0) * ( \
                   UL[2] + WaveA[0]*(RoeEtha-RoeVel*RoeSS) \
                   - 0.5*Rho_StarL*U_Star*U_Star );

/*            // PVRS Approximation
            double Rho_Bar = 0.5 * (WL[0]+WR[0]), \
                   SS_Bar = 0.5 * (SS_L+SS_R);
            double P_Star = fmax(0.0, \
                   0.5*(WL[2]+WR[2]) + 0.5*(WL[1]-WR[1])*Rho_Bar*SS_Bar \
                   ); // To avoid negative value of P_Star
            // P_Star = 0.5*(WL[2]+WR[2]) + 0.5*(WL[1]-WR[1])*Rho_Bar*SS_Bar;
            double U_Star = 0.5*(WL[1]+WR[1]) \
                          + 0.5*(WL[2]-WR[2])/(Rho_Bar*SS_Bar);
            double Rho_StarL = WL[0] + (WL[1]-U_Star)*Rho_Bar/SS_Bar, \
                   Rho_StarR = WR[0] + (U_Star-WR[1])*Rho_Bar/SS_Bar;
*/
/*            // TSRS
            double Rho_Bar = 0.5 * (WL[0]+WR[0]), \
                   SS_Bar = 0.5 * (SS_L+SS_R);
            double P0 = fmax(0.0, \
                   0.5*(WL[2]+WR[2]) + 0.5*(WL[1]-WR[1])*Rho_Bar*SS_Bar \
                   );
            double A_L = 2.0 / ( (Gamma+1.0)*WL[0] ), \
                   A_R = 2.0 / ( (Gamma+1.0)*WR[0] );
            double B_L = Gamma_Ext1 * WL[2], \
                   B_R = Gamma_Ext1 * WR[2];
            double G_L = sqrt(A_L / (P0+B_L)), \
                   G_R = sqrt(A_R / (P0+B_R));
            double P_Star = ( G_L*WL[2] + G_R*WR[2] - (WR[1]-WL[1]) ) \
                   / ( G_L + G_R );
            double U_Star = 0.5*(WL[1]+WR[1]) + 0.5*( \
                   (P_Star-WR[2])*G_R - (P_Star-WL[2])*G_L );
            double Rho_StarL = WL[0] \
                   * (P_Star/WL[2]+Gamma_Ext1) \
                   / (Gamma_Ext1*P_Star/WL[2]+1.0), \
                   Rho_StarR = WR[0] \
                   * (P_Star/WR[2]+Gamma_Ext1) \
                   / (Gamma_Ext1*P_Star/WR[2]+1.0);
*/
            //
            double SS_StarL = sqrt(Gamma*P_Star/Rho_StarL), \
                   SS_StarR = sqrt(Gamma*P_Star/Rho_StarR);
            double Lambda1_L = WL[1] - SS_L, \
                   Lambda1_R = U_Star - SS_StarL;
            if(Lambda1_L < 0.0 && Lambda1_R > 0.0){
                EigVal[0] = Lambda1_L * (Lambda1_R-EigVal[0]) / (Lambda1_R-Lambda1_L);
            }
            double Lambda5_L = U_Star + SS_StarR, \
                   Lambda5_R = WR[1] + SS_R;
            if(Lambda5_L < 0.0 && Lambda5_R > 0.0){
                EigVal[2] = Lambda5_R * (EigVal[2]-Lambda5_L) / (Lambda5_R-Lambda5_L);
            }
            //
            calcFlux(WL, UL, FluxL);
            calcFlux(WR, UR, FluxR);
            //
            for(int i = 0; i < 3; i++){
                *(RoeFlux+i) = 0.5*(FluxL[i] + FluxR[i]);
                for(int j = 0; j < 3; j++){
                    *(RoeFlux+i) = *(RoeFlux+i) - 0.5*WaveA[j]*abs(EigVal[j])*EigVec[i*3+j];
                }
            } 
        }
        void calcTimeStep(){
            static double WaveSpeedMax = -1.0;
            for(int i = 0; i < SolN; i++){
                WaveSpeedMax = fmax( WaveSpeedMax, \
                    abs(*(SolW_ptr+i*3+1)) + sqrt(Gamma* *(SolW_ptr+i*3+2) / *(SolW_ptr+i*3)) );
            }
            TimeStep = CFL * PtsStep / WaveSpeedMax;
        }
        void marchTime(){
            static double GhostLSolW[3], GhostRSolW[3], GhostLSolU[3], GhostRSolU[3];
            for(int i = 0; i < 3; i++){
                GhostLSolW[i] = *(SolW_ptr+i);
                GhostRSolW[i] = *(SolW_ptr+(SolN-1)*3+i); 
            }
            W2U(GhostLSolW, GhostLSolU);
            W2U(GhostRSolW, GhostRSolU);
            // Getting RoeFlux at end points of flux points
            calcRoeFlux(GhostLSolW, SolW_ptr, GhostLSolU, SolU_ptr, Flu_ptr);
            calcRoeFlux(SolW_ptr+(SolN-1)*3, GhostRSolW, SolU_ptr+(SolN-1)*3, GhostRSolU, Flu_ptr+(FluN-1)*3);
            for(int i = 0; i < SolN-1; i++){
                calcRoeFlux(SolW_ptr+i*3, SolW_ptr+(i+1)*3, SolU_ptr+i*3, SolU_ptr+(i+1)*3, Flu_ptr+(i+1)*3);
            }
            // Updating all solution points 0...SolN-1 requires flux points 0...FluN-1
            // So solution points -1...SolN are required.
            for(int i = 0; i < SolN; i++){
                for(int j = 0; j < 3; j++){
                    *(SolU_ptr+i*3+j) = *(SolU_ptr+i*3+j) + TimeStep/PtsStep \
                        * ( *(Flu_ptr+i*3+j) - *(Flu_ptr+(i+1)*3+j) );
                }
            }
            U2W_EntireField();
        }
};
