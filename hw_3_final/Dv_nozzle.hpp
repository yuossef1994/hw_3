//
//  Dv_nozzle.hpp
//  hw_1
//
//  Created by Youssef Z on 2/11/23.
//

#ifndef Dv_nozzle_hpp
#define Dv_nozzle_hpp

#include <stdio.h>
#include <fstream>


class Dv_nozzle {
    
    
private:
    
    static int const _imax_f=61, _imax_c=60, _imin=0;
    int im_f, im_c;
    double A, B, C, D, E, F;  // just variables for calculations
    double k_2=0.2, k_4=0.02;
    double x_f[_imax_f];
    double avg_area;
    double avg_d_area;
    double x_c[_imax_c];
    double delta_t[_imax_c];
    double delta_t_gl=999;
    double CFL=0.01;
    double area_f[_imax_f];
    double _area, _d_area;
    double mu[_imax_c];
    double epsilon_2;
    double _R[3]={0,0,0};
    double _R_iter[3]={0,0,0};
    double eigen_v_avg;
    double _delta_x;
    double _P_0,  _T_0, mach_in, mach_out, T_in, u_in, P_in, rho_in, et_in, T_out, u_out, P_out, rho_out, et_out;
    double _M_0;
    double _mach[_imax_c];
    double T[_imax_c];
    double et[_imax_c];
    double ht[_imax_c];
    double _V[3][_imax_c];
    double _U[3][_imax_c];
    double _F[3][_imax_f];
    double _U_ghost_outflow[3];
    double _V_ghost_outflow[5];
    double _U_ghost_inflow[3];
    double _V_ghost_inflow[5];
    //upwind schemes variables
    
    double kappa=-1;
    double epsilon_upwind=1;
    double _F_C[3][_imax_f];
    double _F_P[3][_imax_f];
    
    double c_plus;
    double c_minus;
    double alpha_plus;
    double alpha_minus;
    double beta_R;
    double beta_L;
    double a[_imax_c];
    double a_L;
    double a_R;
    
    double D_plus;
    double D_minus;
    double P_plus;
    double P_minus;
    double V_L[3];
    double V_R[3];
    double T_L;
    double T_R;
    double r_plus;
    double r_minus;
    double epsi_plus[3][_imax_f];
    double epsi_minus[3][_imax_f];
    double ht_L;
    double ht_R;
    double mach_knight;
    double M_plus;
    double M_minus;
    double M_L;
    double M_R;
    
    
    
    double _D[3][_imax_f]; //JST dissipation scheme
    double gamma = 1.4, R_u=  8314, M_air= 28.96; //gamma air, universal gas const, molecular weight of air
    double R_air = R_u/M_air;
    double epsi_1;
    double rL2norm[3];
    double _rL2initial[3]={0,0,0};
    std::fstream L_vs_Iter; //creation of file to write the results
    
    std::fstream rho_vs_x; //creation of file to write the results
    std::fstream u_vs_x; //creation of file to write the results
    std::fstream p_vs_x; //creation of file to write the results
    
public:
    //constructor to open the file for write
    Dv_nozzle(){
        
     
        
        rho_vs_x.open("/Users/cringedaddy/CFD class/hw_3_back_up/hw_3/hw_3/hw_3/hw_3_final/hw_3_final/rho.csv", std::ios::trunc | std::ios::out);
        
        u_vs_x.open("/Users/cringedaddy/CFD class/hw_3_back_up/hw_3/hw_3/hw_3/hw_3_final/hw_3_final/u.csv", std::ios::trunc | std::ios::out);
        p_vs_x.open("/Users/cringedaddy/CFD class/hw_3_back_up/hw_3/hw_3/hw_3/hw_3_final/hw_3_final/p.csv", std::ios::trunc | std::ios::out);
        L_vs_Iter.open("/Users/cringedaddy/CFD class/hw_3_back_up/hw_3/hw_3/hw_3/hw_3_final/hw_3_final/L2.csv", std::ios::trunc | std::ios::out);
        L_vs_Iter<<"#iter#"<<"##L2_1#"<<"##L2_2#"<<"##L2_3#"<<std::endl;
        
    };
    // destructor to close the file
    ~Dv_nozzle()
    {
        
        rho_vs_x.close();
        u_vs_x.close();
        p_vs_x.close();
        L_vs_Iter.close();
        
    };

   
    void set_geometry(int,double,double);
    void initialization(double,double);
    void set_boundary_cond();
    void euler_explicit();
    double area(double);
    double d_area(double);
    void time_step();
    void print_res();
    void rL2initial();
    double L2norm(int);
};




#endif /* Dv_nozzle_hpp */
