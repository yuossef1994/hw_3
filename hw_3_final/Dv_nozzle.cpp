//
//  Dv_nozzle.cpp
//  hw_1
//
//  Created by Youssef Z on 2/11/23.
//

#include "Dv_nozzle.hpp"
#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;



void Dv_nozzle::set_geometry (int imax, double x_max, double x_min)
{
    
    //imax is number of cells
    //_imax is number of faces
    
     im_f = _imax_f-1;
     im_c = _imax_c-1;
    
    
    
    
    for( int i =_imin; i <= im_f; i++)
    {
      x_f[i]=   x_min + (float(i)/float(im_f))*(x_max- x_min) ;
        
    }
    
  
    for( int j =_imin; j <= im_c; j++)
    {
        x_c[j] =   (x_f[j] + x_f[j+1])/2;
        
    }
   
  
}
    void Dv_nozzle::initialization (double P_0, double T_0)
{
        
        
        
        _P_0= P_0;
        _T_0= T_0;
        
        
        
        
        
        
        
        
     
        // 0 is rho, 1 is velocity,  2 is pressure
        //values calculated at cell as cell average
        for (int i= _imin; i<= im_c; i++)
        {
            
            _M_0 = 0.85*x_c[i] + 1;
            epsi_1= 1.0 + ((gamma-1.0)*_M_0*_M_0)/2.0;
            
             T[i] = _T_0/epsi_1;
      
            
            _V[1][i] = _M_0 * sqrt(gamma*R_air*abs(T[i]));
            _V[2][i]= _P_0/pow(epsi_1,gamma/(gamma-1));
            _V[0][i]= _V[2][i]/(R_air*T[i]);
            et[i] = (R_air/(gamma-1))*T[i] + 0.5* _V[1][i]*_V[1][i];
            _U[0][i]= _V[2][i]/(R_air*T[i]);
            _U[1][i]= _V[0][i]*_V[1][i];
            _U[2][i]=_V[0][i]*et[i];
            
            
            
            _mach[i]=_M_0;
            
        }
     
       
        
        //calculate epsi
        
        for (int i=2;i<im_c;i++)
        {
            //for density
            r_plus= (_V[0][i+1]-_V[0][i])/(_V[0][i]-_V[0][i-1]);
            r_minus= (_V[0][i-1]-_V[0][i-2])/(_V[0][i]-_V[0][i-1]);
            epsi_plus[0][i]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus[0][i]= (r_minus+abs(r_minus))/(1+r_minus);
            //for velocity
            r_plus= (_V[1][i+1]-_V[1][i])/(_V[1][i]-_V[1][i-1]);
            r_minus= (_V[1][i-1]-_V[1][i-2])/(_V[1][i]-_V[1][i-1]);
            epsi_plus[1][i]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus[1][i]= (r_minus+abs(r_minus))/(1+r_minus);
            //for pressure
            r_plus= (_V[2][i+1]-_V[2][i])/(_V[2][i]-_V[2][i-1]);
            r_minus= (_V[2][i-1]-_V[2][i-2])/(_V[2][i]-_V[2][i-1]);
            epsi_plus[2][i]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus[2][i]= (r_minus+abs(r_minus))/(1+r_minus);
            
            
        }
        
        
        
        
        
        //calculate the flux
        for (int i=2;i<=im_c-1;i++)
        {
           
            
            //density
            V_L[0]= _V[0][i-1]+(epsilon_upwind/4)*((1-kappa)*epsi_plus[0][i-1]*(_V[0][i-1]-_V[0][i-2]) + (1+kappa)*epsi_plus[0][i]*(_V[0][i]-_V[0][i-1])        );
            
            V_R[0]= _V[0][i]-(epsilon_upwind/4)*((1-kappa)*epsi_minus[0][i+1]*(_V[0][i+1]-_V[0][i]) + (1+kappa)*epsi_minus[0][i]*(_V[0][i]-_V[0][i-1])        );
            //velocity
            V_L[1]= _V[1][i-1]+(epsilon_upwind/4)*((1-kappa)*epsi_plus[1][i-1]*(_V[1][i-1]-_V[1][i-2]) + (1+kappa)*epsi_plus[1][i]*(_V[1][i]-_V[1][i-1])        );
            
            V_R[1]= _V[1][i]-(epsilon_upwind/4)*((1-kappa)*epsi_minus[1][i+1]*(_V[1][i+1]-_V[1][i]) + (1+kappa)*epsi_minus[1][i]*(_V[1][i]-_V[1][i-1])        );
            //pressure
            V_L[2]= _V[2][i-1]+(epsilon_upwind/4)*((1-kappa)*epsi_plus[2][i-1]*(_V[2][i-1]-_V[2][i-2]) + (1+kappa)*epsi_plus[2][i]*(_V[2][i]-_V[2][i-1])        );
            
            V_R[2]= _V[2][i]-(epsilon_upwind/4)*((1-kappa)*epsi_minus[2][i+1]*(_V[2][i+1]-_V[2][i]) + (1+kappa)*epsi_minus[2][i]*(_V[2][i]-_V[2][i-1])        );
           
            
            a_L=sqrt(gamma*abs(V_L[2]/V_L[0]));
            a_R=sqrt(gamma*abs(V_R[2]/V_R[0]));
            
            M_L=V_L[1]/a_L;
            
            M_R=V_R[1]/a_R;
            mach_knight=(V_L[1]+V_R[1])/(a_L+a_R);
            M_plus= 0.25*(mach_knight+1)*(mach_knight+1);
            M_minus= -0.25*(mach_knight-1)*(mach_knight-1);
            
            
            
            
            
            beta_L=-max(0,1-int(abs(M_L)));
            
            beta_R=-max(0,1-int(abs(M_R)));
            
            
            alpha_plus=0.5*(1+copysign(1,M_L));
            alpha_minus=0.5*(1-copysign(1,M_R));
            c_plus = alpha_plus*(1+beta_L)*M_L -beta_L*M_plus;
            c_minus = alpha_minus*(1+beta_R)*M_R -beta_R*M_minus;
            
            
            ht_L=(gamma/(gamma-1))*(V_L[2]/V_L[0]) + V_L[1]*V_L[1]*0.5;
            ht_R=(gamma/(gamma-1))*(V_R[2]/V_R[0]) + V_R[1]*V_R[1]*0.5;
            
            
            P_plus = M_plus*(2-M_L);
            P_minus = M_minus*(-2-M_R);
            
            D_plus=alpha_plus*(1+beta_L)-beta_L*P_plus;
            
            D_minus=alpha_minus*(1+beta_R)-beta_R*P_minus;
            
            
            //calculate F_P
            _F_P[0][i]=0;
            
            _F_P[1][i]=D_plus*V_L[2]+D_minus*V_R[2];
            
            _F_P[2][i]=0;
            
            //calculate F_c
            
            _F_C[0][i]= V_L[0]*a_L*c_plus*1 + V_R[0]*a_R*c_minus*1;
            
            _F_C[1][i]= V_L[0]*a_L*c_plus*V_L[1] + V_R[0]*a_R*c_minus*V_R[1];
            
            _F_C[2][i]= V_L[0]*a_L*c_plus*ht_L + V_R[0]*a_R*c_minus*ht_R;
            
            _F[0][i]=_F_C[0][i]+_F_P[0][i];
            _F[1][i]=_F_C[1][i]+_F_P[1][i];
            _F[2][i]=_F_C[2][i]+_F_P[2][i];
            
            
        }
        
        
      
        
    }



void Dv_nozzle::set_boundary_cond()
{
    //inflow boundary conditions
    mach_in = 0.5*(3.0*_mach[0] - _mach[1] );
    epsi_1= 1 + ((gamma-1)*mach_in*mach_in)/2;
    T_in = _T_0/epsi_1;
    
    
    u_in =  mach_in * sqrt(gamma*R_air*abs(T_in));
    
    P_in= _P_0/pow(epsi_1,gamma/(gamma-1));
    
    rho_in= P_in/(R_air*T_in);
    
   // et_in = (R_air/(gamma-1))*T_in + 0.5*u_in*u_in;
   
    //extrapolate U and V for 1st ghost cell
    /*
    _U_ghost_inflow[0]=2*_U[0][_imin]-_U[0][_imin];
    _U_ghost_inflow[1]=2*_U[1][_imin]-_U[1][_imin];
    _U_ghost_inflow[2]=2*_U[2][_imin]-_U[2][_imin];
    */
    _V_ghost_inflow[0]=2*u_in-_V[0][_imin];
    _V_ghost_inflow[1]=2*rho_in-_V[1][_imin];
    
    _V_ghost_inflow[2]=2*P_in-_V[2][_imin];
    
    //_V_ghost_inflow[3]=2*et[_imin]-et[_imin];  //et
    _V_ghost_inflow[4]=2*T_in-T[_imin];    //T
  
   
    
    
    
    
    
    
    
    _F[0][0] = rho_in*u_in ;
    _F[1][0]= rho_in*u_in*u_in + P_in;
    _F[2][0]= (gamma/(gamma-1))*u_in*P_in + (rho_in*u_in*u_in*u_in)/2;
  

    
    
    //calculate epsi
    
   
    
        //for density
        r_plus= (_V[0][2]-_V[0][1])/(_V[0][1]-_V[0][0]);
        r_minus= (_V[0][0]-_V_ghost_inflow[0])/(_V[0][1]-_V[0][0]);
        epsi_plus[0][1]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus[0][1]= (r_minus+abs(r_minus))/(1+r_minus);
        //for velocity
        r_plus= (_V[1][2]-_V[1][1])/(_V[1][1]-_V[1][0]);
        r_minus= (_V[1][0]-_V_ghost_inflow[1])/(_V[1][1]-_V[1][0]);
        epsi_plus[1][1]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus[1][1]= (r_minus+abs(r_minus))/(1+r_minus);
        //for pressure
        r_plus= (_V[2][2]-_V[2][1])/(_V[2][1]-_V[2][0]);
        r_minus= (_V[2][0]-_V_ghost_inflow[2])/(_V[2][1]-_V[2][0]);
        epsi_plus[2][1]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus[2][1]= (r_minus+abs(r_minus))/(1+r_minus);
        
        
   
    
    
    
    
    
    
    
       
        
        //density
        V_L[0]= _V[0][0]+(epsilon_upwind/4)*((1-kappa)*epsi_plus[0][0]*(_V[0][0]-_V_ghost_inflow[0]) + (1+kappa)*epsi_plus[0][1]*(_V[0][1]-_V[0][0])        );
        
        V_R[0]= _V[0][1]-(epsilon_upwind/4)*((1-kappa)*epsi_minus[0][2]*(_V[0][2]-_V[0][1]) + (1+kappa)*epsi_minus[0][1]*(_V[0][1]-_V_ghost_inflow[0])        );
        //velocity
        V_L[1]= _V[1][0]+(epsilon_upwind/4)*((1-kappa)*epsi_plus[1][0]*(_V[1][0]-_V_ghost_inflow[1]) + (1+kappa)*epsi_plus[1][1]*(_V[1][1]-_V[1][0])        );
    
        V_R[1]= _V[1][1]-(epsilon_upwind/4)*((1-kappa)*epsi_minus[1][2]*(_V[1][2]-_V[1][1]) + (1+kappa)*epsi_minus[1][1]*(_V[1][1]-_V_ghost_inflow[1])        );
        //pressure
        V_L[2]= _V[2][0]+(epsilon_upwind/4)*((1-kappa)*epsi_plus[2][0]*(_V[2][0]-_V_ghost_inflow[2]) + (1+kappa)*epsi_plus[2][1]*(_V[2][1]-_V[2][0])        );

        V_R[2]= _V[2][1]-(epsilon_upwind/4)*((1-kappa)*epsi_minus[2][2]*(_V[2][2]-_V[2][1]) + (1+kappa)*epsi_minus[2][1]*(_V[2][1]-_V_ghost_inflow[2])        );
       
        
        a_L=sqrt(gamma*abs(V_L[2]/V_L[0]));
        a_R=sqrt(gamma*abs(V_R[2]/V_R[0]));
        
        M_L=V_L[1]/a_L;
        
        M_R=V_R[1]/a_R;
        mach_knight=(V_L[1]+V_R[1])/(a_L+a_R);
        M_plus= 0.25*(mach_knight+1)*(mach_knight+1);
        M_minus= -0.25*(mach_knight-1)*(mach_knight-1);
        
        
        
        
        
        beta_L=-max(0,1-int(abs(M_L)));
        
        beta_R=-max(0,1-int(abs(M_R)));
        
        
        alpha_plus=0.5*(1+copysign(1,M_L));
        alpha_minus=0.5*(1-copysign(1,M_R));
        c_plus = alpha_plus*(1+beta_L)*M_L -beta_L*M_plus;
        c_minus = alpha_minus*(1+beta_R)*M_R -beta_R*M_minus;
        
        
        ht_L=(gamma/(gamma-1))*(V_L[2]/V_L[0]) + V_L[1]*V_L[1]*0.5;
        ht_R=(gamma/(gamma-1))*(V_R[2]/V_R[0]) + V_R[1]*V_R[1]*0.5;
        
        
        P_plus = M_plus*(2-M_L);
        P_minus = M_minus*(-2-M_R);
        
        D_plus=alpha_plus*(1+beta_L)-beta_L*P_plus;
        
        D_minus=alpha_minus*(1+beta_R)-beta_R*P_minus;
        
        
        //calculate F_P
        _F_P[0][1]=0;
        
        _F_P[1][1]=D_plus*V_L[2]+D_minus*V_R[2];
        
        _F_P[2][1]=0;
        
        //calculate F_c
        
        _F_C[0][1]= V_L[0]*a_L*c_plus*1 + V_R[0]*a_R*c_minus*1;
        
        _F_C[1][1]= V_L[0]*a_L*c_plus*V_L[1] + V_R[0]*a_R*c_minus*V_R[1];
        
        _F_C[2][1]= V_L[0]*a_L*c_plus*ht_L + V_R[0]*a_R*c_minus*ht_R;
        
        _F[0][1]=_F_C[0][1]+_F_P[0][1];
        _F[1][1]=_F_C[1][1]+_F_P[1][1];
        _F[2][1]=_F_C[2][1]+_F_P[2][1];
        
        
  
    
    
    
    
    
    
    
    
    
    
    
    
    
    /*
    _F[0][1] = 0.5*(_F[0][0] + _F[0][2]);
    _F[1][1]= 0.5*(_F[1][0] + _F[1][2]);
    _F[2][1]= 0.5*(_F[2][0] + _F[2][2]);
   */
   
//******************************************
  // outflow boundary conditions
    
   
  
    P_out = 125000;
   // P_out = 0.5*(3*_V[2][im_c]-_V[2][im_c-1]);
  // std::cout<<"back pressure is "<<P_out<<std::endl;
  

    rho_out = 0.5*(3*_U[0][im_c]-_U[0][im_c-1]);
    
  
   
    u_out = (0.5*(3*_U[1][im_c]-_U[1][im_c-1]))/rho_out;
  
   _F[0][im_f] = rho_out*u_out ;
   _F[1][im_f]= rho_out*u_out*u_out + P_out;
   _F[2][im_f]= (gamma/(gamma-1))*u_out*P_out + (rho_out*u_out*u_out*u_out)/2;
  
    
    
   
    
   
  
    
    _V_ghost_outflow[0]=2*u_out-_V[0][im_c];
    _V_ghost_outflow[1]=2*rho_out-_V[1][im_c];
    
    _V_ghost_outflow[2]=2*P_out-_V[2][im_c];
    
    //_V_ghost_inflow[3]=2*et[_imin]-et[_imin];  //et
    _V_ghost_outflow[4]=2*T_out-T[im_c];    //T
    
    
    //calculate epsi
    
   
    //for density
    r_plus= (_V_ghost_outflow[0]-_V[0][im_c])/(_V[0][im_c]-_V[0][im_c-1]);
    r_minus= (_V[0][im_c-1]-_V[0][im_c-2])/(_V[0][im_c]-_V[0][im_c-1]);
    epsi_plus[0][im_c]= (r_plus+abs(r_plus))/(1+r_plus);
    epsi_minus[0][im_c]= (r_minus+abs(r_minus))/(1+r_minus);
    //for velocity
    r_plus= (_V_ghost_outflow[1]-_V[1][im_c])/(_V[1][im_c]-_V[1][im_c-1]);
    r_minus= (_V[1][im_c-1]-_V[1][im_c-2])/(_V[1][im_c]-_V[1][im_c-1]);
    epsi_plus[1][im_c]= (r_plus+abs(r_plus))/(1+r_plus);
    epsi_minus[1][im_c]= (r_minus+abs(r_minus))/(1+r_minus);
    //for pressure
    r_plus= (_V_ghost_outflow[2]-_V[2][im_c])/(_V[2][im_c]-_V[2][im_c-1]);
    r_minus= (_V[2][im_c-1]-_V[2][im_c-2])/(_V[2][im_c]-_V[2][im_c-1]);
    epsi_plus[2][im_c]= (r_plus+abs(r_plus))/(1+r_plus);
    epsi_minus[2][im_c]= (r_minus+abs(r_minus))/(1+r_minus);
        
        
   
    
    
    
    
    
    
    
       
        
    //density
    V_L[0]= _V[0][im_c-1]+(epsilon_upwind/4)*((1-kappa)*epsi_plus[0][im_c-1]*(_V[0][im_c-1]-_V[0][im_c-2]) + (1+kappa)*epsi_plus[0][im_c]*(_V[0][im_c]-_V[0][im_c-1])        );
    
    V_R[0]= _V[0][im_c]-(epsilon_upwind/4)*((1-kappa)*epsi_minus[0][im_c+1]*(_V_ghost_outflow[0]-_V[0][im_c]) + (1+kappa)*epsi_minus[0][im_c]*(_V[0][im_c]-_V[0][im_c-1])        );
    //velocity
    V_L[1]= _V[1][im_c-1]+(epsilon_upwind/4)*((1-kappa)*epsi_plus[1][im_c-1]*(_V[1][im_c-1]-_V[1][im_c-2]) + (1+kappa)*epsi_plus[1][im_c]*(_V[1][im_c]-_V[1][im_c-1])        );
    
    V_R[1]=_V[1][im_c]-(epsilon_upwind/4)*((1-kappa)*epsi_minus[1][im_c+1]*(_V_ghost_outflow[1]-_V[1][im_c]) + (1+kappa)*epsi_minus[1][im_c]*(_V[1][im_c]-_V[1][im_c-1])        );
    //pressure
    V_L[2]= _V[2][im_c-1]+(epsilon_upwind/4)*((1-kappa)*epsi_plus[2][im_c-1]*(_V[2][im_c-1]-_V[2][im_c-2]) + (1+kappa)*epsi_plus[2][im_c]*(_V[2][im_c]-_V[2][im_c-1])        );
    
    V_R[2]=_V[2][im_c]-(epsilon_upwind/4)*((1-kappa)*epsi_minus[2][im_c+1]*(_V_ghost_outflow[2]-_V[2][im_c]) + (1+kappa)*epsi_minus[2][im_c]*(_V[2][im_c]-_V[2][im_c-1])        );
       
        
        a_L=sqrt(gamma*abs(V_L[2]/V_L[0]));
        a_R=sqrt(gamma*abs(V_R[2]/V_R[0]));
        
        M_L=V_L[1]/a_L;
        
        M_R=V_R[1]/a_R;
        mach_knight=(V_L[1]+V_R[1])/(a_L+a_R);
        M_plus= 0.25*(mach_knight+1)*(mach_knight+1);
        M_minus= -0.25*(mach_knight-1)*(mach_knight-1);
        
        
        
        
        
        beta_L=-max(0,1-int(abs(M_L)));
        
        beta_R=-max(0,1-int(abs(M_R)));
        
        
        alpha_plus=0.5*(1+copysign(1,M_L));
        alpha_minus=0.5*(1-copysign(1,M_R));
        c_plus = alpha_plus*(1+beta_L)*M_L -beta_L*M_plus;
        c_minus = alpha_minus*(1+beta_R)*M_R -beta_R*M_minus;
        
        
        ht_L=(gamma/(gamma-1))*(V_L[2]/V_L[0]) + V_L[1]*V_L[1]*0.5;
        ht_R=(gamma/(gamma-1))*(V_R[2]/V_R[0]) + V_R[1]*V_R[1]*0.5;
        
        
        P_plus = M_plus*(2-M_L);
        P_minus = M_minus*(-2-M_R);
        
        D_plus=alpha_plus*(1+beta_L)-beta_L*P_plus;
        
        D_minus=alpha_minus*(1+beta_R)-beta_R*P_minus;
        
        
        //calculate F_P
        _F_P[0][im_f-1]=0;
        
        _F_P[1][im_f-1]=D_plus*V_L[2]+D_minus*V_R[2];
        
        _F_P[2][im_f-1]=0;
        
        //calculate F_c
        
        _F_C[0][im_f-1]= V_L[0]*a_L*c_plus*1 + V_R[0]*a_R*c_minus*1;
        
        _F_C[1][im_f-1]= V_L[0]*a_L*c_plus*V_L[1] + V_R[0]*a_R*c_minus*V_R[1];
        
        _F_C[2][im_f-1]= V_L[0]*a_L*c_plus*ht_L + V_R[0]*a_R*c_minus*ht_R;
        
        _F[0][im_f-1]=_F_C[0][im_f-1]+_F_P[0][im_f-1];
        _F[1][im_f-1]=_F_C[1][im_f-1]+_F_P[1][im_f-1];
        _F[2][im_f-1]=_F_C[2][im_f-1]+_F_P[2][im_f-1];
        
  
    
  
    
}

void Dv_nozzle::euler_explicit()
{
    
    for (int i=0;i<=im_c;i++)
    {
        _delta_x=x_f[i+1]-x_f[i];
        avg_area = ((area(x_f[i+1])+area(x_f[i]))/2);
        avg_d_area = ((d_area(x_f[i+1])+d_area(x_f[i]))/2);
        
    _U[0][i]= ( (-_F[0][i+1]   ) *area(x_f[i+1]) + (  _F[0][i]   )*area(x_f[i])
        )*(delta_t_gl/(avg_area *_delta_x)) + _U[0][i];
           
        
        
       
        
    _U[1][i]= (  (-_F[1][i+1]  )*area(x_f[i+1]) + (  _F[1][i]   )*area(x_f[i])
        )*(delta_t_gl/(avg_area* _delta_x)) + _U[1][i];
        
        
       
        
        
        
    _U[2][i]= (  (-_F[2][i+1]  )*area(x_f[i+1]) + (  _F[2][i]   )*area(x_f[i])
        )*(delta_t_gl/(avg_area *_delta_x)) + _U[2][i];
           
        
           
    }
    
    for (int i=0;i<=im_c;i++)
    {
       
        _V[0][i]=_U[0][i];
        _V[1][i]= _U[1][i]/_V[0][i];
        
        et[i] = _U[2][i]/ _V[0][i];
        T[i] = (et[i] - 0.5* _V[1][i]*_V[1][i])*((gamma-1)/R_air);
        
        
        
        _V[2][i]=_V[0][i]*R_air*T[i];
        
        ht[i]=(gamma/(gamma-1))*(_V[2][i]/_V[0][i]) + _V[1][i]*_V[1][i]*0.5;
        _mach[i]=abs(_V[1][i])/sqrt(gamma*R_air*abs(T[i]));
       // cout<< "velocity is    "<<_V[0][i]<<endl;
    }
    
    //calculate epsi
   
    for (int i=2;i<im_c;i++)
    {
        if (max(rL2norm[0],std::max(rL2norm[1],rL2norm[2])) <1e-4)break; // to freeze the limiter
        //for density
        r_plus= (_V[0][i+1]-_V[0][i])/(_V[0][i]-_V[0][i-1]);
        r_minus= (_V[0][i-1]-_V[0][i-2])/(_V[0][i]-_V[0][i-1]);
        epsi_plus[0][i]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus[0][i]= (r_minus+abs(r_minus))/(1+r_minus);
        //for velocity
        r_plus= (_V[1][i+1]-_V[1][i])/(_V[1][i]-_V[1][i-1]);
        r_minus= (_V[1][i-1]-_V[1][i-2])/(_V[1][i]-_V[1][i-1]);
        epsi_plus[1][i]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus[1][i]= (r_minus+abs(r_minus))/(1+r_minus);
        //for pressure
        r_plus= (_V[2][i+1]-_V[2][i])/(_V[2][i]-_V[2][i-1]);
        r_minus= (_V[2][i-1]-_V[2][i-2])/(_V[2][i]-_V[2][i-1]);
        epsi_plus[2][i]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus[2][i]= (r_minus+abs(r_minus))/(1+r_minus);
        
        
        
    }
    
    
    
    
    
    //calculate the flux
    for (int i=2;i<=im_c-1;i++)
    {
       
        
        //density
        V_L[0]= _V[0][i-1]+(epsilon_upwind/4)*((1-kappa)*epsi_plus[0][i-1]*(_V[0][i-1]-_V[0][i-2]) + (1+kappa)*epsi_plus[0][i]*(_V[0][i]-_V[0][i-1])        );
        
        V_R[0]= _V[0][i]-(epsilon_upwind/4)*((1-kappa)*epsi_minus[0][i+1]*(_V[0][i+1]-_V[0][i]) + (1+kappa)*epsi_minus[0][i]*(_V[0][i]-_V[0][i-1])        );
        //velocity
        V_L[1]= _V[1][i-1]+(epsilon_upwind/4)*((1-kappa)*epsi_plus[1][i-1]*(_V[1][i-1]-_V[1][i-2]) + (1+kappa)*epsi_plus[1][i]*(_V[1][i]-_V[1][i-1])        );
        
        V_R[1]= _V[1][i]-(epsilon_upwind/4)*((1-kappa)*epsi_minus[1][i+1]*(_V[1][i+1]-_V[1][i]) + (1+kappa)*epsi_minus[1][i]*(_V[1][i]-_V[1][i-1])        );
        //pressure
        V_L[2]= _V[2][i-1]+(epsilon_upwind/4)*((1-kappa)*epsi_plus[2][i-1]*(_V[2][i-1]-_V[2][i-2]) + (1+kappa)*epsi_plus[2][i]*(_V[2][i]-_V[2][i-1])        );
        
        V_R[2]= _V[2][i]-(epsilon_upwind/4)*((1-kappa)*epsi_minus[2][i+1]*(_V[2][i+1]-_V[2][i]) + (1+kappa)*epsi_minus[2][i]*(_V[2][i]-_V[2][i-1])        );
       
        
        a_L=sqrt(gamma*abs(V_L[2]/V_L[0]));
        a_R=sqrt(gamma*abs(V_R[2]/V_R[0]));
        
        M_L=V_L[1]/a_L;
        
        M_R=V_R[1]/a_R;
        mach_knight=(V_L[1]+V_R[1])/(a_L+a_R);
        M_plus= 0.25*(mach_knight+1)*(mach_knight+1);
        M_minus= -0.25*(mach_knight-1)*(mach_knight-1);
        
        
        
        
        
        beta_L=-max(0,1-int(abs(M_L)));
        
        beta_R=-max(0,1-int(abs(M_R)));
        
        
        alpha_plus=0.5*(1+copysign(1,M_L));
        alpha_minus=0.5*(1-copysign(1,M_R));
        c_plus = alpha_plus*(1+beta_L)*M_L -beta_L*M_plus;
        c_minus = alpha_minus*(1+beta_R)*M_R -beta_R*M_minus;
        
        
        ht_L=(gamma/(gamma-1))*(V_L[2]/V_L[0]) + V_L[1]*V_L[1]*0.5;
        ht_R=(gamma/(gamma-1))*(V_R[2]/V_R[0]) + V_R[1]*V_R[1]*0.5;
        
        
        P_plus = M_plus*(2-M_L);
        P_minus = M_minus*(-2-M_R);
        
        D_plus=alpha_plus*(1+beta_L)-beta_L*P_plus;
        
        D_minus=alpha_minus*(1+beta_R)-beta_R*P_minus;
        
        
        //calculate F_P
        _F_P[0][i]=0;
        
        _F_P[1][i]=D_plus*V_L[2]+D_minus*V_R[2];
        
        _F_P[2][i]=0;
        
        //calculate F_c
        
        _F_C[0][i]= V_L[0]*a_L*c_plus*1 + V_R[0]*a_R*c_minus*1;
        
        _F_C[1][i]= V_L[0]*a_L*c_plus*V_L[1] + V_R[0]*a_R*c_minus*V_R[1];
        
        _F_C[2][i]= V_L[0]*a_L*c_plus*ht_L + V_R[0]*a_R*c_minus*ht_R;
        
        _F[0][i]=_F_C[0][i]+_F_P[0][i];
        _F[1][i]=_F_C[1][i]+_F_P[1][i];
        _F[2][i]=_F_C[2][i]+_F_P[2][i];
        
        
    }
    
    
   
    
    
}
double Dv_nozzle::area(double x)
{
    _area = 0.2 + 0.4*(1 + sin(3.14159265359*(x-0.5)));
    
    
    return _area;
    
}
double Dv_nozzle::d_area(double x)
{
    _d_area = 0.4*3.14159265359*cos(3.14159265359*(x-0.5));
    
    
    return _d_area;
    
}
void Dv_nozzle::time_step()
{
   
    for (int i =0;i<=im_c;i++)
    {
        _delta_x=x_f[i+1]-x_f[i];
        delta_t[i] = (CFL *_delta_x)/(abs(_V[1][i])  + sqrt(gamma*R_air*abs(T[i]))) ;
        delta_t_gl=min(delta_t[i],delta_t_gl);
    }
    
  
    
}

void Dv_nozzle::print_res(){
    rho_vs_x<<"######### x ######   "<<"#########rho######   "<<std::endl;
    u_vs_x<<"######### x ######   "<<"#########ux######   "<<std::endl;
    p_vs_x<<"######### x ######   "<<"#########P######   "<<std::endl;
    for (int i= 0;i<=im_c;i++)
    {
        
        rho_vs_x<<std::setprecision(5)<<x_c[i]<<","<<std::setprecision(5)<< _V[0][i]<<std::endl;
        
        
        
        u_vs_x<<std::setprecision(5)<<x_c[i]<<","<<std::setprecision(14)<< _V[1][i]<<std::endl;
        
        
        
        p_vs_x<<std::setprecision(5)<<x_c[i]<<","<<std::setprecision(5)<< _V[2][i]<<std::endl;
        
        
        
    }
    
    
    
}
double Dv_nozzle::L2norm(int iter){
    
    
    rL2norm[0]=0;
    rL2norm[1]=0;
    rL2norm[2]=0;
    
    for (int i=0;i<=im_c;i++)
    {
        _delta_x=x_f[i+1]-x_f[i];
        _R_iter[0]= (_F[0][i+1]) *area(x_f[i+1]) -(_F[0][i]   )*area(x_f[i])  ;
        
        rL2norm[0] = rL2norm[0] +_R_iter[0]*_R_iter[0];
        
        
        _R_iter[1]= ( _F[1][i+1]   )*area(x_f[i+1]) -(_F[1][i]   )*area(x_f[i]) ;
        
        rL2norm[1] = rL2norm[1] +  _R_iter[1]* _R_iter[1];
        
        _R_iter[2]= ( _F[2][i+1]    )*area(x_f[i+1])- (_F[2][i]  )*area(x_f[i]);
        
        rL2norm[2] = rL2norm[2] +  _R_iter[2]* _R_iter[2];
    }
    double imax_c = _imax_c*1.0;
    
    rL2norm[0]=sqrt(rL2norm[0]/imax_c)/_rL2initial[0];
    rL2norm[1]=sqrt(rL2norm[1]/imax_c)/_rL2initial[1];
    rL2norm[2]=sqrt(rL2norm[2]/imax_c)/_rL2initial[2];
    
/*
    std::cout<< " R eqn1 "<< _R_iter[0]<<std::endl;
    std::cout<< " R eqn2 "<< _R_iter[1]<<std::endl;
    std::cout<< " R eqn3 "<< _R_iter[2]<<std::endl;
    std::cout<< " L eqn1 "<< rL2norm[0]<<std::endl;
    std::cout<< " L eqn2 "<< rL2norm[1]<<std::endl;
    std::cout<< " L eqn3 "<< rL2norm[2]<<std::endl;
    std::cout<< " L inti eqn1 "<< _rL2initial[0]<<std::endl;
    std::cout<< " L inti eqn2 "<< _rL2initial[1]<<std::endl;
    std::cout<< " L inti eqn3 "<< _rL2initial[2]<<std::endl;
    std::cout<< " u is  "<< _V[1][9]<<std::endl;
*/
    L_vs_Iter<<std::setprecision(5)<<iter<<","<<std::setprecision(5)<<rL2norm[0]<<","<<std::setprecision(5)<<rL2norm[1]<<","<<std::setprecision(5)<<rL2norm[2]<<std::endl;
    
    return std::max(rL2norm[0],std::max(rL2norm[1],rL2norm[2]));
    
    
}
void Dv_nozzle::rL2initial()


{
    _rL2initial[0]=0;
    _rL2initial[1]=0;
    _rL2initial[2]=0;
    for (int i=0;i<=im_c;i++)
    {
        _delta_x=x_f[i+1]-x_f[i];
        _R[0]= (_F[0][i+1]) *area(x_f[i+1]) - (_F[0][i])*area(x_f[i])  ;
        
        
        _rL2initial[0] = _rL2initial[0] + _R[0]*_R[0];
        
        
        _R[1]= ( _F[1][i+1]   )*area(x_f[i+1]) - (_F[1][i] )*area(x_f[i]) ;
        
        _rL2initial[1] = _rL2initial[1] + _R[1]*_R[1];
        
        _R[2]= ( _F[2][i+1]   )*area(x_f[i+1]) -(_F[2][i])*area(x_f[i]);
        
        _rL2initial[2] = _rL2initial[2] + _R[2]*_R[2];
    }
    double imax_c = _imax_c*1.0;
    _rL2initial[0] = sqrt(_rL2initial[0]/imax_c);
    _rL2initial[1] = sqrt(_rL2initial[1]/imax_c);
    _rL2initial[2] = sqrt(_rL2initial[2]/imax_c);
    
  
    
    
    
}
