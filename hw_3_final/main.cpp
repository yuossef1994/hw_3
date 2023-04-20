//
//  main.cpp
//  homework_3
//
//  Created by Youssef Z on 2/26/23.
//

#include <iostream>
#include "Dv_nozzle.hpp"
#include "nozzle_exact.hpp"
#include <math.h>

int main(int argc, const char * argv[]) {
    
    
    double A_x, T_0, P_0, Th_area, M_initial;
    T_0=600;
    P_0=300000;
    int max_faces=161;
    
    double x_min =-1;
    double x_max =1;
    int im_f_exact = max_faces-1;
   
    
    double x_f[max_faces];
    double x_c[max_faces-1];
   
   
   
   
   for( int i =0; i <= im_f_exact; i++)
   {
     x_f[i]=   x_min + (float(i)/float(im_f_exact))*(x_max- x_min) ;
       
   }
   
    for( int i =0; i <= im_f_exact-1; i++)
    {
        x_c[i]=   (x_f[i+1]+x_f[i])/2;
        
    }
   
    
    nozzle_exact exact;
    
    for (int i=0;i<=im_f_exact-1;i++)
    {
        A_x = 0.2 + 0.4*(1 + sin(3.14*(x_c[i]-0.5)));
        Th_area = 0.2 + 0.4*(1 + sin(-0.5*3.14));
        
        if ( x_c[i]<0)M_initial=0.5;
          
        else if (x_c[i]==0) M_initial=1;
        else M_initial = 5;
        
            
        exact.initialization(A_x, P_0, T_0, Th_area);
        exact.run_simulation(M_initial);
        exact.print_(x_c[i]);
        
    }
    
    
   
    
    
    
    
    
    int max_iter=100000000;
    double L2norm;
    Dv_nozzle nozzle;
    nozzle.set_geometry(200, 1, -1);
    nozzle.initialization(300000, 600);
    nozzle.set_boundary_cond();
    nozzle.time_step();
    nozzle.euler_explicit();
    nozzle.rL2initial();
    nozzle.L2norm(1);
    
    for(int i=1;i<=max_iter;i++)
    {
        
        nozzle.set_boundary_cond();
        nozzle.time_step();
        nozzle.euler_explicit();
       
        L2norm= nozzle.L2norm(i);
        
        
       if(i%10==0) std::cout<< " iter number  "<<i <<"  L2 norm is  " <<L2norm<<std::endl;
        
        
        if(L2norm < 1e-8){
            
            std::cout<< "solution converged"<<std::endl;
            break;
        }
       
    }
    
    
    nozzle.print_res();
    std::cout<<"h is = "<<abs((x_f[1]-x_f[0]))<<std::endl;
    
    return 0;
}
