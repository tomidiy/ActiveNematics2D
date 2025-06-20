#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdbool.h>
#include <time.h>


int const D = 2;  //dimensionality of space
int const ncomps = 2;   //number of elements of Q vector/tensor
int const Lx = 64;
int const Ly = 64; 



double K = 256.0 * 256.0;
double nematic_coherence_length = 0.5;
double gamma = 10.0 * 256;
double lambda = 0.1;
double Re = 1e-1;
double active_length_scale = 2.0;
int jit_loops = 100;
double dt = 2e-4;
double dt_pseudo = 2e-4;    // dt_pseudo = dt
double p_target_rel_change = 1e-4;   // p_target_rel_change = 1e-4 * dt_pseudo / dt
int max_p_iters = 0;
int udiff_thresh = -1;
int max_steps = (int) 1e7; //1e9;
double max_t = 1e7; //1e9;
double zeta = 100.0;//256.0 * 256.0 / (2.0 * 2.0);   // zeta = K / (active_length_scale * active_length_scale)
double nu = 809.5430810031052 * 1.0;   // nu = sqrt(K / Re);
double C = 256.0 * 256.0 / (0.5 * 0.5);   // C = K / (nematic_coherence_length * nematic_coherence_length);
double A = - 256.0 * 256.0/(0.5 * 0.5);   // A = -C
double S0 = 0.5 * 1.0;           //S0 = sqrt(-A / (4 * C))
//int seed = 12345;
int save_every_n_steps = 100;//10 * 100;  // save_every_n_steps = 10 * jit_loops


double coeff = -1;

char resultspath[50];
char imgpath[50];

struct consts_dict 
{
    double dt;
    double dt_pseudo;
    double nu;
    double A;
    double C;
    double K;
    double lambda;
    double zeta;
    double gamma;
    double p_target_rel_change;
    int max_p_iters;
};


struct consts_dict consts;


void make_directory_if_needed(char* pathname) 
{
    if (access(pathname, F_OK) == -1) 
    {
        if (mkdir(pathname, 0777) == -1) 
        {
            perror("mkdir error");
            exit(EXIT_FAILURE);
        }
    }
}


                                    
void Laplacian(double *arr, double *out, double coeff, int Lx, int Ly)
{
    //Laplacian of a scalar array. This function finds the laplacian of a scalar for each lattice site.
    int xup, xdn, yup, ydn;

    coeff /= 6;
    for(int x=0; x<Lx; x++)
    {
        xup = (x + 1) % Lx;       
        xdn = x < 1 ? (Lx - 1) : (x - 1);   //go to the last index if x is the first index else xdn is the immediate past index
        for(int y=0; y<Ly; y++)
        {
            yup = (y + 1) % Ly;     
            ydn = y < 1 ? (Ly - 1) : (y - 1);
            out[x*Ly + y] = coeff * (-20*arr[x*Ly + y] + 4 * (
                    arr[xup*Ly + y] + arr[x*Ly + ydn] +  arr[x*Ly + yup]
                        +  arr[xdn*Ly + y]) +  arr[xup*Ly + ydn] +  arr[xup*Ly + yup]
                           +  arr[xdn*Ly + ydn] +  arr[xdn*Ly + yup]);
        
           
        }
    }
}


void Laplacian_vector(double *arr, double *out, double coeff_, int Lx, int Ly, int ncomps)
{
    //This function finds the Laplacian of a vector array for each lattice site 
    int xup, xdn, yup, ydn;

    coeff_ /= 6;
    for(int x=0; x<Lx; x++)
    {
        xup = (x + 1) % Lx;
        xdn = x < 1 ? (Lx - 1) : (x - 1);
        for(int y=0; y<Ly; y++)
        {
            yup = (y + 1) % Ly;
            ydn = y < 1 ? (Ly - 1) : (y - 1);
            for(int i=0; i<ncomps; i++)
            {
                out[(x*Ly + y)*ncomps + i] = coeff_ * (-20*arr[(x*Ly + y)*ncomps + i] + 4 * (
                    arr[(xup*Ly + y)*ncomps + i] + arr[(x*Ly + ydn)*ncomps + i] +  arr[(x*Ly + yup)*ncomps + i]
                        +  arr[(xdn*Ly + y)*ncomps + i]) +  arr[(xup*Ly + ydn)*ncomps + i] +  arr[(xup*Ly + yup)*ncomps + i]
                           +  arr[(xdn*Ly + ydn)*ncomps + i] +  arr[(xdn*Ly + yup)*ncomps + i]);
            }
           
        }
    }

}


void div_vector(double *arr, double *out, int Lx, int Ly)
{
    //Calculate the divergence of a vector field for each lattice site.
    int xup, xdn, yup, ydn;

    for(int x=0; x<Lx; x++)
    {
        xup = (x + 1) % Lx;
        xdn = x < 1 ? (Lx - 1) : (x - 1);
        for(int y=0; y<Ly; y++)
        {
            yup = (y + 1) % Ly;
            ydn = y < 1 ? (Ly - 1) : (y - 1);
            out[x*Ly + y] = 0.5 * (arr[(xup*Ly + y)*D + 0] - arr[(xdn*Ly + y)*D + 0] 
                + arr[(x*Ly + yup)*D + 1] - arr[(x*Ly + ydn)*D + 1]);
        }
    }
}


void apply_Q_boundary_conditions(double *Q){}

void apply_u_boundary_conditions(double *u){}

void apply_p_boundary_conditions(double *p){}


void initialize_Q_from_theta(double *Q, double *theta_initial, double S0, int Lx, int Ly, int ncomps) 
{ 
    //The Q tensor, used as a vector, for each lattic site is calculated from the orientation of the lattice site 
    double nx_initial, ny_initial;
    for (int x=0; x<Lx; x++) 
    {
        for (int y=0; y<Ly; y++) 
        {
            nx_initial = cos(theta_initial[x*Ly + y]);
            ny_initial = sin(theta_initial[x*Ly + y]);
            Q[(x*Ly + y)*ncomps + 0] = 2*S0*(nx_initial*nx_initial - 0.5);
            Q[(x*Ly + y)*ncomps + 1] = 2*S0*(nx_initial*ny_initial);
        }
    }
}



void H_S_from_Q(double *u, double *Q, double *H, double *S, double A, double C, double K, double lambda, int D)
{ 
    //Calculations for molecular field and co-rotation tensors for each lattice site. """ 
    int xup, xdn, yup, ydn;
    double dxux, dxuy, dyux, omegaxy, trQsq, lambdaS;

    Laplacian_vector(Q, H, K, Lx, Ly, ncomps); // H = K * ∇²Q
    
    //for(int j=0; j<Lx*Ly*D; j++) {printf("H[%d]: %f\n ", j, H[j]);}
    
    for(int x=0; x<Lx; x++)
    {
        xup = (x + 1) % Lx;
        xdn = x < 1 ? (Lx - 1) : (x - 1);
        for(int y=0; y<Ly; y++)
        {
            yup = (y + 1) % Ly;
            ydn = y < 1 ? (Ly - 1) : (y - 1); 

            dxux = 0.5 * (u[(xup*Ly + y)*D + 0] - u[(xdn*Ly + y)*D + 0]);
            dxuy = 0.5 * (u[(xup*Ly + y)*D + 1] - u[(xdn*Ly + y)*D + 1]);
            dyux = 0.5 * (u[(x*Ly + yup)*D + 0] - u[(x*Ly + ydn)*D + 0]);
            omegaxy = 0.5 * (dxuy - dyux);
            trQsq = 2 * (Q[(x*Ly + y)*ncomps + 0]*Q[(x*Ly + y)*ncomps + 0] + Q[(x*Ly + y)*ncomps + 1]*Q[(x*Ly + y)*ncomps + 1]);
            lambdaS = lambda * sqrt(2*trQsq);
            H[(x*Ly + y)*ncomps + 0] -= (A + C*trQsq) * Q[(x*Ly + y)*ncomps + 0];
            H[(x*Ly + y)*ncomps + 1] -= (A + C*trQsq) * Q[(x*Ly + y)*ncomps + 1];
            
            // using E_xx = du/dx:
            S[(x*Ly + y)*ncomps + 0] = lambdaS * dxux - 2*omegaxy*Q[(x*Ly + y)*ncomps + 1];
            
            // using E_xy = (du/dy + dv/dx)/2:
            S[(x*Ly + y)*ncomps + 1] = lambdaS * 0.5*(dxuy + dyux) + 2*omegaxy*Q[(x*Ly + y)*ncomps + 0];
        }
    }
}



void calculate_Pi(double *Pi_S, double *Pi_A, double *H, double *Q, double lambda, double zeta, int Lx, int Ly, int ncomps)
{
    //Calculation of stress tensor (elastic + active contributions) for each lattice site. The stress tensor is saved as a vector
    
    //symmetric traceless component (two elements)
    for(int i=0; i<Lx*Ly*ncomps; i++){Pi_S[i] = -lambda * H[i] - zeta * Q[i];}
    
    //antisymmetric component (one element)
    for(int x=0; x<Lx; x++)
    { 
        for(int y=0; y<Ly; y++)
        {
            Pi_A[x*Ly + y] = 2 * (Q[(x*Ly + y)*ncomps+0] * H[(x*Ly + y)*ncomps+ 1] - H[(x*Ly + y)*ncomps + 
            0] * Q[(x*Ly + y)*ncomps + 1]);
        }
    }
}



void calculate_pressure_terms(double *u, double *Pi_S, double *pressure_poisson_RHS, int Lx, int Ly)
{ 
    //Calculates the right hand side of the Poisson equation for each lattice site
    int xup, xdn, yup, ydn;
    double dudx, dvdy;

    for(int x=0; x<Lx; x++)
    {
         xup = (x + 1) % Lx;
         xdn = x < 1 ? (Lx - 1) : (x - 1);

         for(int y=0; y<Ly; y++)
         {
            yup = (y + 1) % Ly;
            ydn = y < 1 ? (Ly - 1) : (y - 1);
            dudx = 0.5 * (u[(xup*Ly + y)*D + 0] - u[(xdn*Ly + y)*D + 0]);
            dvdy = 0.5 * (u[(x*Ly + yup)*D + 1] - u[(x*Ly + ydn)*D + 1]);

            pressure_poisson_RHS[x*Ly + y] += (
                //∇•F = ∂ᵢ∂ⱼΠᵢⱼ = (∂x^2 - ∂y^2) Π_{xx} + 2 ∂x ∂y Π_{xy}
                (Pi_S[(xup*Ly + y)*D + 0] + Pi_S[(xdn*Ly + y)*D + 0] - Pi_S[(x*Ly + yup)*D + 0] - Pi_S[(x*Ly + ydn)*D + 0])
                + 0.5 * (Pi_S[(xup*Ly + yup)*D + 1] - Pi_S[(xup*Ly + ydn)*D + 1] - Pi_S[(xdn*Ly + yup)*D + 1] + 
                Pi_S[(xdn*Ly + ydn)*D + 1])
                // - d_i u_j d_j u_i
                - (
                    // (d_x u_x)^2 + (d_y u_y)^2:
                    dudx*dudx + dvdy*dvdy
                    // + 2 d_y u_x d_x u_y:
                    + 0.5 * (u[(x*Ly + yup)*D + 0]-u[(x*Ly + ydn)*D + 0]) * (u[(xup*Ly + y)*D + 1]-u[(xdn*Ly + y)*D + 1])));
        }
    }

}



double relax_pressure_inner_loop(double *p, double *p_aux, double *pressure_poisson_RHS, int Lx, int Ly)
{  
    // Evolves the pressure field towards satisfying the pressure Poisson equation. This is done for each lattice site
    int xup, xdn, yup, ydn;
    
    double sum_abs_p = 0.0;
    for(int x=0; x<Lx; x++)
    {
        xup = (x + 1) % Lx;
        xdn = x < 1 ? (Lx - 1) : (x - 1);

        for(int y=0; y<Ly; y++)
        {
            yup = (y + 1) % Ly;
            ydn = y < 1 ? (Ly - 1) : (y - 1);

            //coefficient depend on choice on Laplacian stencil
            //solved for the central term
            p[x*Ly + y] = 0.05 * (-6 * pressure_poisson_RHS[x*Ly + y] + 4 * (p_aux[xup*Ly + y] + p_aux[x*Ly + yup]
                + p_aux[x*Ly + ydn] + p_aux[xdn*Ly + y]) + p_aux[xup*Ly + yup] + p_aux[xup*Ly + ydn] 
                    + p_aux[xdn*Ly + yup] + p_aux[xdn*Ly + ydn]);

            sum_abs_p += fabs(p[x*Ly + y]);
        }
    }
    return sum_abs_p;
}



int relax_pressure(double *u, double *p, double *Pi_S, double *p_aux, double *pressure_poisson_RHS, 
                                double dt_pseudo, double target_rel_change, int max_p_iters, int Lx, int Ly)
{
    //Find pressure field for each lattic site that maintains incompressibility of flow field
    double sum_abs_p;

    div_vector(u, pressure_poisson_RHS, Lx, Ly);   //RHS = ∇•u
    for(int i=0; i<Lx*Ly; i++)
    {pressure_poisson_RHS[i] *= 1/dt_pseudo;}      //RHS = (∇•u)/∆t
    // RHS += -d_i u_j d_j u_i + ∇•F :
    calculate_pressure_terms(u, Pi_S, pressure_poisson_RHS, Lx, Ly);
    int p_iters = 0;

    bool pressure_has_not_relaxed = 1;

    while(pressure_has_not_relaxed)
    {
        for(int i=0; i<Lx*Ly; i++)
        {
            p_aux[i] = p[i];
        }
        // the pseudotimestep update of the pressure field:
        sum_abs_p = relax_pressure_inner_loop(p, p_aux, pressure_poisson_RHS, Lx, Ly);
        

        apply_p_boundary_conditions(p);
        p_iters += 1;

        //end pressure relaxation if number of pseudotimesteps exceeds max
        if(p_iters >= max_p_iters && max_p_iters > 0){printf("p_iters triggered\n");break;}
        
        //end pressure relaxation if relative change per step falls below 
        //threshold:
        if(sum_abs_p != 0)
        {
            double rel_change = 0.0;
            for(int i=0; i<Lx*Ly; i++)
            {
                rel_change += fabs(p_aux[i] - p[i]);
            }
            rel_change  /= sum_abs_p;

            if(rel_change <= target_rel_change){printf("rel_change triggered\n"); break;}
        }
    }
    return p_iters;
}




void upwind_advective_term(double *u, double *arr, double *out, double coeff, int Lx, int Ly, int ncomps)
{
    //calculate second-order upwind advective derivative -(u•∇)[arr] and add this to 'out'. This is done for each lattice site.
    int xup, xdn, xupup, xdndn, ydn, ydndn, yup, yupup;
    double tmp;

    coeff *= 0.5;

    for(int x=0; x<Lx; x++)
    {   
        xup = (x + 1) % Lx;   //periodic boundary in x
        xdn = x < 1 ? (Lx - 1) : (x - 1);
        xupup = (x + 2) % Lx;
        xdndn = x < 2 ? (Lx - 2 + x) : (x - 2);   //go to the (second last - i) index if x is i position among first 2 indices else xdn is the immediate past 2 index

        for(int y=0; y<Ly; y++)
        {
            // -ux dx A_k
            tmp = coeff * u[(x*Ly + y)*D + 0];
            if(u[(x*Ly + y)*D + 0] > 0)
            {
               for(int q=0; q<ncomps; q++)
               {
                    out[(x*Ly + y)*ncomps + q] += tmp * (3*arr[(x*Ly + y)*ncomps + q] - 4*arr[(xdn*Ly + y)*ncomps + q] 
                        + arr[(xdndn*Ly + y)*ncomps + q]);
               }
            }

            else
            {
                for(int q=0; q<ncomps; q++)
                {
                     out[(x*Ly + y)*ncomps + q] += tmp * (-3*arr[(x*Ly + y)*ncomps + q] + 4*arr[(xup*Ly + y)*ncomps + q]
                         - arr[(xupup*Ly + y)*ncomps + q]);
                }
                   
            } 
            
            // -uy dy A_k
            tmp = coeff * u[(x*Ly + y)*D + 1];
            if(u[(x*Ly + y)*D + 1] > 0)
            {
                ydn = y < 1 ? (Ly - 1) : (y - 1);
                ydndn =  y < 2 ? (Ly - 2 + y) : (y - 2);
                for(int q=0; q<ncomps; q++)
                {
                    out[(x*Ly + y)*ncomps + q] += tmp * (3*arr[(x*Ly + y)*ncomps + q] - 4*arr[(x*Ly + ydn)*ncomps + q] 
                        + arr[(x*Ly + ydndn)*ncomps + q]);
                }
            }
            else
            {
                yup = (y + 1) % Ly;    //periodic boundary in y
                yupup = (y + 2) % Ly;
                for(int q=0; q<ncomps; q++)
                {
                    out[(x*Ly + y)*ncomps + q] += tmp * (-3*arr[(x*Ly + y)*ncomps + q] + 4*arr[(x*Ly + yup)*ncomps + q] 
                        - arr[(x*Ly + yupup)*ncomps+ q]);
                }

            }
        }
    }

}


void u_update_p_Pi_terms(double *dudt, double *u, double *p, double *Pi_S, double *Pi_A, int Lx, int Ly)
{
    //For a lattice site, add velocity update terms from pressure and active+elastic forces
    int xup, xdn, yup, ydn;

    for(int x=0; x<Lx; x++)
    {
        xup = (x + 1) % Lx;
        xdn = x < 1 ? (Lx - 1) : (x - 1);
        for(int y=0; y<Ly; y++)
        {
            yup = (y + 1) % Ly;
            ydn = y < 1 ? (Ly - 1) : (y - 1);
            //pressure + force density from elastic and active stresses
            dudt[(x*Ly + y)*D + 0] += 0.5 * (-(p[xup*Ly + y] - p[xdn*Ly + y]) //[-grad(p)]_x
                //F_x = dx Πxx + dy Πxy: 
                + (Pi_S[(xup*Ly + y)*D + 0] - Pi_S[(xdn*Ly + y)*D + 0])   //dx Πxx
                + ((Pi_S[(x*Ly + yup)*D + 1] + Pi_A[x*Ly + yup])-(Pi_S[(x*Ly + ydn)*D + 1] + Pi_A[x*Ly + ydn]))); //dy Πxy

            dudt[(x*Ly + y)*D + 1] += 0.5 * (-(p[x*Ly + yup] - p[x*Ly + ydn]) //[-grad(p)]_y
                //F_x = dx Πyx + dy Πyy = dx Πyx - dy Πxx: 
                - (Pi_S[(x*Ly + yup)*D + 0] - Pi_S[(x*Ly + ydn)*D + 0])   //-dy Πxx
                + ((Pi_S[(xup*Ly + y)*D + 1] - Pi_A[xup*Ly + y])-(Pi_S[(xdn*Ly + y)*D + 1] - Pi_A[xdn*Ly + y])));
        }
    }

}



void get_Q_update(double *dQ, double *Q, double *H, double *S, double *u, double gamma, int Lx, int Ly, int ncomps)
{
    //Calculate update to Q tensor for each lattice site
    for(int i=0; i<Lx*Ly*ncomps; i++) {dQ[i] = (1/gamma)*H[i] + S[i];}
    upwind_advective_term(u, Q, dQ, coeff=-1, Lx, Ly, ncomps);   //subtract (u•∇)Q
}


void get_u_update(double *dudt, double *u, double *p, double *Pi_S, double *Pi_A, double nu, int Lx, int Ly, int D)
{
    //update flow field for the lattice sites after solving for pressure field
    Laplacian_vector(u, dudt, nu, Lx, Ly, D); //viscous term
    upwind_advective_term(u, u, dudt, coeff=-1, Lx, Ly, D); //convective term
    //pressure and stress tensor (active + elastic) contributions 
    u_update_p_Pi_terms(dudt, u, p, Pi_S, Pi_A, Lx, Ly);  
   
}



int update_step_inner(double *u, double *dudt , double *Q, double *dQdt, double *p, double  *p_aux, 
                            double *pressure_poisson_RHS, double *S, double *H, double *Pi_S, double *Pi_A,
                                    struct consts_dict consts, int Lx, int Ly, int D, int ncomps)
{
    // Euler update for the lattice sites 
    
    double dt = consts.dt;
    double dt_pseudo = consts.dt_pseudo;
    double nu = consts.nu;
    double A = consts.A;
    double C = consts.C;
    double K = consts.K;
    double lambda = consts.lambda;
    double zeta = consts.zeta;
    double gamma = consts.gamma;
    double p_target_rel_change = consts.p_target_rel_change;
    double max_p_iters = consts.max_p_iters;

     // update H and S  
    H_S_from_Q(u, Q, H, S, A, C, K, lambda, D);  
    apply_Q_boundary_conditions(H);

    // update Pi
    calculate_Pi(Pi_S, Pi_A, H, Q, lambda, zeta, Lx, Ly, ncomps);  

    // relax pressure to ensure incompressibility 
    int p_iters = relax_pressure(u, p, Pi_S, p_aux, pressure_poisson_RHS, 
                        dt_pseudo, p_target_rel_change, max_p_iters, Lx, Ly);

    //calculate dQdt 
    get_Q_update(dQdt, Q, H, S, u, gamma, Lx, Ly, ncomps);
    apply_Q_boundary_conditions(dQdt);

    //calculate dudt 
    get_u_update(dudt, u, p, Pi_S, Pi_A, nu, Lx, Ly, D);
    apply_u_boundary_conditions(dudt);

    //update Q, u 
    for (int i=0; i<Lx*Ly*ncomps; i++) {Q[i] += dt * dQdt[i];}
    for (int j=0; j<Lx*Ly*D; j++) {u[j] += dt * dudt[j];}

    return p_iters;
}

double update_step(double *u, double *dudt , double *Q, double *dQdt, double *p, double *p_aux, 
                            double *pressure_poisson_RHS, double *S, double *H, double *Pi_S, double *Pi_A,
                                    struct consts_dict consts, int *stepcount, double *t, double *udiff, int n_steps,
                                         int Lx, int Ly, int D, int ncomps)
{
    //Loop over Euler updates to velocity u, pressure p, and the Q-tensor
    double dt = consts.dt;
     

    int p_iters = 0;
    for (int i = 0; i < n_steps; i++) 
    {
        p_iters += update_step_inner(u, dudt , Q, dQdt, p, p_aux, pressure_poisson_RHS, 
                                                S, H, Pi_S, Pi_A, consts, Lx, Ly, D, ncomps);
        *stepcount += 1;
        *t += dt;
    }
    

    double du_tot_sum = 0;
    double udiff_denom = 0;
    for (int i = 0; i < Lx*Ly*D; i++) 
    {
        du_tot_sum += pow(dudt[i], 2);
        udiff_denom += pow(u[i], 2);
    }
    
    //printf("du_tot_sum: %f  udiff_denom: %f\n", du_tot_sum, udiff_denom);

    if (udiff_denom != 0){*udiff = sqrt(du_tot_sum) / sqrt(udiff_denom);} 
    else {*udiff = 0.0;}

    double avg_p_iters = (double)p_iters / n_steps;
    return avg_p_iters;
}



void get_saddle_splay(double *Q, double *out)
{
    //calculate saddle-splay energy density (for plotting defects)
    for (int x = 0; x < Lx; x++)
    {
        int xup = (x + 1) % Lx;
        int xdn = x < 1 ? (Lx - 1) : (x - 1);
        for (int y = 0; y < Ly; y++)
        {
            int yup = (y + 1) % Ly;
            int ydn = y < 1 ? (Ly - 1) : (y - 1); 
            double twice_dxQxx = Q[(xup*Ly + y)*ncomps + 0] - Q[(xdn*Ly + y)*ncomps + 0];
            double twice_dxQxy = Q[(xup*Ly + y)*ncomps + 1] - Q[(xdn*Ly + y)*ncomps + 1];
            double twice_dyQxx = Q[(x*Ly + yup)*ncomps + 0] - Q[(x*Ly + ydn)*ncomps + 0];
            double twice_dyQxy = Q[(x*Ly +yup)*ncomps + 1] - Q[(x*Ly + ydn)*ncomps + 1];
            out[x*Ly + y] = twice_dxQxy*twice_dyQxx - twice_dxQxx*twice_dyQxy;
        }
    }
}


void E_omega_from_u(double *u, double *E, double *omega)
{
    //Compute the strain rate tensor E and vorticity ω from the flow field u  for each lattice site. (Only used for plots)
    
    for(int x = 0; x < Lx; x++)
    {
        int xup = (x + 1) % Lx;
        int xdn = x < 1 ? (Lx - 1) : (x - 1);
        for (int y = 0; y < Ly; y++)
        {
            int yup = (y + 1) % Ly;
            int ydn = y < 1 ? (Ly - 1) : (y - 1);
            double dxuy = u[(xup*Ly + y)*D + 1] - u[(xdn*Ly + y)*D + 1];
            double dyux = u[(x*Ly + yup)*D + 0] - u[(x*Ly + ydn)*D + 0];
            E[(x*Ly + y)*D + 0] = 0.5 * (u[(xup*Ly + y)*D + 0] - u[(xdn*Ly + y)*D + 0]);
            E[(x*Ly + y)*D + 1] = 0.25 * (dxuy + dyux);
            omega[x*Ly + y] = 0.25 * (dxuy - dyux);
        }
    }
}



void n_dor_from_Q(double *Q, double *nx, double *ny, double *dor) 
{
    //For each lattice site, extract nematic director and degree of order from Q-tensor
    double trQsq, a, b, denom;

    for(int x=0; x<Lx; x++)
    {
        for(int y=0; y<Ly; y++) 
        {
            trQsq = 2 * (pow(Q[(x*Ly + y)*ncomps + 0],2) + pow(Q[(x*Ly + y)*ncomps+ 1], 2));
            dor[x*Ly + y] = 0.5 * sqrt(trQsq);

            a = Q[(x*Ly + y)*ncomps + 0];  //Qxx value
            b = Q[(x*Ly + y)*ncomps + 1];   //Qxy value
            denom = sqrt(b * b + pow(a + sqrt(a * a + b * b), 2));  //normalization for n

            if (denom == 0) {
                nx[x*Ly + y] = 0;
                ny[x*Ly + y] = 0;
                continue;
            }

            nx[x*Ly + y] = (a + sqrt(a * a + b * b)) / denom;
            ny[x*Ly + y] = b / denom; 
        }
    }
}


void run_active_nematic_sim(double *u, double *Q, double *p, char *runlabel, struct consts_dict consts, int Lx, int Ly)
{   
    char run_results_path[50], Q_results_path[50], u_results_path[50], run_img_path[50], p_mean_results_path[50];
    char ss_results_path[50], omega_results_path[50], E_results_path[50], n_dor_results_path[50], p_results_path[50];

    sprintf(run_results_path, "%s%s/",resultspath, runlabel);
    make_directory_if_needed(run_results_path);

    sprintf(Q_results_path, "%s%s/",run_results_path, "Q");
    make_directory_if_needed(Q_results_path);

    sprintf(u_results_path, "%s%s/",run_results_path, "u");
    make_directory_if_needed(u_results_path);

    sprintf(run_img_path, "%s%s/", imgpath, "/");
    make_directory_if_needed(run_img_path);
    
    sprintf(ss_results_path, "%s%s/",run_results_path, "ss");
    make_directory_if_needed(ss_results_path);

    sprintf(omega_results_path, "%s%s/",run_results_path, "omega");
    make_directory_if_needed(omega_results_path);

    sprintf(E_results_path, "%s%s/",run_results_path, "E");
    make_directory_if_needed(E_results_path);

    sprintf(n_dor_results_path, "%s%s/",run_results_path, "n_dor");
    make_directory_if_needed(n_dor_results_path);

    sprintf(p_results_path, "%s%s/",run_results_path, "p");
    make_directory_if_needed(p_results_path);

    sprintf(p_mean_results_path, "%s%s/",run_results_path, "p_mean");
    make_directory_if_needed(p_mean_results_path);

    // save parameters to text file
    FILE *f;
    f = fopen(strcat(run_results_path, strcat(runlabel, "_consts.txt")), "w");
    fprintf(f, "Lx=%d Ly=%d", Lx, Ly);

    // save other constants to file
    fprintf(f, " %s=%.6f", "dt", consts.dt);
    fprintf(f, " %s=%.6f", "dt_pseudo", consts.dt_pseudo);
    fprintf(f, " %s=%.6f", "nu", consts.nu);
    fprintf(f, " %s=%.6f", "A", consts.A);
    fprintf(f, " %s=%.6f", "C", consts.C);
    fprintf(f, " %s=%.6f", "K", consts.K);
    fprintf(f, " %s=%.6f", "Lambda", consts.lambda);
    fprintf(f, " %s=%.6f", "zeta", consts.zeta);
    fprintf(f, " %s=%.6f", "gamma", consts.gamma);
    fprintf(f, " %s=%.6f", "p_target_rel_change", consts.p_target_rel_change);
    fprintf(f, " %s=%.6d", "max_p_ters", consts.max_p_iters);
    
    fclose(f);


    double dudt[Lx][Ly][D];       // velocity update
    double dQdt[Lx][Ly][ncomps];       // Q update 
    double p_aux[Lx][Ly];         // auxiliary array for pressure updates
    double H[Lx][Ly][ncomps];       // molecular field 
    double S[Lx][Ly][ncomps];       // rotational terms for Q
    double Pi_S[Lx][Ly][D];    // stress tensor traceless-symmetric component
    double Pi_A[Lx][Ly];       // stress tensor antisymmetric component
    double div_u[Lx][Ly];       // holds divergence of velocity field
    double pressure_poisson_RHS[Lx][Ly]; // holds right-hand side of pressure-Poisson equation
    
    double E[Lx][Ly][D];
    double omega[Lx][Ly];
    double nx[Lx][Ly];
    double ny[Lx][Ly];
    double dor[Lx][Ly];
    double ss[Lx][Ly];


    // set all arrays to zero
    //memset(dudt, 0, sizeof(u));

    for(int i=0; i<Lx; i++)
    {
        for(int j=0; j<Ly; j++) 
        {
             dudt[i][j][0] = 0.0; 
             dudt[i][j][1] = 0.0;
        }
    }
    memset(dQdt, 0, sizeof(Q));
    memset(p_aux, 0, sizeof(p));
    memset(H, 0, sizeof(Q));
    memset(S, 0, sizeof(Q));
    memset(Pi_S, 0, sizeof(u));
    memset(Pi_A, 0, sizeof(p));
    memset(pressure_poisson_RHS, 0, sizeof(p));
    memset(div_u, 0, sizeof(p));
    
    memset(E, 0, sizeof(u));
    memset(omega, 0, sizeof(p));
    memset(nx, 0, sizeof(p));
    memset(ny, 0, sizeof(p));
    memset(dor, 0, sizeof(p));
    memset(ss, 0, sizeof(p));


    // other initializations
    double udiff = udiff_thresh + 1;
    int stepcount = 0;
    double t = 0.;

    int stepcount_last_save = - save_every_n_steps;
    printf("stepcount last saved: %d\n", stepcount_last_save);
    

    apply_Q_boundary_conditions(Q);
    apply_u_boundary_conditions(u);
    apply_p_boundary_conditions(p);
    
    char filename[50];
    sprintf(filename, "AverageP_%f.txt", zeta);
    FILE *file;
    file = fopen(strcat(p_mean_results_path, filename), "a");
    while(udiff > udiff_thresh & stepcount < max_steps & t < max_t)
    {
        //printf("While loop entered \n");
        clock_t start, end;
        double time_diff=0.;
        int numberOfUpdatesCalled = 0;

        
   
        if(stepcount - stepcount_last_save >= save_every_n_steps)
        {
            for (int i = 0; i < 2; ++i) 
            {
                char label[2];
                char file_name[100];
                strcpy(label, i == 0 ? "Q" : "u");
                sprintf(file_name, "%s%s_%010d.txt", i == 0 ? Q_results_path : u_results_path , label, stepcount);
                FILE* f = fopen(file_name, "w");
                for (int x = 0; x < Lx; x++) 
                {
                    for (int y = 0; y < Ly; y++) 
                    {
                        fprintf(f, "%f %f\n", i == 0 ? Q[(x*Ly + y)*ncomps + 0] : u[(x*Ly + y)*D + 0], 
                                                     i == 0 ? Q[(x*Ly + y)*ncomps + 1] : u[(x*Ly + y)*D + 1]);

                    }
                }
                fclose(f);
            }


            E_omega_from_u(u, &E[0][0][0], &omega[0][0]);
            n_dor_from_Q(Q, &nx[0][0], &ny[0][0], &dor[0][0]);
            get_saddle_splay(Q, &ss[0][0]);

            char label0[2];
            char file_name0[100];
            strcpy(label0, "E");
            sprintf(file_name0, "%s%s_%010d.txt", E_results_path, label0, stepcount);
            FILE* g = fopen(file_name0, "w");
            for (int x = 0; x < Lx; x++) 
                {
                    for (int y = 0; y < Ly; y++) 
                    {
                        fprintf(g, "%f %f\n", E[x][y][0] , E[x][y][1]); 
                    }
                }
            fclose(g);
            

            char label1[6];
            char file_name1[100];
            strcpy(label1, "n_dor");
            sprintf(file_name1, "%s%s_%010d.txt", n_dor_results_path , label1, stepcount);
            FILE* p1 = fopen(file_name1, "w");
            for (int x = 0; x < Lx; x++) 
                {
                    for (int y = 0; y < Ly; y++) 
                    {
                        fprintf(p1, "%f %f %f\n", nx[x][y], ny[x][y], dor[x][y]); 
                    }
                }
            fclose(p1);


            char label2[6];
            char file_name2[100];
            strcpy(label2, "omega");
            sprintf(file_name2, "%s%s_%010d.txt", omega_results_path , label2, stepcount);
            FILE* p2 = fopen(file_name2, "w");
            for (int x = 0; x < Lx; x++) 
                {
                    for (int y = 0; y < Ly; y++) 
                    {
                        fprintf(p2, "%f\n", omega[x][y]); 
                    }
                }            
            fclose(p2);

            char label3[3];
            char file_name3[100];
            strcpy(label3, "ss");
            sprintf(file_name3, "%s%s_%010d.txt", ss_results_path , label3, stepcount);
            FILE* p3 = fopen(file_name3, "w");
            for (int x = 0; x < Lx; x++) 
                {
                    for (int y = 0; y < Ly; y++) 
                    {
                        fprintf(p3, "%f\n", ss[x][y]); 
                    }
                }
            fclose(p3);


            char label4[2];
            char file_name4[100];
            strcpy(label4, "p");
            sprintf(file_name4, "%s%s_%010d.txt", p_results_path , label4, stepcount);
            FILE* p4 = fopen(file_name4, "w");
            for (int x = 0; x < Lx*Ly; x++) 
                {
                    fprintf(p4, "%f\n", p[x]);       
                }   
            fclose(p4);


            stepcount_last_save = stepcount;
        }
        

        else
        {
            int n_steps = t>0 ? jit_loops : 1 ;
            bool save_plot = 0;
            // do one step in first loop to catch errors and get printout
            start =  clock();
            update_step(u, &dudt[0][0][0] , Q, &dQdt[0][0][0], p, &p_aux[0][0], 
                                &pressure_poisson_RHS[0][0], &S[0][0][0], &H[0][0][0], &Pi_S[0][0][0], &Pi_A[0][0],
                                    consts, &stepcount, &t, &udiff, n_steps, Lx, Ly, D, ncomps);

            numberOfUpdatesCalled+=1;
            end = clock();
            printf("dudt: %f\n", dudt[0][0][0]);
        }
        
        time_diff += (double)(end - start) / CLOCKS_PER_SEC; 
        printf("average time difference is: %f\n", time_diff/(numberOfUpdatesCalled));
        
        div_vector(u, &div_u[0][0], Lx, Ly);
  

        double pressureSum = 0.0;
        for(int i=0; i<Lx*Ly; i++){pressureSum += p[i];}
        fprintf(file, "%f  %d\n",(double) pressureSum/(Lx*Ly), stepcount);
        

        printf("stepcount: %d\n", stepcount);
       
        printf("u: %f\n", *u);
    }
    fclose(file);
}




void ss_plot_function(double *ss, double vmax, double *result) 
{
    //Threshold on aboslute value of saddle-splay to show defect
    for (int x = 0; x < Lx; x++) 
    {
        for (int y = 0; y < Ly; y++) 
        {
            if (fabs(ss[x*Ly + y]) < vmax) {result[x*Ly + y] = NAN;} 
            else {result[x*Ly + y] = ss[x*Ly + y];}
        }
    }
}



    


int main(void)
{
    char resultspath[] = "./Results/";
    char runlabel[] = "test ";
    char imgpath[] = "./Images/";

    char bc_label[] = "Periodic boundary";
    char run_label[] = "my_test_name";

    //package parameters as dictionary
    struct consts_dict consts = 
    {
        .dt = dt,
        .dt_pseudo = dt_pseudo,
        .nu = nu,
        .A = A,
        .C = C,
        .K = K,
        .lambda = lambda,
        .zeta = zeta,
        .gamma = gamma,
        .p_target_rel_change = p_target_rel_change,
        .max_p_iters = max_p_iters
    };

    double theta_initial[Lx*Ly];
    srand(time(NULL));
    for(int x = 0; x < Lx; x++) {
        for(int y = 0; y < Ly; y++) {
            theta_initial[x*Ly + y] = M_PI * rand() / RAND_MAX;
            theta_initial[x*Ly + y] += 1.0 * M_PI * rand() / RAND_MAX;
            //theta_initial[x*Ly + y] = M_PI/2;
        }
    }
    


    /*  
    char file_holder[100];       
    sprintf(file_holder, "Theta_initial.txt");
    FILE* pf = fopen(file_holder, "w");
    for (int i = 0; i < Lx*Ly; i++) 
        {
            fprintf(pf, "%f\n", theta_initial[i]); 
        }
    fclose(pf);
    */

    /*
    FILE *file = fopen("Theta_initial.txt", "r");
    if(file == NULL){perror("Error opening file"); return 1;}
    char buffer[1024];
    int valueCount = 0;
    while (fgets(buffer, sizeof(buffer), file) != NULL)
    {
        buffer[strcspn(buffer, "\n")] = '\0';
        double value = strtod(buffer, NULL);
        theta_initial[valueCount] = value;
        //printf("Theta initial: %f\n", theta_initial[valueCount]);

        valueCount++;
    }
    fclose(file);
    printf("valueCount: %d\n", valueCount);
    */
    

    // initialize the velocity field, pressure field, and Q-tensor
    double u[Lx][Ly][D];
    for (int x = 0; x < Lx; x++)
    {
        for(int y = 0; y < Ly; y++)
        {
            u[x][y][0] = 0.0; 
            u[x][y][1] = 0.0;
        }
    }


    double p[Lx][Ly];
    for (int x = 0; x < Lx; x++) 
    {
        for(int y = 0; y < Ly; y++){p[x][y] = 0.0;}
    }

    double Q[Lx][Ly][ncomps];
    initialize_Q_from_theta(&Q[0][0][0], theta_initial, S0, Lx, Ly, ncomps);

    printf("Boundary conditions: %s\n", bc_label);


    run_active_nematic_sim(&u[0][0][0], &Q[0][0][0], &p[0][0], run_label, consts, Lx, Ly);

    return 0;
}