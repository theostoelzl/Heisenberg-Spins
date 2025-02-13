/*
COMPILE THE CODE WITH: g++ -std=c++11 -g -O3 -mcmodel=medium IL_MC_heatbath_split_cubic.cpp -o IL_MC_cubic.o
 (change -mcmodel=large if you want bigger Nh grid)
RUN THE CODE WITH: ./MC.o [spin input file] [temp] [rng_seed]
 ([spin input file]=none if you want to start from a random arrangement)
 
"system.txt" should exist in the directory with information as follows:
jx jy jz (?)
u0 umin mmin
eqsweeps avsweeps
 (j1 should be +/-1)
*/
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <random>
#include <string>
#include <vector>
#include <iomanip>

using namespace std;

//array size for spin lattice
#define n1 20
#define n2 20
#define n3 20
double spins[n1][n2][n3];
double J[3];

//grid of points for CDF integral
#define Nh 100001
#define Nm 1024 //ALWAYS 2^N (where N is some integer)
double CDF[Nh][Nm]; //cumulative distribution function

//a fudge factor that should be ignored if possible
double x;

//Monte carlo functions
int PBC(int n, int nmax);
double local_field(double arr[n1][n2][n3], double js[3], int i, int j, int k);
double total_X(double arr[n1][n2][n3], double js[3], double u0, double umin, double mmin);
double total_U(double arr[n1][n2][n3], double js[3], double u0, double umin, double mmin);
double total_energy(double arr[n1][n2][n3], double js[3], double u0, double umin, double mmin);

//Functions for sampling prob dist
double boltz(double h, double m, double u0, double umin, double mmin, double kT); //bolztmann factor
double part_func(double h, double u0, double umin, double mmin, double kT); //1site partition function
void gen_CDF(double hmin, double hmax, double u0, double umin, double mmin, double kT); //does integrals
double int_CDF(double h, double hmin, double hmax, double m); //interpolates CDF
//new functions for the mapping function with an irregular grid
//Non-uniform grid in R space. Points are precisley where they lie on CDF axis. So no information is lost
int R_to_array(int harr, double R);
double int_M_adaptive(double h, double hmin, double hmax, double R);

int main(int argc, char *argv[]){
    //debug
    ofstream debug;
    debug.open("debug.txt");
    
    //command line input arguments
    string spinin_file=argv[1]; //source of spin structure
    double kT=stod(argv[2]);
    double rngseed=stod(argv[3]);
    
    //system parameters to be read from system file
    double Jtot;
    double u0, umin, mmin; //umin is a relic, always zero
    int eqsweeps, avsweeps, sweeps; //equilibriation, averaging and total num of swps
    
    //reading system information
    ifstream systin_str("system.txt");
    systin_str >> J[0] >> J[1] >> J[2]
               >> u0 >> umin >> mmin
               >> eqsweeps >> avsweeps
               >> x;
    sweeps=eqsweeps+avsweeps;
    // overwrite u0
    u0 = stod(argv[4]);
    // total J is sum of nearest neighbours all aligned with spin 1
    Jtot=4*abs(J[0]) + 4*abs(J[1]) + 4*abs(J[2]);
    
    //scaling units for m* in calculation
    //calculation is done in units where m*=1 (historical reasons)
    u0*=(mmin*mmin);
    kT*=(mmin*mmin);
    
    //random number things
    mt19937 rng; //mersenne twister generator
    rng.seed(time(NULL)+100000*rngseed);
    uniform_int_distribution<int> site_picker1(0,n1-1);
    uniform_int_distribution<int> site_picker2(0,n2-1);
    uniform_int_distribution<int> site_picker3(0,n3-1);
    uniform_real_distribution<double> uni_dist(0,1);
    
    //assigning initial spins
    cout << "(*------------*)\n";
    if(spinin_file=="none"){
        cout << "(*) Randomising initial spin config... ";
        for(int i=0; i<n1; i++){
            for(int j=0; j<n2; j++){
                for(int k=0; k<n3; k++){
                    spins[i][j][k] = 2*uni_dist(rng)-1;
                }
            }
        }
    }
    else{
        ifstream spin_in;
        spin_in.open(spinin_file);
        cout << "(*) Reading initial spin config from "+spinin_file+"... ";
        for(int i=0; i<n1; i++){
            for(int j=0; j<n2; j++){
                for(int k=0; k<n3; k++){
                    spin_in >> spins[i][j][k];
                }
            }
        }
        spin_in.close();
    }
    cout << "DONE\n";
    
    //setting up probability distribution sampler
    double hmin=-Jtot, hmax=Jtot;
    cout << "(*) Generating cumulative distribution function... " << flush;
    gen_CDF(hmin, hmax, u0, umin, mmin, kT);
    ofstream cdf_file;
    cdf_file.open("cdf.txt", fstream::app);
    for (int i=0; i<Nm; i++) {
	cdf_file << CDF[70000][i] << endl;
    }
    cdf_file.close();
    cout << "DONE\n";

    //Debug for testing probability distribution
//    ofstream foutp;
//    ofstream foutd ("ha.txt");
//    for(int i=1; i<=20; i++){
//        double h = (hmax-hmin)*(2*uni_dist(rng)-1)/2;
//
//
//        foutd << i << "\t" << h << endl;
//        cout << i << "\t" << h << endl;
//
//        foutp.open("p_sample_"+to_string(i)+".txt");
//        for(int j=0; j<50000; j++){
//            foutp << int_M(h, hmin, hmax, uni_dist(rng)) << endl;
//        }
//        foutp.close();
//    }
    
    //Initialising averages. In order:
    // energy, squared energy, total spin, total spin squared
    // X energy, squared X energy, U energy, squared U energy
    // spin per site, squared spin per site
    // correlation function (currently out of action)
    double en_avg=0, en2_avg=0, s_avg=0, s2_avg=0;
    double ex_avg=0, ex2_avg=0, eu_avg=0, eu2_avg=0;
    double si_avg[n1][n2][n3]={ }, si2_avg[n1][n2][n3]={ };
    double s_corr[n3]={ };
    //Running energy totals
    double enX=total_X(spins, J, u0, umin, mmin);
    double enU=total_U(spins, J, u0, umin, mmin);
    double toten=enX+enU;
    // absolute total magnetisation per sweep
    double abs_magn[sweeps] = {};
    int avsamp=0; //counts the number of data points in average
    
    cout << "(*) Starting Monte Carlo simulation at temperature " << kT/(mmin*mmin) << "... " << flush;
    for(int i=0; i<sweeps; i++){
        for(int j=0; j<n1*n2*n3; j++){
            //pick random site and calculate field
            int a=site_picker1(rng);
            int b=site_picker2(rng);
            int c=site_picker3(rng);
            double s_old=spins[a][b][c];
            double h=local_field(spins, J, a, b, c);
   
            //computing initial local energy
            double ex_before=-s_old*h;
            double eu_before=s_old*(u0+umin)/(mmin*mmin)*(-2*s_old + s_old*s_old*s_old/(mmin*mmin));
            double len_before=ex_before+eu_before;
          
            //pick new spin from distribution
            double r=uni_dist(rng);
            double s_new=int_M_adaptive(h, hmin, hmax, r);
            
            //compute change in energy and assign new spin
            double ex_after=-s_new*h;
            double eu_after=s_new*(u0+umin)/(mmin*mmin)*(-2*s_new + s_new*s_new*s_new/(mmin*mmin));
            double len_after=ex_after+eu_after;
            toten+=len_after-len_before;
            enX+=ex_after-ex_before;
            enU+=eu_after-eu_before;
            spins[a][b][c]=s_new;
        }
       // debug << toten << endl;
        if(i>=eqsweeps && i%1==0){
            avsamp++;
            en_avg+=toten;
            ex_avg+=enX;
            eu_avg+=enU;
            en2_avg+=toten*toten;
            ex2_avg+=enX*enX;
            eu2_avg+=enU*enU;
	    abs_magn[i] = 0;
            for(int a=0; a<n1; a++){
                for(int b=0; b<n2; b++){
                    for(int c=0; c<n3; c++){
                        double s=spins[a][b][c];
                        si_avg[a][b][c]+=s;
			abs_magn[i] += s;
                        si2_avg[a][b][c]+=s*s;
//                        if(a==0 && b==0){ //correlation function
//                            for(int m=0; m<n3; m++){
//                                s_corr[c]+=spins[a][b][c]*spins[a][b][PBC(m+c,n3)];
//                            }
//                        }
                    }
                }
            }
	    abs_magn[i] = abs(abs_magn[i]);
        }
        
    }
    cout << "DONE\n(*------------*)\n";
   
    //output streams
    ofstream enavg_out;
    ofstream en2avg_out;
    ofstream eSavg_out;
    ofstream eS2avg_out;
    ofstream savg_out;
    ofstream abs_magn_out;
    ofstream s2avg_out;
    ofstream siavg_out;
    ofstream si2avg_out;
    ofstream scorr_out;
    ofstream s_out;
    ofstream a_out;
    
    //output files
    enavg_out.open("energy.txt", fstream::app);
    en2avg_out.open("energy2.txt", fstream::app);
    eSavg_out.open("energysplit.txt", fstream::app);
    eS2avg_out.open("energysplit2.txt", fstream::app);
    savg_out.open("spin_total.txt", fstream::app);
    abs_magn_out.open("abs_magnetisation.txt", fstream::app);
    s2avg_out.open("spin2_total.txt", fstream::app);
    siavg_out.open("spins.txt", fstream::app);
    si2avg_out.open("spins2.txt", fstream::app);
    scorr_out.open("spin_corr.txt", fstream::app);
    s_out.open("spins_after.txt"); //snapshot of spins at the end for resuming
    
    avsweeps=avsamp;
    enavg_out << kT/(mmin*mmin) << "\t" << setprecision(12) << en_avg/avsweeps/(mmin*mmin) << endl;
    en2avg_out << kT/(mmin*mmin) << "\t" << setprecision(12) << en2_avg/avsweeps/(mmin*mmin*mmin*mmin) << endl;
    eSavg_out << kT/(mmin*mmin) << "\t" << setprecision(12) << ex_avg/avsweeps/(mmin*mmin) << "\t" << eu_avg/avsweeps/(mmin*mmin) << endl;
    eS2avg_out << kT/(mmin*mmin) << "\t" << setprecision(12) << ex2_avg/avsweeps/(mmin*mmin) << "\t" << eu2_avg/avsweeps/(mmin*mmin) << endl;    siavg_out << kT/(mmin*mmin) << "\t";
    si2avg_out << kT/(mmin*mmin) << "\t";
    for(int i=0; i<n1; i++){
        for(int j=0; j<n2; j++){
            for(int k=0; k<n3; k++){
                siavg_out << si_avg[i][j][k]/avsweeps << "\t";
                si2avg_out << si2_avg[i][j][k]/avsweeps << "\t";
                s_avg+=si_avg[i][j][k];
                s2_avg+=si2_avg[i][j][k];
                s_out << spins[i][j][k] << endl;
                if(i==0 && j==0){
                    scorr_out << s_corr[k]/avsweeps/(double)n3 << "\t";
                }
            }
        }
    }
    double abs_magn_sum;
    abs_magn_sum = 0;
    for (int i=0; i<sweeps; i++) {
	if (i >= eqsweeps) {
		abs_magn_sum += abs_magn[i];
	}
    }
    abs_magn_out << kT/(mmin*mmin) << "\t" << abs_magn_sum/avsweeps << endl;
    savg_out << kT/(mmin*mmin) << "\t" << s_avg/avsweeps << endl;
    s2avg_out << kT/(mmin*mmin) << "\t" << s2_avg/avsweeps << endl;
    
    enavg_out.close();
    en2avg_out.close();
    savg_out.close();
    abs_magn_out.close();
    s2avg_out.close();
    siavg_out.close();
    si2avg_out.close();
    scorr_out.close();
    s_out.close();

    return 0;
}

//***********************************//
//-----------------------------------//
//       Monte Carlo FUNCTIONS       //
//-----------------------------------//
//***********************************//

//PBC=Periodic Boundary Conditions
int PBC(int n, int nmax){
    while(n>=nmax){ n-=nmax; }
    while(n<0){ n+=nmax; };
    return n;
}

//Computes local field. This is where lattice information enters
double local_field(double arr[n1][n2][n3], double js[3], int i, int j, int k){
    i=PBC(i,n1);
    j=PBC(j,n2);
    k=PBC(k,n3);

    return (js[0]*(arr[i][j][PBC(k+1,n3)] + arr[i][j][PBC(k-1,n3)]) +
            js[1]*(arr[i][PBC(j+1,n2)][k] + arr[i][PBC(j-1,n2)][k]) +
            js[2]*(arr[PBC(i+1,n1)][j][k] + arr[PBC(i-1,n1)][j][k]));
}

//Computes TOTAL energy, computed using local_field function
double total_energy(double arr[n1][n2][n3], double js[3], double u0, double umin, double mmin){
    double e=0;
    for(int i=0; i<n1; i++){
        for(int j=0; j<n2; j++){
            for(int k=0; k<n3; k++){
                double m=arr[i][j][k];
                e+=u0+m*(-0.5*local_field(arr, js, i, j, k) +
                                    (u0+umin)/(mmin*mmin)*(-2*m + m*m*m/(mmin*mmin))
                                    );
            }
        }
    }
    
    return e;
}

//Computes EXCHANGE energy
double total_X(double arr[n1][n2][n3], double js[3], double u0, double umin, double mmin){
    double ex=0;
    for(int i=0; i<n1; i++){
        for(int j=0; j<n2; j++){
            for(int k=0; k<n3; k++){
                double m=arr[i][j][k];
                ex+=-0.5*m*local_field(arr, js, i, j, k);
            }
        }
    }
    
    return ex;
}

//Compute ON-SITE energy
double total_U(double arr[n1][n2][n3], double js[3], double u0, double umin, double mmin){
    double eu=0;
    for(int i=0; i<n1; i++){
        for(int j=0; j<n2; j++){
            for(int k=0; k<n3; k++){
                double m=arr[i][j][k];
                eu+=u0+m*((u0+umin)/(mmin*mmin)*(-2*m + m*m*m/(mmin*mmin)));
            }
        }
    }
    
    return eu;
}

//***********************************//
//-----------------------------------//
//        M sampler FUNCTIONS        //
//-----------------------------------//
//***********************************//

//boltzmann factor
double boltz(double h, double m, double u0, double umin, double mmin, double kT){
    return exp( -(x+u0+m*(-h + (u0+umin)*(-2*m + m*m*m/(mmin*mmin))/(mmin*mmin)))/kT );
}

//partition function integral over m
// uses trapezium method
double part_func(double h, double u0, double umin, double mmin, double kT){
    double dm=4.0/(Nm-1);
    double Zmid=0;
    for(int i=0; i<Nm; i++){
        Zmid+=boltz(h, -2+i*dm, u0, umin, mmin, kT);
    }
    double z=0.5*dm*(boltz(h, -2, u0, umin, mmin, kT)+2*Zmid+boltz(h, 2, u0, umin, mmin, kT));
    return z;
}

//integrates probability (boltz/Z) to get cumulative dist function
void gen_CDF(double hmin, double hmax, double u0, double umin, double mmin, double kT){
    double h, dh=(hmax-hmin)/(Nh-1); //because we need to data points at either end!
    double dm=4.0/(Nm-1);
    
    double Z, cdf, cdf_add;
    
    for (int i=0; i<Nh; i++){
        h=hmin+i*dh;
        //calculate partition function for this h
        Z=part_func(h, u0, umin, mmin, kT);
        
        //calculate CDF
        // Lots of safety features in here, best to leave them be
        cdf=0;
        CDF[i][0]=0; //force the first point to be 0
        for(int j=1; j<Nm; j++){
            cdf_add = 0.5*dm*(boltz(h, -2+(j-1)*dm, u0, umin, mmin, kT)
                         + boltz(h, -2+j*dm, u0, umin, mmin, kT));
            cdf += cdf_add;

	    if(cdf==0){
                CDF[i][j]=0;
            }
            else if(std::isinf(cdf)==1){
                CDF[i][j]=1;
            }
            else{
                if(std::isnan(Z)==1){
                    cout << "NOT DONE\n Z is nan, pick better X\n";
                    exit(-1);
                }
                else if(std::isinf(Z)==1){
                    cout << "NOT DONE\n Z is infinite, pick better X\n";
                    exit(-1);
                }
                else if(Z==0){
                    cout << "NOT DONE\n Z is zero, pick better X\n";
                    exit(-1);
                }
                CDF[i][j]=cdf/Z;
            }
        }
        for(int j=0; j<Nm; j++){
            CDF[i][j] = CDF[i][j] / CDF[i][Nm-1]; //force the last point to be 1 by normalising
	}
    }
}

//Bilinear interpolation of the CDF array
double int_CDF(double h, double hmin, double hmax, double m){
    // find grid points either side of h coord
    double dh=(hmax-hmin)/(Nh-1);
    double harr=(h-hmin)/dh;
    int h1arr=floor(harr), h2arr=ceil(harr);
    
    //find grid points either side of m coord
    double dm=4.0/(Nm-1);
    double marr=(m+1.0)/dm;
    int m1arr=floor(marr), m2arr=ceil(marr);
    
    double outp;
    
    if(h1arr!=h2arr){
        if(m1arr!=m2arr){
            double m1=-1+dm*m1arr, m2=-1+dm*m2arr;
            double h1=hmin+dh*h1arr, h2=hmin+dh*h2arr;
            outp = (CDF[h1arr][m1arr]*(h2-h)*(m2-m) +
                    CDF[h2arr][m1arr]*(h-h1)*(m2-m) +
                    CDF[h1arr][m2arr]*(h2-h)*(m-m1) +
                    CDF[h2arr][m2arr]*(h-h1)*(m-m1))/(h2-h1)/(m2-m1);
        }
        else if(m1arr==m2arr){
            double h1=hmin+dh*h1arr, h2=hmin+dh*h2arr;
            outp = CDF[h1arr][(int)marr] + (CDF[h2arr][(int)marr]-CDF[h1arr][(int)marr])*(h-h1)/(h2-h1);
        }
    }
    else{
        if(m1arr!=m2arr){
            double m1=-1+dm*m1arr, m2=-1+dm*m2arr;
            outp = CDF[(int)harr][m1arr] + (CDF[(int)harr][m2arr]-CDF[(int)harr][m1arr])*(m-m1)/(m2-m1);
        }
        else if(m1arr==m2arr){
            outp = CDF[(int)harr][(int)marr];
        }
    }
    
    return outp;
}

//Converts a given value of R into array units. Gives the array value below it.
// Given a random number 0<R<1, this scans through the CDF and finds the
// grid point in the CDF array below where this R corresponds too
int R_to_array(int harr, double R){
    int lowi=0, highi=Nm-1;
    int midhighi=(highi+1)/2, midlowi=midhighi-1;
    double lowR, midlowR, midhighR, highR;
    //Uses a bisection type algorithm, which is why
    // Nm should always be a power of 2
    while(highi-lowi!=1){
        lowR=CDF[harr][lowi];
        highR=CDF[harr][highi];
        midlowR=CDF[harr][midlowi];
        midhighR=CDF[harr][midhighi];
        
        if(R>lowR && R<midlowR){ //if it's in the bottom half
            highi=midlowi;
            midhighi=lowi+(highi+1-lowi)/2;
            midlowi=midhighi-1;
        }
        else if(R>midhighR && R<highR){ //if it's in the top half
            lowi=midhighi;
            midhighi=lowi+(highi+1-lowi)/2;
            midlowi=midhighi-1;
        }
        else if(R>midlowR && R<midhighR){ //if it's in the middle
            lowi=midlowi;
            highi=midhighi;
        }
    }
    
    return lowi;
}

//Mapping function. Returns m for any R and h by linearly interpolating array
// This inverts the CDF *on the fly* using R_to_array. Four points to be
// interpolated between are now *non-rectangular*, so need to use something fancy
double int_M_adaptive(double h, double hmin, double hmax, double R){
    double dm=4.0/(Nm-1.0);
    double dh=(hmax-hmin)/(Nh-1);
    double harr=(h-hmin)/dh;
    int h1arr=floor(harr), h2arr=ceil(harr);
    
    double outp;
    
    if(h1arr!=h2arr){
        //converting R into array units
        int R11arr=R_to_array(h1arr, R);
        int R12arr=R11arr+1;
        int R21arr=R_to_array(h2arr, R);
        int R22arr=R21arr+1;

        //all the coords we need
        double h1=hmin+dh*h1arr, h2=hmin+dh*h2arr;
        double R11=CDF[h1arr][R11arr], R12=CDF[h1arr][R12arr];
        double R21=CDF[h2arr][R21arr], R22=CDF[h2arr][R22arr];
        
        //now do some algebra with ^^those^^ things to interpolate
        //convenient definitions
        double Dh=h2-h1;
        double DR1=R12-R11, DR2=R22-R21;
        //quadratic coefficients
        double a=Dh*(DR1-DR2);
        double b=(h1-h)*(DR1-DR2)-Dh*DR1;
        double g=(h-h1)*DR1;
        
        //finding s and t
        double s, s1, s2, t;
        if(a==0){ s=-g/b;} //don't need to solve quadratic (often true to double prec)
        else{ //need to solve quadratic
            double ss=sqrt(b*b-4*a*g);
            s1=(-b+ss)/(2*a), s2=(-b-ss)/(2*a);
            if(s1>0 && s1<1){ s=s1; }
            else{ s=s2; }
        }
        t=(R-R12-s*(R22-R12))/(R11-R12+s*(DR1-DR2));
        
        //bilinear interpolation
        outp = (-2+dm*R12arr)*(1-s)*(1-t) +
        (-2+dm*R22arr)*s*(1-t) +
        (-2+dm*R11arr)*(1-s)*t +
        (-2+dm*R21arr)*s*t;
    }
    else{//(if we lie directly on a grid point) NEVER EVEN GETS USED SO IS INCOMPLETE
        cout << "Should have implemented this, huh";
	    
	int R1arr=R_to_array(h1arr, R);
        int R2arr=R1arr+1;
        
        double M1=-2+dm*R1arr;
        double M2=-2+dm*R2arr;
        double R1=CDF[h1arr][R1arr];
        double R2=CDF[h1arr][R2arr];
        outp = M1 + (R-R1)*(M2-M1)/(R2-R1);
    }
    
    return outp;
}

