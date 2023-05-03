#include<iostream>
#include<fstream>
#include<vector>
#include<cstdlib>
#include<math.h>
#include<complex>
#include<sstream>
#include<iomanip>
#include<numeric>
using namespace std ;

#define PI (3.141592653589793238462643)
const int q=3; // q spin states
const int L=32; // linear system size
//const double T=0.25;// temperature in units of J, wants 0.25 and 0.5
const double MAX_TEMP = 2.0;
const double TEMP_STEP = 0.02;

constexpr int power(int base, int exp) {
    int result = 1;
    for (int i = 0; i < exp; i++) {
        result *= base;
    }
    return result;
}

const int NDIM=2;
const int N = power(L, NDIM);
// total number of spins
//const double pconnect = 1-exp(-1/T); // connection probability

const int NCLUSTERS=1; // # of cluster builds in one MC step.
const int NESTEPS=10000; // # of equilibrium MC step.
const int NMSTEPS=10000; // # of measurement MC step .
const int NBINS=1 ; // #of measurement bins: "experiments"

vector<int> S(N); // the spin array: the spin state at each site, now 2D
vector<int>M(q); // number of spins in different states: q=3 for us
vector<complex<double> > W(q); // order parameter weights

void printLattice(int q, int L, double T, int dim) {
    
    cout << "q = " << q << endl;
    cout << "L = " << L << endl;
    cout << "T = " << T << "J" << endl;
    cout << NDIM << " dimension\n" << endl;
}


int writeCorrelation(string prefix, complex<double> *m0cmr, complex<double> *mr, complex<double> m0c);
//int writeMag(string filename, double mean_mag, double T);
int writeMoment(string filename, double T, complex<double> m, double m1, double m2, double m4);

// ERROR CODE
#define ERR_NBR -1
#define ERR_FILE -2
// lattice handling:
enum dirs{RIGHT,LEFT,UP,DOWN};
int indx(int x, int y) {return L*y+x;}; //make an indx on every site
int xpos(int i) {return i%L;}
int ypos(int i) {return i/L;}

int Nbr(int i, int dir)
{
    int x=xpos(i);
    int y=ypos(i);
    switch(dir) {
        case RIGHT: return indx((x+1)%L,y);
        case LEFT: return indx((x-1+L)%L,y);
        case UP: return indx(x,(y+1)%L);
        case DOWN: return indx(x,(y-1+L)%L);
    }
    return ERR_NBR;
}

void FlipandBuildFrom (int s, double pconnect) {
    int oldstate(S[s]) ,newstate((S[s]+1)%q);
    S[s]=newstate; // flip spin
    M[oldstate]--;
    M[newstate]++;
    
    for ( int dir =0; dir < NDIM*2; dir++) // go thru neighbours
    {
        int j=Nbr(s,dir);
        if(S[j] == oldstate)
            if(rand()/(RAND_MAX+1.) <pconnect)
               {FlipandBuildFrom(j, pconnect);}
    }
}

string makePrefix(int q, int N, double T) {
    string filename;
    ostringstream T_2d;
    T_2d << std::fixed << std::setprecision(2) << T;
    std::string T_str = T_2d.str();
    filename = to_string(q) + "," + to_string(L) + "," + T_str;
    return filename;
}

int main()
{
    // initialize order parameter weights
    
    string m_filename = to_string(L) + " magnetisations_vs_T.csv";
    // writes header
    ofstream m_file(m_filename);
    m_file << "T," << "m," << "m1," << "m2," <<"m4" <<endl;
    m_file.close();
    
    double T = 0.0;
    while (T < MAX_TEMP) {
        double pconnect = 1-exp(-1/T);
        printLattice(q,L,T,NDIM);
        for(int s=0; s<q; s++) W[s] = complex<double>(cos(2*PI*s/q), sin(2*PI*s/q)); // W is the array of all possible magnetisation
        for(int i=0; i<N; i++) S[i]=0; // initialize to the spin=0 state, S is the list of all spins
        for(int s=1; s<q; s++) M[s]=0; // initialize counters, M is how many of are in each state
        M[0]=N;
        srand((unsigned) time(0)); // initialize random number gen.
        
        // equilibriate
        for(int t=0; t<NESTEPS; t++)
            for(int c=0; c<NCLUSTERS; c++)
            {
                FlipandBuildFrom(rand()%N, pconnect);
            }
        
        // measure
        for(int n=0; n<NBINS; n++)
        {
            complex<double> m(0. ,0.);
            double m1=0, m2=0, m4=0; //
            complex<double> m0c = 0.;
            complex<double> mr[N], m0cmr[N] = {0.}; // container for all of r's
            
            
            
            for(int t=0; t<NMSTEPS; t++) {
                
                // build and flip
                for(int c=0; c<NCLUSTERS; c++) FlipandBuildFrom(rand()%N, pconnect);
                complex<double> tm(0. ,0.);
                
                
                for(int s=0; s<q; s++){tm+=W[s]*double(M[s]);}
                tm/=N;
                double tm1=abs(tm); double tm2=tm1*tm1;
                m+=tm ; m1+=tm1 ; m2+=tm2 ; m4+=tm2*tm2;
                // add all MC steps for m_r
                
                
                m0c += conj(W[S[0]]);
                
                for (int j=0; j<N; j++) {
                    mr[j] += W[S[j]];
                    m0cmr[j] += conj(W[S[0]]) * W[S[j]];
                }
                
            }
            
            m/=NMSTEPS ; m1/=NMSTEPS ; m2/=NMSTEPS ; m4/=NMSTEPS;
            m0c /= NMSTEPS;
            for (int j=0; j<N; j++) {
                mr[j]/=NMSTEPS;
                m0cmr[j]/=NMSTEPS;
            }
            
            
            // write to file
            
            // file prefix
            string prefix = makePrefix(q, L, T);
            
            
            // Correlation
            // compute C(r), just the first term or both? Both!

            if (NDIM == 1 && (T == 0.25 || T == 0.5)) {
                if (writeCorrelation(prefix, m0cmr, mr, m0c))
                    return ERR_FILE;
            }
            
            
            // moments
            if (writeMoment(m_filename, T, m, m1, m2, m4))
                return ERR_FILE;
            
        }
        T += TEMP_STEP;
        cout << "average magnetisation per site written to " + m_filename << endl;

    }
}

int writeCorrelation(string prefix, complex<double> *m0cmr, complex<double> *mr, complex<double> m0c) {
    complex<double> Cr[N];
    string corr_filename = "correlation," + prefix + ".csv";
    ofstream corr_file(corr_filename);
    
    if (!corr_file.is_open())
        return ERR_FILE;
    for (int j=0; j<N; j++) {
        Cr[j] = m0cmr[j] - m0c * mr[j];
        //        cout << j << "," << Cr[j].real() << endl;
        corr_file << j << "," << Cr[j].real() << endl;
    }
    corr_file.close();
    cout << "Correlation written to file: " + corr_filename << endl;
    return 0;
}
    


int writeMoment(string m_filename, double T, complex<double> m, double m1, double m2, double m4) {
    ofstream moment_file(m_filename, ios_base::app);

    if (!moment_file.is_open())
        return ERR_FILE;

    moment_file << T << "," << m.real() << "," << m1 << "," << m2 << "," << m4 << endl ;
    moment_file.close();
    cout << "moment written to file: " + m_filename << endl;
    return 0;
}
