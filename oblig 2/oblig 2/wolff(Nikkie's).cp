#include <iostream>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <math.h>
#include <complex>
#include <iomanip>

using namespace std;

#define PI 3.14159265358979323846
//k_B = 1;
const int q = 3; //number of spin states; we have 3: 0, 1 and 2
const int L = 16; //number of sites; L = N
const double T = 0.25; //temperature in units J
//changing temperature has to be done manually, a bit hacky but it will have to do for now

int dim = 1; //dimension of the lattice; 1D chain case

const int N = L;
const double pconnect = 1. - exp(-1. / T); //probability of connecting two sites; this expression works because our temperature is in units of J
//might need to make a function for connecting probability so that it can be computed for different temperatures

const int Nclusters = 1; //number of cluster builds in one MC step
const int NEsteps = 100000; //number of equilibrium MC steps - for equilibration/burn-in i assume
const int NMCsteps = 100000; //number of MC steps for measurement
const int Nbins = 10; //number of measurement bins/experiments

vector<int> Spins(N); //array of spins
vector<int> M(q); //number of spins in different states - will be used in flipping the spins
vector<complex<double>> W(q); //oder parameter weights
vector<complex<double>> m0r(N); //order parameter for each site
vector<complex<double>> mr(N); //order parameter for each site
vector<complex<double>> Cr(N); //order parameter for each site


//function to initialize the lattice with the spins
enum directions {right, left}; //1D chain case
// enum directions2D {right, left, up, down}; //2D square lattice case
// int indx(int x, int y) {return y * L + x;}; //make an indx on every site


int indx(int x) {return x;}
//https://stackoverflow.com/questions/49124741/implement-2d-array-coordinates-in-1d-array-in-c
//the source above explains the intuition behind converting 2D lattice coordinates to 1D array
int xpos(int i) {return i%L;}
// int ypos(int i) {return i / L;}

int Nbr(int i, int dir)
{
    int x = xpos(i);
    // int y = ypos(i);

    switch(dir)
    {
        // case ::right: return indx((x + 1)%L, y); //return index of the right neighbour
        // case ::left: return indx((x - 1 + L)%L, y); //return index of the left neighbour
        // case ::up: return indx(x, (y + 1)%L);
        // case ::down: return indx(x, (y - 1 + L)%L);

        case ::right: return indx((x + 1)%L); //return index of the right neighbour
        case ::left: return indx((x - 1 + L)%L); //return index of the left neighbour
        default: return -1; //if the direct ion is not right or left, return -1
    
    }

}

//implementing cluster-flipping algorithm
void FlipAndBuild_from(int s) //might need to implement temperature as an argument later as well when looking for phase transitions
{
    int oldstate(Spins[s]);
    int newstate((Spins[s] + 1) % q); //save the old state of the site
    //newstate is modulated so that it is a state between 0,1 and 2 in our case of Potts model with q = 3
    Spins[s] = newstate; //flip the spin of the site
    M[oldstate]--; //decrease the number of spins in the old state by 1
    M[newstate]++; //increase the number of spins in the new state by 1

    for (int dir = 0; dir < 2; dir++) //loop over all neighbours s
    //2 neighbours in 1D case, 4 neighbours in 2D case
        {
            int j = Nbr(s, dir); //get the index of the neighbour
            if (Spins[j] == oldstate && std::rand() / (RAND_MAX + 1.) < pconnect) //if the neighbour has the same spin as the site and the random number is less than the probability of connecting the two sites
            {
                FlipAndBuild_from(j); //call the function again
            }
        }
}


int main()
{
    for(int s = 0; s < q; s++)
    {
        W[s] = complex<double>(cos(2 * PI * s / q), sin(2 * PI * s / q));
    }
    
  
    for(int i = 0; i < N; i++)
    {
        Spins[i] = 0; //initialise all spins to be in the 0 state
        m0r[i] = 0; //initialise all m0r to be 0
        mr[i] = 0; //initialise all mr to be 0
    }

    for(int s = 1; s < q; s++)
    {
        M[s] = 0; //initialise counters
    }
    M[0] = N;


    srand((unsigned) time(0)); //seed the random number generator


    //equilibration/burn-in
    for(int t = 0; t < NEsteps; t++)
    {
        for(int c = 0; c < Nclusters; c++)
        {
            FlipAndBuild_from(rand() % N); //flip the spin of a random site and build the cluster
        }
    }

    vector<complex<double>> m(Nbins);
    double m1[Nbins];
    double m2[Nbins];
    double m4[Nbins];

    //do your measurements
    for(int n = 0; n < Nbins; n++)
    {
        // vector<double> m1(N, 0.);
        // vector<double> m2(N,0);
        // vector<double> m4(0.);

        // complex<double> m(0., 0.);
        // double m1 = 0; double m2 = 0; double m4 = 0;


        complex<double> m0(0., 0.);

        m[n] = 0.0;
        m1[n] = 0.0;
        m2[n] = 0.0;
        m4[n] = 0.0;

        for(int t = 0; t < NMCsteps; t++)
        {
            for(int c = 0; c < Nclusters; c++)
            {
                FlipAndBuild_from(rand() % N); //flip the spin of the site and build the cluster
            }

            complex<double> tm(0., 0.);
            for(int s = 0; s < q; s++)
            {
                tm += W[s] * double(M[s]); //compute the order parameter
            }

            tm/=N; //normalise the order parameter
            double tm1 = abs(tm); //compute the first moment of the order parameter
            double tm2 = tm1 * tm1; //compute the second moment of the order parameter

            // m += tm; m1 += tm1; m2 += tm2; m4 += tm2 * tm2; //compute the moments of the order parameter

            complex<double> m0conj = conj(W[Spins[0]]); //compute the order parameter for each site
            m0 += m0conj; //compute the order parameter for each site
            
            for(int j = 0; j < N; j++)
            {
                m0r[j] += m0conj * W[Spins[j]]; //compute the order parameter for each site
                mr[j] += W[Spins[j]]; //compute the order parameter for each site
            }

            m[n]+= tm; m1[n] += tm1, m2[n] += tm2, m4[n] += tm2 * tm2; //compute the magnetic moments
        }
        
        //normalise everything by number of MC steps
        m[n]/=NMCsteps; m1[n]/=NMCsteps; m2[n]/=NMCsteps; m4[n]/=NMCsteps;
        m0 /= NMCsteps; //normalise the order parameter for each site

        // cout << m << ' ' << m1 << ' ' << m2 << ' ' << m4 << endl ;

        for(int j = 0; j < N; j++)
        {
            m0r[j] /= NMCsteps;
            mr[j] /= NMCsteps;

        //print length of m0r for a sanity check
        // cout << m0r.size() << endl;

        //compute the correlation function: i think that this is computing the 'experimental' average
        for(int j = 0; j < N; j++)
        {
            Cr[j] += m0r[j] - m0 * mr[j]; //do not understand why there is a summation here: it could be for the experimental average
        }

        // cout << Cr.size() << endl;
        // cout << m1[n] << endl;

    }

    //  //print size of Cr for each j

    // for(int j = 0; j < N; j++)
    // {
    //     cout << Cr[j] << endl;
    // }

    for (int j = 0; j < N; j++)
    {
        Cr[j] /= Nbins;
        // std::cout << Cr[j].real() << endl;
    }


    //compute average magnetic moments for each experiment
    //wha we have, are magnetic moments of each site, for each individual experiment.
    //now, we want to average over these magnetic moments for each site, and basically obtain
    //the average magnetic moments for each experiment.
    //doing this for many experiments will let us compute the statistics, such as standard deviations from
    //each experiment etc.
    //basically, u will get the zeroth, first, second and forth magnetic moments for each temperature that
    //u run this programme for
    //and this will be computed for as many experiments as u wish - here, we choose Nbins = 10, so 10 exps!

    double m_avg = 0.0;
    double m1_avg = 0.0;
    double m2_avg = 0.0;
    double m4_avg = 0.0;

    // here, we sum over the magnetic moments at each site to get a running total
    for (int n = 0; n < Nbins; n++)
    {
        m_avg += m[n].real();
        m1_avg += m1[n];
        m2_avg += m2[n];
        m4_avg += m4[n];
    }

    //and here, we normalise it - so, we obtain one value for a magnetic moment per experiment!
    m_avg /= Nbins;
    m1_avg /= Nbins;
    m2_avg /= Nbins;
    m4_avg /= Nbins;

    // cout << m1_avg << endl;

    // ofstream file("correlation_0.25.txt");

    //     if (file.is_open()) {
    //     // Write vector values to the file

    //     for (int j = 0; j < N; j++)
    //     {
    //         file << Cr[j].real() << " ";
    //     }

    //     // Close the file
    //     file.close();

    //     cout << "Vector values saved to file successfully." << endl;
    // } else {
    //     cerr << "Failed to open file for writing." << endl;
    // }


    //not ready yet

    ofstream file("correlation_0.25.txt");

        if (file.is_open()) {
        // Write vector values to the file

        for (int j = 0; j < N; j++)
        {
            file << Cr[j].real() << " ";
        }

        // Close the file
        file.close();

    //     cout << "Magnetic moments saved to file successfully." << endl;
    // } else {
    //     cerr << "Failed to open file for writing." << endl;
    // }

    //use temperature to name the file
    // ofstream file("magnetic_moments_T15.txt");
    // ostringstream oss; // Create ostringstream object
    // oss << fixed << std::setprecision(2) << T; // Set precision to 1 and convert float to string
    // ofstream file("magnetic_moments_L" + to_string(L) + "_T"  + to_string(T) + ".txt");

    //     if (file.is_open()) {
    // //     // Write vector values to the file

    //     for (int j = 0; j < N; j++)
    //     {
    //         file << m_avg << " " << m1_avg << " " << m2_avg << " " << m4_avg << "\n";
    //     }

    //     // Close the file
    //     file.close();

    //     cout << "Magnetic moments saved to file successfully." << endl;
    // } else {
    //     cerr << "Failed to open file for writing." << endl;
    // }


    }
    }
}



















