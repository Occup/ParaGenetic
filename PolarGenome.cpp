#include "PolarGenome.hpp"

using namespace std;

extern unsigned int NUM_POLs;
random_device PolarGenome::randev;
mt19937_64 PolarGenome::rangen;
normal_distribution<double> PolarGenome::dist_phi {0, STD_PHI};
normal_distribution<double> PolarGenome::dist_theta {0, STD_THT};

PolarGenome::PolarGenome() :
    GAGenome(Init, Mutate, Compare), Polar()
{
    rangen.seed(randev());
    crossover(Cross);
    evaluator(Objective);
}

PolarGenome & PolarGenome::operator = (const PolarGenome & orig)
{
    if (& orig != this) copy(dynamic_cast<const GAGenome &>(orig));
    return * this;
}

GAGenome * PolarGenome::clone(CloneMethod flag) const
{
    if (flag == ATTRIBUTES) {
        return dynamic_cast<GAGenome *>(new PolarGenome());
    }
    return dynamic_cast<GAGenome *>(new PolarGenome(* this));
}

void PolarGenome::copy(const GAGenome & orig)
{ 
    GAGenome::copy(orig);
    const PolarGenome & genome = reinterpret_cast<const PolarGenome &>(orig);
    Phi_vec = genome.Phi_vec;
    Theta_vec = genome.Theta_vec;
    D_vec = genome.D_vec;
    CrossCos = genome.CrossCos;
    len = genome.len;
    Polar::score = genome.Polar::score;
    score_err = genome.score_err;
}

void PolarGenome::Init(GAGenome & orig)
{
    PolarGenome & genome = reinterpret_cast<PolarGenome &>(orig);
    double phi[] = {0.0000000e+00,   
                    1.0471976e+00,   
                    2.0943951e+00,   
                    3.1415927e+00,   
                    4.1887902e+00,   
                    5.2359878e+00,   
                    0.0000000e+00,   
                    1.5707963e+00,   
                    3.1415927e+00,   
                    4.7123890e+00,   
                    0.0000000e+00,
                    0.0000000e+00,   
                    1.0471976e+00,   
                    2.0943951e+00,   
                    3.1415927e+00,   
                    4.1887902e+00,   
                    5.2359878e+00,   
                    0.0000000e+00,   
                    1.5707963e+00,   
                    3.1415927e+00,   
                    4.7123890e+00,   
                    0.0000000e+00};
    double theta[] = {7.8539816e-01,
                      1.0471976e+00,
                      5.2359878e-01,
                      7.8539816e-01,
                      1.0471976e+00,
                      5.2359878e-01,
                      1.0471976e+00,
                      7.8539816e-01,
                      5.2359878e-01,
                      5.2359878e-01,
                      1.5707963e+00,
                      7.8539816e-01,
                      1.0471976e+00,
                      5.2359878e-01,
                      7.8539816e-01,
                      1.0471976e+00,
                      5.2359878e-01,
                      1.0471976e+00,
                      7.8539816e-01,
                      5.2359878e-01,
                      5.2359878e-01,
                      1.5707963e+00};
    for (int k = 0; k < NUM_POLs; ++k) {
        phi[k] += dist_phi(rangen);
        theta[k] += dist_theta(rangen);
    }
    genome.read_params(phi, theta, NUM_POLs);
}

int PolarGenome::Mutate(GAGenome & orig, double pMu)
{
    PolarGenome & genome = reinterpret_cast<PolarGenome &>(orig);
    if (pMu <= 0 || GARandomFloat() > pMu) {
        return 0;
    }
    rangen.seed(randev());
    for (int ind = 0; ind < genome.len; ++ind) {
        genome.Phi_vec[ind] += dist_phi(rangen);
        genome.Theta_vec[ind] += dist_theta(rangen);
    }
    genome.Euler2Vec();
    genome.Polar::score = numeric_limits<double>::min();
    return 1;
}

double PolarGenome::Compare(const GAGenome & orig0, const GAGenome & orig1)
{
    const PolarGenome & genome0 = reinterpret_cast<const PolarGenome &>(orig0);
    const PolarGenome & genome1 = reinterpret_cast<const PolarGenome &>(orig1);

    double difference = 0;
    for (int k = 0; k < NUM_POLs; ++k) {
        difference += fabs(genome0.Phi_vec[k] - genome1.Phi_vec[k]);
        difference += fabs(genome0.Theta_vec[k] - genome1.Theta_vec[k]);
    }

    return difference;
}

double PolarGenome::Objective(GAGenome & orig)
{
    PolarGenome & genome = reinterpret_cast<PolarGenome &>(orig);
    // cout << "calculating Error of genome..." << endl;
    // cout << "calculating following genome:" << endl
    //      << genome;
    return genome.ErrorEstimation();
}

int PolarGenome::Cross(const GAGenome & p0, const GAGenome & p1,
                       GAGenome * c0, GAGenome * c1)
{
    const PolarGenome & mom = reinterpret_cast<const PolarGenome &>(p0);
    const PolarGenome & dad = reinterpret_cast<const PolarGenome &>(p1);
    PolarGenome * sis = reinterpret_cast<PolarGenome *>(c0);
    PolarGenome * bro = reinterpret_cast<PolarGenome *>(c1);
    int nc = 0;
    const int cross = GARandomInt(1, NUM_POLs - 2);
    if (c0) {
        ++nc;
        sis->copy(dad);
        for (int k = cross; k < NUM_POLs; ++k) {
            sis->Phi_vec[k] = mom.Phi_vec[k];
            sis->Theta_vec[k] = mom.Theta_vec[k];
        }
        sis->Polar::score = numeric_limits<double>::min();
        sis->Euler2Vec();
    }

    if (c1) {
        ++nc;
        bro->copy(mom);
        for (int k = cross; k < NUM_POLs; ++k) {
            bro->Phi_vec[k] = dad.Phi_vec[k];
            bro->Theta_vec[k] = dad.Theta_vec[k];
        }
        bro->Polar::score = numeric_limits<double>::min();
        bro->Euler2Vec();
    }
    return nc;
}
