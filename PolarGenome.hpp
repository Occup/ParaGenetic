#include <ga/ga.h>
#include <ga/garandom.h>
#include "Polar.hpp"
#include <limits>
#include <cmath>
#include <random>

#define STD_PHI 1.5
#define STD_THT 1.2



class PolarGenome : public GAGenome, public Polar
{
private:
    static std::mt19937_64 rangen;
    static std::random_device randev;
    static std::normal_distribution<double> dist_phi;
    static std::normal_distribution<double> dist_theta;
public:
    GADefineIdentity("PolarGenome", 233);
    static void Init(GAGenome&);
    static int Mutate(GAGenome&, double);
    static double Compare(const GAGenome &, const GAGenome &);
    static double Objective(GAGenome &);
    static int Cross(const GAGenome &, const GAGenome &,
                     GAGenome *, GAGenome *);

    
    PolarGenome();
    PolarGenome(const PolarGenome & orig) { copy(dynamic_cast<const GAGenome &>(orig)); };
    virtual ~PolarGenome() {};
    PolarGenome & operator = (const PolarGenome &);
    PolarGenome(const GAGenome & gag) { copy(gag); };
    // PolarGenome & operator = (const GAGenome &);
    virtual GAGenome * clone(CloneMethod flag = CONTENTS) const;
    virtual void copy(const GAGenome &);
    double operator () (int i, int j) { return i == 0 ? Phi_vec[j] : Theta_vec[j]; }
    friend std::ostream & operator << (std::ostream & os, PolarGenome & PG) { os << dynamic_cast<Polar &>(PG); return os; };
};
