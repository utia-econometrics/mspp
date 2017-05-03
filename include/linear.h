#ifndef LINEAR_H
#define LINEAR_H

#include "commons.h"

class linearobjective: public objective
{
        linearobjective(const std::vector<double>& coefs) : fcoefs(coefs) {}
        double getcoef(unsigned int i) { return fcoefs[i]; }
    private:
        std::vector<double> fcoefs;
};

class linearconstraint: public constraint
{
    linearconstraint(const std::vector<double>& lhs, double rhs, constraint::type t);
};

template<class Xi>
class linearproblem: public problem<linearobjective,linearconstraint,Xi>
{
    public:
        linearproblem(std::vector<unsigned int>& soldim) :
           problem<linearobjective,linearconstraint,Xi>(soldim) {}
        virtual ~linearproblem();
    protected:

    private:
};

struct biglpresult
{
    scenario<stagesolution> solution;
    double optimalvalue;
};

template<class Xi>
class biglpsolver
{
public:
    static void solve(const linearproblem<Xi>& p, const scenario<Xi>& s,
                      biglpresult& r);
};


#endif // LINEAR_H
