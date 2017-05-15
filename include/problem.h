#ifndef PROBLEM_H
#define PROBLEM_H

#include "commons.h"
#include "probability.h"
#include <limits>
#include <string>

static_assert(std::numeric_limits<double>::is_iec559, "IEEE 754 required");

const double inf = std::numeric_limits<double>::infinity();
const double minf = -inf;

struct varinfo
{
public:
    enum type { R, Rplus, Rminus };
    double l;
    double h;

    varinfo(double al=minf, double ah=inf)
      : l(al),h(ah)
    {}
    varinfo(type at)
    {
        switch(at)
        {
        case R: l=minf; h=inf; break;
        case Rplus: l=0; h=inf; break;
        case Rminus: l=minf; h=0; break;
        }
    }
};

typedef std::vector<varinfo> varinfo_list;
typedef std::shared_ptr<varinfo_list> varinfo_list_ptr;

template<typename C>
using constraint_list=std::vector<C>;

template<typename C>
using constraint_list_ptr = std::shared_ptr<constraint_list<C>>;

template<typename O>
using objective_ptr = std::shared_ptr<O>;

template<class O, class C, class Xi>
class problem : public object
{
    public:

        /** Default constructor */
        problem(const std::vector<unsigned int>& d) :
             fd(d)
        {
            fsumd=0;
            for(int i=0;i < fd.size();i++)
                fsumd += d[i];
        }

        unsigned int varsuntilstage(unsigned int k)
        {
            assert(k <= T());
            unsigned int s;
            for(int i=0; i<=k; i++)
                s+=fd[i];
            return s;
        }
        unsigned int T() { return fd.size(); }

        const unsigned int d(unsigned int k) const
                      { return fd[k]; }
        unsigned int sumd() const { return fsumd; }

        virtual void f(
                const history<Xi>& xih,
                unsigned int k,
                objective_ptr<O>& f) const = 0;
        virtual void constraints(
                const history<Xi>& xih,
                unsigned int k,
                varinfo_list_ptr& xs,
                constraint_list_ptr<C>& constraints
                ) const = 0;

    private:
        std::vector<unsigned int> fd;
        unsigned int fsumd;
};


#endif // PROBLEM_H
