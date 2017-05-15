#ifndef LINEAR_H
#define LINEAR_H

#include "problem.h"
#include "probability.h"
#include <sstream>
#include <iostream>

class linearfunction: public object
{
public:
    linearfunction(const std::vector<double>& acoefs)
       : coefs(acoefs) {}

    linearfunction(unsigned int dim=0)
             : coefs(dim,0.0) {}
    std::vector<double> coefs;
    unsigned int dim() const { return coefs.size(); }
};


typedef std::shared_ptr<linearfunction> linearfunction_ptr;

class linearconstraint
{
public:
    enum type {eq, geq, leq};
    linearconstraint(const std::vector<double>& alhs,
                       type at, double arhs) :
         lhs(alhs), t(at), rhs(arhs) {}

    linearconstraint(unsigned int lhsdim, type at = eq)
        : lhs(lhsdim,0), t(at), rhs(0.0) {}
    std::vector<double> lhs;
    type t;
    double rhs;
    unsigned int dim() { return lhs.size(); }
private:
};

typedef std::vector<linearconstraint> linearconstraint_list;
typedef std::shared_ptr<linearconstraint_list> linearconstraint_list_ptr;

template<class Xi>
class linearproblem: public problem<linearfunction,linearconstraint,Xi>
{
    public:
        linearproblem(const std::vector<unsigned int>& soldim) :
           problem<linearfunction,linearconstraint,Xi>(soldim) {}
        virtual void gamma(
                const scenario<Xi>& s,
                const path& pth,
                const history<Xi>& xih,
                unsigned int k,
                linearfunction_ptr& r) const
        {
            assert(k<pth.size());
            assert(k<xih.size());

            if(s.horizon() < k+2)
                throw "Insufficient horizon of a scenarion in call of linearproblem::gamma";

            path np(pth);
            if(np.size() == k+1)
                np.resize(k+1);

            history<Xi> nh(xih);
            if(nh.size() == k+1)
                nh.resize(k+1);

            linearfunction* lfp = new linearfunction(this->d(k));
            r.reset(lfp);

            unsigned int n=s.ts().numbranches(np,k+1);

            std::cout << "k: " << k << " n: "             << n << std::endl;

            for(np[k+1]=0; np[k+1]<n; np[k+1]++)
            {
                nh[k+1] = s.x()(np,k+1);
                probability pr = s.p()(np,k+1);
                linearfunction_ptr fp;
                this->f(nh,k+1,fp);
                assert(fp->coefs.size()==this->d(k));
std::cout << "  pr: " << pr << " c1: " << fp->coefs[0] << std::endl;
                for(int i=0; i<this->d(k); i++)
                    r->coefs[i] += pr*(*fp).coefs[i];
            }
            assert(r.get()!=nullptr);
        }
};

template<class Xi>
using linearproblem_ptr = std::shared_ptr<linearproblem<Xi>>;

class lpsolver
{
public:
    virtual void solve(const varinfo_list& vars,
            const linearfunction& objective,
            const constraint_list<linearconstraint>& constraints,
            std::vector<double>& sol,
            double& objvalue) const = 0;
};

template<class Xi>
class biglpsolver : public scenariocallback<Xi>
{
public:
    biglpsolver(linearproblem_ptr<Xi>& pp, const scenario_ptr<Xi>& sp)
     : fp(pp), fs(sp)
    {
        if(fs->ts().numbranches(path{0},0) != 1)
            throw "scemarop must have a single root";
        if(fp->T()>fs->horizon())
            throw "infufficicient scenario horizon for a problem";
    }

private:
    /// state variables
    linearfunction fobj;
    varinfo_list fvars;
    linearconstraint_list fconstraints;
    unsigned int fdim;
public:
    virtual void callback(
                      const scenario<Xi>& s,
                      const path& pth,
                      const history<double>& h,
                      const probhistory& ph,
                      unsigned int k)
    {
        varinfo_list_ptr vars;
        constraint_list_ptr<linearconstraint> constraints;

        linearfunction_ptr objective;
        fp->gamma(s,pth,h,k,objective);

        fp->constraints(h,k,vars,constraints);

        assert(vars != nullptr);
        assert(objective != nullptr);

        unsigned int thisstagedim = fp->d(k);
        assert(vars->size()==thisstagedim);
        assert(objective->dim()==thisstagedim);

        std::vector<unsigned int> ofs(k+1);

        unsigned int o = 0;
        for(int i=0;i<=k;i++)
        {
            ofs[i] = o+fs->ts().relnodeindex(pth,i) * fp->d(i);
            o += fs->ts().numnodes(i) * fp->d(i);
        }

        probability up = fs->p().uncprob(pth,k);

        for(unsigned int i=0; i<thisstagedim; i++)
        {
            assert(fobj.coefs[ofs[k]+i]==0.0);
            fobj.coefs[ofs[k]+i] = up * objective->coefs[i];

            varinfo& v=fvars[ofs[k]+i]=(*vars)[i];

        }
        if(constraints)
            for(unsigned int j=0; j<constraints->size(); j++)
            {
                linearconstraint& s = (*constraints)[j];

                fconstraints.push_back(linearconstraint(fdim));
                linearconstraint& d = fconstraints[fconstraints.size()-1];

                unsigned int src=0;
                for(unsigned int i=0; i<=k; i++)
                {
                    unsigned int dst=ofs[i];
                    for(unsigned int r=0; r<fp->d(i); r++)
                        d.lhs[dst++] = s.lhs[src++];
                    d.rhs = s.rhs;
                    d.t = s.t;
                }
//std::cout << "src " << src << " cs " << s.lhs.size() << std::endl;
                assert(src==s.lhs.size());
            }
    }

    void solve(const lpsolver& lps)
    {
        fdim = 0;
        for(int i=0; i<fs->ts().horizon(); i++)
           fdim += fs->ts().numnodes(i) * fp->d(i);

        fobj= linearfunction(fdim);
        fvars.resize(fdim);
        fconstraints.clear();

        fs->foreachnode(this, fs->horizon()-fp->T());

        std::vector<double> x(fdim);
        double ov;
        lps.solve(fvars,fobj,fconstraints,x,ov);
        for(int i=0; i<x.size(); i++)
            std::cout << "x" << i << "=" << x[i] << std::endl;
        std::cout << "optimal:" << ov << std::endl;
    }

private:
    linearproblem_ptr<Xi> fp;
    scenario_ptr<Xi> fs;
};

#endif // LINEAR_H
