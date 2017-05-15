#ifndef SCENARIOS_H
#define SCENARIOS_H

#include "commons.h"
#include "probability.h"


class nxtreestructure: public treestructure
{
public:

    nxtreestructure(const std::vector<unsigned int>& ns)
           : fns(ns), fstageoffsets(ns.size()), fstagesizes(ns.size())
    {
        fnumnodes = 0;
        for(unsigned int i=0, s = 1; i<ns.size(); i++)
        {
            s *= ns[i];
            fstageoffsets[i] = fnumnodes;
            fstagesizes[i] = s;
            fnumnodes += s;
        }
    }

    virtual unsigned int numbranches(const path& p, unsigned int k) const
    {
        assert(k < horizon());
        return fns[k];
    }

    virtual unsigned int horizon() const
    {
        return fns.size();
    }

    virtual unsigned int numnodes(unsigned int k) const
    {
        assert(k < horizon());
        return fstagesizes[k];
    }

    virtual unsigned int stageoffset(unsigned int k) const
    {
        assert(k < horizon());

        return fstageoffsets[k];
    }

    virtual unsigned int relnodeindex(const path& p, unsigned int k) const
    {
        assert(p.size() > 0);
        assert(p.size() <= horizon());
        assert(k < horizon());

        unsigned int a=0;

        unsigned int n=1;

        for(int i=k; i>=0; i--)
        {
            a += p[i] * n;
            n*=fns[i];
        }

        assert(a < fnumnodes);

        return a;
    }
    virtual unsigned int totalnumnodes() const
    {
        return fnumnodes;
    }

private:
    std::vector<unsigned int> fns;
    std::vector<unsigned int> fstageoffsets;
    std::vector<unsigned int> fstagesizes;
    unsigned int fnumnodes;
};


class uniformtreeprobs: public treeprobs
{
public:
    /** Default constructor */
    uniformtreeprobs(const treestructure_ptr& ts) : treeprobs(ts) {}
    /** Default destructor */
    virtual double operator()(const path& p, unsigned int k) const
    {
        return 1.0 / ts().numbranches(p,k);
    }

};


#endif // SCENARIOS_H
