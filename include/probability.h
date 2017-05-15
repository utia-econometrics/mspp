#ifndef PROBABILITY_H
#define PROBABILITY_H

#include "commons.h"

typedef std::vector<unsigned int> path;

class treestructure : public object
{
public:
    treestructure() {}
    virtual ~treestructure() {}

    virtual unsigned int horizon() const = 0;
    virtual unsigned int totalnumnodes() const = 0;

    virtual unsigned int numbranches
       (const path& p, unsigned int k) const = 0;
    unsigned int numbranches(const path& p) const { return numbranches(p,p.size()-1); }

    virtual unsigned int numnodes(unsigned int k) const = 0;
    virtual unsigned int stageoffset(unsigned int k) const = 0;
    virtual unsigned int relnodeindex
      (const path& p, unsigned int k) const = 0;

    unsigned int nodeindex(const path& p, unsigned int k) const
    {
        return stageoffset(k)+relnodeindex(p,k);
    }

    unsigned int nodeindex(const path& p) const {return nodeindex(p,p.size()-1); }
};

typedef std::shared_ptr<treestructure> treestructure_ptr;

template <class T>
using history = std::vector<T>;

template<class T>
class tree
{
public:
    /** Default constructor */
    tree(const treestructure_ptr& ts) :  fts(ts) {}
    /** Default destructor */
    virtual ~tree() {}

    virtual T operator()(const path& p, unsigned int k) const = 0;

    const treestructure& ts() const { return *fts; }
    const treestructure_ptr& ts_ptr() const { return fts; }
    const unsigned int horizon() { return fts->horizon(); }

private:
    treestructure_ptr fts;
};

template <class T>
using tree_ptr = std::shared_ptr<tree<T>>;

template <class T>
class generaltree: public tree<T>
{
public:
    generaltree(const treestructure_ptr& ts) : tree<T>(ts), fdata(ts->totalnumnodes()) {}

    virtual T operator()(const path& p, unsigned int k) const
    {
        return fdata[this->ts().nodeindex(p,k)];
    }
    virtual T& operator()(const path& p, unsigned int k)
    {
        return fdata[this->ts().nodeindex(p,k)];
    }

    void setnode(const path& p, unsigned int k, const T& newval )
    {
        fdata[this->ts().nodeindex(p,k)] = newval;
    }
protected:
    std::vector<T> fdata;
};


typedef double probability;

class treeprobs: public tree<probability>
{
public:
    treeprobs(const treestructure_ptr& ts) : tree<probability>(ts) {}
    virtual probability uncprob(const path& pth, unsigned int k) const
    {
        probability p=1.0;
        for(unsigned int i=0; i<=k; i++)
            p*=(*this)(pth,i);
    }
};

using treeprobs_ptr = std::shared_ptr<treeprobs>;

using probhistory = history<probability>;

template <class Xi>
class scenario;

template <class Xi>
class scenariocallback
{
public:
    virtual void callback(
                   const scenario<Xi>& s,
                   const path& p,
                   const history<Xi>& h,
                   const probhistory& ph,
                   unsigned int k) = 0;
};

template <class Xi>
class scenario
{
public:
    scenario(const tree_ptr<Xi>& ax, const treeprobs_ptr& ap)
      : fp(ap), fx(ax)
    {
        assert(ax->ts_ptr().get()==ap->ts_ptr().get());
    }
    const treestructure& ts() const { return fx->ts();}
    const treeprobs& p() const {return *fp; }
    const tree<Xi>& x() const { return *fx;}
    unsigned int horizon() const { return fx->ts().horizon(); }
private:
    void doforeachnode(
             scenariocallback<Xi> *callee,
             path& pth, history<Xi>& xih, probhistory& ph, unsigned int k,
             unsigned int omit
            )
    {
        bool cb = k+1 < ts().horizon()-omit;
        unsigned int n=ts().numbranches(pth,k);
        for(unsigned int i=0; i<n; i++)
        {
            pth[k]=i;
            xih[k]=(*fx)(pth,k);
            ph[k]=(*fp)(pth,k);
            callee->callback(*this, pth,xih,ph,k);
            if(cb)
               doforeachnode(callee, pth, xih, ph, k+1,omit) ;
        }
    }
public:
    void foreachnode(scenariocallback<Xi> *callee, unsigned int omit = 0)
    {
         unsigned int h = ts().horizon();
         if(h>omit)
         {
             path pth(h);
             history<Xi> xih(h);
             probhistory ph(h);
             doforeachnode(callee, pth, xih ,ph, 0, omit);
         }
    }

private:
    treeprobs_ptr fp;
    tree_ptr<Xi> fx;
};

template<class Xi>
using scenario_ptr=std::shared_ptr<scenario<Xi>>;

#endif // PROBABILITY_H
