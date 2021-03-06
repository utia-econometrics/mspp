#ifndef LPSOLVERS_H
#define LPSOLVERS_H

#include "linear.h"


class cplexlpsolver : public lpsolver
{
public:
    virtual void solve(const varinfo_list& vars,
            const linearfunction& objective,
            const constraint_list<linearconstraint>& constraints,
            std::vector<double>& sol,
            double& objvalue) const;
};


#endif // LPSOLVERS_H
