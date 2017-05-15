#include <iostream>
#include "linear.h"
#include "scenarios.h"
#include "lpsolvers.h"
#include <fstream>

using namespace std;

class csvlpsolver : public lpsolver
{
public:
    virtual void solve(const varinfo_list& vars,
            const linearfunction& objective,
            const constraint_list<linearconstraint>& constraints,
            std::vector<double>& x) const
    {
        using namespace std;
        ofstream f("problem.csv");
        f << "variables" << endl;
        for(int i=0;i<vars.size(); i++)
            f << "x" << i << ",";
        f << endl;
        f << endl;

        f << "objective function" << endl;
        for(int i=0;i<objective.coefs.size(); i++)
           f << objective.coefs[i] << ",";
        f << endl;

        f << "lower consttraints" << endl;
        for(int i=0;i<vars.size(); i++)
            if(vars[i].l == minf)
                f << ",";
            else
                f << vars[i].l << ",";
        f << endl;

        f << "upper consttraints" << endl;
        for(int i=0;i<vars.size(); i++)
            if(vars[i].h == inf)
                f << ",";
            else
                f << vars[i].h << ",";

        f << endl;

        for(int y=0; y<3; y++)
        {
            linearconstraint::type t = (linearconstraint::type) y;
            f << "consttraints ";
            switch(t)
            {
                case linearconstraint::eq:
                    f << "=";
                break;
                case linearconstraint::geq:
                    f << ">=";
                break;
                case linearconstraint::leq:
                    f << "<=";
                break;
            }
            f << endl;
            for(int i=0; i<constraints.size(); i++)
            {
                const linearconstraint& c = constraints[i];
                if(c.t == t)
                {
                    for(int j=0; j<c.lhs.size(); j++)
                        f << c.lhs[j] << ",";
                    f << c.rhs << endl;
                }
             }
        }
        throw "Test solver only produces a csv file";
    }
};


/**
 * @brief The testproblem class

http://homepages.cae.wisc.edu/~linderot/classes/ie495/lecture4.pdf

*/

class testproblem2: public linearproblem<double>
{

public:
    testproblem2() : linearproblem({2,2})
    {

    }

    virtual void f(
            const history<double>& xih,
            unsigned int k,
            linearfunction_ptr& f) const
    {
       f.reset(new linearfunction({1.0,1.0}));
    }

    virtual void constraints(
            const history<double>& h,
            unsigned int stage,
            varinfo_list_ptr& vars,
            constraint_list_ptr<linearconstraint>& constraints
            ) const
    {
        if(stage==0)
        {
            vars.reset(new varinfo_list);
            vars->push_back(varinfo(varinfo::Rplus));
            vars->push_back(varinfo(varinfo::Rplus));
        }
        else if(stage==1)
        {

            vars.reset(new varinfo_list);
            vars->push_back(varinfo(varinfo::Rplus));
            vars->push_back(varinfo(varinfo::Rplus));

            double omega=h[1];

            constraints.reset(new linearconstraint_list);
            constraints->push_back(linearconstraint({omega,1,1,0}, linearconstraint::geq, 7.0));
            constraints->push_back(linearconstraint({omega,1,0,1}, linearconstraint::geq, 4.0));
        }
        else if(stage==2)
        {

            vars.reset(new varinfo_list);
            vars->push_back(varinfo(varinfo::Rplus));
            vars->push_back(varinfo(varinfo::Rplus));

            double omega=h[2];

            constraints.reset(new linearconstraint_list);
            constraints->push_back(linearconstraint({0,0,omega,1,1,0}, linearconstraint::geq, 8.0));
            constraints->push_back(linearconstraint({0,0,omega,1,0,1}, linearconstraint::geq, 12.0));
        }
        else
            assert(1);

    }


};





class scenariolister : public scenariocallback<double>
{
public:
    virtual void callback(
                      const path& p,
                      const history<double>& h,
                      const probhistory& ph,
                      unsigned int k)
    {
        std::cout << k << " " << h[k] << " " << ph[k] << std::endl;
    }
};

/*
tbd

zrusit kacka u operatoru
zrusit setnode
zvazit slouceni treestructure a tree
zvazit spojeni scenare a probabiity
prochazeni stromu jednotne
*/

int main(int argc, char *argv[])
{
    const int numleaves = 5;
    treestructure_ptr fp(new nxtreestructure({1,numleaves,1}));
    generaltree<double>* g=new generaltree<double>(fp);
    tree_ptr<double> xp(g);
    treeprobs_ptr tp(new uniformtreeprobs(fp));
    (*g)(path{0},0)=1;
    for(unsigned int i=0; i<numleaves; i++)
    {
        g->setnode(path{0,i},1,(double)i / numleaves);
//        for(unsigned int j=0; j<numleaves; j++)
//            g->setnode(path{0,i,j},2,(double) i / numleaves + (double)j / numleaves);
    }
    scenario_ptr<double> sp(new scenario<double>(xp,tp));
    linearproblem_ptr<double> pp(new testproblem2());

//    scenariolister l;
//    s.foreachnode(&l);

    biglpsolver<double> b(pp,sp);
    cplexlpsolver cps;
    b.solve(cps);
}
