#ifndef ASIDE_H
#define ASIDE_H

template<class Xi>
struct atom
{
    Xi x;
    double p;
    atom(): x(0), p(0) {}
    atom(const Xi& ax, double ap): x(ax), p(ap) {}
};

template <class Xi>
using scenario = adapted<atom<Xi>>;


#endif // ASIDE_H
