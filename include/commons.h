#ifndef COMMONS_H
#define COMMONS_H

#include <vector>
#include <memory>

// stages index from zero to T-1

typedef std::vector<unsigned int> pathindex;

typedef std::vector<double> stagesolution;


class objective
{
};

struct range
{
public:
    enum type {infinf,linf,infh,lh};
    double l;
    double h;
};

class constraint
{
public:
    enum type {equal, lessequal};
};

class filtration
{
public:
    filtration(unsigned int T): fT(T) {};
    virtual ~filtration();
    virtual void numatoms(const pathindex& aindex) const = 0;
    virtual unsigned int numatoms() const = 0;
    virtual unsigned int absindex(const pathindex& aindex) const = 0;
private:
    unsigned int fT;
};

typedef std::shared_ptr<filtration> filtrationptr;

template<class Xi>
class scenario
{
    public:
        /** Default constructor */
        scenario(const filtrationptr& f) : ff(f), data(f->numatoms()) {}
        /** Default destructor */
        virtual ~scenario() {};

        virtual Xi& operator()(const pathindex& aindex);

        virtual std::vector<Xi&> operator[](const pathindex& aindex);

        const filtration& getfiltration() { return *ff; }
    protected:

    private:
        std::vector<Xi> data;
        filtrationptr ff; //!
};



template<class O, class C, class Xi>
class problem
{
    public:
        class stageinfo
        {
            O* fO;
            std::vector<C*> fC;
            std::vector<range> fvars;
        public:
            stageinfo(): fO(0) {}
            ~stageinfo()
            {
                delete fO;
                for(int i=0; i< fC.size(); i++)
                    delete fC[i];
            }
            void set_O(O* o) { if(o) delete o; fO = o; }
            void add_C(C* c) { fC.push_back(c); }
            void add_var(const range& r) { fvars.push_back(r); }

        };


        /** Default constructor */
        problem(std::vector<unsigned int>& soldims);
        /** Default destructor */
        virtual ~problem();

        virtual const std::vector<unsigned int>& getsoldims() { return fsoldims; }
        virtual void getstage(const std::vector<Xi>& path, stageinfo& result) const = 0;

    protected:

    private:
        std::vector<unsigned int> fsoldims;
};




template<class X>
class distribution
{
    public:
        /** Default constructor */
        distribution(unsigned int adim): dim(adim)
        {}

        distribution(const distribution& adist): dim(adist.dim) {};
        /** Default destructor */
        virtual ~distribution() = 0;

    protected:

    private:
        unsigned int dim;
};


#endif // COMMONS_H
