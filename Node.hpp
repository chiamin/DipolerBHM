#ifndef Node_HPP
#define Node_HPP
#include <list>
#include <sstream>
#include <cstdlib>

template <typename Temp>
struct Node
{
    typedef typename std::list<Node>::iterator Iter;
    const static typename std::list<Node<Temp>>::iterator NULL_ITER;

    Node () {}
    Node (int s, double t, bool e=false) : site(s), time(t), end(e), assoc(6,NULL_ITER) {}

    double  time;
    int     site;
    bool    end;               // indicate whether this is an end-node
    bool    fix=false;
    Temp    z;
    std::vector<Iter> assoc;          // link to the neighbor nodes just upper or at the same time

    // Conj
    Iter conj=NULL_ITER;
    bool has_conj () const { return this->conj != NULL_ITER; }
    friend void conj (Iter& it1, Iter& it2) { it1->conj = it2; it2->conj = it1; }
    // Conj for dipole
    Iter conj2=NULL_ITER;
    bool has_conj2 () const { return this->conj2 != NULL_ITER; }
    friend void conj2 (Iter& it1, Iter& it2) { it1->conj2 = it2; it2->conj2 = it1; }

    // Info
    template <typename Latt>
    const std::string info (const Latt& latt) const;
    const std::string info () const;
};
template <typename Temp>
const typename std::list<Node<Temp>>::iterator Node<Temp>::NULL_ITER;

template <typename Temp>
template <typename Latt>
const std::string Node<Temp> :: info (const Latt& latt) const
{
    std::stringstream sstr;
    sstr << "self:\n";
    sstr << "  x, y, z = " << latt(this->site).x << ", " << latt(this->site).y << ", " << latt(this->site).z << "\n";
    sstr << this->info();
    for(int i : latt(this->site).nbis)
    {
        sstr << "assoc " << i << ":\n";
        sstr << assoc.at(i)->info();
    }
    return sstr.str();
}

template <typename Temp>
const std::string Node<Temp> :: info () const
{
    std::stringstream sstr;
    sstr << "  address = " << this << "\n";
    sstr << "  site = " << this->site << "\n";
    sstr << "  time = " << this->time << "\n";
    if (fix)
        sstr << "  fix = true" << "\n";
    if (this->conj == NULL_ITER)
        sstr << "  conj = NULL" << "\n";
    else
        sstr << "  conj = " << &*conj << "\n";
    sstr << this->z << "\n";
    return sstr.str();
}
#endif
