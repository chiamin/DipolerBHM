#ifndef __WORM_HPP_CMC__
#define __WORM_HPP_CMC__
#include <cassert>
#include "WLDiagram.hpp"
using namespace std;

struct Worm
{
    Iter   it;
    bool   up;
};

tuple<Worm,Iter> insert (WLDiagram& wld, Iter it, double time, bool create, bool up)
{
    auto [bef, aft] = wld.insert2 (it, time, create);

    Worm worm;
    worm.up = up;
    if (up)
    {
        worm.it = aft;
        bef->fix = true;
        return {worm, bef};
    }
    else
    {
        worm.it = bef;
        aft->fix = true;
        return {worm, aft};
    }
}

void walk (WLDiagram& wld, Worm& worm, double time)
{
    assert (time <= (++Iter(worm.it))->time);
    assert (time >= (--Iter(worm.it))->time);

    auto z = worm.it->z;
    Iter aft = ++Iter(worm.it);
    wld.remove (worm.it);
    worm.it = wld.insert (aft, time);
    worm.it->z = z;
}

void halt (WLDiagram& wld, Worm& worm)
{
    Iter aft = ++Iter(worm.it);
    Iter it = (worm.up ? aft : --Iter(worm.it));
    double time = it->time;

    auto z = worm.it->z;
    wld.remove (worm.it);
    worm.it = wld.insert (aft, time);
    worm.it->z = z;
}

tuple<Iter,Iter> hop (WLDiagram& wld, Worm& worm, const Iter it_to)
{
    double time = worm.it->time;
    bool create = (worm.up != worm.it->z.create);

    auto [it_bef, it_aft] = wld.insert2 (it_to, time, create);

    Iter it_from = worm.it;
    if (worm.up)
    {
        conj (worm.it, it_bef);
        worm.it = it_aft;
        return {it_from, it_bef};   // return the nodes for the link
    }
    else
    {
        conj (worm.it, it_aft);
        worm.it = it_bef;
        return {it_from, it_aft};   // return the nodes for the link
    }
}

void del_hop (WLDiagram& wld, Worm& worm)
{
    Iter it = (worm.up ? ++Iter(worm.it) : --Iter(worm.it));
    if (!it->has_conj())
    {
        cout << "Ahead is not a link" << endl;
        throw;
    }
    auto it_conj = it->conj;

    wld.remove (worm.it);
    wld.remove (it);
    worm.it = it_conj;
    worm.it->conj = NULL_ITER;
    worm.it->conj2 = NULL_ITER;
}

void cross_end (WLDiagram& wld, Worm& worm)
{
    int site = worm.it->site;
    Iter it_beg = wld.beg (site),
         it_end = wld.end (site);

    if (worm.up)
    {
        it_end->z = worm.it->z;
        wld.remove (worm.it);
        Iter it_to = ++Iter(it_beg);
        worm.it = wld.insert (it_to, 0.);
        worm.it->z = it_end->z;
    }
    else
    {
        Iter it_beg2 = ++Iter(worm.it);
        it_end->z = it_beg2->z;
        auto z = worm.it->z;
        wld.remove (worm.it);
        worm.it = wld.insert (it_end, wld.beta());
        worm.it->z = z;
    }
}

void cross_node (WLDiagram& wld, Worm& worm)
{
    Iter it = worm.it;
    Iter it_to;
    if (worm.up)
    {
        it++;
        it_to = ++Iter(it);
    }
    else
    {
        it--;
        it_to = it;
    }

    auto z = it->z;
    it->z = worm.it->z;
    wld.remove (worm.it);
    worm.it = wld.insert (it_to, it->time);
    worm.it->z = z;
}

void remove_worms (WLDiagram& wld, Worm& worm)
{
    Iter it_fix = (worm.up ? ++Iter(worm.it) : --Iter(worm.it));
    assert (it_fix->fix);
    wld.remove (worm.it);
    wld.remove (it_fix);
}
#endif
