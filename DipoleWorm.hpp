#ifndef __DIPOLEWORM_HPP_CMC__
#define __DIPOLEWORM_HPP_CMC__
#include <iostream>
#include "Worm.hpp"
using namespace std;

using Dir = LatticeSite::Direction;

struct DipoleWorm
{
    Worm   worm1, worm2;
    Dir    dir;
    bool   up () const { return worm1.up; }
};

DipoleWorm insert (WLDiagram& wld, int site, double time, bool create, bool up, Dir dir)
{
    DipoleWorm worms;
    worms.dir = dir;

    // Assume PBC
    int nbi = get_neighbor_index (dir, true);
    int site2 = wld.latt(site).nb.at(nbi);

    auto it1 = wld.find_aft (site, time);
    auto it2 = wld.find_aft (site2, time);
    auto [worm1, fix1] = insert (wld, it1, time, create, up);
    auto [worm2, fix2] = insert (wld, it2, time, !create, up);
    worms.worm1 = worm1;
    worms.worm2 = worm2;
    return worms;
}

tuple<Iter,Iter> get_dipole_nb_aft (const Iter it1, const Iter it2, bool up, bool go_positive, LatticeSite::Direction dir)
{
    assert (it1->site < it2->site);

    int nbi = get_neighbor_index (dir, go_positive);

    Iter it_nb1, it_nb2;
    if (go_positive)
    {
        it_nb1 = it2;
        it_nb2 = it2->assoc.at(nbi);
        if (up)
            it_nb1++;
    }
    else
    {
        it_nb1 = it1->assoc.at(nbi);
        it_nb2 = it1;
        if (up)
            it_nb2++;
    }
    return {it_nb1, it_nb2};
}

// worm1 is always at the left of worm2
//
// Ex: Hop to the right
//
//              |    |
//              |    |
//         x    x    |
//         |         |
//         |         |
//
// After hopping:
//
//        up                   down
//
//        |    |                |    |
//        |    |                |    |
//        x    |                o----o      --
//   o----o    x           o----o    x        }  same time
//   |    o----o           |    x    |      --
//   |         |           |         |
//   |         |           |         |
//
tuple<Iter,Iter,Iter,Iter> hop (WLDiagram& wld, DipoleWorm& worms, bool go_positive)
{
    auto [it_nb1, it_nb2] = get_dipole_nb_aft (worms.worm1.it, worms.worm2.it, worms.up(), go_positive, worms.dir);

    auto [it_from1, it_to1] = hop (wld, worms.worm1, it_nb1);
    auto [it_from2, it_to2] = hop (wld, worms.worm2, it_nb2);

    // Second link for the dipole hoppings
    conj2 (it_from1, it_from2);
    conj2 (it_to1, it_to2);
    return {it_from1, it_to1, it_from2, it_to2};
}

void del_hop_conj2 (WLDiagram& wld, Worm& worm)
{
    Iter it = (worm.up ? ++Iter(worm.it) : --Iter(worm.it));
    if (!it->has_conj2())
    {
        cout << "Ahead is not a link" << endl;
        throw;
    }
    auto it_conj = it->conj2;

    wld.remove (worm.it);
    wld.remove (it);
    worm.it = it_conj;
    worm.it->conj = NULL_ITER;
    worm.it->conj2 = NULL_ITER;
}

void delete_hop (WLDiagram& wld, DipoleWorm& worms)
{
    Iter it1 = worms.worm1.it,
         it2 = worms.worm2.it;
    if (worms.up())
    {
        it1++;
        it2++;
    }
    else
    {
        it1--;
        it2--;
    }
    if (it1->conj == it2)
    {
        del_hop_conj2 (wld, worms.worm1);
        del_hop_conj2 (wld, worms.worm2);
    }
    else
    {
        del_hop (wld, worms.worm1);
        del_hop (wld, worms.worm2);
    }
}
#endif
