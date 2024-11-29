#ifndef __LATTICESITE_HPP_CMC__
#define __LATTICESITE_HPP_CMC__
#include <vector>

struct LatticeSite
{
    enum Direction { X=0, Y=1, Z=2 };

    int         x, y, z;        // coordinate index
    bool        s;              // sublattice index
    std::vector<int> nb;        // neighboring sites
    std::vector<int> opp;       // opposite neighboring sites
    std::vector<int> nbis;      // neighboring indices
    std::vector<bool> has_nbi;  // neighboring indices exist or not

    LatticeSite () : nb(6), opp(6), has_nbi(6,false) {}
};
#endif
