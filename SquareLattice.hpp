#ifndef __SQUARELATTICE_HPP_CMC__
#define __SQUARELATTICE_HPP_CMC__
#include "LatticeSite.hpp"

class SquareLattice
{
    public:
        SquareLattice (int lx, int ly, int lz, bool xpb, bool ypbc, bool zpbc);

        const LatticeSite& operator() (int i) const { return _site.at(i); }

        int  size () const { return _site.size(); }
        bool xpbc () const { return _xpbc; }
        bool ypbc () const { return _ypbc; }
        bool zpbc () const { return _zpbc; }
        int  lx   () const { return _lx; }
        int  ly   () const { return _ly; }
        int  lz   () const { return _lz; }

        int  get_index (int x, int y, int z) const;

        void save (std::ostream& s) const;
        void load (std::istream& s);

    private:
        void set_nbs     ();
        void set_opp_nbs ();

        int  _lx, _ly, _lz;
        bool _xpbc, _ypbc, _zpbc;
        int  _size;
        int  _wx, _wy, _wz; // weights of x, y, z from coordinate index to site index

        std::vector<LatticeSite> _site;
};

// Get neighboring index along a certain axis for positive or negative direction
int get_neighbor_index (LatticeSite::Direction dir,  bool pos)
{
    int nbi;
    if (dir == LatticeSite::Direction::X)
    {
        if (pos)
            nbi = 1;
        else
            nbi = 0;
    }
    else if (dir == LatticeSite::Direction::Y)
    {
        if (pos)
            nbi = 3;
        else
            nbi = 2;
    }
    else if (dir == LatticeSite::Direction::Z)
    {
        if (pos)
            nbi = 5;
        else
            nbi = 4;
    }
    return nbi;
}

SquareLattice::SquareLattice (int lx, int ly, int lz, bool xpbc, bool ypbc, bool zpbc)
: _lx (lx)
, _ly (ly)
, _lz (lz)
, _wx (1)
, _wy (lx)
, _wz (lx*ly)
, _xpbc (xpbc)
, _ypbc (ypbc)
, _zpbc (zpbc)
{
    assert(lx > 0 and ly > 0 and lz > 0);
    _site.resize (lx*ly*lz);
    this->set_nbs();
    this->set_opp_nbs();
}

int pbc_coor (int x, int L)
{
    while (x < 1)
        x += L;
    while (x > L)
        x -= L;
    return x;
}

int SquareLattice::get_index (int x, int y, int z) const
{
    x = pbc_coor (x, _lx)-1;
    y = pbc_coor (y, _ly)-1;
    z = pbc_coor (z, _lz)-1;
    return x + y*_wy + z*_wz;
}

// Set x, y, z, nb, and nbs for each site
void SquareLattice::set_nbs ()
{
    for(int x = 1; x <= _lx; x++)
        for(int y = 1; y <= _ly; y++)
            for(int z = 1; z <= _lz; z++)
            {
                int index = get_index (x, y, z);
                auto& site = _site[index];

                site.x = x;
                site.y = y;
                site.z = z;

                site.nb[0] = get_index (x-1, y, z);
                site.nb[1] = get_index (x+1, y, z);
                site.nb[2] = get_index (x, y-1, z);
                site.nb[3] = get_index (x, y+1, z);
                site.nb[4] = get_index (x, y, z-1);
                site.nb[5] = get_index (x, y, z+1);

                // x
                if (_lx == 2)
                {
                    if (x == 1)
                    {
                        site.nbis.push_back (1);
                        site.has_nbi[1] = true;
                    }
                    else
                    {
                        site.nbis.push_back (0);
                        site.has_nbi[0] = true;
                    }
                }
                else if (_lx != 1)
                {
                    if (x > 1 or _xpbc)
                    {
                        site.nbis.push_back (0);
                        site.has_nbi[0] = true;
                    }
                    if (x < _lx or _xpbc)
                    {
                        site.nbis.push_back (1);
                        site.has_nbi[1] = true;
                    }
                }
                // y
                if (_ly == 2)
                {
                    if (y == 1)
                    {
                        site.nbis.push_back (3);
                        site.has_nbi[3] = true;
                    }
                    else
                    {
                        site.nbis.push_back (2);
                        site.has_nbi[2] = true;
                    }
                }
                else if (_ly != 1)
                {
                    if (y > 1 or _ypbc)
                    {
                        site.nbis.push_back (2);
                        site.has_nbi[2] = true;
                    }
                    if (y < _ly or _ypbc)
                    {
                        site.nbis.push_back (3);
                        site.has_nbi[3] = true;
                    }
                }
                // z
                if (_lz == 2)
                {
                    if (z == 1)
                    {
                        site.nbis.push_back (5);
                        site.has_nbi[5] = true;
                    }
                    else
                    {
                        site.nbis.push_back (4);
                        site.has_nbi[4] = true;
                    }
                }
                else if (_lz != 1)
                {
                    if (z > 1 or _zpbc)
                    {
                        site.nbis.push_back (4);
                        site.has_nbi[4] = true;
                    }
                    if (z < _lz or _zpbc)
                    {
                        site.nbis.push_back (5);
                        site.has_nbi[5] = true;
                    }
                }
            }
}

void SquareLattice::set_opp_nbs ()
{
    for(int i = 0; i < this->size(); i++)
        for(int nbi : _site[i].nbis)
        {
            int nb = _site[i].nb[nbi];
            int opp_nb = 0;
            while (_site[nb].nb.at(opp_nb) != i)
                opp_nb++;
            _site[i].opp.at(nbi) = opp_nb;
        }
}

void SquareLattice::save (std::ostream& s) const
{
    s.write((char*)&_lx, sizeof(_lx));
    s.write((char*)&_ly, sizeof(_ly));
    s.write((char*)&_lz, sizeof(_lz));
    s.write((char*)&_xpbc, sizeof(_xpbc));
    s.write((char*)&_ypbc, sizeof(_ypbc));
    s.write((char*)&_wx, sizeof(_wx));
    s.write((char*)&_wy, sizeof(_wy));
    s.write((char*)&_wz, sizeof(_wz));
}

void SquareLattice::load (std::istream& s)
{
    s.read((char*)&_lx, sizeof(_lx));
    s.read((char*)&_ly, sizeof(_ly));
    s.read((char*)&_lz, sizeof(_lz));
    s.read((char*)&_xpbc, sizeof(_xpbc));
    s.read((char*)&_ypbc, sizeof(_ypbc));
    s.read((char*)&_wx, sizeof(_wx));
    s.read((char*)&_wy, sizeof(_wy));
    s.read((char*)&_wz, sizeof(_wz));
    _site.clear();
    _site.resize(_lx*_ly*_lz);
    this->set_nbs();
    this->set_opp_nbs();
}
#endif
