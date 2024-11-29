#ifndef __WLDIAGRAM_HPP_CMC__
#define __WLDIAGRAM_HPP_CMC__
#include "Diagram.hpp"
#include "SquareLattice.hpp"
using namespace std;

class Operator
{
    public:
        int  nbef;
        bool create;

        double Ebef=std::nan("");
        double Eaft=std::nan("");

        int naft () const { return create ? nbef+1 : nbef-1; }

        friend std::ostream& operator<< (std::ostream& os, const Operator& op)
        {
            os << "  nbef = " << op.nbef << "\n";
            os << "  creat = " << op.create << "\n";
            return os;
        }
};

using Latt = SquareLattice;
using Iter = Diagram<Operator,Latt>::Iter;
const auto NULL_ITER = Node<Operator>::NULL_ITER;


class WLDiagram : public Diagram<Operator,Latt>
{
  public:
    typedef Diagram<Operator,Latt>::const_Iter const_Iter;
    typedef Diagram<Operator,Latt>             Diag;

    WLDiagram (int Lx, int Ly, int Lz, double Beta, bool xpbc, bool ypbc, bool zpbc, int n_init=0);

    std::tuple<Iter,Iter> insert2 (const Iter it_to, double time, bool create);

    void info_plot_lines (std::ofstream& ofs) const;
    void check           (int site) const;
};

std::tuple<Iter,Iter> WLDiagram::insert2 (const Iter it_to, double time, bool create)
{
    auto [it_bef, it_aft] = Diagram::insert2 (it_to, time);

    int n = it_to->z.nbef;
    if (create)
    {
        it_aft->z.nbef = n+1;
        it_bef->z.nbef = n;
        it_aft->z.create = false;
        it_bef->z.create = true;
    }
    else
    {
        assert (n != 0);
        it_aft->z.nbef = n-1;
        it_bef->z.nbef = n;
        it_aft->z.create = true;
        it_bef->z.create = false;
    }

#ifdef DEBUG_MODE
    check (site1, __func__);
    check (site2, __func__);
#endif
#ifdef DEBUG_TRACK_MODE
    __DEBUG_TRACK
#endif
    return {it_bef, it_aft};
}

WLDiagram :: WLDiagram (int Lx, int Ly, int Lz, double Beta, bool xpbc, bool ypbc, bool zpbc, int n_init)
: Diagram<Operator,Latt> (Lx, Ly, Lz, Beta, xpbc, ypbc, zpbc)
{
    for(int i = 0; i < Diag::SITEN; ++i)
    {
        Diag::it_end[i]->z.nbef = n_init;
    }
#ifdef DEBUG_MODE
    for(int site = 0; site < SITEN; ++site)
        check (site, __func__);
#endif
#ifdef DEBUG_TRACK_MODE
__DEBUG_TRACK
#endif
}

void WLDiagram :: info_plot_lines (std::ofstream& ofs) const
{
    auto symb = [] (auto it)
    {
        if (it->fix)
            return "x ";
        else if (!it->end and !it->has_conj())
            return "* ";
        else
            return ". ";
    };
    for(const auto& line : nodeList)
    {
        auto it = line.begin();
        ofs << symb(it) << it->site << " " << it->time << endl;
        for(it++; it != line.end(); it++)
        {
            ofs << symb(it) << it->site << " " << it->time << endl;
            auto it_bef = it;
            it_bef--;
            ofs << "- " << it->site << " " << it->time << " " << it_bef->time << " " << it->z.nbef << endl;
        }
    }
}

void WLDiagram :: check (int site) const
{
  std::ofstream ofs ("Check.info");
  if (!ofs) {
    std::cout << __FILE__ << " :: " << __func__ << " : cannot open file 'Check.error'\n";
    throw;
  }

  // Check End-nodes n match
  const_Iter it_1st = ++(nodeList[site].begin());
  if (it_1st->z.nbef != it_end[site]->z.nbef) {
    std::cout << "User check error\n";
    std::cout << "\n*** end-node n not match ***\n\n";
    ofs << "\n*** end-node n not match ***\n\n";
    ofs << "first node info:\n----\n" << it_1st->info(latt) << "\n";
    ofs << "end node info:\n----\n" << it_end[site]->info(latt) << "\n";
    ofs.close();
    throw;
  }

  for(const_Iter it = ++(nodeList[site].begin()); it != it_end[site]; ++it) {
    const_Iter it_aft = ++const_Iter(it);
    // Check creation/annihilation
    if (it != it_end[site]) {
      if (it->z.create) {
        if (it->z.nbef != (it_aft->z.nbef - 1)) {
          std::cout << "User check error\n";
          std::cout << "\n*** creation/annihilation error ***\n\n";
          ofs << "\n*** creation/annihilation error ***\n\n";
          ofs << "self info:\n----\n" << it->info(latt) << "\n";
          ofs << "aft info:\n----\n" << it_aft->info(latt) << "\n";
          ofs.close();
          throw;
        }
      }
      else {
        if (it->z.nbef != (it_aft->z.nbef + 1)) {
          std::cout << "User check error\n";
          std::cout << "\n*** creation/annihilation error ***\n\n";
          ofs << "\n*** creation/annihilation error ***\n\n";
          ofs << "self info:\n----\n" << it->info(latt) << "\n";
          ofs << "aft info:\n----\n" << it_aft->info(latt) << "\n";
          ofs.close();
          throw;
        }
      }
    }
    // Check n > 0
    if (it->z.nbef < 0) {
      std::cout << "User check error\n";
      std::cout << "\n*** n < 0 error ***\n\n";
      ofs << "\n*** n < 0 error ***\n\n";
      ofs << "self info:\n----\n" << it->info(latt) << "\n";
      ofs.close();
      throw;
    }
  }
  ofs.close();
}
#endif
