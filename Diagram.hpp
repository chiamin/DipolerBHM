#ifndef DIAGRAM_HPP_CMC
#define DIAGRAM_HPP_CMC
#include <cmath>
#include <fstream>
#include <sstream>
#include <limits>
#include <typeinfo>
#include <vector>
#include <tuple>
#include "Node.hpp"
using namespace std;

template <typename Info, typename Latt>
class Diagram
{
    public:
        typedef typename      std::list<Node<Info>>::iterator       Iter;
        typedef typename      std::list<Node<Info>>::const_iterator const_Iter;
        const static typename std::list<Node<Info>>                 NULL_LIST;

        Diagram (int Lx, int Ly, int Lz, double Beta, bool xpbc, bool ypbc, bool zpbc);

        Iter                  insert   (const Iter goal, double t);
        std::tuple<Iter,Iter> insert2  (const Iter goal, double t);
        void                  remove   (const Iter aft);
        Iter                  find_aft (int site, double time);

        // Information
        double beta   () const         { return BETA; }
        int    size   () const         { return SITEN; }
        Iter   beg    (int site) const { return it_beg[site]; }
        Iter   end    (int site) const { return it_end[site]; }

        // Debug
        std::string  info            () const;
        std::string  info            (const Node<Info>& node) const;
        void         diag_info_plot  (std::string filename, std::string ext="diag") const;
        void         check           (int site, std::string funcname="unknown") const;
        void         check_neighbors (int site, std::string funcname) const;
        void         do_check        (int site, std::ofstream& ofs, const std::string& funcname) const;

                void info_plot_links  (std::ofstream& ofs) const;
                void info_plot_assoc  (std::ofstream& ofs) const;
        virtual void info_plot_line   (const const_Iter& it, std::ofstream& ofs) const;
        virtual void write            (std::string name);
        virtual void read             (std::string name);

        Latt latt;
        std::vector<std::list<Node<Info>>>  nodeList;

    protected:
        int    SITEN;
        double BETA;

        std::vector<Iter>                   it_beg, it_end;

        void link     (const Iter it, const Iter it_aft, int nbi);
        void unlink   (const Iter it, const Iter it_aft, int nbi);
};
template <typename Info, typename Latt>
const typename std::list<Node<Info>> Diagram<Info,Latt>::NULL_LIST;

template <typename Info, typename Latt>
Diagram<Info,Latt> :: Diagram (int Lx, int Ly, int Lz, double Beta, bool xpbc, bool ypbc, bool zpbc)
: SITEN (Lx*Ly*Lz)
, BETA (Beta)
, latt (Lx, Ly, Lz, xpbc, ypbc, zpbc)
{
    nodeList.resize (SITEN);
    it_beg.resize (SITEN);
    it_end.resize (SITEN);
    // Insert terminal nodeList
    for(int site = 0; site < SITEN; ++site)
    {
        Node<Info> node_head (site, 0., true);
        Node<Info> node_tail (site, BETA, true);
        nodeList.at(site).push_back (node_head);
        nodeList.at(site).push_back (node_tail);
        it_beg[site] = nodeList.at(site).begin();
        it_end[site] = it_beg[site];
        ++it_end[site];
    }
    // Links the nodeList
    Iter it, it_nb;
    for(int site = 0; site < SITEN; ++site)
    {
        for(int nbi : latt(site).nbis)
        {
            // find the neighbor site
            int neighbor = latt(site).nb.at(nbi);
            it = nodeList.at(site).begin();
            it_nb = nodeList[neighbor].begin();
            // link the head nodeList
            it->assoc.at(nbi) = it_nb;
            // link the tail nodeList
            (++it)->assoc.at(nbi) = (++it_nb);
        }
    }
#ifdef DEBUG_MODE
    for(int site = 0; site < SITEN; ++site)
        check (site, __func__);
#endif
}

template <typename Info, typename Latt>
typename Diagram<Info,Latt>::Iter Diagram<Info,Latt> :: insert (const Iter goal, double t)
{
    int site = goal->site;
    // Insert nodes:
    auto it = nodeList.at(site).insert (goal, Node<Info>(site, t));
    // Set associate:
    for(int nbi : latt(site).nbis)
    {
        link (it, goal, nbi);
    }
    return it;
}

template <typename Info, typename Latt>
std::tuple <typename Diagram<Info,Latt>::Iter, typename Diagram<Info,Latt>::Iter>
Diagram<Info,Latt> :: insert2 (const Iter goal, double t)
{
    int site = goal->site;
    // Insert nodes:
    auto bef = nodeList.at(site).insert (goal, Node<Info>(site, t));
    auto aft = nodeList.at(site).insert (goal, Node<Info>(site, t));

    // Set associate:
    for(int nbi : latt(site).nbis)
    {
        link (bef, goal, nbi);
        aft->assoc.at(nbi) = bef->assoc.at(nbi);
    }
    return {bef, aft};
}

template <typename Info, typename Latt>
void Diagram<Info,Latt> :: remove (const Iter it)
{
    Iter aft = ++Iter(it);
    int site = aft->site;

    // Unlink the node:
    for(int nbi : latt(site).nbis)
    {
        unlink (it, aft, nbi);
    }
    // Remove
    nodeList.at(site).erase (it);
}

template <typename Info, typename Latt>
inline typename Diagram<Info,Latt>::Iter Diagram<Info,Latt> :: find_aft (int site, double t)
{
    Iter it = nodeList.at(site).begin();
    while (it->time < t) ++it;
    return it;
}

template <typename Info, typename Latt>
void Diagram<Info,Latt> :: link (const Iter it, const Iter it_aft, int nbi)
{
    Iter it_bef = --Iter(it);
    Iter ass = it_bef->assoc.at(nbi);
    int opp = latt(it->site).opp.at(nbi);

    if (ass->end and ass->time == 0)
        ass++;

    while (ass->time < it->time)// or (ass->time == it->time and ass->assoc.at(opp) == it_aft))
    {
        ass->assoc.at(opp) = it;
        ass++;
    }
    while (ass->time == it->time and ass->assoc.at(opp) == it_aft)
    {
        ass->assoc.at(opp) = it;
        ass++;
    }
    it->assoc.at(nbi) = ass;
}

template <typename Info, typename Latt>
void Diagram<Info,Latt> :: unlink (const Iter it, const Iter it_aft, int nbi)
{
    Iter ass = it->assoc.at(nbi);
    ass--;
    int opp = latt(it->site).opp.at(nbi);
    while (ass->assoc.at(opp) == it)
    {
        ass->assoc.at(opp) = it_aft;
        ass--;
    }
}

// -------------------------------------------------------------
template <typename Info, typename Latt>
void Diagram<Info,Latt> :: info_plot_links (std::ofstream& ofs) const
{
    for(const auto& line : nodeList)
    {
        for(auto it = line.begin(); it != line.end(); it++)
        {
            if (it->has_conj())
            {
                auto it2 = it->conj;
                if (it->site < it2->site)
                    ofs << "<-> " << it->site << " " << it->time << " " << it2->site << " " << it->time << std::endl;
            }
        }
    }
}

template <typename Info, typename Latt>
void Diagram<Info,Latt> :: info_plot_assoc (std::ofstream& ofs) const
{
    for(const auto& line : nodeList)
    {
        for(auto it = line.begin(); it != line.end(); it++)
        {
            int site = it->site;
            for(int nbi : latt(site).nbis)
            {
                auto ass = it->assoc.at(nbi);
                ofs << "-> " << it->site << " " << it->time << " " << ass->site << " " << ass->time << std::endl;
            }
        }
    }
}

/*
template <typename Info, typename Latt>
int Diagram<Info,Latt> :: hopdir (const const_Iter& it) const
{
  if (it->_conji == -1) return 0;
  int dir, di, dj, L;
  int i = it->site, j = it->conj()->site;
  if (latt(i).xi() != latt(j).xi()) {
    // x direction
    dir = 1;
    di = latt(i).xi();
    dj = latt(j).xi();
    L = latt.Lx();
  }
  else if (latt(i).yi() != latt(j).yi()) {
    // y direction
    dir = 2;
    di = latt(i).yi();
    dj = latt(j).yi();
    L = latt.Ly();
  }
  else if (latt(i).zi() != latt(j).zi()) {
    // z direction
    dir = 3;
    di = latt(i).zi();
    dj = latt(j).zi();
    L = latt.Lz();
  }
  if (di < dj)
    dir *= -1;
  if ((di == 0 && dj == L) || (di == L && dj == 0))
    // crossing boundary
    dir *= -1;
  return dir;
}
*/
template <typename Info, typename Latt>
std::string Diagram<Info,Latt> :: info (const Node<Info>& node) const
{
  std::stringstream sstr;
  sstr << "location: " << &node << "\n"
       << "site = " << node._site << "\n"
       << "time = " << node._time << "\n";
       //<< "n bef = " << node.nbef() << "\n";
  for(int i = 0; i < latt(node._site).nbs(); ++i) {
    sstr << "assoc_memory[" << i << "]: " << &(node.assoc.at(i)) << "\n";
    sstr << "assoc_site [" << i << "] = " << node.assoc.at(i)->site << "\n";
    sstr << "assoc_time [" << i << "] = " << node.assoc.at(i)->time << "\n";
  }
  return sstr.str();
}

template <typename Info, typename Latt>
std::string Diagram<Info,Latt> :: info () const
{
  std::stringstream sstr;
  Iter it;
  sstr << " --- Information of Worldline Diagram ---\n\n";
  for(int is = 0; is < SITEN; ++is) {
    sstr << "site " << is << ":\n\n";
    for(it = nodeList[is].begin(); it != nodeList[is].end(); ++it) {
      sstr << info(*it) << "\n";
      }
    }
  sstr << " --- End of Worldline Information ---\n\n";
  return sstr.str();
}

template <typename Info, typename Latt>
void Diagram<Info,Latt> :: diag_info_plot (std::string filename, std::string ext) const
{
  std::ofstream ofs ((filename+".line."+ext).c_str());
  std::ofstream ofs2 ((filename+".hop."+ext).c_str());
  if (!ofs) {
    std::cout << __FILE__ << " :: " << __func__ << " : cannot open file '" << (filename+".line."+ext) << "'\n";
    exit(100);
  }
  if (!ofs2) {
    std::cout << __FILE__ << " :: " << __func__ << " : cannot open file '" << (filename+".hop."+ext) << "'\n";
    exit(100);
  }
  for(int i = 0; i < SITEN; ++i) {
    ofs << latt(i).x() << " " << latt(i).y() << " " << latt(i).z();
    for(const_Iter it = nodeList[i].begin(); it != nodeList[i].end(); ++it) {
      // Write line information
      info_plot_line (it, ofs);
      // Write hopping information
      //info_plot_hop (it, ofs2);
    }
    ofs << "\n";
  }
  ofs.close();
  ofs2.close();
}

template <typename Info, typename Latt>
void Diagram<Info,Latt> :: info_plot_line (const const_Iter& it, std::ofstream& ofs) const
{
  ofs << " " << it->time;
}
/*
template <typename Info, typename Latt>
void Diagram<Info,Latt> :: info_plot_hop (const const_Iter& it, std::ofstream& ofs) const
{
  int i = it->site;
  if (!(it->is_hopping())) {
    if (it != nodeList[i].begin() && it != it_end[i])
      // worms
      ofs << 0 << " " << latt(i).xi() << " " << latt(i).yi() << " " << latt(i).zi() << " " << it->time << "\n";
  }
  else if (it->site > it->conj()->site) {
    int hp = hopdir (it);
    int j = it->conj()->site;
    int dsite = hp > 0 ? 1 : -1;
    if (abs(hp) == 1) {
      // x direction
      if ((latt(i).xi() - latt(j).xi()) > 1)
        // cross boundary
        ofs << -1 << " " << latt(i).xi() << " " << latt(j).xi();
      else if ((latt(i).xi() - latt(j).xi()) < -1)
        // cross boundary
        ofs << -1 << " " << latt(j).xi() << " " << latt(i).xi();
      else
        ofs << 1 << " " << latt(i).xi() << " " << (latt(i).xi() - dsite);
      ofs << " " << latt(i).yi() << " " << latt(i).zi() << " " << it->time << "\n";
    }
    else if (abs(hp) == 2) {
      // y direction
      if ((latt(i).yi() - latt(j).yi()) > 1)
        ofs << -2 << " " << latt(i).yi() << " " << latt(j).yi();
      else if ((latt(i).yi() - latt(j).yi()) < -1)
        ofs << -2 << " " << latt(j).yi() << " " << latt(i).yi();
      else
        ofs << 2 << " " << latt(i).yi() << " " << (latt(i).yi() - dsite);
      ofs << " " << latt(i).xi() << " " << latt(i).zi() << " " << it->time << "\n";
    }
    else if (abs(hp) == 3) {
      // z direction
      if ((latt(i).zi() - latt(j).zi()) > 1)
        ofs << -3 << " " << latt(i).zi() << " " << latt(j).zi();
      else if ((latt(i).zi() - latt(j).zi()) < -1)
        ofs << -3 << " " << latt(j).zi() << " " << latt(i).zi();
      else
        ofs << 3 << " " << latt(i).zi() << " " << (latt(i).zi() - dsite);
      ofs << " " << latt(i).xi() << " " << latt(i).yi() << " " << it->time << "\n";
    }
  }
}
*/
template <typename Info, typename Latt>
void Diagram<Info,Latt> :: check_neighbors (int site, std::string funcname) const
{
  check (site, funcname);
  for(int nbi : latt(site).nbis) {
    int nbsite = latt.nb (site, nbi);
    check (nbsite, funcname);
  }
}

template <typename Info, typename Latt>
void Diagram<Info,Latt> :: write (std::string name)
{
  ofstream ofs2 ((name+".latt").c_str());
  latt.save (ofs2);
  // Rename the old bakcup-file if exist
  std::string fullname = name+".diagram";
  std::ifstream ifs (fullname.c_str());
  bool exist = ifs.good();
  ifs.close();
  std::string oldname = fullname + ".old";
  if (exist) rename (fullname.c_str(), oldname.c_str());
  std::ofstream ofs (fullname.c_str(), std::ios::binary);
  if (!ofs) {
    std::cout << "Error open file: DoubleDiagram: write: " << fullname << "\n";
    exit(100);
  }
  // Write parameters
  ofs.write((char*)&SITEN, sizeof(SITEN));
  // Write nodeList
  for(int i = 0; i < SITEN; ++i) {
    int size = nodeList[i].size();
    ofs.write((char*)&size, sizeof(size));
    for(Iter it = nodeList[i].begin(); it != nodeList[i].end(); ++it)
      ofs.write((char*)&(*it),sizeof(*it));
  }
  // Write associate
  for(int i = 0; i < SITEN; ++i)
    for(Iter it = nodeList[i].begin(); it != nodeList[i].end(); ++it) {
      for(int nbi : latt(i).nbis) {
        Iter it_ass = it->assoc.at(nbi);
        int nbsite = latt(i).nb.at(nbi);
        int dist = distance (nodeList[nbsite].begin(), it_ass);
        ofs.write((char*)&dist, sizeof(dist));
      }
    }
  ofs.close();
  // remove old backup file
  if (exist) std::remove (oldname.c_str());
}

template <typename Info, typename Latt>
void Diagram<Info,Latt> :: read (std::string name)
{
  ifstream ifs2 ((name+".latt").c_str());
  latt.load (ifs2);
  // Check the old backup file
  std::string fullname = name+".diagram";
  std::string oldname = fullname + ".old";
  std::ifstream ifs_old (oldname.c_str());
  bool old_exist = ifs_old.good();
  ifs_old.close();
  std::ifstream ifs;
  if (old_exist) ifs.open (oldname.c_str(), std::ios::binary);
  else ifs.open (fullname.c_str(), std::ios::binary);
  if (!ifs) {
    std::cout << "Error open file: DoubleDiagram: read: " << fullname << "\n";
    exit(100);
  }
  // Read parameters
  ifs.read((char*)&SITEN, sizeof(SITEN));
  nodeList.resize (SITEN);
  it_beg.resize (SITEN);
  it_end.resize (SITEN);
  // Read nodeList
  for(int i = 0; i < SITEN; ++i) {
    int size;  ifs.read((char*)&size, sizeof(size));
    nodeList[i].resize(size, Node<Info>());
    for(Iter it = nodeList[i].begin(); it != nodeList[i].end(); ++it)
      ifs.read((char*)&(*it),sizeof(*it));
  }
  // Read associate
  for(int i = 0; i < SITEN; ++i) {
    for(Iter it = nodeList[i].begin(); it != nodeList[i].end(); ++it) {
      for(int nbi: latt(i).nbis) {
        int nbsite = latt(i).nb.at(nbi);
        int dist;  ifs.read((char*)&dist, sizeof(dist));
        Iter it_ass = nodeList[nbsite].begin();
        advance (it_ass, dist);
        it->assoc.at(nbi) = it_ass;
      }
    }
    // Set it_beg[], it_end[]
    it_beg[i] = nodeList[i].begin();
    it_end[i] = --Iter(nodeList[i].end());
  }
  ifs.close();
}

template <typename Info, typename Latt>
void Diagram<Info,Latt> :: check (int site, std::string funcname) const
{
  std::ofstream ofs ("Check.info");
  if (!ofs) {
    std::cout << __FILE__ << " :: " << __func__ << " : cannot open file 'Check.error'\n";
    exit(100);
  }

  // Check end-nodes
  for(int nbi : latt(site).nbis)
    if (!(nodeList.at(site).begin()->assoc.at(nbi)->end()) || !(it_end[site]->assoc.at(nbi)->end())) {
      std::cout << "User check error\n";
      std::cout << "\n In function  " << funcname << "\n\n";
      std::cout << "\n*** End-node associate error ***\n\n";
      ofs << "\n In function  " << funcname << "\n\n";
      ofs << "\n*** End-node associate error ***\n\n";
      ofs << "self info:\n----\n" << nodeList.at(site).begin()->info(latt) << "\n";
      ofs.close();
      exit(100);
    }

  for(const_Iter it = ++(nodeList.at(site).begin()); it != it_end[site]; ++it) {
    // Check conjugate
    /*
    if (it->_conji != -1) {
      if (it->conj()->conj() != it) {
        std::cout << "User check error 1\n";
        std::cout << "\n In function  " << funcname << "\n\n";
        std::cout << "\n*** Conjugate error ***\n\n";
        ofs << "\n In function  " << funcname << "\n\n";
        ofs << "\n*** Conjugate error ***\n\n";
        ofs << "self info:\n----\n" << it->info(latt) << "\n";
        ofs << "conjugate info:\n----\n" << it->conj()->info(latt) << "\n";
        ofs.close();
        exit(100);
      }
      if (it->conj()->time != it->time ) {
        std::cout << "User check error 2\n";
        std::cout << "\n In function  " << funcname << "\n\n";
        std::cout << "\n*** Conjugate time error ***\n\n";
        ofs << "\n In function  " << funcname << "\n\n";
        ofs << "\n*** Conjugate time error ***\n\n";
        ofs << "self info:\n----\n" << it->info(latt) << "\n";
        ofs << "conjugate info:\n----\n" << it->conj()->info(latt) << "\n";
        ofs.close();
        exit(100);
      }
    }
    */
    // Check site
    if (it->site < 0) {
      std::cout << "User check error 9\n";
      std::cout << "\n In function  " << funcname << "\n\n";
      std::cout << "\n*** Site < 0 error ***\n\n";
      ofs << "\n In function  " << funcname << "\n\n";
      ofs << "\n*** Site < 0 error ***\n\n";
      ofs.close();
      exit(100);
    }
    const_Iter itpre = it;
    --itpre;
    if (it->site != itpre->site ) {
      std::cout << "User check error 3\n";
      std::cout << "\n In function  " << funcname << "\n\n";
      std::cout << "\n*** Site error ***\n\n";
      ofs << "\n In function  " << funcname << "\n\n";
      ofs << "\n*** Site error ***\n\n";
      ofs.close();
      exit(100);
    }
    // Check time
    if (it->time < itpre->time) {
      std::cout << "User check error 4\n";
      std::cout << "\n In function  " << funcname << "\n\n";
      std::cout << "\n*** Time error ***\n\n";
      ofs << "\n In function  " << funcname << "\n\n";
      ofs << "\n*** Time error ***\n\n";
      ofs << "self info:\n----\n" << it->info(latt) << "\n";
      ofs << "pre-node info:\n----\n" << itpre->info(latt) << "\n";
      ofs.close();
      exit(100);
    }
    // For associate:
    for(int nbi : latt(site).nbis) {
      Iter assoc = it->assoc.at(nbi);
      Iter assbef = assoc;
      --assbef;
      int opp = latt(site).opp(nbi);
      // Check recursive associate
      /*
      if (nbi != it->_conji && assoc->assoc.at(opp) == it) {
        std::cout << "User check error 5\n";
        std::cout << "\n In function  " << funcname << "\n\n";
        std::cout << "\n*** Recursive associate error ***\n\n";
        ofs << "\n In function  " << funcname << "\n\n";
        ofs << "\n*** Recursive associate error ***\n\n";
        ofs << "self info:\n----\n" << it->info(latt) << "\n";
        ofs << "associate info:\n----\n" << assoc->info(latt) << "\n";
        ofs.close();
        exit(100);
      }
      */
      // Check associate time
      if (assoc->time < it->time ) {
        std::cout << "User check error 6\n";
        std::cout << "\n In function  " << funcname << "\n\n";
        std::cout << "\n*** Assoc time error ***\n\n";
        ofs << "\n In function  " << funcname << "\n\n";
        ofs << "\n*** Assoc time error ***\n\n";
        ofs << "self info:\n----\n" << it->info(latt) << "\n";
        ofs << "associate info:\n----\n" << assoc->info(latt) << "\n";
        ofs.close();
        exit(100);
      }
      // Check associate site
      if (assoc->site != latt(site).nb.at(nbi) ) {
        std::cout << "User check error 7\n";
        std::cout << "\n In function  " << funcname << "\n\n";
        std::cout << "\n*** assoc site error ***\n\n";
        ofs << "\n In function  " << funcname << "\n\n";
        ofs << "\n*** assoc site error ***\n\n";
        ofs << "self info:\n----\n" << it->info(latt) << "\n";
        ofs << "associate-before info:\n----\n" << assbef->info(latt) << "\n";
        ofs.close();
        exit(100);
      }
      // Check associate-before time
      if (assbef->time > it->time ) {
        std::cout << "User check error 8\n";
        std::cout << "\n In function  " << funcname << "\n\n";
        std::cout << "\n*** assoc between error ***\n\n";
        ofs << "\n In function  " << funcname << "\n\n";
        ofs << "\n*** assoc between error ***\n\n";
        ofs << "self info:\n----\n" << it->info(latt) << "\n";
        ofs << "ass info:\n----\n" << assoc->info(latt) << "\n";
        ofs << "assbef info:\n----\n" << assbef->info(latt) << "\n";
        ofs.close();
        exit(100);
      }
    }
  }
  ofs.close();
}
#endif
