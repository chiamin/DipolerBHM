#ifndef __LIVEDIPOLEWORM_HPP_CMC__
#define __LIVEDIPOLEWORM_HPP_CMC__
#include <random>
#include "DipoleWorm.hpp"
#include "ParaBox.hpp"
#include "ProbTable.hpp"

int tmp_i = 0;

enum Status { UNKNOWN, FREE, HALTED, BOUNCED, REMOVED, PASS_NEIGHBOR, CROSS_NODE, CROSS_END };

ostream& operator<< (ostream& os, Status status)
{
    if (status == UNKNOWN)
        os << "UNKNOWN";
    else if (status == FREE)
        os << "FREE";
    else if (status == HALTED)
        os << "HALTED";
    else if (status == BOUNCED)
        os << "BOUNCED";
    else if (status == REMOVED)
        os << "REMOVED";
    else if (status == PASS_NEIGHBOR)
        os << "PASS_NEIGHBOR";
    else if (status == CROSS_NODE)
        os << "CROSS_NODE";
    else if (status == CROSS_END)
        os << "CROSS_END";
    return os;
}

struct Para
{
    double UHALF, t, EOFFSET;
    std::vector<double> MU;
    int nMIN=0;
    int nMAX=std::numeric_limits<int>::max();

    // For measurement
    int N_hop;
};

double onsite_energy (int site, int n, const Para& para)
{
    return (para.UHALF * (n-1) - para.MU.at(site)) * n;
}

void update_onsite_energy (Iter& it, const Para& para)
{
    int site = it->site;
    it->z.Ebef = onsite_energy (site, it->z.nbef, para);
    it->z.Eaft = onsite_energy (site, it->z.naft(), para);
}

void check_onsite_energy (const WLDiagram& wld, const Para& para)
{
    for(const auto& line : wld.nodeList)
    {
        for(auto it = ++line.begin(); it != line.end(); it++)
        {
            int site = it->site;
            auto Ebef = onsite_energy (site, it->z.nbef, para);
            if (it->z.Ebef != Ebef)
                throw;
            assert (it->z.Ebef == Ebef);

            if (it != wld.end(site))
            {
                auto Eaft = onsite_energy (site, it->z.naft(), para);
                if (it->z.Eaft != Eaft)
                    throw;
                assert (it->z.Eaft == Eaft);
            }
        }
    }
}

inline tuple<Iter,Iter> get_blocks (const DipoleWorm& worms)
{
    Iter it_block1 = worms.worm1.it,
         it_block2 = worms.worm2.it;
    if (worms.up())
    {
        it_block1++;
        it_block2++;
    }
    else
    {
        it_block1--;
        it_block2--;
    }
    return {it_block1, it_block2};
}

Status liveWalk (WLDiagram& wld, DipoleWorm& worms, double dtime)
{
    assert (worms.worm1.it->time == worms.worm2.it->time);

    // Find block
    auto [it_block1, it_block2] = get_blocks (worms);
    double dt_block1 = std::abs (it_block1->time - worms.worm1.it->time),
           dt_block2 = std::abs (it_block2->time - worms.worm2.it->time);

    // No block
    if (dtime < dt_block1 and dtime < dt_block2)
    {
        // Free walk
        double time = worms.worm1.it->time + (worms.up() ? dtime : -dtime);
        walk (wld, worms.worm1, time);
        walk (wld, worms.worm2, time);

#ifdef PRINT_DETAIL
cout << "FREE" << endl;
#endif

        return FREE;
    }
    // Two blocks at the same time
    else if (dtime > dt_block1 and dtime > dt_block2 and it_block1->time == it_block2->time)
    {
        if (it_block1->end and it_block2->end)
        {
            // Cross end
            cross_end (wld, worms.worm1);
            cross_end (wld, worms.worm2);
#ifdef PRINT_DETAIL
cout << "CROSS_END" << endl;
#endif
            return CROSS_END;
        }
        else if (it_block1->fix and it_block2->fix)
        {
            // Remove worms
            remove_worms (wld, worms.worm1);
            remove_worms (wld, worms.worm2);
#ifdef PRINT_DETAIL
cout << "REMOVED" << endl;
#endif
            return REMOVED;
        }
        // Both blocks can be crossed
        else if (it_block1->z.create == worms.worm1.it->z.create and
                 it_block2->z.create == worms.worm2.it->z.create)
        {
            // Cross two nodes
            cross_node (wld, worms.worm1);
            cross_node (wld, worms.worm2);
#ifdef PRINT_DETAIL
cout << "CROSS_NODE" << endl;
#endif
            return CROSS_NODE;
        }
        // Cannot cross
        else
        {
            // Halt
            halt (wld, worms.worm1);
            halt (wld, worms.worm2);
#ifdef PRINT_DETAIL
cout << "HALTED" << endl;
#endif
            return HALTED;
        }
    }
    // One block
    else
    {
        // Block 1
        if (dt_block1 < dt_block2)
        {
            assert (dtime > dt_block1);
            walk (wld, worms.worm2, it_block1->time);
            if (it_block1->z.create == worms.worm1.it->z.create)
            {
                // Cross node
                cross_node (wld, worms.worm1);
#ifdef PRINT_DETAIL
cout << "CROSS_NODE 1" << endl;
#endif
                return CROSS_NODE;
            }
            else
            {
                // Halt
                halt (wld, worms.worm1);
                // Bounce
                worms.worm1.up = !worms.worm1.up;
                worms.worm2.up = !worms.worm2.up;
#ifdef PRINT_DETAIL
cout << "BOUNCED 1" << endl;
#endif
                return BOUNCED;
            }
        }
        // Block 2
        else
        {
            assert (dtime > dt_block2);
            walk (wld, worms.worm1, it_block2->time);
            if (it_block2->z.create == worms.worm2.it->z.create)
            {
                // Cross node
                cross_node (wld, worms.worm2);
#ifdef PRINT_DETAIL
cout << "CROSS_NODE 2" << endl;
#endif
                return CROSS_NODE;
            }
            else
            {
                // Halt
                halt (wld, worms.worm2);
                // Bounce
                worms.worm1.up = !worms.worm1.up;
                worms.worm2.up = !worms.worm2.up;
#ifdef PRINT_DETAIL
cout << "BOUNCED 2" << endl;
#endif
                return BOUNCED;
            }
        }
    }
    cout << "Unknown status" << endl;
    throw;
    return UNKNOWN;
}

bool check_nlimit (int n, bool creat, int nmax, int nmin)
{
    if (n == nmax and creat)
        return false;
    else if (n == nmin and !creat)
        return false;
    return true;
}

tuple<double,double> shift_energy (double E1, double E2, double EOFFSET)
{
    double Esmall = (E1 <= E2 ? E1 : E2);
    E1 += -Esmall + EOFFSET;
    E2 += -Esmall + EOFFSET;
    return {E1, E2};
}

tuple<Iter,Iter> find_conj (const Iter it1, const Iter it2, bool up)
{
    Iter conj1 = it1->conj,
         conj2 = it2->conj;
    if (conj1 == it2 and conj2 == it1)
    {
        // Search for the next one
        Iter conj_next1 = conj1,
             conj_next2 = conj2;
        if (up)
        {
            conj_next1++;
            conj_next2++;
        }
        else
        {
            conj_next1--;
            conj_next2--;
        }

        if (conj_next1->time == conj1->time)
        {
            if (conj_next2->time == conj1->time)
            {
                cout << "Ambigous conjugate node" << endl;
                throw;
            }
        }
    }
}

tuple<Iter,Iter> find_it_to (const Iter it1, const Iter it2)
{
    if (it1->conj == it2)
        return {it1->conj2, it2->conj2};
    else
        return {it1->conj, it2->conj};
}

void act (WLDiagram& wld, DipoleWorm& worms, std::mt19937& rand, Status status, double Ewalk, Para& para)
{
    std::uniform_real_distribution<double> uniform_dis (0., 1.);
    ProbTable table;

    int site1 = worms.worm1.it->site,
        site2 = worms.worm2.it->site;
    // 1. Free
    if (status == FREE)
    {
        // 1.1 Bounce
        table.add (Ewalk);
        // 1.2 Insert hopping
        for(bool go_positive : {false, true})
        {
            int nbi = get_neighbor_index (worms.dir, go_positive);
            // Check if the neighbor exsit for both worms
            if (wld.latt(site1).has_nbi.at(nbi)
            and wld.latt(site2).has_nbi.at(nbi))
            {
                auto [it_nb1, it_nb2] = get_dipole_nb_aft (worms.worm1.it, worms.worm2.it, worms.up(), go_positive, worms.dir);
                int n1 = it_nb1->z.nbef,
                    n2 = it_nb2->z.nbef;
                bool creat1 = (worms.worm1.up != worms.worm1.it->z.create),
                     creat2 = (worms.worm2.up != worms.worm2.it->z.create);
                // Check if particle number can be created or annihilated
                if (check_nlimit (n1, creat1, para.nMAX, para.nMIN)
                and check_nlimit (n2, creat2, para.nMAX, para.nMIN))
                {
                    double weight1 = para.t * (creat1 ? n1+1 : n1),
                           weight2 = para.t * (creat2 ? n2+1 : n2);
                    double weight = weight1 * weight2;
                    table.add (weight, go_positive);
                }
            }
        }
        // 1.3 choose
#ifdef LOCAL_OPTIMAL
        table.locally_optimal();
#endif
        int choice = table.get_choice (uniform_dis(rand));
        if (choice == 0)
        {
#ifdef PRINT_DETAIL
cout << "Free bounce" << endl;
#endif
            // Bounce
            worms.worm1.up = !worms.worm1.up;
            worms.worm2.up = !worms.worm2.up;
        }
        else
        {
#ifdef PRINT_DETAIL
cout << "Hop" << endl;
#endif
            // Hop
            auto [it_from1, it_to1, it_from2, it_to2] = hop (wld, worms, table.nbi(choice));
            // Update onsite energy
            update_onsite_energy (it_to1, para);
            update_onsite_energy (worms.worm1.it, para);
            update_onsite_energy (it_to2, para);
            update_onsite_energy (worms.worm2.it, para);
            // Measure hopping
            para.N_hop += 2;
        }
    }
    // 2. Halt
    else if (status == HALTED)
    {
        // 2.1 bounce
        int n1 = (worms.worm1.it->z.create ? worms.worm1.it->z.naft() : worms.worm1.it->z.nbef);
        int n2 = (worms.worm2.it->z.create ? worms.worm2.it->z.naft() : worms.worm2.it->z.nbef);
        double weight = para.t * n1 * para.t * n2;
        table.add (weight);

        // 2.2 remove interaction
        auto [it_block1, it_block2] = get_blocks (worms);
        auto [it_to1, it_to2] = find_it_to (it_block1, it_block2);

        double Ebef = it_to1->z.Ebef + it_to2->z.Ebef,
               Eaft = it_to1->z.Eaft + it_to2->z.Eaft;
        tie(Ebef,Eaft) = shift_energy (Ebef, Eaft, para.EOFFSET);
        weight = (worms.up() ? Eaft : Ebef);
        table.add (weight);

        // 2.3 relink interaction
        for(bool go_positive : {false, true})
        {
            int nbi = get_neighbor_index (worms.dir, go_positive);

            // If neighbor does not exist, skip
            if (!wld.latt(it_to1->site).has_nbi.at(nbi) or !wld.latt(it_to2->site).has_nbi.at(nbi))
                continue;

            auto [it_nb1, it_nb2] = get_dipole_nb_aft (it_to1, it_to2, !worms.up(), go_positive, worms.dir);

            // If hop to the original site, skip
            if (it_nb1->site == worms.worm1.it->site)
                continue;

            int n1 = it_nb1->z.nbef,
                n2 = it_nb2->z.nbef;
            bool creat1 = (worms.worm1.up == worms.worm1.it->z.create),
                 creat2 = (worms.worm2.up == worms.worm2.it->z.create);
            // Check whether cannot creat or annihilate particle anymore
            if (check_nlimit (n1, creat1, para.nMAX, para.nMIN)
            and check_nlimit (n2, creat2, para.nMAX, para.nMIN))
            {
                double weight1 = para.t * (creat1 ? n1+1 : n1),
                       weight2 = para.t * (creat2 ? n2+1 : n2),
                weight = weight1 * weight2;
                table.add (weight, go_positive);
            }
        }

        // 2.4 Choose
#ifdef LOCAL_OPTIMAL
        table.locally_optimal();
#endif
        int choice = table.get_choice (uniform_dis(rand));
        if (choice == 0)        // bounce
        {
#ifdef PRINT_DETAIL
cout << "Blocked bounce" << endl;
#endif
            worms.worm1.up = !worms.worm1.up;
            worms.worm2.up = !worms.worm2.up;
        }
        else if (choice == 1)   // delete hopping
        {
#ifdef PRINT_DETAIL
cout << "Delete hop" << endl;
#endif
            delete_hop (wld, worms);
            // Measure hopping
            para.N_hop -= 2;
        }
        else                    // relink hopping
        {
#ifdef PRINT_DETAIL
cout << "Relink hop" << endl;
#endif
            delete_hop (wld, worms);
            worms.worm1.up = !worms.worm1.up;
            worms.worm2.up = !worms.worm2.up;
            auto [it_from1, it_to1, it_from2, it_to2] = hop (wld, worms, table.nbi(choice));
            // Update onsite energy
            update_onsite_energy (it_to1, para);
            update_onsite_energy (worms.worm1.it, para);
            update_onsite_energy (it_to2, para);
            update_onsite_energy (worms.worm2.it, para);
        }
    }
}

void updateZ (WLDiagram& wld, std::mt19937& rand, Para& para)
{
    std::uniform_real_distribution<double> uniform_dis (0., 1.);
    std::uniform_int_distribution<int> choose2 (0, 1);

    // Insert dipole worms
    // 1. Set direction
    DipoleWorm worms;
    worms.dir = LatticeSite::Direction::X;
    // 2. Choose sites
    int site_range = (wld.latt.xpbc() ? wld.latt.lx()-1 : wld.latt.lx()-2);
    std::uniform_int_distribution<int> distrib (0, site_range);
    int site1 = distrib (rand);
    int site2 = site1 + 1;
    // 3. Choose time, up, create
    double time = wld.beta() * uniform_dis(rand);
    bool up = choose2 (rand);
    bool create = choose2 (rand);
    // 4. Find the iterators to insert
    auto it1 = wld.find_aft (site1, time);
    auto it2 = wld.find_aft (site2, time);
    // 5. If particle cannot be inserted, finish the update
    int n1 = it1->z.nbef,
        n2 = it2->z.nbef;
    if (!check_nlimit (n1, create, para.nMAX, para.nMIN)
    or  !check_nlimit (n2, !create, para.nMAX, para.nMIN))
        return;
    // 6. Insert worms
    auto [worm1, fix1] = insert (wld, it1, time, create, up);
    auto [worm2, fix2] = insert (wld, it2, time, !create, up);
    worms.worm1 = worm1;
    worms.worm2 = worm2;
    //    Update onsite energy
    update_onsite_energy (worm1.it, para);
    update_onsite_energy (fix1, para);
    update_onsite_energy (worm2.it, para);
    update_onsite_energy (fix2, para);

    // Update the world-line diagram
    Status status = UNKNOWN;
    do
    {
check_onsite_energy (wld, para);
        double t0 = worms.worm1.it->time;
        // Compute the diagonal energy
        double Ebef = worms.worm1.it->z.Ebef + worms.worm2.it->z.Ebef,
               Eaft = worms.worm1.it->z.Eaft + worms.worm2.it->z.Eaft;
        // Get the time to walk
        auto [Ebef_shift, Eaft_shift] = shift_energy (Ebef, Eaft, para.EOFFSET);
        double Ewalk = (worms.up() ? Ebef_shift : Eaft_shift);
        double dtime = -log (uniform_dis(rand)) / Ewalk;

        // Walk
        status = liveWalk (wld, worms, dtime);

        // Update onsite energy
        /*double t1;
        if (status == CROSS_END)
            t1 = (worms.up() ? wld.beta() : 0.);
        else if (status == REMOVED)
            t1 = fix1->time;
        else
            t1 = worms.worm1.it->time;*/

        // Act
        act (wld, worms, rand, status, Ewalk, para);
    }
    while (status != REMOVED);
}

void set_onsiteE_for_end_nodes (WLDiagram& wld, Para& para)
{
    for(int i = 0; i < wld.latt.size(); i++)
    {
        Iter it = wld.end(i);
        update_onsite_energy (it, para);
    }
    para.N_hop = 0;
}
#endif
