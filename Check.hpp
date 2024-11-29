#ifndef __CHECK_HPP_CMC__
#define __CHECK_HPP_CMC__
#include "LiveDipoleWorm.hpp"

void check_onsiteE_node (const Iter it, const Para& para)
{
    double Ebef = onsite_energy (it->site, it->z.nbef, para);
    double Eaft = onsite_energy (it->site, it->z.naft(), para);
    assert (Ebef == it->z.Ebef);
    assert (Eaft == it->z.Eaft);
}

void check_Nhop (WLDiagram& wld, const Para& para)
{
    double E = 0.;
    int Nhop = 0;
    for(int i = 0; i < wld.nodeList.size(); i++)
    {
        double time_pre = 0.;
        for(Iter it = ++wld.nodeList[i].begin(); it != wld.nodeList[i].end(); it++)
        {
            if (!it->end)
                Nhop += 1;
            if (!(it == wld.beg(i)))
            {
                check_onsiteE_node (it, para);
                time_pre = it->time;
            }
        }
    }
    assert (Nhop/2 == para.N_hop);
}
#endif
