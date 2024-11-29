#ifndef __MEASURE_HPP_CMC__
#define __MEASURE_HPP_CMC__
#include "WLDiagram.hpp"

void get_onsite_obs (const WLDiagram& wld, double& E, vector<double>& ns)
{
    E = 0.;
    ns.resize (wld.nodeList.size());
    for(int i = 0; i < wld.nodeList.size(); i++)
    {
        const auto& line = wld.nodeList[i];
        ns[i] = 0.;
        //double tsum = 0.;
        double t0 = 0.;
        for(auto it = ++line.begin(); it != line.end(); it++)
        {
            double t = it->time;
            E += (t - t0) * it->z.Ebef;
            ns[i] += (t - t0) * it->z.nbef;
            //tsum += (t - t0);
            t0 = t;
        }
        //assert (std::abs(tsum - wld.beta()) < 1e-12);
    }
}

int get_Np (const WLDiagram& wld)
{
    int N = 0;
    for(int i = 0; i < wld.size(); i++)
    {
        N += wld.end(i)->z.nbef;
    }
    return N;
}
#endif
