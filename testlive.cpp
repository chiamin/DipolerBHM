#include <iostream>
#include <fstream>
#include <map>
//#define PRINT_DETAIL
#define LOCAL_OPTIMAL
#include "LiveDipoleWorm.hpp"
#include "Check.hpp"
#include "ReadInput.h"
#include "NUpdate.hpp"
#include "Measure.hpp"

class Observables
{
    public:
        void measure (const std::string& name, double value)
        {
            if (obs.count(name) == 0)
            {
                obs[name].first = value;
                obs[name].second = 1;
            }
            else
            {
                auto& [val, count] = obs[name];
                val += value;
                count += 1;
            }
        }
        double get (const std::string& name)
        {
            auto [val, count] = obs[name];
            return val / double(count);
        }

    private:
        std::map<std::string, std::pair<double,int>> obs;
};

class VectorObservables
{
    public:
        void measure (const std::string& name, const vector<double>& value)
        {
            if (obs.count(name) == 0)
            {
                obs[name].first = value;
                obs[name].second = 1;
            }
            else
            {
                auto& [val, count] = obs[name];
                for(int i = 0; i < val.size(); i++)
                    val.at(i) += value.at(i);
                count += 1;
            }
        }
        vector<double> get (const std::string& name)
        {
            auto [val, count] = obs[name];
            for(auto& vi : val)
                vi /= double(count);
            return val;
        }

    private:
        std::map<std::string, std::pair<vector<double>,int>> obs;
};

int main (int argc, char *argv[])
{
    ifstream input = open_file (argv[1]);

    auto L        = read_value<int> (input,"L");
    auto mu       = read_value<double> (input,"mu");
    auto xpbc     = read_value<bool> (input,"xpbc");
    auto init_n   = read_value<int> (input,"init_n");
    auto beta     = read_value<double> (input,"beta");
    auto nMAX     = read_value<int> (input,"nMAX");
    auto MC_steps = read_value<int> (input,"MC_steps");
    auto seed     = read_value<unsigned long> (input,"seed");
    auto t        = read_value<double> (input,"t");
    auto U        = read_value<double> (input,"U");

    Para para;
    para.UHALF = 0.5 * U;
    para.t = std::sqrt(t);
    para.EOFFSET = read_value<double> (input,"EOffset");
    para.MU.resize (L, mu);

cout << "t = " << para.t << endl;

    WLDiagram wld (L, 1, 1, beta, xpbc, false, false, init_n);
    set_onsiteE_for_end_nodes (wld, para);

//check_onsiteE_Nhop (wld, para);

    Observables obs;
    VectorObservables obs_n;

    std::mt19937 rand (seed);
    double Edi;
    vector<double> ns;
    for(int i = 0; i < MC_steps; i++)
    {
//        ParNumUpdate (wld, rand, para, 1);
        updateZ (wld, rand, para);

assert (para.N_hop % 2 == 0);
int N_hop2 = int(para.N_hop) / 2;

        get_onsite_obs (wld, Edi, ns);
        int Npi = get_Np (wld);
        obs.measure ("Ed", Edi);
        obs.measure ("Np", Npi);
        obs.measure ("N_hop", N_hop2);//para.N_hop);
        obs_n.measure ("n", ns);


        double Ed = obs.get("Ed")/wld.beta();
        double Ek = -obs.get("N_hop")/wld.beta();
        double E = Ed + Ek;
        double Np = obs.get("Np");
        auto ns_mean = obs_n.get("n");

/*if (std::abs(Np - init_n*L) > 1e-8)
{
ofstream ofs ("test.wld");
wld.info_plot_lines (ofs);
wld.info_plot_links (ofs);
wld.info_plot_assoc (ofs);
cout << "Np = " << Np << endl;
exit(0);
}*/

        cout << i << " Ed, Ek, E, n = " << Ed << " " << Ek << " " << E << " " << Np/L << endl;
        for(int j = 0; j < ns_mean.size(); j++)
            cout << "n " << j << " " << ns_mean.at(j)/wld.beta() << endl;
//check_onsiteE_Nhop (wld, para);
//int gg; cin>>gg;

    }

ofstream ofs ("test.wld");
wld.info_plot_lines (ofs);
wld.info_plot_links (ofs);
wld.info_plot_assoc (ofs);

    return 0;
}
