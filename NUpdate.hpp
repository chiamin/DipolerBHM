#ifndef __NUPDATE_H_CMC__
#define __NUPDATE_H_CMC__

void shift_energy (vector<double>& Eds, vector<double>& Ens, const vector<bool>& valid)
{
    // Diagonal energy
    double Emin = std::numeric_limits<double>::infinity();//*std::min_element(std::begin(Eds), std::end(Eds));
    for(int i = 0; i < Eds.size(); i++)
    {
        if (valid.at(i) and Eds.at(i) < Emin)
            Emin = Eds.at(i);
    }
cout << "Emin = " << Emin << endl;
    for(auto& E : Eds)
        E -= Emin;
    // Operator energy (from c or cdag)
    Emin = *std::min_element(std::begin(Ens), std::end(Ens));
    for(auto& E : Ens)
        E /= Emin;
}

void ParNumUpdate (WLDiagram& wld, std::mt19937& rand, const Para& para, int nrange)
{
    int N = wld.size();

    // Define random number distributions
    std::uniform_int_distribution<int> distrib (0, N-1);
    std::uniform_real_distribution<double> uniform_dis (0., 1.);

    // Choose a random site
    int site = distrib (rand);

    // Get the worldline for the site
    auto& line = wld.nodeList.at(site);

    // Define variables
    vector<int> dns (2*nrange+1);                           // the amount of particle number to be changed
    std::iota (std::begin(dns), std::end(dns), -nrange);    // dns = [-nrange, -nrange+1, ..., nrange]

    vector<double> Eds (dns.size(), 0.);                    // diagonal energy
    vector<double> Ens (dns.size(), 1.);                    // energy from c or cdag operator
    vector<bool> valid (dns.size(), true);                  // whether it is possible to change the particle number

    vector<vector<double>> Eds_each (dns.size());           // diagonal energy for each segment
    for(auto& ei : Eds_each)
        ei.reserve (line.size()-1);

    // Iterate over the worldline segments and store their energies
    double t0 = 0.;
    for(auto it = ++line.begin(); it != line.end(); it++)
    {
        double dt = it->time - t0;
        // For each change of n
        for(int i = 0; i < dns.size(); i++)
        {
            int dn = dns[i];
            int n = it->z.nbef;
            int n_larger = (it->z.create ? n+1 : n);
            int n_new = n + dn;

            // Check if it is possible changing to the target particle number
            if (n_new > para.nMAX or n_new < para.nMIN)
            {
                valid.at(i) = false;
                Eds.at(i) = std::numeric_limits<double>::max();
            }
            if (valid.at(i))
            {
                // Diagonal energy
                double Ed;
                if (dn == 0)
                    Ed = it->z.Ebef;
                else
                    Ed = onsite_energy (site, n_new, para);
                Eds.at(i) += dt * Ed;
                Eds_each.at(i).push_back (Ed);
                // Operator energy
                if (!it->end)
                {
//cout << "++ " << i << " " << n_larger+dn << " " << std::sqrt(n_larger+dn) << " " << Ens.at(i) << endl;
                    Ens.at(i) *= std::sqrt(n_larger+dn);
                }
            }
        }
        dt = it->time;
    }

//for(int i = 0; i < Eds.size(); i++)
//cout << "!! " << Eds.at(i) << " " << Ens.at(i) << endl;

    //shift_energy (Eds, Ens, valid);

    // Choose the particle number according to their weights
    ProbTable table;
    for(int i = 0; i < dns.size(); i++)
    {
//cout << "* i " << i << endl;
        if (valid.at(i))
        {
//cout << (Ens.at(i) * exp(-Eds.at(i))) << " " << Ens.at(i) << " " << Eds.at(i) << endl;
            table.add (Ens.at(i) * exp(-Eds.at(i)), i);
        }
    }
    int choice = table.get_choice (uniform_dis(rand));
    int ii = table.nbi(choice);

    // Update the worldline
    if (dns.at(ii) != 0)
    {
        auto it_E = Eds_each.at(ii).begin();
        for(auto it = ++line.begin(); it != line.end(); it++)
        {
            auto it_bef = --Iter(it);
            it->z.nbef += dns.at(ii);
            it->z.Ebef = *it_E;
            it_bef->z.Eaft = *it_E;
            it_E++;
        }
    }
}

#endif
