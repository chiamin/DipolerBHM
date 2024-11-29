#include <iostream>
#include <fstream>
#include "DipoleWorm.hpp"
#include "ParaBox.hpp"

double get_time (const DipoleWorm& worms, double dt)
{
    double t = worms.worm1.it->time;
    if (worms.up())
        t += dt;
    else
        t -= dt;
    return t;
}

int main (int argc, char *argv[])
{
    WLDiagram wld (6, 1, 1, 2, true, true, true, 1);

    bool up = false;
    int site = 3;
    double time = 1;
    bool create = true;
    cout << "Insert" << endl;
    auto worms = insert (wld, site, time, create, up, LatticeSite::Direction::X);
    cout << "Walk" << endl;
    walk (wld, worms.worm1, get_time(worms,0.5));
    walk (wld, worms.worm2, get_time(worms,0.5);
    cout << "Hop" << endl;
    hop (wld, worms, false);
    cout << "Walk" << endl;
    walk (wld, worms.worm1, get_time(worms,0.2));
    walk (wld, worms.worm2, get_time(worms,0.2));
    cout << "Hop" << endl;
    hop (wld, worms, true);
    cout << "Walk" << endl;
    walk (wld, worms.worm1, get_time(worms,0.2));
    walk (wld, worms.worm2, get_time(worms,0.2));
    // bounce
    cout << "Bounce" << endl;
    worms.worm1.up = !worms.worm1.up;
    worms.worm2.up = !worms.worm2.up;
    cout << "Delete hop" << endl;
    del_hop (wld, worms.worm1);
    del_hop (wld, worms.worm2);
    cout << "Bounce" << endl;
    worms.worm1.up = !worms.worm1.up;
    worms.worm2.up = !worms.worm2.up;
    cout << "Hop" << endl;
    hop (wld, worms, false);
    cout << "Walk" << endl;
    walk (wld, worms.worm1, get_time(worms,0.2));
    walk (wld, worms.worm2, get_time(worms,0.2));
    // Cross end
    cout << "Cross end" << endl;
    cross_end (wld, worms.worm1);
    cross_end (wld, worms.worm2);
    cout << "Walk" << endl;
    walk (wld, worms.worm1, get_time(worms,0.2));
    walk (wld, worms.worm2, get_time(worms,0.2));
    cout << "Hop" << endl;
    hop (wld, worms, true);
    cout << "Walk" << endl;
    walk (wld, worms.worm1, get_time(worms,0.2));
    walk (wld, worms.worm2, get_time(worms,0.2));
    // Cross node
    cout << "Cross node" << endl;
    walk (wld, worms.worm1, std::abs(1.-worms.worm1.it->time));
    cross_node (wld, worms.worm2);
    cout << "Walk" << endl;
    walk (wld, worms.worm1, get_time(worms,0.2));
    walk (wld, worms.worm2, get_time(worms,0.2));
    cout << "Bounce" << endl;
    worms.worm1.up = !worms.worm1.up;
    worms.worm2.up = !worms.worm2.up;
    cout << "Hop" << endl;
    hop (wld, worms, true);
    cout << "Walk" << endl;
    walk (wld, worms.worm1, get_time(worms,0.1));
    walk (wld, worms.worm2, get_time(worms,0.1));
    cout << "Remove worms" << endl;
    remove_worms (wld, worms.worm1);
    remove_worms (wld, worms.worm2);

    ofstream ofs ("test.wld");
    wld.info_plot_lines (ofs);
    wld.info_plot_links (ofs);
    wld.info_plot_assoc (ofs);
    return 0;
}
