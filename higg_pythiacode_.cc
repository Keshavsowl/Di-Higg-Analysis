#include "Pythia8/Pythia.h"
#include <iostream>

using namespace Pythia8;
using namespace std;

int main() {
    Pythia pythia;
    pythia.readString("Higgs:useBSM = on");
    pythia.readString("HiggsBSM:ffbar2A3H2 = on");
    pythia.readString("25:m0 = 125.0");
    pythia.readString("25:onMode = off");
    pythia.readString("25:onIfMatch = 22 22"); 
    pythia.readString("35:m0 = 125.0"); 
    pythia.readString("35:onMode = off");
    pythia.readString("35:onIfMatch = 5 -5");

    pythia.init();
    for (int iEvent = 0; iEvent < 100; ++iEvent) { 
        pythia.next();
        bool hasPhotonDecay = false;
        bool hasBottomDecay = false;
        int photonIdx1 = -1;
        int photonIdx2 = -1;
        int bottomIdx1 = -1;
        int bottomIdx2 = -1;

        for (int iH = 0; iH < pythia.event.size(); ++iH) {
            if (pythia.event[iH].id() == 25) {
                int iDau1 = pythia.event[iH].daughter1();
                int daughterID1 = pythia.event[iDau1].id();

                int iDau2 = pythia.event[iH].daughter2();
                int daughterID2 = pythia.event[iDau2].id();

                if (daughterID1 == 22 && daughterID2 == 22) {
                    hasPhotonDecay = true;
                    photonIdx1 = iDau1; 
                    photonIdx2 = iDau2;
                }
            }
            if (pythia.event[iH].id() == 35) {
                int iDau1 = pythia.event[iH].daughter1();
                int daughterID1 = pythia.event[iDau1].id();

                int iDau2 = pythia.event[iH].daughter2();
                int daughterID2 = pythia.event[iDau2].id();

                if ((daughterID1 == 5 && daughterID2 == -5) || (daughterID1 == -5 && daughterID2 == 5)) {
                    hasBottomDecay = true;
                    if (bottomIdx1 == -1)
                        bottomIdx1 = iDau1;
                    else
                        bottomIdx2 = iDau1; 
                }
            }
        }

        if (hasPhotonDecay && hasBottomDecay) {
            cout << "Event " << iEvent << " contains both photon and bottom decays of the Higgs." << endl;
            cout << "Photon 1: px = " << pythia.event[photonIdx1].px() << ", py = " << pythia.event[photonIdx1].py()
                 << ", pz = " << pythia.event[photonIdx1].pz() << ", pseudorapidity = " << pythia.event[photonIdx1].eta()
                 << ", energy = " << pythia.event[photonIdx1].e() << endl;
            cout << "Photon 2: px = " << pythia.event[photonIdx2].px() << ", py = " << pythia.event[photonIdx2].py()
                 << ", pz = " << pythia.event[photonIdx2].pz() << ", pseudorapidity = " << pythia.event[photonIdx2].eta()
                 << ", energy = " << pythia.event[photonIdx2].e() << endl;
            cout << "Bottom Quark 1: px = " << pythia.event[bottomIdx1].px() << ", py = " << pythia.event[bottomIdx1].py()
                 << ", pz = " << pythia.event[bottomIdx1].pz() << ", pseudorapidity = " << pythia.event[bottomIdx1].eta()
                 << ", energy = " << pythia.event[bottomIdx1].e() << endl;
            if (bottomIdx2 != -1) {
                cout << "Bottom Quark 2: px = " << pythia.event[bottomIdx2].px() << ", py = " << pythia.event[bottomIdx2].py()
                     << ", pz = " << pythia.event[bottomIdx2].pz() << ", pseudorapidity = " << pythia.event[bottomIdx2].eta()
                     << ", energy = " << pythia.event[bottomIdx2].e() << endl;
            }
        } else if (hasPhotonDecay) {
            cout << "Event " << iEvent << " contains photon decay of one Higgs and bottom decay of the other Higgs." << endl;
            cout << "Photon 1: px = " << pythia.event[photonIdx1].px() << ", py = " << pythia.event[photonIdx1].py()
                 << ", pz = " << pythia.event[photonIdx1].pz() << ", pseudorapidity = " << pythia.event[photonIdx1].eta()
                 << ", energy = " << pythia.event[photonIdx1].e() << endl;
            cout << "Photon 2: px = " << pythia.event[photonIdx2].px() << ", py = " << pythia.event[photonIdx2].py()
                 << ", pz = " << pythia.event[photonIdx2].pz() << ", pseudorapidity = " << pythia.event[photonIdx2].eta()
                 << ", energy = " << pythia.event[photonIdx2].e() << endl;
            cout << "Bottom Quark 1: px = " << pythia.event[bottomIdx1].px() << ", py = " << pythia.event[bottomIdx1].py()
                 << ", pz = " << pythia.event[bottomIdx1].pz() << ", pseudorapidity = " << pythia.event[bottomIdx1].eta()
                 << ", energy = " << pythia.event[bottomIdx1].e() << endl;
        } else if (hasBottomDecay) {
            cout << "Event " << iEvent << " contains bottom decay of one Higgs and photon decay of the other Higgs." << endl;
            cout << "Photon 1: px = " << pythia.event[photonIdx1].px() << ", py = " << pythia.event[photonIdx1].py()
                 << ", pz = " << pythia.event[photonIdx1].pz() << ", pseudorapidity = " << pythia.event[photonIdx1].eta()
                 << ", energy = " << pythia.event[photonIdx1].e() << endl;
            cout << "Photon 2: px = " << pythia.event[photonIdx2].px() << ", py = " << pythia.event[photonIdx2].py()
                 << ", pz = " << pythia.event[photonIdx2].pz() << ", pseudorapidity = " << pythia.event[photonIdx2].eta()
                 << ", energy = " << pythia.event[photonIdx2].e() << endl;
            cout << "Bottom Quark 1: px = " << pythia.event[bottomIdx1].px() << ", py = " << pythia.event[bottomIdx1].py()
                 << ", pz = " << pythia.event[bottomIdx1].pz() << ", pseudorapidity = " << pythia.event[bottomIdx1].eta()
                 << ", energy = " << pythia.event[bottomIdx1].e() << endl;
            if (bottomIdx2 != -1) {
                cout << "Bottom Quark 2: px = " << pythia.event[bottomIdx2].px() << ", py = " << pythia.event[bottomIdx2].py()
                     << ", pz = " << pythia.event[bottomIdx2].pz() << ", pseudorapidity = " << pythia.event[bottomIdx2].eta()
                     << ", energy = " << pythia.event[bottomIdx2].e() << endl;
            }
        }
    }
    pythia.stat();
    return 0;
}

