#include "profileLib/profile.hpp"



int main() {

    config inputs;
    inputs.delta_main = 0.00124;
    inputs.Ny = 1000;
    inputs.rhoe = 0.044;
    inputs.Te = 55.2;
    inputs.Tw = 97.5;
    inputs.pe = inputs.rhoe * gcon * inputs.Te;
    inputs.ue = 869.1;
    
    Profile prof(inputs); 
    prof.get_full_profile();
    cout << endl << "-- Program finished." << endl;
}