#include<iostream>
#include<vector>
#include<cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <algorithm> // std::clamp

using namespace std;

constexpr double gcon = 287.0;
constexpr double gam = 1.4;
constexpr double pi = 3.141592653;
constexpr double Pr = 0.72;

inline double smoothstep(double t) {
    t = std::clamp(t, 0.0, 1.0);
    return t*t*(3.0 - 2.0*t);                 // 3t^2 - 2t^3
}

struct config {
    double ue, pe, rhoe, Te, Tw, Retau, delta, Ny;
};


class Profile {

    private:

    double k, B, PI;
    double ue, pe, rhoe, rhow, Te, Tw, Retau, delta, Ny;

    double mu_w, nu_w, u_tau;


    vector<double> u_plus, u_plus_old, up_guess, y_plus, u_vd;
    vector<double> u, y, T, rho, v, w;
    vector<double> S, G, g, x;
    double G_top;
    vector<double> u_ue;


    double a0, a1, a2;

    public:

    Profile(config inputs); 
    void get_profile();

    void solve_uplus(); 
    void create_yplus();
    void create_u_guess();
    void writeCsvProfile(const std::string& filename);
    void writeTecplotProfile(const std::string& filename);
    void sutherlands();
    void createG();
    double invert_G(double S);
    void find_u();
    void find_T_rho();
    void add_wake();
    double Newton_uplus(double& yp, double& guess); 
    void smooth_u(int span);
};