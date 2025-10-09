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

struct config {
    double ue, pe, rhoe, Te, Tw, delta_main, Ny;
};


class Profile {

    private:

    double k, B, PI;
    double ue, pe, rhoe, rhow, Te, Tw, delta_main, Ny, Me, r, Taw, mu_e;

    double mu_w, nu_w;


    vector<double> u_plus, u_plus_old, up_guess, y_plus, u_vd;
    vector<double> u, ul, ur, yl, y, yr, rhol, rho, rhor, Tl, T, Tr, v, w;
    vector<double> S, G, g, x;
    double G_top;
    vector<double> u_ue;


    double delta;               // delta currently being used
    vector<double> deltas;      // Vector of deltas

    double u_tau;               // u_tau currently being used
    vector<double> u_taus;      // Vector of u_taus

    double Re_tau;              // Re_tau currently being used
    vector<double> Re_taus;     // Vector of Re_taus

    vector<double> xs;          // Vector of x-positions

    double a0, a1, a2;

    public:

 
    void get_full_profile();
    tuple<vector<double>, vector<double>, vector<double>> get_single_uy_profile(int num);



    pair<vector<double>, vector<double>> find_T_rho(vector<double>& ui);
    void add_wake(vector<double>& yi);





    
    Profile(config inputs);
    void compute_all_BL_vars();
    void create_u_guess();
    void createG();
    void create_yplus();
    vector<double> get_y_from_yplus();
    vector<double> find_u();
    double invert_G(double S);

    void solve_uplus(); 
    double Newton_uplus(double& yp, double& guess); 

    inline double sutherlands(double T);
    void writeCsvProfile(const std::string& filename);
    void writeTecplotProfile(const std::string& filename);
    void generateLine(const std::string& filename);

    void solve_v();
    
    vector<double> linear_interpolate(
    const vector<double>& y_data,
    const vector<double>& f_data,
    const vector<double>& y_new);


};

vector<double> thomas_solve(const std::vector<double>& a,  // subdiag: a[0]=0
             const std::vector<double>& b,  // diag
             const std::vector<double>& c,  // superdiag: c[N-1]=0
             const std::vector<double>& d);  // RHS
