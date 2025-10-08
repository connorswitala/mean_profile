#include "profile.hpp"

// Class constructor
Profile::Profile(config inputs) {

    cout << endl << "=== Flow conditions ===" << endl;
    ue = inputs.ue;
    cout << "-- ue [m/s] = " << ue << endl;

    pe = inputs.pe;
    cout << "-- pe [Pa] = " << pe << endl;

    Te = inputs.Te;
    cout << "-- Te [K] = " << Te << endl;

    Tw = inputs.Tw;
    cout << "-- Tw [K] = " << Tw << endl;

    Me  = ue / sqrt(gam * gcon * Te);
    cout << "-- Me = " << Me << endl;

    r   = cbrt(Pr);                     
    Taw = Te * (1.0 + r * 0.5 * (gam - 1.0) * Me * Me);
    cout << "-- Taw [K] = " << Taw << endl;

    rhoe = inputs.rhoe;
    cout << "-- rhoe [kg/m^3] = " << rhoe << endl;

    rhow = pe/(gcon * Tw);
    cout << "-- rhow [kg/m^3] = " << rhow << endl;

    delta_main = inputs.delta_main; 
    cout << "-- delta [m] = " << delta_main << endl;

    Ny = inputs.Ny;
    cout << "-- Ny = " << Ny << endl;

    mu_w = sutherlands(Tw);
    nu_w = mu_w / rhow;
    cout << "-- mu_w [m^2/s] = " << mu_w << endl;

    mu_e = sutherlands(Te);
    cout << "-- mu_e [m^2/s] = " << mu_e << endl;
    
    k = 0.41;
    B = 5.2;

    u_plus_old = vector<double>(Ny);
    u_plus = vector<double>(Ny);
    y_plus = vector<double>(Ny);
    u = vector<double>(Ny);
    ul = vector<double>(Ny);
    ur = vector<double>(Ny);
    rho = vector<double>(Ny);
    rhol = vector<double>(Ny);
    rhor = vector<double>(Ny);
    y = vector<double>(Ny);
    yl = vector<double>(Ny);
    yr = vector<double>(Ny);    
    S = vector<double>(Ny);
    v = vector<double>(Ny);
}

// This is the main function that does everthing.
void Profile::get_full_profile() {

    create_u_guess();
    create_yplus();
    createG();
    solve_uplus();
    compute_all_BL_vars();

    cout << endl << "=== Left state ===" << endl;
    tie(ul, rhol, yl) = get_single_uy_profile(0);
    cout << endl << "=== Center state ===" << endl;
    tie(u, rho, y) = get_single_uy_profile(1);
    cout << endl << "=== Right state ===" << endl;
    tie(ur, rhor, yr) = get_single_uy_profile(2);
    
    cout << endl << "-- Beginning interpolation" << endl;

    ul = linear_interpolate(yl, ul, y);
    ur = linear_interpolate(yr, ur, y);
    rhol = linear_interpolate(yl, rhol, y);
    rhor = linear_interpolate(yr, rhor, y);

    solve_v();

    writeCsvProfile("u_plus_plot.csv");
    generateLine("line1.dat");
    // writeTecplotProfile("u_plus_plot.dat");

    if (u[Ny - 10] == 869.1) cout << "-- h_cutoff = " << y[Ny - 10] << endl;

    
}

tuple<vector<double>, vector<double>, vector<double>> Profile::get_single_uy_profile(int num) {

    u_tau = u_taus[num];
    Re_tau = Re_taus[num];
    delta = deltas[num];

    vector<double> yi(Ny);
    vector<double> ui(Ny);       

    yi = get_y_from_yplus(); 
    add_wake(yi);
    ui = find_u();
    auto[Ti, rhoi] = find_T_rho(ui);
    if (num == 1) T = Ti;
    return {ui, rhoi, yi};
}


// This function solves the u+ vector given a y+ vector.
void Profile::solve_uplus() {

    for (int j = 0; j < Ny; ++j) 
        u_plus[j] = Newton_uplus(y_plus[j], up_guess[j]);
}

void Profile::compute_all_BL_vars() {
    // ----- Inputs assumed set on the class -----
    // Te, Taw, Tw, rhoe, rhow, ue, mu_e, nu_w, delta_main

    // Reynolds number per unit length: Re_x = (rho_e * U_e / mu_e) * x
    const double Re_unit = rhoe * ue / mu_e;           // [1/m]
    if (Re_unit <= 0.0) throw std::runtime_error("Re_unit must be > 0");
    if (delta_main <= 0.0) throw std::runtime_error("delta_main must be > 0");

    // --- Compressibility correction pieces (your F) ---
    const double a = Taw / Te;
    const double b = Tw  / Te;

    const double disc = (a + b) * (a + b) - 4.0 * b;   // inside sqrt
    auto clamp01 = [](double x){ return x < -1.0 ? -1.0 : (x > 1.0 ? 1.0 : x); };

    double F = 1.0;  // fallback if the formula gets invalid
    if (disc > 0.0) {
        const double root = std::sqrt(disc);
        const double kappa = (a + b - 2.0) / root;
        const double vpar  = (a - b)       / root;     // renamed from 'v' to avoid shadowing
        const double s1 = std::asin(clamp01(kappa));
        const double s2 = std::asin(clamp01(vpar));
        const double denom = (s1 + s2) * (s1 + s2);
        if (denom > 0.0) F = (a - 1.0) / denom;
    }
    if (!std::isfinite(F) || F <= 0.0) F = 1.0;        // keep things sane

    // --- Helpers for turbulent ZPG correlations (δ/x and Cf laws) ---
    auto delta_of_x = [&](double x){
        // δ/x = 0.16 Re_x^{-1/7}  =>  δ = 0.16 * x * (Re_unit * x)^{-1/7}
        return 0.16 * x * std::pow(Re_unit * x, -1.0/7.0);
    };
    auto Cf_of_x = [&](double x){
        // Cf = 0.027 Re_x^{-1/7}, with your compressibility divider F
        return (0.027 * std::pow(Re_unit * x, -1.0/7.0)) / F;
    };

    // Invert δ(x) to get the center location x_c for the given δ_c = delta_main:
    // δ = 0.16 * Re_unit^{-1/7} * x^{6/7}  =>  x = [δ/0.16 * Re_unit^{1/7}]^{7/6}
    const double delta_c = delta_main;
    const double xc = std::pow( (delta_c / 0.16) * std::pow(Re_unit, 1.0/7.0), 7.0/6.0 );
    if (!std::isfinite(xc) || xc <= 0.0) throw std::runtime_error("xc invalid from δ inversion");

    // Choose your left/right spacing around xc
    const double dx = 0.001;     // meters (adjust as needed, but keep > 0)
    const double xl = xc - dx;
    const double xr = xc + dx;
    if (xl <= 0.0) throw std::runtime_error("xl must be > 0 (move xc or shrink dx)");

    // --- Compute δ, Cf, τw, uτ, Reτ at each station ---
    auto station = [&](double x){
        const double delta = delta_of_x(x);
        const double Cf    = Cf_of_x(x);
        const double tau_w = 0.5 * Cf * rhoe * ue * ue;     // τw uses edge density for Cf law
        const double u_tau = std::sqrt(tau_w / rhow);        // friction velocity based on ρ_w
        const double Re_tau = u_tau * delta / nu_w;
        return std::tuple<double,double,double,double>(delta, u_tau, Re_tau, Cf);
    };

    double delta_l, u_tau_l, Re_tau_l, Cf_l;
    double delta_m, u_tau_m, Re_tau_m, Cf_m;
    double delta_r, u_tau_r, Re_tau_r, Cf_r;

    std::tie(delta_l, u_tau_l, Re_tau_l, Cf_l) = station(xl);
    std::tie(delta_m, u_tau_m, Re_tau_m, Cf_m) = station(xc);
    std::tie(delta_r, u_tau_r, Re_tau_r, Cf_r) = station(xr);

    // --- Store outputs on the class ---
    xs      = {xl, xc, xr};
    deltas  = {delta_l, delta_m, delta_r};
    u_taus  = {u_tau_l, u_tau_m, u_tau_r};
    Re_taus = {Re_tau_l, Re_tau_m, Re_tau_r};

    // --- Report ---
    std::cout << "\n=== Boundary layer values ===\n";
    std::cout << "-- Left:   delta = " << delta_l << "\t u_tau = " << u_tau_l
              << "\t Re_tau = " << Re_tau_l << "\t Cf = " << Cf_l << "\t x = " << xl << "\n";
    std::cout << "-- Center: delta = " << delta_m << "\t u_tau = " << u_tau_m
              << "\t Re_tau = " << Re_tau_m << "\t Cf = " << Cf_m << "\t x = " << xc << "\n";
    std::cout << "-- Right:  delta = " << delta_r << "\t u_tau = " << u_tau_r
              << "\t Re_tau = " << Re_tau_r << "\t Cf = " << Cf_r << "\t x = " << xr << "\n";
}


//T This function uses a Newton solver to solve from u+ from y+ using Spalding Law of the wall
double Profile::Newton_uplus(double& yp, double& guess) {

    double residual = 1;
    double f, fp, upold, upnew;
    upold = guess;

    while (residual > 1e-10) {
        f = upold
            + exp(-k*B) * (exp(k*upold) - 1 
            - k * upold
            - 0.5 * (k * upold) * (k * upold)
            - 0.333333 * (k * upold) * (k * upold) * (k * upold) ) - yp;

        fp = 1
            + exp(-k*B) * ( - 0.5 * k * k * k * upold * upold 
            - k * k * upold 
            + k * exp(k * upold)
            + k);

        upnew = upold  - f/fp;

        residual = fabs(upnew - upold);
        upold = upnew;
    }

    return upnew;
}

// This function adds Cole's wake term to u+
void Profile::add_wake(vector<double>& yi) {

    double eta;
    double ydelta = delta * rhow * u_tau / mu_w;
    double u_int = Newton_uplus(ydelta, up_guess[Ny / 2]);

    PI = 0.5 * (ue/u_tau * G_top - u_int);
    cout << "-- PI = " << PI << endl;
    
    for (int j = 0; j < Ny; ++j) {
        eta = min(1.0, yi[j] / delta);
        S[j] = u_tau/ue * (u_plus[j] + PI * 2 * sin(0.5 * pi * eta) * sin(0.5 * pi * eta));
        // cout << "S[" << j << "] = " << S[j] << endl;
    }
}

// This function finds u from u+_VD
vector<double> Profile::find_u() {

    vector<double> ui(Ny);

    for (int j = 0; j < Ny; ++j) {
       double u_ue = invert_G(S[j]);
       ui[j] = u_ue * ue;
    }

    return ui;
    cout << "-- u[1] = " << ui[0] << "\tu[Ny] = " << ui[Ny - 1] << endl;
}

// This function finds temperature and density from u(y)
pair<vector<double>, vector<double>> Profile::find_T_rho(vector<double>& ui) {

    vector<double> Ti(Ny);
    vector<double> rhoi(Ny);

    for (int j = 0; j < Ny; ++j) {
        Ti[j] = Te * (a0 + a1 * ui[j]/ue + a2 * ui[j] * ui[j] / (ue * ue) );
        rhoi[j] = pe / (gcon * Ti[j]);
    }

    cout << "-- T[1] = " << Ti[0] << "\tT[Ny] = " << Ti[Ny - 1] << endl;
    cout << "-- rho[1] = " << rhoi[0] << "\trho[Ny] = " << rhoi[Ny - 1] << endl;    

    return make_pair(Ti, rhoi);
}

// This functions computes viscosity at the wall using Sutherland's law
inline double Profile::sutherlands(double T) {
    return 1.716e-5 * pow(T/273.5, 1.5) * (273.15 + 110.4)/(T + 110.4);
}

// This function creates a guess for u+ for the Newton solver
void Profile::create_u_guess() {

    up_guess = vector<double>(Ny);

    for (int j = 0; j < Ny; ++j) {

        double g = 1/k * log(1 + k * y_plus[j]) + B;
        up_guess[j] = min(y_plus[j], g);        
    }
}

// This function creates a vector of y+ from 0.1 to 400
void Profile::create_yplus() {

    double y_plus_min = 0.0;
    double y_plus_max = 275.0;

    for (int i = 0; i < Ny; ++i) {
        y_plus[i] = y_plus_min + i * y_plus_max / (Ny - 1);
    }

}

vector<double> Profile::get_y_from_yplus() {
    for (int i = 0; i < Ny; ++i) {
        y[i] = y_plus[i] * nu_w / u_tau;
    }
    return y;
}


// This function takes in the S value and finds the unique G value for it (uses binary search and interpolates)
double Profile::invert_G(double S) {
    // Clamp to domain
    if (S <= G.front()) return 0.0;
    if (S >= G.back())  return 1.0;  // consider nudging Pi or utau so S<=G(1)

    // Binary search for k: G[k-1] <= S <= G[k]
    int lo = 1, hi = (int)G.size() - 1;
    while (lo <= hi) {
        int mid = (lo + hi) >> 1;
        if (G[mid] < S) lo = mid + 1;
        else            hi = mid - 1;
    }
    int k = lo; // now G[k-1] <= S <= G[k]

    // Linear interpolation in (G,x)
    double t = (S - G[k-1]) / (G[k] - G[k-1]);
    return x[k-1] + t * (x[k] - x[k-1]);
}

// This functions creates a vector of cumulative Van Driest integrals for U / U_e from 0 to 1. It uses Walz's quadratic
// temperature variation to find rho(y)
void Profile::createG() {
    // Number of intervals; we will have N+1 points in [0,1]
    const int N = 700;
    x.resize(N + 1);
    g.resize(N + 1);
    G.resize(N + 1);

    // Walz quadratic coefficients: T(x)/Te = a0 + a1 x + a2 x^2, x = U/Ue
    a0 = Tw / Te;
    a1 = (Taw - Tw) / Te;
    a2 = (Te  - Taw) / Te;   // typically ≤ 0 for adiabatic/isothermal

    double dx = 1.0 / static_cast<double>(N);

    // Cumulative integral G(x) = ∫_0^x sqrt( rho/ρw ) dξ
    // with rho/ρw = a0 / (a0 + a1 x + a2 x^2)
    for (int k = 0; k <= N; ++k) {
        x[k] = k * dx; // 0 ... 1

        double denom = a0 + a1 * x[k] + a2 * x[k] * x[k];
        // tiny floor to avoid sqrt of negative due to FP error
        double ratio = a0 / max(denom, 1e-300);
        g[k] = std::sqrt(ratio);

        if (k == 0) {
            G[k] = 0.0;
        } else {
            G[k] = G[k - 1] + 0.5 * (g[k - 1] + g[k]) * dx; // trapezoid
        }

        // cout << "G[" << x[k] << "] = " << G[k] << endl;
    }

    G_top = G[N - 1];

    // Optional: sanity checks
    // assert(G.back() > 0);           // integral at x=1 must be positive
    // assert(std::is_sorted(x.begin(), x.end()));
}

// Write a line profile for Tecplot (.dat)
void Profile::writeTecplotProfile(const std::string& filename) {

    const std::string& zone_title = "profile";

    std::ofstream ofs(filename);
    if (!ofs) throw std::runtime_error("writeTecplotProfile: cannot open file: " + filename);

    ofs << "TITLE = \"Mean Profile Solver\"\n";
    ofs << "VARIABLES = \"y_plus\", \"y_me\" , \"u_plus\" , \"u_me\" , \"T_me\" , \"rho_me\"" << endl;
    ofs << "ZONE T=\"" << zone_title << "\", I=" << Ny << ", F=POINT\n";

    for (int j = 0; j < Ny; ++j)
        ofs << y_plus[j] << " " << y[j] << " " << u_plus[j] << " " << u[j] << " " << T[j] << " " << rho[j] << "\n";
}

#include <fstream>
#include <iomanip>

void Profile::generateLine(const std::string& filename) {
    std::ofstream ofs(filename);
    if (!ofs) throw std::runtime_error("generateLine: cannot open file: " + filename);

    ofs << "VARIABLES=x,y,z,n,rho,u,v,w,t,p\n";
    ofs << "ZONE f=point,i=   " << Ny << "\n";

    // scientific notation, 15 digits after the decimal, uppercase 'E'
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs << std::setprecision(15) << std::uppercase;

    for (int j = 0; j < Ny; ++j) {
        ofs << 0.0      << "   "
            << y[j]     << "   "
            << 0.0      << "   "
            << y[j]     << "   "
            << rho[j]   << "   "
            << u[j]     << "   "
            << v[j]     << "   "
            << 0.0      << "   "
            << T[j]     << "   "
            << pe       << "\n";
    }
}


void Profile::solve_v() {

    double dx = xs[2] - xs[0]; 
    vector<double> rhov(Ny);
    rhov[0] = 0;

    for (int j = 1; j < Ny; ++j) {
        rhov[j] = (( -(y[j] - y[j-1]) * (rhor[j] * ur[j] - rhol[j] * ul[j]) / dx ) + rho[j-1] * v[j-1]); 
    }


    double R_max = 0.0;
    double R;

    for (int j = 1; j < Ny; ++j) {
        R = (rho[j] * v[j] - rho[j-1] * v[j-1])/(y[j] - y[j-1]) + (rhor[j] * ur[j] - rhol[j] * ul[j])/(dx);
        if (R > R_max) R_max = R;
    }

    for (int j = 0; j < Ny; ++j) {
        v[j] = rhov[j] / rho[j];  
    }

    cout << "-- Max Residual = " << R_max << endl;

}

void Profile::writeCsvProfile(const std::string& filename) {
    std::ofstream ofs(filename);
    if (!ofs) throw runtime_error("writeCsvProfile: cannot open file: " + filename);

    ofs << "y, ul, u, ur, v, T, rho\n";

    for (int i = 0; i < Ny; ++i)
        ofs << y[i] << "," << ul[i] << "," << u[i] << "," << ur[i] << "," << v[i] << "," << T[i] << "," << rho[i] << "\n";
}


vector<double> Profile::linear_interpolate(
    const vector<double>& y_data,
    const vector<double>& f_data,
    const vector<double>& y_new)
{
    if (y_data.size() != f_data.size()) {
        throw invalid_argument("y_data and f_data must be the same size.");
    }
    if (y_data.size() < 2) {
        throw invalid_argument("Need at least two data points to interpolate.");
    }

    vector<double> f_new(y_new.size());

    for (size_t j = 0; j < y_new.size(); ++j) {
        double y = y_new[j];

        // Handle out-of-bounds with clamping (or could extrapolate if desired)
        if (y <= y_data.front()) {
            f_new[j] = f_data.front();
            continue;
        }
        if (y >= y_data.back()) {
            f_new[j] = f_data.back();
            continue;
        }

        // Find interval [y_data[i], y_data[i+1]]
        size_t i = 0;
        while (i + 1 < y_data.size() && y > y_data[i+1]) {
            ++i;
        }

        double x0 = y_data[i];
        double x1 = y_data[i+1];
        double f0 = f_data[i];
        double f1 = f_data[i+1];

        // Linear interpolation formula
        f_new[j] = f0 + (f1 - f0) * ( (y - x0) / (x1 - x0) );
    }

    return f_new;
}




