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
    cout << "-- delta [m] = " << delta << endl;


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
    

    // solve_uplus();
    // createG();
    // add_wake();
    // find_u(); 
    // // smooth_u(50);
    // find_T_rho();


    // writeCsvProfile("u_plus_plot.csv");
    // writeTecplotProfile("u_plus_plot.dat");
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

    return {ui, rhoi, yi};
}


// This function solves the u+ vector given a y+ vector.
void Profile::solve_uplus() {

    for (int j = 0; j < Ny; ++j) 
        u_plus[j] = Newton_uplus(y_plus[j], up_guess[j]);
}

void Profile::compute_all_BL_vars() {

    double a = Taw / Te;
    double b = Tw / Te;
    double Re = rhoe * ue / mu_e;

    double kk = (a + b - 2) / sqrt(((a + b) * (a + b) - 4 * b));
    double v = (a - b) / sqrt((a + b) * (a + b) - 4 * b);
    double F = (a - 1)/((asin(kk) + asin(v)) * (asin(kk) + asin(v)));


    // Center
    double delta_c = delta_main;
    double xc = pow(delta_c / 0.16 * pow(Re, 1.0/7.0), 7.0/6.0);
    double Cf = 0.027 * pow(Re * xc, -1.0/7.0);
    Cf = Cf / F;

    double tau_w_c = 0.5 * Cf * rhoe * ue * ue;
    double u_tau_c = sqrt(tau_w_c / rhow);
    double Re_tau_c = u_tau_c * delta_c / nu_w;

    double dx = 0.0005;

    // Left
    double xl = xc - dx;
    double delta_l = 0.16 * pow(Re * xl, -1.0/7.0) * xl;
    Cf = 0.027 * pow(Re * xl, -1.0/7.0);
    Cf = Cf / F;

    double tau_w_l = 0.5 * Cf * rhoe * ue * ue;
    double u_tau_l = sqrt(tau_w_l/rhow);
    double Re_tau_l = u_tau_l * delta_l / nu_w;

    // Right
    double xr = xc + dx;
    double delta_r = 0.16 * pow(Re * xr, -1.0/7.0) * xr;
    Cf = 0.027 * pow(Re * xr, -1.0/7.0);
    Cf = Cf / F;

    double tau_w_r = 0.5 * Cf * rhoe * ue * ue;
    double u_tau_r = sqrt(tau_w_r/rhow);
    double Re_tau_r = u_tau_r * delta_r / nu_w;

    xs = vector<double>(3);
    xs = {xl, xc, xr};

    u_taus = vector<double>(3);
    u_taus = {u_tau_l, u_tau_c, u_tau_r};

    Re_taus = vector<double>(3);
    Re_taus = {Re_tau_l, Re_tau_c, Re_tau_r};

    deltas = vector<double>(3);
    deltas = {delta_l, delta_c, delta_r};

    cout << endl << "=== Boundary layer values ===" << endl;
    cout << "-- Left: delta = " << deltas[0] << "\t u_tau = " << u_taus[0] << "\t Re_tau = " << Re_taus[0] << "\t x = " << xl << endl;
    cout << "-- Center: delta = " << deltas[1] << "\t u_tau = " << u_taus[1] << "\t Re_tau = " << Re_taus[1] << "\t x = " << xc  << endl;
    cout << "-- Right: delta = " << deltas[2] << "\t u_tau = " << u_taus[2] << "\t Re_tau = " << Re_taus[2] << "\t x = " << xr  << endl;
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

    double y_plus_min = 0.1;
    double y_plus_max = 400.0;

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

// Write a profile for ParaView (.csv)
void Profile::writeCsvProfile(const std::string& filename) {

    int precision = 10;

    if (y_plus.size() != u_plus.size())
        throw std::runtime_error("writeCsvProfile: y_plus and u_plus must have same length.");

    const size_t N = y_plus.size();
    if (N == 0) throw runtime_error("writeCsvProfile: empty vectors.");

    std::ofstream ofs(filename);
    if (!ofs) throw runtime_error("writeCsvProfile: cannot open file: " + filename);

    ofs << "y_plus, y, u_plus, u, T, rho\n";
    ofs << std::setprecision(precision) << std::scientific;

    for (std::size_t i = 0; i < N; ++i)
        ofs << y_plus[i] << "," << y[i] << "," << u_plus[i] << "," << u[i] << "," << T[i] << "," << rho[i] << "\n";
}








