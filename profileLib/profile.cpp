#include "profile.hpp"


Profile::Profile(config inputs) {

    ue = inputs.ue;
    pe = inputs.pe;
    rhoe = inputs.rhoe;
    Te = inputs.Te;
    Tw = inputs.Tw;
    Retau = inputs.Retau;
    delta = inputs.delta; 
    Ny = inputs.Ny;

    rhow = pe/(gcon * Tw);
    cout << "-- rho_w = " << rhow << endl;
    sutherlands();
    u_tau = Retau * nu_w / delta;
    cout << "-- u_tau = " << u_tau << endl;

    k = 0.41;
    B = 5.2;

    u_plus_old = vector<double>(Ny);
    u_plus = vector<double>(Ny);
    y_plus = vector<double>(Ny);
    u = vector<double>(Ny);
    y = vector<double>(Ny);
    S = vector<double>(Ny);
}

void Profile::get_profile() {
    solve_uplus();
    createG();
    add_wake();
    find_u(); 
    // smooth_u(50);
    find_T_rho();


    writeCsvProfile("u_plus_plot.csv");
    writeTecplotProfile("u_plus_plot.dat");
}

void Profile::solve_uplus() {

    create_yplus();
    create_u_guess();

    for (int j = 0; j < Ny; ++j) 
        u_plus[j] = Newton_uplus(y_plus[j], up_guess[j]);
    

    cout << "-- u+ solved for" << endl;
}

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

void Profile::add_wake() {

    double eta;
    double ydelta = delta * rhow * u_tau / mu_w;
    double u_int = Newton_uplus(ydelta, up_guess[Ny / 2]);

    PI = 0.5 * (ue/u_tau * G_top - u_int);
    cout << "-- PI = " << PI << endl;

    for (int j = 0; j < Ny; ++j) {
        eta = min(1.0, y[j] / delta);
        u_plus[j] += PI * 2 * sin(0.5 * pi * eta) * sin(0.5 * pi * eta);
        S[j] = u_tau/ue * u_plus[j];
        // cout << "S[" << j << "] = " << S[j] << endl;
    }

}



void Profile::find_u() {

    u = vector<double>(Ny);

    for (int j = 0; j < Ny; ++j) {
       double u_ue = invert_G(S[j]);
       u[j] = u_ue * ue;
    }

    cout << "-- u[1] = " << u[0] << "\tu[Ny] = " << u[Ny - 1] << endl;
}
// Blend u to Ue over the next `span` points, starting at i0 (where y[i0] > delta)
void Profile::smooth_u(int span) {

    // find first index above delta
    int i0 = -1;
    for (int i = 0; i < Ny; ++i) {
        if (y[i] > delta) { i0 = i; break; }
    }
    if (i0 < 0) return; // nothing to do (top not reached)

    // cap the span so we don't run past the array
    int i1 = i0 + span;

    // linearly blend from u[i0] to Ue across [i0 .. i1]
    for (int k = i0; k <= i1; ++k) {
        double t = (y[k] - i0) / (i1 - i0);          // raw 0..1
        double s = smoothstep(t);                   // or smoothstep(t)
        u[k] = (1.0 - s) * u[k] + s * ue;            // monotone ease to Ue

    }
    // hold constant above the ramp
    for (int k = i1 + 1; k < Ny; ++k) u[k] = ue;
}


void Profile::find_T_rho() {

    T = vector<double>(Ny);
    rho = vector<double>(Ny);


    for (int j = 0; j < Ny; ++j) {
        T[j] = Te * (a0 + a1 * u[j]/ue + a2 * u[j] * u[j] / (ue * ue) );
        rho[j] = pe / (gcon * T[j]);
    }

    cout << "-- T[1] = " << T[0] << "\tT[Ny] = " << T[Ny - 1] << endl;
    cout << "-- rho[1] = " << rho[0] << "\trho[Ny] = " << rho[Ny - 1] << endl;    
}



void Profile::sutherlands() {
    mu_w = 1.716e-5 * pow(Tw/273.5, 1.5) * (273.15 + 110.4)/(Tw + 110.4);
    nu_w = mu_w/rhow;

    cout << "-- mu_w = " << mu_w << endl;
    cout << "-- nu_w = " << nu_w << endl;
}
void Profile::create_u_guess() {

    up_guess = vector<double>(Ny);

    for (int j = 0; j < Ny; ++j) {

        double g = 1/k * log(1 + k * y_plus[j]) + B;
        up_guess[j] = min(y_plus[j], g);        
    }
    cout << "-- guess for u+ created." << endl;
}
void Profile::create_yplus() {

    double y_plus_min = 0.1;
    double y_plus_max = 400.0;

    for (int i = 0; i < Ny; ++i) {
        y_plus[i] = y_plus_min + i * y_plus_max / (Ny - 1);
        y[i] = y_plus[i] * nu_w / u_tau;
    }

    cout << "-- y and y+ created." << endl;


}

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
void Profile::createG() {
    // Number of intervals; we will have N+1 points in [0,1]
    const int N = 700;
    x.resize(N + 1);
    g.resize(N + 1);
    G.resize(N + 1);

    // Edge Mach and recovery temperature (Walz)
    const double Me  = ue / sqrt(gam * gcon * Te);    // gcon = R (J/kg/K)
    const double r   = cbrt(Pr);                      // ≈ Pr^(1/3)
    const double Taw = Te * (1.0 + r * 0.5 * (gam - 1.0) * Me * Me); // K

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

    cout << "-- G created." << endl;
}


// Write a line profile for Tecplot (.dat)
// Columns: y_plus, u_plus (POINT format)
void Profile::writeTecplotProfile(const std::string& filename) {

    const std::string& zone_title = "profile";
    int precision = 10;

    if (y_plus.size() != u_plus.size())
        throw std::runtime_error("writeTecplotProfile: y_plus and u_plus must have same length.");
    const std::size_t N = y_plus.size();
    if (N == 0) throw std::runtime_error("writeTecplotProfile: empty vectors.");

    std::ofstream ofs(filename);
    if (!ofs) throw std::runtime_error("writeTecplotProfile: cannot open file: " + filename);

    ofs << "TITLE = \"u_plus vs y_plus\"\n";
    ofs << "VARIABLES = \"y_plus\", \"u_plus\"\n";
    ofs << "ZONE T=\"" << zone_title << "\", I=" << N << ", F=POINT\n";
    ofs << std::setprecision(precision) << std::scientific;

    for (std::size_t i = 0; i < N; ++i)
        ofs << y_plus[i] << " " << u_plus[i] << "\n";
}

// Write a profile for ParaView (.csv)
// Header row then data rows: y_plus,u_plus
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








