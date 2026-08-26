#include "header.hpp"
#include "models/confluent-wet.hpp"
#include "error_msg.hpp"
#include "random.hpp"
#include "lb.hpp"
#include "tools.hpp"

using namespace std;
namespace opt = boost::program_options;

// from main.cpp:
extern unsigned nthreads;

ConfluentWet::ConfluentWet(unsigned LX_, unsigned LY_, unsigned BC_)
  : Nematic(LX_, LY_, BC_)
{
  // Nematic declares Gamma, xi, tau, friction, LL, CC and zeta with no in-class
  // initialiser, and the constructor runs before the runcard is parsed, so every one of
  // them has to be given a value here or an omitted option silently reads uninitialised
  // memory. LL is set for xi_N = sqrt(LL/(2 CC)) = 1.41 at CC = 0.1.
  Gamma = 0.05; xi = 0.4; tau = 1.; friction = 0.;
  CC = 0.1; LL = 0.4; zeta = 0.;
  angle_deg = 0.; noise = 0.05;
}

void ConfluentWet::Initialize()
{
  Nematic::Initialize();

  // ---- parameters that cannot possibly be meant, rejected loudly -----------
  if(tau <= .5)
    throw error_msg("confluent-wet: tau must exceed 1/2 (the LB viscosity is "
                    "rho*cs^2*(tau-1/2) and would be negative), got tau=", tau, ".");
  if(rho <= 0.)
    throw error_msg("confluent-wet: rho must be positive, got ", rho, ".");
  if(CC <= 0.)
    throw error_msg("confluent-wet: CC must be positive, got ", CC, ".");
  if(LL <= 0.)
    throw error_msg("confluent-wet: LL must be positive (Frank elasticity), got ", LL, ".");
  if(friction < 0.)
    throw error_msg("confluent-wet: friction must be non-negative, got ", friction, ".");
  if(director_config != "uniform" && director_config != "defect-pair")
    throw error_msg("confluent-wet: director-config must be 'uniform' or 'defect-pair', "
                    "got '", director_config, "'.");
  if(director_config == "defect-pair" && defect_sep <= 0.)
    throw error_msg("confluent-wet: director-config=defect-pair needs defect-sep > 0.");
  if(init_order < 0.)
    throw error_msg("confluent-wet: initial-order must be non-negative.");
  if(chi_config != "uniform" && chi_config != "noise")
    throw error_msg("confluent-wet: chi-config must be 'uniform' or 'noise', got '",
                    chi_config, "'.");
  if(chi0 < 0. || chi0 > 1.)
    throw error_msg("confluent-wet: chi0 must lie in [0,1].");
  if(m0 < 0. || m0 > 1.)
    throw error_msg("confluent-wet: m0 must lie in [0,1].");
  if(tau_m <= 0.)
    throw error_msg("confluent-wet: tau-m must be positive.");
  if(chi_width <= 0.)
    throw error_msg("confluent-wet: chi-width must be positive (a hard switch is "
                    "chi-width -> 0, not 0).");
  if(switch_sign != 1 && switch_sign != -1)
    throw error_msg("confluent-wet: switch-sign must be +1 or -1.");
  if(pmem_width < 0.)
    throw error_msg("confluent-wet: pmem-width must be non-negative.");
  if(Dbio < 0.)
    throw error_msg("confluent-wet: Dbio must be non-negative.");
  // Explicit five-point Laplacian with dt = dx = 1 is stable only for D <= 1/4 in 2D.
  // Rejected rather than warned about: an unstable Dbio corrupts chi and m silently,
  // because the [0,1] clamp hides the blow-up and leaves a plausible-looking binary field.
  if(Dbio > 0.2)
    throw error_msg("confluent-wet: Dbio=", Dbio, " exceeds the explicit-diffusion "
                    "stability limit (D <= 1/4 in 2D at dt = dx = 1; 0.2 leaves margin "
                    "for the advection term). If the flow needs more than this to keep "
                    "chi and m resolved -- the requirement is Dbio ~ 9*u_rms/l_a for a "
                    "3-cell Batchelor scale -- then zeta is too large for the lattice, "
                    "not the other way round.");

  // tau_chi <= 0 is the frozen-chi control, not an error.
  use_switching = (tau_chi > 0.);

  // NOTE ON tau_m AND THE ACOUSTIC TIME -- deliberately NOT a check.
  //
  // The pressure field is established by sound waves, so a memory faster than the acoustic
  // time integrates the transient rather than the load. But the relevant acoustic time is
  // set by the PRESSURE CORRELATION LENGTH L_P, not by the box: P at a point is a Poisson
  // solution dominated by sources within ~L_P, and the contributions of distant sources
  // cancel for a statistically isotropic source field. The floor is therefore
  //
  //     tau_m >~ sqrt(3) * L_P ,   NOT   sqrt(3) * L
  //
  // and L_P is measurable only after the fact -- at L = 200 the velocity and chi structure
  // scales came out at 25-42 lattice units against a box of 200, so a box-based bound is
  // roughly an order of magnitude too strict. An earlier version of this file enforced
  // sqrt(3)*L as a hard startup error; that was wrong, and rejected legitimate runs. Check
  // tau_m against the measured L_P of a frozen-chi calibration run instead.

  chi.SetSize(LX, LY, Type);
  chiN.SetSize(LX, LY, Type);
  chi_tmp.SetSize(LX, LY, Type);
  m.SetSize(LX, LY, Type);
  mN.SetSize(LX, LY, Type);
  m_tmp.SetSize(LX, LY, Type);
  ux_mat.SetSize(LX, LY, Type);
  uy_mat.SetSize(LX, LY, Type);
  sigma_bulk.SetSize(LX, LY, Type);
  sigma_active_xx.SetSize(LX, LY, Type);
  sigma_active_yx.SetSize(LX, LY, Type);
  pressure.SetSize(LX, LY, Type);
  pressure_lb.SetSize(LX, LY, Type);
}

void ConfluentWet::ConfigureAtNode(unsigned k)
{
  const double x = GetXPosition(k);
  const double y = GetYPosition(k);

  double theta;
  if(director_config == "defect-pair")
  {
    // Analytic +1/2 / -1/2 winding, textbook single pair:
    //   theta = angle + .5*atan2(y-y+, x-x+) - .5*atan2(y-y-, x-x-),
    // +1/2 on the left, -1/2 on the right, separated along x by defect_sep. Q is
    // single-valued off the two cores because the atan2 branch cuts are 2pi jumps in
    // 2*theta. Use PLAIN (non-wrapped) displacements: minimum-imaging them puts a jump on
    // the half-box line that injects a spurious extra pair.
    const double xpl = LX/2.0 - defect_sep/2.0, ypl = LY/2.0;
    const double xmi = LX/2.0 + defect_sep/2.0, ymi = LY/2.0;
    theta = angle + .5*atan2(y-ypl, x-xpl) - .5*atan2(y-ymi, x-xmi)
          + noise*M_PI*(random_real()-.5);
  }
  else
    theta = angle + noise*M_PI*(random_real()-.5);

  QQxx[k] = init_order*cos(2*theta);
  QQyx[k] = init_order*sin(2*theta);

  // uniform density at rest: the layer is confluent, so this is the whole domain
  ux[k] = uy[k] = ux_mat[k] = uy_mat[k] = 0;
  n[k]  = rho;
  ff[k] = GetEquilibriumDistribution(ux[k], uy[k], n[k]);
  // total for the conservation check
  ftot  = accumulate(begin(ff[k]), end(ff[k]), ftot);
}

void ConfluentWet::Configure()
{
  Nematic::Configure();
  ConfigurePhenotype();
}

void ConfluentWet::ConfigurePhenotype()
{
  for(unsigned k=0; k<DomainSize; ++k)
    m[k] = m0;

  if(chi_config == "uniform" || chi_noise == 0.)
  {
    for(unsigned k=0; k<DomainSize; ++k)
      chi[k] = chi0;
    return;
  }

  // White noise smoothed by diffusion into a correlated field of the requested length,
  // then rescaled to the requested standard deviation. chi_tmp is borrowed as scratch; it
  // holds no state before the first update.
  for(unsigned k=0; k<DomainSize; ++k)
    chi[k] = random_normal();
  // Go through the virtual boundary handler, never Apply*() directly: bc=0 has no boundary
  // layer at all, so ApplyPBC() would index past the field.
  BoundaryConditionsFields();

  constexpr double smooth_rate = 0.2;
  const unsigned smooth_steps =
    chi_length > 0 ? static_cast<unsigned>(ceil(chi_length*chi_length/(2*smooth_rate))) : 0u;

  for(unsigned s=0; s<smooth_steps; ++s)
  {
    for(unsigned k=0; k<DomainSize; ++k)
      chi_tmp[k] = chi[k] + smooth_rate*laplacian(chi, get_neighbours(k), sD);
    swap(chi.get_data(), chi_tmp.get_data());
    BoundaryConditionsFields();
  }

  double mean = 0., var = 0.;
  for(unsigned k=0; k<DomainSize; ++k) mean += chi[k];
  mean /= DomainSize;
  for(unsigned k=0; k<DomainSize; ++k) var += (chi[k]-mean)*(chi[k]-mean);
  var = sqrt(var/DomainSize);

  const double scale = var > 0 ? chi_noise/var : 0.;
  for(unsigned k=0; k<DomainSize; ++k)
  {
    double v = chi0 + scale*(chi[k]-mean);
    chi[k] = v < 0 ? 0. : (v > 1 ? 1. : v);
  }
}

double ConfluentWet::ChiStar(double mm) const
{
  return .5*(1. + switch_sign*tanh((mm - mc)/chi_width));
}

double ConfluentWet::MemoryTarget(unsigned k) const
{
  const double P = pressure[k];
  return pmem_width > 0
    ? .5*(1. + tanh((P - pmem)/pmem_width))
    : (P > pmem ? 1. : 0.);
}

void ConfluentWet::ClampUnit(ScalarField& fld, unsigned k)
{
  if(fld[k] < 0.)      fld[k] = 0.;
  else if(fld[k] > 1.) fld[k] = 1.;
}

void ConfluentWet::UpdateNodeQuantities(unsigned k)
{
  const auto& d = get_neighbours(k);
  const auto& f = ff[k];
  const double Qxx = QQxx[k];
  const double Qyx = QQyx[k];

  // ---- moments of the distribution ----------------------------------------
  const double nn = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];
  const double vx = (f[1] - f[2] + f[5] - f[6] - f[7] + f[8])/nn;
  const double vy = (f[3] - f[4] + f[5] - f[6] + f[7] - f[8])/nn;

  const double del2Qxx = laplacian(QQxx, d, sD);
  const double dxQxx   = derivX   (QQxx, d, sB);
  const double dyQxx   = derivY   (QQxx, d, sB);
  const double del2Qyx = laplacian(QQyx, d, sD);
  const double dxQyx   = derivX   (QQyx, d, sB);
  const double dyQyx   = derivY   (QQyx, d, sB);

  // ---- free energy and molecular field ------------------------------------
  // f_bulk = .5*CC*term^2  =>  H = -df/dQ = 2*CC*term*Q. The factor 2 is what makes the
  // molecular field consistent with the bulk stress below (and matches confluent-memory).
  const double term = 1. - Qxx*Qxx - Qyx*Qyx;
  const double Hxx  = 2*CC*term*Qxx + LL*del2Qxx;
  const double Hyx  = 2*CC*term*Qyx + LL*del2Qyx;

  // ---- stress -------------------------------------------------------------
  // sigmaB carries the entire trace; sigmaF/S/A are exactly traceless. With no phase field
  // there is no Gibbs-Duhem -mu*phi term, so the isotropic thermodynamic stress reduces to
  // sigma_iso = f - n df/dn = f_bulk itself.
  const double sigmaB = .5*CC*term*term;
  const double active = zeta*(1. - chi[k]);
  const double sigmaF = 2*xi*( (Qxx*Qxx-1.)*Hxx + Qxx*Qyx*Hyx )
    - active*Qxx
    + LL*(dyQxx*dyQxx + dyQyx*dyQyx - dxQxx*dxQxx - dxQyx*dxQyx);
  const double sigmaS = 2*xi*( Qyx*Qxx*Hxx + (Qyx*Qyx-1.)*Hyx )
    - active*Qyx
    - 2*LL*(dxQxx*dyQxx + dxQyx*dyQyx);
  const double sigmaA = 2*(Qxx*Hyx - Qyx*Hxx);

  n[k]       =  nn;
  ux[k]      =  vx;
  uy[k]      =  vy;
  HHxx[k]    =  Hxx;
  HHyx[k]    =  Hyx;
  dxQQxx[k]  =  dxQxx;
  dxQQyx[k]  =  dxQyx;
  dyQQxx[k]  =  dyQxx;
  dyQQyx[k]  =  dyQyx;
  sigmaXX[k] =  sigmaF + sigmaB;
  sigmaYY[k] = -sigmaF + sigmaB;
  sigmaXY[k] =  sigmaS + sigmaA;
  sigmaYX[k] =  sigmaS - sigmaA;

  // ---- diagnostics --------------------------------------------------------
  sigma_bulk[k]      = sigmaB;
  sigma_active_xx[k] = -active*Qxx;
  sigma_active_yx[k] = -active*Qyx;
  // Pi_ij = -p_LB delta_ij + sigma_ij + viscous, so the mechanical pressure is
  // -1/2 Tr(Pi) = p_LB - sigma_bulk. Referenced to the resting state (n = rho, term = 0).
  pressure[k]    = (nn - rho)/3. - sigmaB;
  pressure_lb[k] = (nn - rho)/3.;
}

void ConfluentWet::UpdateQuantities()
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    UpdateNodeQuantities(k);
}

void ConfluentWet::ComputeNodeMaterialVelocity(unsigned k)
{
  const auto& d = get_neighbours(k);
  const double nn = n[k];

  const double divSx = derivX(sigmaXX, d, sB) + derivY(sigmaXY, d, sB);
  const double divSy = derivX(sigmaYX, d, sB) + derivY(sigmaYY, d, sB);

  // A forced LB carries the force trapezoidally in time, so the physical velocity is
  // u = (Sum_v f_v c_v + F/2)/n with F = div(sigma) - friction*u. Friction makes that
  // implicit, but it closes in one line:
  //     u*(1 + friction/(2n)) = u_code + div(sigma)/(2n).
  const double denom = 1. + friction/(2*nn);
  ux_mat[k] = (ux[k] + divSx/(2*nn))/denom;
  uy_mat[k] = (uy[k] + divSy/(2*nn))/denom;
}

void ConfluentWet::ComputeMaterialVelocity()
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    ComputeNodeMaterialVelocity(k);
}

void ConfluentWet::UpdateNodeFields(unsigned k, bool first)
{
  const auto& d = get_neighbours(k);

  const double nn  = n[k];
  const double Qxx = QQxx[k];
  const double Qyx = QQyx[k];
  const double Hxx = HHxx[k];
  const double Hyx = HHyx[k];

  // Bare LB moment: used ONLY for the collision (the simple forcing scheme is built
  // around f^eq(u_code)) and for the friction term of the force.
  const double vx = ux[k];
  const double vy = uy[k];
  // Material velocity: used for every physical transport below.
  const double ax = ux_mat[k];
  const double ay = uy_mat[k];

  const double dxux = derivX(ux_mat, d, sB);
  const double dyux = derivY(ux_mat, d, sB);
  const double dxuy = derivX(uy_mat, d, sB);
  const double dyuy = derivY(uy_mat, d, sB);

  const double dxSxx = derivX(sigmaXX, d, sB);
  const double dySxy = derivY(sigmaXY, d, sB);
  const double dxSyx = derivX(sigmaYX, d, sB);
  const double dySyy = derivY(sigmaYY, d, sB);

  const double expansion = dxux + dyuy;
  const double shear     = .5*(dxuy + dyux);
  const double vorticity = .5*(dxuy - dyux);
  const double traceQL   = Qxx*(dxux - dyuy) + 2*Qyx*shear;

  // ---- Beris-Edwards, advected with the MATERIAL velocity -----------------
  const double Dxx = Gamma*Hxx - ax*dxQQxx[k] - ay*dyQQxx[k] - 2*vorticity*Qyx
    + xi*((Qxx+1)*(2*dxux - traceQL) + 2*Qyx*shear - expansion);
  const double Dyx = Gamma*Hyx - ax*dxQQyx[k] - ay*dyQQyx[k] + 2*vorticity*Qxx
    + xi*(Qyx*(expansion - traceQL) + 2*shear);

  // ---- phenotype and memory: intensive scalars, material derivative -------
  // D_t chi = (chi*(m) - chi)/tau_chi + Dbio lap chi
  // D_t m   = (g(P)    - m)/tau_m     + Dbio lap m
  // Centred differences: CFL ~ 4e-3 makes the weak instability of centred advection under
  // Heun stepping (growth ~ CFL^4/8) irrelevant, while upwinding would inject a numerical
  // diffusivity |u|dx/2 that would swamp Dbio. The [0,1] clamp below catches any
  // dispersive overshoot.
  const double chi_k = chi[k];
  const double m_k   = m[k];

  // Dbio is the SAME for both: diffusion here is cell motility, and a cell carries its
  // phenotype and its memory together. The two transport operators are therefore identical
  // by construction and the fields differ only through their sources.
  const double chi_rhs =
    - ax*derivX(chi, d, sB) - ay*derivY(chi, d, sB)
    + Dbio*laplacian(chi, d, sD)
    + (use_switching ? (ChiStar(m_k) - chi_k)/tau_chi : 0.);

  const double m_rhs =
    - ax*derivX(m, d, sB) - ay*derivY(m, d, sB)
    + Dbio*laplacian(m, d, sD)
    + (MemoryTarget(k) - m_k)/tau_m;

  // ---- lattice Boltzmann ---------------------------------------------------
  const double Fx = dxSxx + dySxy - friction*vx;
  const double Fy = dxSyx + dySyy - friction*vy;
  const auto fe = GetEquilibriumDistribution(vx, vy, nn);

  if(first)
  {
    double Qxx_noise = 0., Qxy_noise = 0.;
    LBNode ff_noise = {{0.}};

    if(Q_fluct)
    {
      const double Q_stren = sqrt(Gamma*Q_kBT);
      Qxx_noise = Q_stren*random_real();
      Qxy_noise = Q_stren*random_real();
    }
    if(u_fluct)
      ff_noise = GenerateNoiseDistribution(sqrt(3.*nn*u_kBT*(2.*tau-1.)/tau/tau));

    QNxx[k] = Qxx + .5*Dxx + Qxx_noise;
    QNyx[k] = Qyx + .5*Dyx + Qxy_noise;
    QQxx[k] = QNxx[k] + .5*Dxx;
    QQyx[k] = QNyx[k] + .5*Dyx;

    chiN[k]    = chi_k + .5*chi_rhs;
    chi_tmp[k] = chi_k + chi_rhs;
    mN[k]      = m_k   + .5*m_rhs;
    m_tmp[k]   = m_k   + m_rhs;

    for(unsigned v=0; v<lbq; ++v)
    {
      fn[k][v] = ff[k][v] + .5*((fe[v]-ff[k][v])/tau
          + w[v]*(Fx*xdir(v) + Fy*ydir(v))) + ff_noise[v];
      ff[k][v] = fn[k][v] + .5*((fe[v]-ff[k][v])/tau
          + w[v]*(Fx*xdir(v) + Fy*ydir(v)));
    }
  }
  else
  {
    QQxx[k] = QNxx[k] + .5*Dxx;
    QQyx[k] = QNyx[k] + .5*Dyx;

    chi_tmp[k] = chiN[k] + .5*chi_rhs;
    m_tmp[k]   = mN[k]   + .5*m_rhs;

    for(unsigned v=0; v<lbq; ++v)
      ff[k][v] = fn[k][v] + .5*((fe[v]-ff[k][v])/tau
          + w[v]*(Fx*xdir(v) + Fy*ydir(v)));
  }

  // chi and m are fractions: the continuum equations preserve [0,1] (both relax towards
  // targets in [0,1] and advection/diffusion obey a maximum principle), so anything
  // outside is discretisation error.
  ClampUnit(chi_tmp, k);
  ClampUnit(m_tmp,   k);
}

void ConfluentWet::UpdateFields(bool first)
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    UpdateNodeFields(k, first);

  // chi and m are READ from their neighbours inside the loop (advection, diffusion), so
  // unlike Q -- whose derivatives were cached in UpdateQuantities -- they cannot be
  // updated in place.
  swap(chi.get_data(), chi_tmp.get_data());
  swap(m.get_data(),   m_tmp.get_data());
}

void ConfluentWet::BoundaryConditionsFields()
{
  Nematic::BoundaryConditionsFields();

  // chi and m ARE differentiated (advection, Laplacian), so they need real boundary values.
  switch(BC)
  {
    case 0:
      break;
    case 1:
    {
      auto apply_bc = [](ScalarField& field) {
        field.ApplyPBC(PBCWall::LeftRight);
        field.ApplyNeumann(Wall::Front);
        field.ApplyNeumann(Wall::Back);
        field.ApplyNeumann(Corner::RightBack, Wall::Back);
        field.ApplyNeumann(Corner::RightFront, Wall::Front);
        field.ApplyNeumann(Corner::LeftBack, Wall::Back);
        field.ApplyNeumann(Corner::LeftFront, Wall::Front);
      };
      apply_bc(chi);
      apply_bc(m);
      break;
    }
    default:
      chi.ApplyPBC();
      m.ApplyPBC();
  }
}

void ConfluentWet::BoundaryConditionsFields2()
{
  Nematic::BoundaryConditionsFields2();

  // ux_mat/uy_mat are differentiated in UpdateNodeFields; the rest are pure diagnostics and
  // are handled here only so the boundary layer is not written out as garbage when bc != 0.
  switch(BC)
  {
    case 0:
      break;
    case 1:
    {
      auto apply_bc = [](ScalarField& field) {
        field.ApplyPBC(PBCWall::LeftRight);
        field.ApplyNeumann(Wall::Front);
        field.ApplyNeumann(Wall::Back);
        field.ApplyNeumann(Corner::RightBack, Wall::Back);
        field.ApplyNeumann(Corner::RightFront, Wall::Front);
        field.ApplyNeumann(Corner::LeftBack, Wall::Back);
        field.ApplyNeumann(Corner::LeftFront, Wall::Front);
      };
      apply_bc(ux_mat);
      apply_bc(uy_mat);
      apply_bc(sigma_bulk);
      apply_bc(sigma_active_xx);
      apply_bc(sigma_active_yx);
      apply_bc(pressure);
      apply_bc(pressure_lb);
      break;
    }
    default:
      ux_mat.ApplyPBC();
      uy_mat.ApplyPBC();
      sigma_bulk.ApplyPBC();
      sigma_active_xx.ApplyPBC();
      sigma_active_yx.ApplyPBC();
      pressure.ApplyPBC();
      pressure_lb.ApplyPBC();
  }
}

void ConfluentWet::Step()
{
  BoundaryConditionsFields();
  UpdateQuantities();
  BoundaryConditionsFields2();
  // Needs the sigma halo, hence after BoundaryConditionsFields2; and its own halo is
  // needed by UpdateFields, hence the second call.
  ComputeMaterialVelocity();
  BoundaryConditionsFields2();
  UpdateFields(true);

  BoundaryConditionsLB();
  Move();

  for(unsigned s=1; s<=npc; ++s)
  {
    BoundaryConditionsFields();
    UpdateQuantities();
    BoundaryConditionsFields2();
    ComputeMaterialVelocity();
    BoundaryConditionsFields2();
    UpdateFields(false);
  }
}

void ConfluentWet::RuntimeChecks()
{
  // Conservation tolerance: RELATIVE to the conserved total, with a floor of 1 so small
  // domains keep the historical absolute-1 behaviour of Nematic. An absolute bound is
  // domain-size blind and trips on the O(1e-6) relative per-step flux imbalance of a large
  // lattice (400x400 at rho=40 carries ftot ~ 6.4e6) even though mass is conserved to
  // roundoff; a real leak drifts past the relative bound cumulatively. Same fix as
  // Lyotropic::RuntimeChecks.
  constexpr double conserve_rel_tol = 1e-4;

  double fcheck = 0;
  for(unsigned k=0; k<DomainSize; ++k)
    fcheck = accumulate(begin(ff[k]), end(ff[k]), fcheck);

  if(abs(ftot-fcheck) > max(1.0, conserve_rel_tol*abs(ftot)))
    throw error_msg("f is not conserved (", ftot, "/", fcheck, ")");
}

option_list ConfluentWet::GetOptions()
{
  auto options = Nematic::GetOptions();

  // options[0] is Nematic's "Model options"
  options[0].add_options()
    ("Dbio", opt::value<double>(&Dbio),
     "biomass (cell-motility) diffusivity; applied identically to chi and m")
    ("tau-chi", opt::value<double>(&tau_chi),
     "phenotype relaxation time (<= 0 freezes chi)")
    ("chi-width", opt::value<double>(&chi_width),
     "width of the tanh in chi*(m)")
    ("mc", opt::value<double>(&mc),
     "memory threshold in chi*(m)")
    ("switch-sign", opt::value<int>(&switch_sign),
     "+1: compression drives chi up (contact inhibition); -1: the opposite branch")
    ("tau-m", opt::value<double>(&tau_m),
     "memory relaxation time")
    ("pmem", opt::value<double>(&pmem),
     "pressure threshold of the memory source g(P)")
    ("pmem-width", opt::value<double>(&pmem_width),
     "smoothing width of g(P); 0 is a sharp step");

  // options[1] is Nematic's "Initial configuration options"
  options[1].add_options()
    ("director-config", opt::value<string>(&director_config),
     "director init: uniform (angle+noise) or defect-pair")
    ("defect-sep", opt::value<double>(&defect_sep),
     "separation of the +/- 1/2 defect pair (lattice units, defect-pair mode)")
    ("initial-order", opt::value<double>(&init_order),
     "initial nematic order amplitude")
    ("chi-config", opt::value<string>(&chi_config),
     "phenotype init: uniform or noise")
    ("chi0", opt::value<double>(&chi0),
     "initial mean phenotype")
    ("chi-noise", opt::value<double>(&chi_noise),
     "standard deviation of the initial phenotype noise")
    ("chi-length", opt::value<double>(&chi_length),
     "correlation length of the initial phenotype noise")
    ("m0", opt::value<double>(&m0),
     "initial memory");

  return options;
}
