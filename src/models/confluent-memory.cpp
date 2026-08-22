#include "header.hpp"
#include "models/confluent-memory.hpp"
#include "error_msg.hpp"
#include "serialization.hpp"
#include "lb.hpp"
#include "random.hpp"

using namespace std;
namespace opt = boost::program_options;

// from main.cpp:
extern unsigned nthreads;

ConfluentMemory::ConfluentMemory(unsigned LX_, unsigned LY_, unsigned BC_)
  : LyotropicWithDivision(LX_, LY_, BC_)
{
  // The inherited Lyotropic parameters have no in-class initialisers, so give
  // them defaults here (the constructor runs before the options are parsed).
  // AA and KK default to zero and MUST stay there: this model has no phase field.
  AA = 0.; KK = 0.;
  CC = 0.1; LL = 0.5; xi = 1.;
  GammaP = 0.05; GammaQ = 0.03; friction = 100.;
  zeta = 0.; zetaI = 0.;
  tauNem = tauIso = 1.;          // unused: there is no lattice Boltzmann step
  conc = 1.1; level = 0.; radius = 0.;
  angle_deg = 0.; noise = 1.;
}

void ConfluentMemory::Initialize()
{
  LyotropicWithDivision::Initialize();

  // ---- the confluent contract, enforced loudly rather than silently ---------
  if(init_config != "full")
    throw error_msg("confluent-memory is a confluent model: it requires "
                    "config=full (no free boundary, no vacuum), got '",
                    init_config, "'.");
  if(BC != 0)
    throw error_msg("confluent-memory requires bc=0 (periodic, no boundary "
                    "layer): a confluent layer has no boundary.");
  if(AA != 0. || KK != 0.)
    throw error_msg("confluent-memory has no phase field: AA and KK must be 0 "
                    "(got AA=", AA, ", KK=", KK, "). The double well would tear "
                    "the layer into holes.");
  if(friction <= 0.)
    throw error_msg("confluent-memory is a dry model: it requires friction > 0.");
  if(B <= 0.)
    throw error_msg("confluent-memory needs B > 0: the crowding modulus is the "
                    "only thing holding the layer together and the only source "
                    "of mechanical pressure.");
  if(tau_m <= 0.)
    throw error_msg("confluent-memory requires tau-m > 0.");
  if(chi_width <= 0.)
    throw error_msg("chi-width must be positive (it is the width of the tanh in "
                    "chi*(m); a hard switch is chi-width -> 0, not 0).");
  if(switch_sign != 1 && switch_sign != -1)
    throw error_msg("switch-sign must be +1 (compression -> passive) or -1.");
  if(pmem_width < 0.)
    throw error_msg("pmem-width must be non-negative.");
  if(m0 < 0. || m0 > 1.)
    throw error_msg("m0 must lie in [0,1].");
  if(chi0 < 0. || chi0 > 1.)
    throw error_msg("chi0 must lie in [0,1].");
  if(phi_floor <= 0.)
    throw error_msg("phi-floor must be positive.");
  // Silent-death trap: SetMasks() weights the division rate by
  // (1 - phi/phi_critical), which is NEGATIVE in a layer denser than
  // phi_critical, so no division event is ever drawn.
  if(division_phi_critical_factor && conc > phi_critical)
    throw error_msg("division-phi-critical-factor=1 with conc=", conc,
                    " > phi-critical=", phi_critical, " silently disables all "
                    "division (the crowding factor 1-phi/phi_critical is "
                    "negative). Set division-phi-critical-factor=0.");

  // tau_chi <= 0 is the frozen-chi control, not an error.
  use_switching = (tau_chi > 0.);

  chi.SetSize(LX, LY, Type);
  phichi.SetSize(LX, LY, Type);
  phichiN.SetSize(LX, LY, Type);
  phichi_tmp.SetSize(LX, LY, Type);
  m.SetSize(LX, LY, Type);
  phim.SetSize(LX, LY, Type);
  phimN.SetSize(LX, LY, Type);
  phim_tmp.SetSize(LX, LY, Type);
}

void ConfluentMemory::Configure()
{
  // Uniform density conc plus the director specified by director-config; the
  // base routine already rejects anything but config=full for defect-pair.
  LyotropicWithDivision::Configure();

  if(init_frame.empty())
    ConfigurePhenotype();
  else
    ConfigureFromFrame();
}

void ConfluentMemory::ConfigurePhenotype()
{
  if(chi_config == "uniform" || chi_noise == 0.)
  {
    for(unsigned k=0; k<DomainSize; ++k)
      chi[k] = chi0;
  }
  else if(chi_config == "noise")
  {
    // White noise smoothed by diffusion into a correlated field of the requested
    // length, then rescaled to the requested standard deviation. phichi_tmp is
    // borrowed as scratch here; it holds no state before the first update.
    for(unsigned k=0; k<DomainSize; ++k)
      chi[k] = random_normal();
    // Go through the virtual boundary handler, never Apply*() directly: bc=0 has
    // no boundary layer at all, so ApplyPBC() would index past the field.
    BoundaryConditionsFields();

    constexpr double smooth_rate = 0.2;
    const unsigned smooth_steps =
      chi_length > 0 ? static_cast<unsigned>(ceil(chi_length*chi_length/(2*smooth_rate))) : 0u;

    for(unsigned s=0; s<smooth_steps; ++s)
    {
      for(unsigned k=0; k<DomainSize; ++k)
      {
        const auto& d = get_neighbours(k);
        phichi_tmp[k] = chi[k] + smooth_rate*laplacian(chi, d, sD);
      }
      swap(chi.get_data(), phichi_tmp.get_data());
      BoundaryConditionsFields();
    }

    double mean = 0.;
    for(unsigned k=0; k<DomainSize; ++k) mean += chi[k];
    mean /= DomainSize;
    double var = 0.;
    for(unsigned k=0; k<DomainSize; ++k) var += (chi[k]-mean)*(chi[k]-mean);
    var = sqrt(var/DomainSize);

    const double scale = var > 0 ? chi_noise/var : 0.;
    for(unsigned k=0; k<DomainSize; ++k)
      chi[k] = min(1., max(0., chi0 + scale*(chi[k]-mean)));
  }
  else
    throw error_msg("error: chi configuration '", chi_config, "' unknown.");

  for(unsigned k=0; k<DomainSize; ++k)
  {
    phichi[k] = phi[k]*chi[k];
    phim[k]   = phi[k]*m0;
  }

  ProjectDensity(phichi, phi);
  ProjectDensity(phim, phi);
  UpdateDerivedFields();
}

void ConfluentMemory::ConfigureFromFrame()
{
  const string content = read_file_to_string(init_frame);

  auto load_into = [&](const string& nm, ScalarField& fld)
  {
    const auto vals = read_frame_field(content, nm);
    if(vals.size() != DomainSize)
      throw error_msg("init-frame field '", nm, "' has ", vals.size(),
                      " values but this run has DomainSize=", DomainSize,
                      " (grid-size mismatch between snapshot and run).");
    auto& data = fld.get_data();
    for(unsigned k=0; k<DomainSize; ++k)
      data[k] = vals[k];
  };

  load_into("QQxx", QQxx);
  load_into("QQyx", QQyx);
  load_into("phi",  phi);

  // The phenotype either carries over, or is reset uniform for the homogeneous
  // control arm; the memory always carries over (it is a property of the cells).
  if(init_frame_uniform_chi)
    for(unsigned k=0; k<DomainSize; ++k)
      phichi[k] = phi[k]*chi0;
  else
    load_into("phichi", phichi);
  load_into("phim", phim);

  ProjectDensity(phichi, phi);
  ProjectDensity(phim, phi);
  UpdateDerivedFields();

  // Retarget the mass reference used by conserve_phi and the diagnostics.
  double total = 0.;
  for(unsigned k=0; k<DomainSize; ++k)
    total += phi[k];
  totalphi = countphi = ptot = total;
}

void ConfluentMemory::UpdateDerivedFields()
{
  for(unsigned k=0; k<DomainSize; ++k)
  {
    // No vacuum in a confluent layer, so no masking or neighbour fallback: the
    // division is unconditional. UpdateQuantities() guards phi against reaching
    // the floor, so this cannot divide by (near) zero unnoticed.
    const double p = phi[k] > phi_floor ? phi[k] : phi_floor;
    chi[k] = phichi[k]/p;
    m[k]   = phim[k]/p;
  }

  BoundaryConditionsFields();
}

double ConfluentMemory::SwitchingSource(double chi_eff, double m_eff) const
{
  // tau_chi * D_t chi = chi*(m) - chi, with the target
  //   chi*(m) = .5*(1 + switch_sign*tanh((m - mc)/chi_width))
  // Because chi* lies strictly inside (0,1) and is an attractor, chi cannot leave
  // [0,1] on its own; the [0,phi] clamp in the update is only a round-off guard,
  // not -- as in the earlier linear form -- what sets the steady state.
  const double chi_star = .5*(1 + switch_sign*tanh((m_eff - mc)/chi_width));
  return (chi_star - chi_eff)/tau_chi;
}

void ConfluentMemory::ProjectDensity(ScalarField& fld, const ScalarField& density)
{
  for(unsigned k=0; k<DomainSize; ++k)
  {
    const double upper = density[k] > 0 ? density[k] : 0.;
    if(fld[k] < 0)
      fld[k] = 0;
    else if(fld[k] > upper)
      fld[k] = upper;
  }
}

double ConfluentMemory::MemorySource(unsigned k, double m_eff) const
{
  // With AA = KK = 0 there is no surface stress, so sigma_bulk is already the
  // surface-free half-trace and P = -sigma_bulk is the crowding pressure.
  const double press = -sigma_bulk[k];
  const double g = pmem_width > 0
    ? .5*(1 + tanh((press - pmem)/pmem_width))
    : (press > pmem ? 1. : 0.);
  return phi[k]*(g - m_eff)/tau_m;
}

double ConfluentMemory::UpwindFace(const ScalarField& fld, unsigned k,
                                   unsigned neighbour, double dphi_k) const
{
  // dphi_k > 0 means material enters node k from the neighbour, so the incoming
  // parcel carries the neighbour's value; otherwise it carries this node's.
  return dphi_k >= 0 ? fld[neighbour] : fld[k];
}

void ConfluentMemory::UpdateNodeQuantities(unsigned k)
{
  const auto& d = get_neighbours(k);
  const double Qxx = QQxx[k];
  const double Qyx = QQyx[k];
  const double p   = phi[k];

  const double del2Qxx = laplacian(QQxx, d, sD);
  const double dxQxx   = derivX   (QQxx, d, sB);
  const double dyQxx   = derivY   (QQxx, d, sB);
  const double del2Qyx = laplacian(QQyx, d, sD);
  const double dxQyx   = derivX   (QQyx, d, sB);
  const double dyQyx   = derivY   (QQyx, d, sB);

  // ---- free energy: crowding (one-sided) + density-coupled nematic bulk -----
  const double q2   = Qxx*Qxx + Qyx*Qyx;
  const double term = Snem*p - q2;
  const double Hxx  = 2*CC*term*Qxx + LL*del2Qxx;
  const double Hyx  = 2*CC*term*Qyx + LL*del2Qyx;

  const double dp        = p - phi_critical;
  const double f_crowd   = dp > 0 ? .5*B*dp*dp : 0.;
  const double mu_crowd  = dp > 0 ? B*dp       : 0.;
  const double mu        = mu_crowd + CC*Snem*term;

  // ---- stress: isotropic (crowding + nematic bulk) and deviatoric ----------
  const double sigmaB = f_crowd + .5*CC*term*term
    - (mu_crowd + CC*Snem*term)*p
    + zetaI*(conc - p);
  const double active = zeta*p*(1 - chi[k]);
  const double sigmaF = 2*xi*( (Qxx*Qxx-1)*Hxx + Qxx*Qyx*Hyx )
    - active*Qxx
    + LL*(dyQxx*dyQxx + dyQyx*dyQyx - dxQxx*dxQxx - dxQyx*dxQyx);
  const double sigmaS = 2*xi*( Qyx*Qxx*Hxx + (Qyx*Qyx-1)*Hyx )
    - active*Qyx
    - 2*LL*(dxQxx*dyQxx + dxQyx*dyQyx);
  const double sigmaA = 2*(Qxx*Hyx - Qyx*Hxx);

  HHxx[k]    = Hxx;
  HHyx[k]    = Hyx;
  MU[k]      = mu;
  dxQQxx[k]  = dxQxx;
  dxQQyx[k]  = dxQyx;
  dyQQxx[k]  = dyQxx;
  dyQQyx[k]  = dyQyx;
  sigmaXX[k] =  sigmaF + sigmaB;
  sigmaYY[k] = -sigmaF + sigmaB;
  sigmaXY[k] =  sigmaS + sigmaA;
  sigmaYX[k] =  sigmaS - sigmaA;
  // sigmaF is traceless, so the half-trace is sigmaB alone; with no phase field
  // this is already the surface-free pressure that drives the memory.
  sigma_bulk[k]      = sigmaB;
  sigma_active_xx[k] = -active*Qxx;
  sigma_active_yx[k] = -active*Qyx;
}

void ConfluentMemory::UpdateQuantities()
{
  double sum = 0.;
  double pmin = phi[0];

  #pragma omp parallel for reduction(+:sum) reduction(min:pmin) \
                           num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
  {
    sum += phi[k];
    if(phi[k] < pmin) pmin = phi[k];
    UpdateNodeQuantities(k);
  }

  countphi = sum;
  phi_min  = pmin;

  // A confluent layer that has opened a hole is no longer the model we solved.
  // Stop rather than let chi = phichi/phi and m = phim/phi silently blow up.
  if(phi_min < phi_floor)
    throw error_msg("confluent-memory: phi fell to ", phi_min, " (below phi-floor=",
                    phi_floor, "). The layer is no longer confluent -- raise B, "
                    "lower zeta or alpha, or reconsider the parameters.");
}

void ConfluentMemory::ComputeNodeVelocity(unsigned k)
{
  const auto& d = get_neighbours(k);

  const double dxSxx = derivX(sigmaXX, d, sB);
  const double dySxy = derivY(sigmaXY, d, sB);
  const double dxSyx = derivX(sigmaYX, d, sB);
  const double dySyy = derivY(sigmaYY, d, sB);

  // Overdamped force balance: friction*u = div(sigma). No lattice Boltzmann.
  const double vx = (dxSxx + dySxy)/friction;
  const double vy = (dxSyx + dySyy)/friction;

  ux[k]     = vx;
  uy[k]     = vy;
  ux_phi[k] = vx*phi[k];
  uy_phi[k] = vy*phi[k];
}

void ConfluentMemory::ComputeVelocity()
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    ComputeNodeVelocity(k);
}

void ConfluentMemory::UpdateNodeFields(unsigned k, bool first)
{
  const auto& d = get_neighbours(k);

  const double Qxx = QQxx[k];
  const double Qyx = QQyx[k];
  const double p   = phi[k];
  const double vx  = ux[k];
  const double vy  = uy[k];

  const double dxux = derivX(ux, d, sB);
  const double dyux = derivY(ux, d, sB);
  const double dxuy = derivX(uy, d, sB);
  const double dyuy = derivY(uy, d, sB);

  const double expansion = dxux + dyuy;
  const double shear     = .5*(dxuy + dyux);
  const double vorticity = .5*(dxuy - dyux);
  const double traceQL   = Qxx*(dxux - dyuy) + 2*Qyx*shear;

  // ---- Beris-Edwards ------------------------------------------------------
  const double Dxx = GammaQ*HHxx[k] - vx*dxQQxx[k] - vy*dyQQxx[k] - 2*vorticity*Qyx
    + xi*((Qxx+1)*(2*dxux - traceQL) + 2*Qyx*shear - expansion);
  const double Dyx = GammaQ*HHyx[k] - vx*dxQQyx[k] - vy*dyQQyx[k] + 2*vorticity*Qxx
    + xi*(Qyx*(expansion - traceQL) + 2*shear);

  // ---- the shared face flux J = phi*u - GammaP*grad(mu) --------------------
  // The face momentum is the pair-symmetric .5*(u_phi[k]+u_phi[n]), which makes
  // each per-face increment antisymmetric (dphi_px@k = -dphi_mx@n). That is what
  // lets the upwind weighting below conserve Sum(phichi) and Sum(phim) exactly;
  // the extra .5*u_phi[k] terms cancel over the four faces, so phi is unaffected.
  const double dphi_px = GammaP*(MU[d[1]] - MU[k]) - .5*(ux_phi[k] + ux_phi[d[1]]);
  const double dphi_mx = GammaP*(MU[d[2]] - MU[k]) + .5*(ux_phi[k] + ux_phi[d[2]]);
  const double dphi_py = GammaP*(MU[d[3]] - MU[k]) - .5*(uy_phi[k] + uy_phi[d[3]]);
  const double dphi_my = GammaP*(MU[d[4]] - MU[k]) + .5*(uy_phi[k] + uy_phi[d[4]]);

  const double drift = conserve_phi ? (countphi - totalphi)/DomainSize : 0.;

  const double chi_eff = chi[k];
  const double m_eff   = m[k];

  const double growth_rate = alpha*(division_mask[k] ? 1. : 0.)
                           - beta *(death_mask[k]    ? 1. : 0.);
  const double R = chi_eff*p*growth_rate;

  const double Dphi = dphi_px + dphi_mx + dphi_py + dphi_my - drift + R;

  // ---- phenotype ----------------------------------------------------------
  const double phichiTransport =
    UpwindFace(chi, k, d[1], dphi_px)*dphi_px
  + UpwindFace(chi, k, d[2], dphi_mx)*dphi_mx
  + UpwindFace(chi, k, d[3], dphi_py)*dphi_py
  + UpwindFace(chi, k, d[4], dphi_my)*dphi_my
  - chi_eff*drift;

  const double phenotypeDiffusion =
    .5*(phi[d[1]] + p)*(chi[d[1]] - chi_eff)
  - .5*(p + phi[d[2]])*(chi_eff - chi[d[2]])
  + .5*(phi[d[3]] + p)*(chi[d[3]] - chi_eff)
  - .5*(p + phi[d[4]])*(chi_eff - chi[d[4]]);

  // Switching is driven by the MEMORY, and is a relaxation towards chi*(m) with the
  // exact timescale tau_chi -- the mirror of the memory's own relaxation towards g(P).
  const double Sswitch = use_switching ? p*SwitchingSource(chi_eff, m_eff) : 0.;

  const double Dphichi = phichiTransport + Dchi*phenotypeDiffusion + Sswitch
                       + (growTogether ? chi_eff*R : R);

  // ---- mechanical memory --------------------------------------------------
  const double phimTransport =
    UpwindFace(m, k, d[1], dphi_px)*dphi_px
  + UpwindFace(m, k, d[2], dphi_mx)*dphi_mx
  + UpwindFace(m, k, d[3], dphi_py)*dphi_py
  + UpwindFace(m, k, d[4], dphi_my)*dphi_my
  - m_eff*drift;

  // m_eff*R: biomass created by division inherits the memory of its mother.
  const double Dphim = phimTransport + MemorySource(k, m_eff) + m_eff*R;

  // ---- predictor / corrector ----------------------------------------------
  if(first)
  {
    QNxx[k] = Qxx + .5*Dxx;
    QNyx[k] = Qyx + .5*Dyx;
    phn[k]  = p   + .5*Dphi;

    QQxx[k]    = Qxx + Dxx;
    QQyx[k]    = Qyx + Dyx;
    phi_tmp[k] = p   + Dphi;

    phichiN[k]   = phichi[k] + .5*Dphichi;
    phichi_tmp[k] = phichi[k] + Dphichi;
    phimN[k]     = phim[k] + .5*Dphim;
    phim_tmp[k]  = phim[k] + Dphim;
  }
  else
  {
    QQxx[k]    = QNxx[k] + .5*Dxx;
    QQyx[k]    = QNyx[k] + .5*Dyx;
    phi_tmp[k] = phn[k]  + .5*Dphi;

    phichi_tmp[k] = phichiN[k] + .5*Dphichi;
    phim_tmp[k]   = phimN[k]   + .5*Dphim;
  }

  // Both carried densities obey 0 <= . <= phi, since chi and m lie in [0,1].
  const double upper = phi_tmp[k] > 0 ? phi_tmp[k] : 0.;
  if(phichi_tmp[k] < 0)          phichi_tmp[k] = 0;
  else if(phichi_tmp[k] > upper) phichi_tmp[k] = upper;
  if(phim_tmp[k] < 0)            phim_tmp[k] = 0;
  else if(phim_tmp[k] > upper)   phim_tmp[k] = upper;
}

void ConfluentMemory::UpdateFields(bool first)
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    UpdateNodeFields(k, first);

  swap(phi.get_data(),    phi_tmp.get_data());
  swap(phichi.get_data(), phichi_tmp.get_data());
  swap(phim.get_data(),   phim_tmp.get_data());
  UpdateDerivedFields();
}

void ConfluentMemory::BoundaryConditionsFields()
{
  LyotropicWithDivision::BoundaryConditionsFields();

  // bc=0 (enforced in Initialize) has no boundary layer, so there is nothing to
  // do; the switch is kept so the model still behaves if bc is ever relaxed.
  switch(BC)
  {
    case 0:
      break;
    case 1:
    case 2:
      phichi.ApplyNeumannChannel(); chi.ApplyNeumannChannel();
      phim  .ApplyNeumannChannel(); m  .ApplyNeumannChannel();
      break;
    case 3:
    case 4:
      phichi.ApplyNeumann(); chi.ApplyNeumann();
      phim  .ApplyNeumann(); m  .ApplyNeumann();
      break;
    default:
      phichi.ApplyPBC(); chi.ApplyPBC();
      phim  .ApplyPBC(); m  .ApplyPBC();
  }
}

void ConfluentMemory::Step()
{
  SetMasks();

  BoundaryConditionsFields();
  UpdateQuantities();
  Lyotropic::BoundaryConditionsFields2();
  ComputeVelocity();
  Lyotropic::BoundaryConditionsFields2();
  UpdateFields(true);

  for(unsigned n=1; n<=npc; ++n)
  {
    BoundaryConditionsFields();
    UpdateQuantities();
    Lyotropic::BoundaryConditionsFields2();
    ComputeVelocity();
    Lyotropic::BoundaryConditionsFields2();
    UpdateFields(false);
  }
}

option_list ConfluentMemory::GetOptions()
{
  auto options = LyotropicWithDivision::GetOptions();

  options[0].add_options()
    ("B", opt::value<double>(&B),
     "crowding modulus (one-sided, active above phi-critical)")
    ("Snem", opt::value<double>(&Snem),
     "preferred value of one half Tr(Q^2) per unit density")
    ("Dchi", opt::value<double>(&Dchi),
     "phenotype diffusion coefficient")
    ("tau-chi", opt::value<double>(&tau_chi),
     "phenotype relaxation time; tau_chi*D_t chi = chi*(m) - chi. "
     "<= 0 disables switching (frozen-chi control)")
    ("chi-width", opt::value<double>(&chi_width),
     "width of the tanh in chi*(m); small = hard switch, large = linear response")
    ("switch-sign", opt::value<int>(&switch_sign),
     "+1: compression drives chi towards passive/grow (contact inhibition); -1: opposite")
    ("mc", opt::value<double>(&mc),
     "memory threshold in chi*(m)")
    ("growTogether", opt::value<int>(&growTogether),
     "phi*chi growth source: 0 uses R, 1 uses chi*R (chi-preserving)")
    ("chi-config", opt::value<string>(&chi_config),
     "phenotype initialization mode: noise or uniform")
    ("chi0", opt::value<double>(&chi0),
     "initial phenotype mean")
    ("chi-noise", opt::value<double>(&chi_noise),
     "initial phenotype standard deviation")
    ("chi-length", opt::value<double>(&chi_length),
     "initial phenotype correlation length")
    ("tau-m", opt::value<double>(&tau_m),
     "relaxation time of the mechanical memory m")
    ("pmem", opt::value<double>(&pmem),
     "pressure threshold of the memory response g(P)")
    ("pmem-width", opt::value<double>(&pmem_width),
     "smoothing width of g(P); 0 gives a sharp step")
    ("m0", opt::value<double>(&m0),
     "initial mechanical memory")
    ("phi-floor", opt::value<double>(&phi_floor),
     "density below which the layer is declared no longer confluent and the run stops")
    ("init-frame", opt::value<string>(&init_frame),
     "seed Q, phi, phichi and phim from this frame*.json snapshot")
    ("init-frame-uniform-chi", opt::value<int>(&init_frame_uniform_chi),
     "with init-frame, set chi uniform (=chi0) instead of loading it");

  return options;
}
