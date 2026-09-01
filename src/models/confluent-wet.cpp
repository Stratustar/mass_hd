#include <algorithm>

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
extern unsigned nsubsteps;

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

ConfluentWet::~ConfluentWet()
{
  // The video streams are appended to for the whole run and are the only file handles this
  // model owns; closing them here means a run that ends by exception still leaves a stream
  // whose last complete frame is readable.
  if(video_open)
  {
    vid_u.close(); vid_p.close(); vid_m.close(); vid_chi.close(); vid_meta.close();
  }
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
  if(chi_config != "uniform" && chi_config != "noise" &&
     chi_config != "stripe"  && chi_config != "blocks" &&
     chi_config != "binary-noise")
    throw error_msg("confluent-wet: chi-config must be 'uniform', 'noise', 'stripe', "
                    "'blocks' or 'binary-noise', got '", chi_config, "'.");
  if(chi_config == "binary-noise")
  {
    if(chi_length <= 0.)
      throw error_msg("confluent-wet: chi-config=binary-noise needs chi-length > 0, the "
                      "correlation length of the pattern; at 0 the median split would be "
                      "applied to white noise and every node would be its own domain.");
    if(DomainSize % 2)
      throw error_msg("confluent-wet: chi-config=binary-noise splits the nodes exactly "
                      "half and half at the median and needs an even node count, got ",
                      DomainSize, ".");
  }
  if(chi0 < 0. || chi0 > 1.)
    throw error_msg("confluent-wet: chi0 must lie in [0,1].");
  if(chi_lo < 0. || chi_lo > 1. || chi_hi < 0. || chi_hi > 1.)
    throw error_msg("confluent-wet: chi-lo and chi-hi must lie in [0,1].");
  if(m_lo < 0. || m_lo > 1. || m_hi < 0. || m_hi > 1.)
    throw error_msg("confluent-wet: m-lo and m-hi must lie in [0,1].");
  if(chi_config == "stripe" && LX % 2)
    throw error_msg("confluent-wet: chi-config=stripe splits the box along x into two "
                    "equal halves and needs an even LX, got LX=", LX, ".");
  if(chi_config == "blocks")
  {
    if(chi_block == 0)
      throw error_msg("confluent-wet: chi-config=blocks needs chi-block > 0.");
    if(LX % chi_block || LY % chi_block)
      throw error_msg("confluent-wet: chi-block=", chi_block, " must divide both LX=", LX,
                      " and LY=", LY, "; a partial block at the edge would break the "
                      "exactly-half area fraction the mixed initial condition relies on.");
    if((LX/chi_block)*(LY/chi_block) % 2)
      throw error_msg("confluent-wet: chi-config=blocks needs an EVEN number of blocks so "
                      "the two phases can be exactly half and half; LX/chi-block * "
                      "LY/chi-block = ", (LX/chi_block)*(LY/chi_block), " is odd.");
  }
  // The activity law and the open-loop override.
  if(zeta0_frac < 0. || zeta0_frac > 1.)
    throw error_msg("confluent-wet: zeta0-frac is the activity floor as a FRACTION of "
                    "zeta and must lie in [0,1], got ", zeta0_frac, ".");
  if(open_loop != 0 && open_loop != 1)
    throw error_msg("confluent-wet: open-loop must be 0 or 1.");
  if(open_loop && zeta_open < 0.)
    throw error_msg("confluent-wet: open-loop needs zeta-open >= 0, got ", zeta_open, ".");
  if(!open_loop && zeta_open != 0.)
    throw error_msg("confluent-wet: zeta-open=", zeta_open, " is set but open-loop=0, so "
                    "it would be silently ignored and the run would be at a different "
                    "activity than its runcard reads. Set open-loop=1 or drop zeta-open.");
  // nstep_done counts Step() calls, which is the runcard time only at nsubsteps = 1.
  if(chi_freeze_steps && nsubsteps != 1)
    throw error_msg("confluent-wet: chi-freeze-steps counts Step() calls and is only equal "
                    "to runcard time at nsubsteps=1, got nsubsteps=", nsubsteps, ".");
  // The video stream.
  if(nvideo)
  {
    if(video_stride == 0 || LX % video_stride || LY % video_stride)
      throw error_msg("confluent-wet: video-stride=", video_stride, " must be positive and "
                      "divide both LX=", LX, " and LY=", LY, ".");
    if(video_p_scale <= 0. || video_u_scale <= 0.)
      throw error_msg("confluent-wet: nvideo>0 needs video-p-scale and video-u-scale > 0. "
                      "They are the CLIPPING limits of the stored byte (P in "
                      "[-scale,+scale], |u| in [0,scale]), not the colour-bar limits, so "
                      "set them WIDE -- the campaign convention is 6*sigma_P and 6*u_rms "
                      "stored, +/-3*sigma_P and 3*u_rms displayed. An auto-scale is "
                      "deliberately not offered: a per-run scale makes two runs "
                      "incomparable by eye, which is the whole point of the video.");
  }
  if(m0 < 0. || m0 > 1.)
    throw error_msg("confluent-wet: m0 must lie in [0,1].");
  if(tau_m <= 0.)
    throw error_msg("confluent-wet: tau-m must be positive.");
  if(chi_width < 0.)
    throw error_msg("confluent-wet: chi-width must be non-negative; 0 is the HARD STEP "
                    "chi* = Theta(switch-sign*(m - mc)), matching pmem-width = 0 for "
                    "g(P), and anything above it is the tanh.");
  if(switch_sign != 1 && switch_sign != -1)
    throw error_msg("confluent-wet: switch-sign must be +1 or -1.");
  if(pmem_width < 0.)
    throw error_msg("confluent-wet: pmem-width must be non-negative.");

  if(ntracer > 0 && tracer_count == 0)
    throw error_msg("confluent-wet: ntracer > 0 needs tracer-count > 0.");
  if(tracer_count > 0 && ntracer == 0)
    throw error_msg("confluent-wet: tracer-count > 0 needs ntracer > 0, the tracer "
                    "sampling interval in steps.");
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
  ConfigureTracers();
}

void ConfluentWet::SmoothToLength(ScalarField& fld, double length)
{
  // Diffusion is the cheapest way to give white noise a prescribed correlation length:
  // after time t the Gaussian kernel has width sqrt(2 D t), so smooth_rate*steps =
  // length^2/2 lands on `length`. chi_tmp is borrowed as scratch -- it holds no state
  // before the first update.
  //
  // ONLY VALID FOR chi. The boundary handler between sweeps is the model's own, which
  // fills chi and m and nothing else; passing any other field would diffuse it against a
  // stale halo. Both call sites pass chi.
  constexpr double smooth_rate = 0.2;
  const unsigned steps =
    length > 0 ? static_cast<unsigned>(ceil(length*length/(2*smooth_rate))) : 0u;

  for(unsigned s=0; s<steps; ++s)
  {
    for(unsigned k=0; k<DomainSize; ++k)
      chi_tmp[k] = fld[k] + smooth_rate*laplacian(fld, get_neighbours(k), sD);
    swap(fld.get_data(), chi_tmp.get_data());
    // Go through the virtual boundary handler, never Apply*() directly: bc=0 has no
    // boundary layer at all, so ApplyPBC() would index past the field.
    BoundaryConditionsFields();
  }
}

void ConfluentWet::ConfigurePhenotype()
{
  for(unsigned k=0; k<DomainSize; ++k)
    m[k] = m0;

  // A SEPARATE RNG STREAM FOR THE PHENOTYPE PATTERN, when asked for.
  //
  // Configure() lays the director down FIRST and spends one draw per node doing it, so two
  // runs differing only in `seed` differ in Q as well as in chi. A campaign that puts
  // several (chi, m) starts on ONE shared flow needs the opposite: hold `seed` fixed, vary
  // `chi-seed`, and the initial Q is bit-identical across the set while the phenotype
  // pattern changes. chi-seed = 0 keeps the historical behaviour -- the global stream,
  // shared with the director -- so every earlier runcard reproduces exactly.
  std::mt19937 chi_gen(chi_seed);
  std::normal_distribution<double> chi_gauss(0., 1.);
  const bool own_stream = chi_seed > 0;
  auto rnd_normal = [&]() { return own_stream ? chi_gauss(chi_gen) : random_normal(); };
  auto rnd_below  = [&](unsigned b) { return own_stream ? chi_gen()%b : randu()%b; };

  // ---- two-phase initial conditions ---------------------------------------
  // All three set chi AND m, because a phase is a (chi, m) pair: preparing a half-box at
  // chi = chi_hi while leaving m at the global m0 starts that half OFF its own fixed
  // point, and it would then relax on tau_m towards the other phase whether or not the
  // physics is bistable -- which is exactly the question being asked.
  if(chi_config == "stripe")
  {
    // Split along x, left half (x < LX/2) is the HIGH phase. Periodic boundaries make this
    // two flat interfaces, not one; the front speed is read off both and they are an
    // internal consistency check on each other.
    for(unsigned k=0; k<DomainSize; ++k)
    {
      const bool hi = GetXPosition(k) < LX/2;
      chi[k] = hi ? chi_hi : chi_lo;
      m[k]   = hi ? m_hi   : m_lo;
    }
    return;
  }

  if(chi_config == "blocks")
  {
    const unsigned nbx = LX/chi_block, nby = LY/chi_block, nb = nbx*nby;
    // EXACTLY half of the blocks are high: build the list, shuffle it, take the first
    // half. An independent coin per block would leave the initial area fraction scattered
    // by +/- 1/(2 sqrt(nb)) -- 1.25% at nb = 1600 -- and the mixed initial condition
    // exists precisely to start on the 50/50 line, so that scatter would be a systematic
    // seed-dependent bias in the very quantity being measured.
    std::vector<unsigned char> hi(nb, 0);
    for(unsigned b=0; b<nb/2; ++b) hi[b] = 1;
    for(unsigned b=nb; b>1; --b)         // Fisher-Yates
      swap(hi[b-1], hi[rnd_below(b)]);

    for(unsigned k=0; k<DomainSize; ++k)
    {
      const unsigned b = (GetXPosition(k)/chi_block)*nby + GetYPosition(k)/chi_block;
      chi[k] = hi[b] ? chi_hi : chi_lo;
      m[k]   = hi[b] ? m_hi   : m_lo;
    }
    return;
  }

  if(chi_config == "binary-noise")
  {
    // The mixed start with a TUNABLE DOMAIN SIZE: a correlated Gaussian field of length
    // chi_length, split at its own median into exactly half (chi_hi, m_hi) and half
    // (chi_lo, m_lo). `blocks` puts the domains on a lattice and `stripe` puts one at the
    // box scale; this is the same initial condition with neither artifact.
    //
    // THE SPLIT IS AT THE MEDIAN, NOT AT A FIXED VALUE, and that is the whole design. A
    // fixed threshold would let the sample mean of the smoothed field -- a genuinely
    // random number, since smoothing to length l leaves only ~(L/l)^2 independent patches
    // -- scatter the initial area fraction from run to run. The loop is positive
    // feedback, so that scatter would be amplified into exactly the quantity being
    // measured. Sorting and cutting in the middle makes the area fraction 0.5 by
    // construction, the same guarantee `blocks` gets from its shuffle.
    for(unsigned k=0; k<DomainSize; ++k)
      chi[k] = rnd_normal();
    BoundaryConditionsFields();
    SmoothToLength(chi, chi_length);

    std::vector<unsigned> idx(DomainSize);
    for(unsigned k=0; k<DomainSize; ++k) idx[k] = k;
    std::sort(idx.begin(), idx.end(),
              [&](unsigned p, unsigned q) { return chi[p] < chi[q]; });
    for(unsigned i=0; i<DomainSize; ++i)
    {
      const bool hi = i >= DomainSize/2;
      chi[idx[i]] = hi ? chi_hi : chi_lo;
      m[idx[i]]   = hi ? m_hi   : m_lo;
    }
    return;
  }

  if(chi_config == "uniform" || chi_noise == 0.)
  {
    for(unsigned k=0; k<DomainSize; ++k)
      chi[k] = chi0;
    return;
  }

  // White noise smoothed into a correlated field of the requested length, then rescaled to
  // the requested standard deviation.
  for(unsigned k=0; k<DomainSize; ++k)
    chi[k] = rnd_normal();
  BoundaryConditionsFields();
  SmoothToLength(chi, chi_length);

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

double ConfluentWet::Activity(unsigned k) const
{
  // OPEN LOOP: the activity is prescribed and chi is not written back into the stress.
  // chi and m still evolve -- their dynamics is what the calibration measures -- but the
  // arrow chi -> activity is cut, which is what makes f(zeta_eff) a function rather than a
  // fixed point condition.
  if(open_loop) return zeta_open;
  // CLOSED LOOP with an activity FLOOR:
  //     zeta_eff = zeta [ z0 + (1 - z0)(1 - chi) ]
  // z0 = 0 recovers the pre-2026-09 law zeta*(1-chi). z0 > 0 keeps the passive phase
  // weakly active, so it still stirs, still generates pressure fluctuations, and can still
  // be driven back -- without it chi = 1 is an absorbing state and coexistence is
  // impossible by construction rather than by physics.
  return zeta*(zeta0_frac + (1. - zeta0_frac)*(1. - chi[k]));
}

double ConfluentWet::ChiStar(double mm) const
{
  // chi_width = 0 is the HARD STEP, chi* = Theta(s*(m - mc)); on the s = -1 branch the
  // campaign runs, that reads chi* = Theta(mc - m). Written as a separate branch for the
  // same reason MemoryTarget has one: tanh((m-mc)/w) has no w -> 0 limit a floating-point
  // division can represent. The tie m == mc resolves to 0, matching g's P == pmem.
  return chi_width > 0
    ? .5*(1. + switch_sign*tanh((mm - mc)/chi_width))
    : (switch_sign*(mm - mc) > 0. ? 1. : 0.);
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
  const double active = Activity(k);
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
    + (use_switching && nstep_done >= chi_freeze_steps
       ? (ChiStar(m_k) - chi_k)/tau_chi : 0.);

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

  // Same stage, same velocity field, same integrator as the two lines above.
  AdvanceTracers(first);
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

  ++nstep_done;
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

// =============================================================================
// The video stream

void ConfluentWet::VideoOpen(const string& dir)
{
  // Truncating (not appending): a rerun into the same output directory must not glue its
  // frames onto the previous run's, which would be invisible in the byte stream and would
  // silently corrupt every time series computed from it.
  const auto mode = ios::out | ios::binary | ios::trunc;
  vid_u  .open((dir + "video_u.u8"  ).c_str(), mode);
  vid_p  .open((dir + "video_P.u8"  ).c_str(), mode);
  vid_m  .open((dir + "video_m.u8"  ).c_str(), mode);
  vid_chi.open((dir + "video_chi.u8").c_str(), mode);
  vid_meta.open((dir + "video_meta.csv").c_str(), ios::out | ios::trunc);
  if(!vid_u || !vid_p || !vid_m || !vid_chi || !vid_meta)
    throw error_msg("confluent-wet: cannot open the video streams in '", dir, "'.");

  // Header: everything the renderer needs that is not in parameters.json, plus the exact
  // per-frame scalars, which are computed in double precision and are therefore NOT
  // affected by the uint8 quantisation of the fields.
  vid_meta << "# confluent-wet video stream\n"
           << "# fields video_{u,P,m,chi}.u8: raw uint8, "
           << LX/video_stride << "x" << LY/video_stride
           << " per frame, x-major (reshape (LX/stride, LY/stride))\n"
           << "# stored range: chi,m in [0,1]; P in [" << -video_p_scale << ","
           << video_p_scale << "]; |u| in [0," << video_u_scale << "]\n"
           << "t,chi_mean,chi_std,m_mean,m_std,u_rms,P_mean,P_std,P_clip,u_clip\n";
  vid_meta.precision(10);
  video_open = true;
}

void ConfluentWet::VideoPack(vector<unsigned char>& out, const ScalarField& fld,
                             double lo, double hi) const
{
  // BLOCK MEAN, not subsampling. Subsampling a turbulent field aliases: the defect cores
  // that dominate |u| and P are a few lattice units across, so every other one would be
  // dropped or doubled depending on where it sits, and the video would flicker with the
  // sampling phase rather than with the flow.
  const unsigned sx = video_stride, nx = LX/sx, ny = LY/sx;
  const double norm = 1./(sx*sx), span = hi - lo;
  for(unsigned bx=0; bx<nx; ++bx)
    for(unsigned by=0; by<ny; ++by)
    {
      double acc = 0.;
      for(unsigned dx=0; dx<sx; ++dx)
        for(unsigned dy=0; dy<sx; ++dy)
          acc += fld[GetDomainIndex(bx*sx+dx, by*sx+dy)];
      const double x = (acc*norm - lo)/span;
      out[by + ny*bx] = x <= 0. ? 0 : (x >= 1. ? 255
                                     : static_cast<unsigned char>(lround(255.*x)));
    }
}

void ConfluentWet::ConfigureTracers()
{
  tr_x.clear(); tr_y.clear(); tr_xN.clear(); tr_yN.clear();
  tracer_nx = tracer_ny = 0;
  if(ntracer == 0 || tracer_count == 0) return;

  // STRATIFIED, not random. The box holds ~(L/L_P)^2 independent pressure patches -- of
  // order 4e4 at L = 1000 -- so a uniform grid samples them all exactly once while a
  // random draw adds placement variance for nothing. Deterministic, hence no seed: two
  // runs differing only in `seed` still carry identical tracer grids, which is what makes
  // the seed-to-seed comparison of a Lagrangian statistic meaningful.
  const double aspect = double(LX)/double(LY);
  tracer_nx = max(1u, unsigned(lround(sqrt(double(tracer_count)*aspect))));
  tracer_ny = max(1u, tracer_count/tracer_nx);
  const double dx = double(LX)/tracer_nx, dy = double(LY)/tracer_ny;
  for(unsigned i=0; i<tracer_nx; ++i)
    for(unsigned j=0; j<tracer_ny; ++j)
    {
      tr_x.push_back((i + .5)*dx);
      tr_y.push_back((j + .5)*dy);
    }
  tr_xN.assign(tr_x.size(), 0.);
  tr_yN.assign(tr_y.size(), 0.);
}

double ConfluentWet::InterpolateField(const ScalarField& f, double x, double y) const
{
  // BILINEAR, and deliberately not better. It is the second-order interpolant whose
  // gradient reduces exactly to the centred difference the field solver uses, so the
  // tracer and the fields carry the same truncation order. A cubic or spline interpolant
  // would make the tracer MORE accurate than the field it is meant to represent, which is
  // the wrong kind of mismatch: the trajectory would then be a path of a flow the fields
  // are not being advected by.
  const double xf = floor(x), yf = floor(y);
  const double fx = x - xf,  fy = y - yf;
  const int iL = int(LX), iM = int(LY);
  const unsigned i0 = unsigned((int(xf) % iL + iL) % iL);
  const unsigned j0 = unsigned((int(yf) % iM + iM) % iM);
  const unsigned i1 = (i0 + 1) % LX;
  const unsigned j1 = (j0 + 1) % LY;
  return (1.-fx)*(1.-fy)*f[GetDomainIndex(i0, j0)]
       +      fx *(1.-fy)*f[GetDomainIndex(i1, j0)]
       + (1.-fx)*     fy *f[GetDomainIndex(i0, j1)]
       +      fx *     fy *f[GetDomainIndex(i1, j1)];
}

void ConfluentWet::AdvanceTracers(bool first)
{
  if(tr_x.empty()) return;

  // Mirrors UpdateNodeFields term for term. There, chi is carried as
  //     predictor : chiN = chi + .5 rhs(state_0) ,  chi = chi + rhs(state_0)
  //     corrector : chi  = chiN + .5 rhs(state_1)
  // which is Heun at npc = 1 and its fixed-point iteration beyond. The tracer runs the
  // same two lines with rhs = u_mat interpolated at the tracer, and it is called from
  // inside UpdateFields, so the velocity it reads is the SAME ux_mat, at the same stage,
  // that the fields are being advected by -- not a re-derived or time-averaged copy.
  const double Lx = double(LX), Ly = double(LY);
  for(unsigned p=0; p<tr_x.size(); ++p)
  {
    const double vx = InterpolateField(ux_mat, tr_x[p], tr_y[p]);
    const double vy = InterpolateField(uy_mat, tr_x[p], tr_y[p]);
    double nx, ny;
    if(first)
    {
      tr_xN[p] = tr_x[p] + .5*vx;
      tr_yN[p] = tr_y[p] + .5*vy;
      nx = tr_x[p] + vx;
      ny = tr_y[p] + vy;
    }
    else
    {
      nx = tr_xN[p] + .5*vx;
      ny = tr_yN[p] + .5*vy;
    }
    // The anchors tr_xN/tr_yN are left UNWRAPPED on purpose: wrapping only the final
    // position keeps predictor and corrector on the same branch across a periodic edge.
    // |u| << 1 (CFL ~ 4e-3), so they never stray more than a lattice unit outside.
    tr_x[p] = nx - Lx*floor(nx/Lx);
    tr_y[p] = ny - Ly*floor(ny/Ly);
  }
}

void ConfluentWet::TracerOpen(const string& dir)
{
  const auto mode = ios::out | ios::binary | ios::trunc;
  trc_dat.open((dir + "tracer.f32").c_str(), mode);
  trc_meta.open((dir + "tracer_meta.csv").c_str(), ios::out | ios::trunc);
  if(!trc_dat || !trc_meta)
    throw error_msg("confluent-wet: cannot open the tracer streams in '", dir, "'.");

  trc_meta << "# confluent-wet Lagrangian tracer stream\n"
           << "# tracer.f32: float32, " << tr_x.size() << " tracers x 4 values per sample,\n"
           << "#   laid out [P, m, x, y] per tracer, tracer-major; reshape (nsample, "
           << tr_x.size() << ", 4)\n"
           << "# grid " << tracer_nx << "x" << tracer_ny << " (requested " << tracer_count
           << "), sampled every " << ntracer << " steps\n"
           << "# P and m are bilinear interpolations at the tracer, the same interpolant\n"
           << "#   used to advect it; x,y are wrapped into [0,LX) x [0,LY)\n"
           << "t,n_tracer,P_mean,P_std\n";
  tracer_open = true;
}

void ConfluentWet::WriteTracers(const string& dir, unsigned t)
{
  if(ntracer == 0 || tr_x.empty() || t % ntracer) return;
  if(!tracer_open) TracerOpen(dir);

  const unsigned np = tr_x.size();
  vector<float> buf(4*np);
  double pm = 0., p2 = 0.;
  for(unsigned p=0; p<np; ++p)
  {
    const double P = InterpolateField(pressure, tr_x[p], tr_y[p]);
    const double M = InterpolateField(m,        tr_x[p], tr_y[p]);
    buf[4*p+0] = float(P);
    buf[4*p+1] = float(M);
    buf[4*p+2] = float(tr_x[p]);
    buf[4*p+3] = float(tr_y[p]);
    pm += P; p2 += P*P;
  }
  pm /= np;
  const double psd = sqrt(max(0., p2/np - pm*pm));

  trc_dat.write((char*)buf.data(), buf.size()*sizeof(float));
  trc_meta << t << ',' << np << ',' << pm << ',' << psd << '\n';

  // Flushed every sample, for the reason spelled out in WriteAuxiliary: nothing ever
  // destructs this model, so an unflushed tail is simply lost -- and for a correlation
  // function the lost tail is the longest lag, i.e. the part being measured.
  trc_dat.flush();
  trc_meta.flush();
}

void ConfluentWet::WriteAuxiliary(const string& dir, unsigned t)
{
  // The tracers run on their own, finer clock, so this is before the video gate.
  WriteTracers(dir, t);

  if(nvideo == 0 || t % nvideo) return;
  if(!video_open) VideoOpen(dir);

  const unsigned nx = LX/video_stride, ny = LY/video_stride;
  vector<unsigned char> buf(nx*ny);

  // |u| has to be built; it is not a stored field. speed_tmp is local rather than a member
  // because this runs once every ~10^2 steps and the allocation is 1% of the write.
  ScalarField speed;
  speed.SetSize(LX, LY, Type);
  double u2 = 0., cm = 0., c2 = 0., mm = 0., m2 = 0., pm = 0., p2 = 0.;
  unsigned pclip = 0, uclip = 0;
  for(unsigned k=0; k<DomainSize; ++k)
  {
    const double sp = sqrt(ux_mat[k]*ux_mat[k] + uy_mat[k]*uy_mat[k]);
    speed[k] = sp;
    u2 += sp*sp;
    cm += chi[k];      c2 += chi[k]*chi[k];
    mm += m[k];        m2 += m[k]*m[k];
    pm += pressure[k]; p2 += pressure[k]*pressure[k];
    if(fabs(pressure[k]) > video_p_scale) ++pclip;
    if(sp > video_u_scale)                ++uclip;
  }
  const double N = DomainSize;
  cm /= N; mm /= N; pm /= N;
  const double csd = sqrt(max(0., c2/N - cm*cm));
  const double msd = sqrt(max(0., m2/N - mm*mm));
  const double psd = sqrt(max(0., p2/N - pm*pm));

  VideoPack(buf, speed,    0.,             video_u_scale); vid_u  .write((char*)buf.data(), buf.size());
  VideoPack(buf, pressure, -video_p_scale, video_p_scale); vid_p  .write((char*)buf.data(), buf.size());
  VideoPack(buf, m,        0.,             1.);            vid_m  .write((char*)buf.data(), buf.size());
  VideoPack(buf, chi,      0.,             1.);            vid_chi.write((char*)buf.data(), buf.size());

  vid_meta << t << ',' << cm << ',' << csd << ',' << mm << ',' << msd << ','
           << sqrt(u2/N) << ',' << pm << ',' << psd << ','
           << pclip/N << ',' << uclip/N << '\n';

  // EVERY stream is flushed every frame, and the binary ones are not an afterthought: the
  // model object is deliberately never deleted (ModelPtr's destructors are commented out),
  // so no destructor ever runs and nothing closes these files at exit. Without an explicit
  // flush the last partially-filled buffer is simply dropped -- measured: a 101-frame run
  // left exactly 100 frames on disk while the CSV, which was already being flushed, listed
  // 101. That is the worst possible failure mode for this stream, because the missing frame
  // is the LAST one and every 'what did it settle to' question reads the end of the run.
  // At ~10^2 steps between frames the cost is not measurable.
  vid_u.flush(); vid_p.flush(); vid_m.flush(); vid_chi.flush(); vid_meta.flush();
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
     "width of the tanh in chi*(m); 0 is the hard step chi* = Theta(switch-sign*(m-mc))")
    ("mc", opt::value<double>(&mc),
     "memory threshold in chi*(m)")
    ("switch-sign", opt::value<int>(&switch_sign),
     "+1: compression drives chi up (contact inhibition); -1: the opposite branch")
    ("tau-m", opt::value<double>(&tau_m),
     "memory relaxation time")
    ("pmem", opt::value<double>(&pmem),
     "pressure threshold of the memory source g(P)")
    ("pmem-width", opt::value<double>(&pmem_width),
     "smoothing width of g(P); 0 is a sharp step")
    ("zeta0-frac", opt::value<double>(&zeta0_frac),
     "activity floor as a fraction of zeta: zeta_eff = zeta*(z0 + (1-z0)*(1-chi)). "
     "0 (default) is the pre-2026-09 law zeta*(1-chi)")
    ("open-loop", opt::value<int>(&open_loop),
     "1: m and chi evolve but the activity is held at zeta-open and chi is not fed back")
    ("zeta-open", opt::value<double>(&zeta_open),
     "the prescribed activity in open-loop mode")
    ("chi-freeze-steps", opt::value<unsigned>(&chi_freeze_steps),
     "hold the phenotype switching source off for this many steps (transport stays on)")
    ("nvideo", opt::value<unsigned>(&nvideo),
     "steps between video-stream frames; 0 disables the stream")
    ("video-stride", opt::value<unsigned>(&video_stride),
     "block-averaging factor of the video lattice (must divide LX and LY)")
    ("video-p-scale", opt::value<double>(&video_p_scale),
     "stored clipping limit of P in the video stream (P in [-scale,+scale]); set it WIDE, "
     "the colour-bar limit is chosen by the renderer")
    ("video-u-scale", opt::value<double>(&video_u_scale),
     "stored clipping limit of |u| in the video stream (|u| in [0,scale])")
    ("frame-light", opt::value<int>(&frame_light),
     "1: write only the seven fields the analysis reads, dropping the LB populations and "
     "the stress decomposition (~4x smaller frames, not restartable)");

  // options[1] is Nematic's "Initial configuration options"
  options[1].add_options()
    ("director-config", opt::value<string>(&director_config),
     "director init: uniform (angle+noise) or defect-pair")
    ("defect-sep", opt::value<double>(&defect_sep),
     "separation of the +/- 1/2 defect pair (lattice units, defect-pair mode)")
    ("initial-order", opt::value<double>(&init_order),
     "initial nematic order amplitude")
    ("chi-config", opt::value<string>(&chi_config),
     "phenotype init: uniform, noise, binary-noise, stripe or blocks")
    ("chi0", opt::value<double>(&chi0),
     "initial mean phenotype")
    ("chi-noise", opt::value<double>(&chi_noise),
     "standard deviation of the initial phenotype noise")
    ("chi-length", opt::value<double>(&chi_length),
     "correlation length of the initial phenotype noise (noise and binary-noise)")
    ("chi-seed", opt::value<unsigned>(&chi_seed),
     "seed of the phenotype pattern alone; 0 (default) shares the run's global generator. "
     "Non-zero decouples the pattern from the director, so a set of runs at one `seed` and "
     "several `chi-seed` share their initial Q exactly")
    ("m0", opt::value<double>(&m0),
     "initial memory")
    ("chi-lo", opt::value<double>(&chi_lo),
     "phenotype of the LOW phase (chi-config = stripe, blocks or binary-noise)")
    ("chi-hi", opt::value<double>(&chi_hi),
     "phenotype of the HIGH phase (chi-config = stripe, blocks or binary-noise)")
    ("m-lo", opt::value<double>(&m_lo),
     "memory of the LOW phase (chi-config = stripe, blocks or binary-noise)")
    ("m-hi", opt::value<double>(&m_hi),
     "memory of the HIGH phase (chi-config = stripe, blocks or binary-noise)")
    ("chi-block", opt::value<unsigned>(&chi_block),
     "block edge in lattice units for chi-config = blocks; must divide LX and LY")
    ("ntracer", opt::value<unsigned>(&ntracer),
     "Lagrangian tracer sampling interval in steps (0 = tracers off)")
    ("tracer-count", opt::value<unsigned>(&tracer_count),
     "requested number of tracers; rounded to a near-square stratified grid");

  return options;
}
