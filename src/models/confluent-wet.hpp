#ifndef MODELS_CONFLUENT_WET_HPP_
#define MODELS_CONFLUENT_WET_HPP_

#include <fstream>

#include "models/nematic.hpp"

/** Confluent, incompressible active nematic with phenotype and mechanical memory,
 *  on a WET (lattice Boltzmann) foundation.
 *
 * Rebuilt without a phase field. Derives from Nematic -- the traditional wet active
 * nematic -- and NOT from the Lyotropic line, so there is no phi at all: no binary
 * order parameter, no double well, no interface term, no Cahn-Hilliard mobility, no
 * crowding modulus B, no proliferation, no division/death masks, no vacuum. The layer
 * is the whole domain and the only density in the problem is the LB density n.
 *
 *   f      = 1/2 CC (1 - q^2)^2  +  Frank(LL)          , q^2 = Qxx^2 + Qyx^2
 *   H      = 2 CC term Q + LL grad^2 Q                 , term = 1 - q^2
 *   D_t Q  = Gamma H + corotation + xi-alignment         (Beris-Edwards)
 *   sigma  = sigma_bulk delta + traceless(elastic + Frank + active)
 *   active stress = -zeta_eff(chi) Q ,  zeta_eff = zeta [z0 + (1 - z0)(1 - chi)]
 *   LB     : f_v <- f_v + (f^eq_v - f_v)/tau + w_v (F.c_v) , F = div(sigma) - friction u
 *
 *   D_t chi = (chi*(m) - chi)/tau_chi + Dbio lap chi , chi*(m)= .5(1+s tanh((m-mc)/chi_width))
 *   D_t m   = (g(P)    - m)/tau_m     + Dbio lap m   , g(P)   = .5(1+  tanh((P-pmem)/pmem_width))
 *
 * Either width may be set to ZERO, which replaces its tanh by the hard step it is the
 * smoothing of -- chi* = Theta(s(m - mc)) and g = Theta(P - pmem). The two are independent
 * switches; the 2026-09 step campaign sets both.
 *
 * LANDAU-DE GENNES CONVENTION. H = 2 CC term Q, matching go-or-grow / confluent-memory,
 * NOT the H = CC term Q of nematic.cpp. That file pairs H = CC term Q with a bulk stress
 * 1/2 CC term^2, but H = -dF/dQ forces f = 1/4 CC term^2 for that H: its isotropic stress
 * is twice its own free energy. Doubling H instead of halving the stress fixes the
 * inconsistency AND makes CC mean the same thing here as in the dry confluent-memory runs,
 * so the two model lines can be compared at equal parameters. The nematic correlation
 * length (defect core radius) is then
 *
 *   xi_N = sqrt(LL/(2 CC))
 *
 * and it must be kept above ~1.5 lattice units or the cores -- where sigma_bulk and the
 * pressure signal live -- are simply not resolved.
 *
 * PRESSURE. The flow is incompressible (weakly compressible in the usual LB sense,
 * div u = O(Ma^2)), so the pressure is the Lagrange multiplier of the constraint, not an
 * equation of state; no crowding modulus is needed to generate it. The isotropic stress
 * has two pieces,
 *
 *   p_LB       = n/3                              (LB equation of state, = the multiplier)
 *   sigma_bulk = 1/2 CC (1 - q^2)^2               (Landau-de Gennes bulk)
 *
 * sigma_bulk is a GENUINE constitutive stress, not bookkeeping: for a free energy that does
 * not depend on the density, sigma_iso = f - n df/dn = f. (With phi present this is the
 * familiar Gibbs-Duhem f - mu phi; we are simply in the df/dn = 0 case.) It is therefore
 * kept inside sigmaXX/sigmaYY and fed to the LB as a force, exactly as the constitutive
 * theory prescribes. It vanishes in the ordered state (term -> 0) and reaches 1/2 CC only
 * where the order melts, i.e. in defect cores.
 *
 * The mechanical pressure -- the isotropic part of the traction, and what the memory
 * responds to -- is the half-trace of the TOTAL momentum flux:
 *
 *   P = -1/2 Tr(Pi) = p_LB - sigma_bulk
 *
 * The subtraction is not the removal of a spurious term: the constraint reaction inside
 * p_LB rises to meet sigma_bulk (an isotropic stress does no work in incompressible flow),
 * and P is what is left over -- the net load. Both fields are written to every frame,
 * referenced to the resting state (n = rho, sigma_bulk = 0):
 *
 *   pressure    = (n - rho)/3 - sigma_bulk    physical; drives m
 *   pressure_lb = (n - rho)/3                 diagnostic; its difference from `pressure`
 *                                             must reproduce 1/2 CC (1-q^2)^2 exactly,
 *                                             which checks the trace/traceless split
 *
 * NOTE that in a periodic box mass conservation pins <n> = rho exactly, so <P> is set by
 * the initial density and NOT by the dynamics: only P - <P> is meaningful, and pmem has to
 * be calibrated per regime rather than carried over from the dry model, where
 * P = 1/2 B (phi^2-1) had an absolute zero at phi = phi_critical.
 *
 * MATERIAL VELOCITY. In a forced lattice Boltzmann the fluid velocity is NOT the bare
 * moment of the distribution: the force is applied trapezoidally in time, so
 *
 *   u_phys = (Sum_v f_v c_v + F/2)/n = u_code + F/(2n)
 *
 * With friction folded in this closes exactly,
 *
 *   u_mat = [u_code + div(sigma)/(2n)] / [1 + friction/(2n)] ,
 *
 * which is what ComputeMaterialVelocity() stores in ux_mat/uy_mat. This is not a detail:
 * near a defect core div(sigma) is dominated by the Frank stress and reaches ~0.14 at
 * LL = 0.4, xi_N = 1.4, so F/(2n) ~ 1.8e-3 against a typical flow speed of ~4e-3 -- a 40%
 * offset, precisely where the pressure signal lives and where chi and m most need to be
 * transported correctly. The LB collision keeps using the bare u_code (the simple forcing
 * scheme is constructed around f^eq(u_code) and must not be half-modified); every PHYSICAL
 * use -- advection of Q, chi and m, the velocity gradients in Beris-Edwards, and the output
 * -- uses u_mat.
 *
 * TRANSPORT OF chi AND m. Both are INTENSIVE (a phenotype fraction and a memory in [0,1]),
 * so the correct equation is the material derivative D_t = d_t + u.grad, not the
 * conservative d_t + div(. u): the two differ by chi div u, which is small but wrong.
 * Because nothing has to be conserved, none of the dry model's phi-weighted face-flux
 * machinery is needed -- plain centred differences suffice. The CFL number is ~4e-3, so
 * the weak instability of centred advection under Heun stepping grows like CFL^4/8 ~ 3e-11
 * per step and is irrelevant, while first-order upwinding would inject a numerical
 * diffusivity |u|dx/2 ~ 2e-3 that would swamp any physical Dbio. A clamp back into [0,1]
 * after each update is the safety net against dispersive overshoot.
 *
 * ONE DIFFUSIVITY FOR BOTH FIELDS. Diffusion here is not a property of the phenotype or of
 * the memory separately: it is the motility of the cells themselves, and a cell carries its
 * phenotype and its memory together. There is therefore a SINGLE coefficient Dbio, applied
 * identically to chi and to m. Formally, if the biomass density c diffuses with Dbio then
 * the extensive densities obey d_t(c X) + div(c X u) = div(Dbio grad(c X)); the layer is
 * confluent and incompressible, so c is uniform and this reduces to D_t X = ... + Dbio lap X
 * for X = chi and for X = m alike.
 *
 * The transport operators of chi and m are thus IDENTICAL BY CONSTRUCTION, at any Dbio. The
 * two fields differ only through their source terms, which makes any discrepancy beyond
 * that a bug and makes the pair a clean two-field diagnostic.
 *
 * Dbio DEFAULTS TO ZERO. At Dbio = 0 the model is NOT a reaction-diffusion system. It is
 *
 *   hyperbolic transport (pure advection, no intrinsic length)
 *     + a local two-stage relaxation CASCADE m -> chi along each material trajectory
 *     + one NONLOCAL elliptic coupling, chi -> P, through lap P = d_i d_j sigma_dev,ij.
 *
 * Two consequences are worth stating because they are easy to get wrong. Turing patterns
 * are impossible (no diffusion, let alone two different diffusivities), so any spatial
 * structure comes from stirring plus the pressure coupling. And the LOCAL system cannot
 * oscillate: m drives chi but chi does not drive m locally, and a cascade of two
 * first-order lags is unconditionally stable -- so any oscillation or travelling wave that
 * appears is necessarily hydrodynamic, closed through chi -> activity -> flow -> P.
 *
 * WHAT REPLACES THE DIFFUSIVE CUTOFF. With Dbio = 0 nothing stops chaotic advection from
 * thinning filaments indefinitely; what saves it is the relaxation itself, which resets a
 * material element towards chi*(m) before the strain has had long enough to compress it:
 *
 *   thinnest structure ~ l_a exp(-lambda tau) ,   lambda ~ |grad u| ~ u/l_a
 *
 * so the field stays resolved only while lambda*tau <~ ln(l_a/3). At the step-1 test
 * parameters (zeta = 1.8e-3, l_a = 15, eta = 6.67, hence u ~ 4e-3 and lambda ~ 2.7e-4 per
 * step) that bound is tau <~ 6000 steps: tau_m = 2000 is comfortable, tau_chi = 6000 sits
 * exactly on the edge. Pushing tau_chi far beyond the eddy turnover time 1/lambda ~ 3700
 * steps therefore needs either more lattice units per l_a or a non-zero Dbio. Switching Dbio
 * on restores a genuine cutoff at the Batchelor scale sqrt(Dbio/lambda), so keeping that at
 * ~3 lattice units would take Dbio ~ 9*lambda ~ 2.4e-3 here -- note in passing that the dry
 * model's Dchi = 2e-4 was NOT a cutoff (its Batchelor scale is 0.86 lattice units); that
 * model also relied on the relaxation. Watch the clamp-hit count and the high-k end of the
 * chi spectrum.
 */
class ConfluentWet : public Nematic
{
protected:
  /** Phenotype in [0,1]: chi = 1 is the passive/grow branch, so the active stress is
   *  -zeta*(1-chi)*Q. Predictor and buffer for the predictor-corrector step. */
  ScalarField chi, chiN, chi_tmp;
  /** Mechanical memory in [0,1]: the exponentially weighted fraction of recent time the
   *  material element spent above the pressure threshold. */
  ScalarField m, mN, m_tmp;

  /** Material (physical) velocity, u_code + F/(2n) with friction resolved. Everything
   *  physical uses this; the LB collision keeps the bare moment in ux/uy. */
  ScalarField ux_mat, uy_mat;

  /** Isotropic part of the thermodynamic stress, 1/2 CC term^2. The code builds
   *  sigmaXX = sigmaF + sigma_bulk, sigmaYY = -sigmaF + sigma_bulk, so sigma_bulk IS the
   *  half-trace of sigma and everything else is exactly traceless. */
  ScalarField sigma_bulk;
  /** Active stress -zeta*(1-chi)*Q, stored for the stress decomposition (traceless) */
  ScalarField sigma_active_xx, sigma_active_yx;
  /** Physical mechanical pressure -1/2 Tr(Pi) = (n - rho)/3 - sigma_bulk */
  ScalarField pressure;
  /** Gauge-dependent hydrodynamic pressure (n - rho)/3; diagnostic only */
  ScalarField pressure_lb;

  /** ACTIVITY FLOOR, as a fraction of zeta:
   *
   *      zeta_eff(chi) = zeta [ zeta0_frac + (1 - zeta0_frac)(1 - chi) ]
   *
   *  so chi = 0 gives the full zeta and chi = 1 gives zeta0_frac*zeta instead of ZERO.
   *  A layer whose passive phenotype is exactly non-active is a singular limit: at chi = 1
   *  the flow dies completely, P collapses to the equilibrium value everywhere, and the
   *  memory can never be re-excited -- chi = 1 is then an absorbing state and no two-phase
   *  coexistence is possible whatever the parameters. The floor is what makes the passive
   *  phase a PHASE rather than a trap.
   *
   *  DEFAULT 0, which is exactly the pre-2026-09 model (active = zeta*(1-chi)), so every
   *  earlier runcard reproduces bit for bit. The 2026-09 memory campaign sets 0.3
   *  EXPLICITLY in every run.dat. */
  double zeta0_frac = 0.;

  /** OPEN LOOP: m and chi evolve exactly as usual, but the activity is held at the
   *  constant zeta_open and chi is NOT written back into the stress.
   *
   *  This is the calibration mode, and it is not the same as freezing chi (tau_chi <= 0).
   *  Freezing chi would also freeze the activity at zeta_eff(chi0), but it would kill the
   *  phenotype dynamics we want to observe; what the activity ladder needs is the response
   *  of the flow AND of m to a PRESCRIBED activity, so that f -- the fraction of time a
   *  material point spends above the pressure threshold -- can be tabulated against
   *  zeta_eff with the feedback cut. Closing the loop then only requires solving
   *  chi_bar = chi*(f(zeta_eff(chi_bar))) on that table. */
  int open_loop = 0;
  double zeta_open = 0.;

  /** Biomass (cell-motility) diffusivity, applied IDENTICALLY to chi and m: a cell
   *  carries its phenotype and its memory together. Defaults to 0. */
  double Dbio = 0.;
  /** Phenotype relaxation time. tau_chi <= 0 disables switching (frozen-chi control). */
  double tau_chi = 0.;
  /** Width and centre of the switch in chi*(m).
   *
   *  chi-width = 0 is the HARD STEP chi* = Theta(switch_sign*(m - mc)), which on the
   *  s = -1 branch reads chi* = Theta(mc - m). It is a separate branch rather than a very
   *  small width, exactly as pmem-width = 0 already is for g(P), because tanh((m-mc)/w)
   *  has no w -> 0 limit a floating-point division can represent.
   *
   *  THE STEP IS NOT AN INFINITELY SHARP SWITCH. What smooths chi* once w_chi is gone is
   *  the spatial spread of m itself: the space-averaged chi* is the CDF of m, whose
   *  central slope corresponds to an effective width sqrt(2 pi)/2 * sigma_m. The model
   *  therefore loses a parameter AND its sharpness follows tau_m automatically, since
   *  sigma_m falls with tau_m -- which is what makes tau_m, and not a width, the single
   *  knob that decides whether the bistability is realised. */
  double chi_width = 0.15, mc = 0.5;
  /** +1: sustained compression drives chi UP, towards passive/grow (contact inhibition,
   *  a self-limiting negative feedback since the activity is zeta*(1-chi)). -1 is the
   *  opposite branch. */
  int switch_sign = 1;
  /** Phenotype initialisation: mode, mean, std and correlation length of the noise.
   *
   *  Modes:
   *    uniform   chi = chi0, m = m0 everywhere
   *    noise     correlated Gaussian noise of width chi_noise and length chi_length
   *    binary-noise
   *              the same correlated Gaussian field, but split at its OWN MEDIAN into
   *              exactly half (chi_hi, m_hi) and half (chi_lo, m_lo). The mixed start
   *              with a tunable domain size: `blocks` on a random lattice, `stripe` at
   *              the box scale, this one at chi_length. Splitting at the median rather
   *              than at a fixed value is what keeps the initial area fraction on the
   *              50/50 line run after run -- a fixed threshold would let the sample mean
   *              of the smoothed field scatter it, and that is the very quantity the
   *              positive feedback amplifies.
   *    stripe    two half-boxes split along x: x < LX/2 carries (chi_hi, m_hi), the rest
   *              (chi_lo, m_lo).  The front-propagation geometry: one flat interface, two
   *              of them under periodic boundaries, and both phases prepared at their own
   *              self-consistent (chi, m) so the interface is the only thing out of
   *              equilibrium.
   *    blocks    a chi_block x chi_block checkerboard of randomly assigned (chi_hi, m_hi)
   *              and (chi_lo, m_lo) cells, EXACTLY half of each (a shuffled assignment,
   *              not an independent coin per block -- with 1600 blocks an independent coin
   *              leaves a +/- 1.25% drift in the initial area fraction, and the whole point
   *              of the mixed initial condition is that it starts on the 50/50 line). */
  std::string chi_config = "uniform";
  double chi0 = 0.5, chi_noise = 0., chi_length = 0.;
  /** SEED OF THE PHENOTYPE PATTERN ALONE, 0 = share the run's global generator.
   *
   *  Configure() lays down the director FIRST and consumes one random draw per node doing
   *  it, so two runs differing only in `seed` differ in Q as well as in chi. A campaign
   *  that compares several (chi, m) starts against ONE shared flow needs the opposite:
   *  fix `seed`, vary `chi-seed`, and the initial Q -- and hence the whole hydrodynamic
   *  initial condition -- is bit-identical across the set while the phenotype pattern
   *  changes. Non-zero draws the pattern from a local generator instead, which also
   *  leaves the global stream untouched for everything downstream. */
  unsigned chi_seed = 0;
  /** The two phases used by chi-config = stripe and blocks */
  double chi_lo = 0., chi_hi = 1., m_lo = 0., m_hi = 1.;
  /** Block edge in lattice units (blocks mode); must divide both LX and LY */
  unsigned chi_block = 0;
  /** Steps during which the phenotype SWITCHING SOURCE is held off.
   *
   *  Transport and diffusion of chi keep running -- freezing those too would hold a
   *  perfectly sharp initial interface against a flow that is meanwhile developing, and
   *  the front would then start from a state the dynamics can never produce. What this
   *  freezes is only the reaction, so the flow and the pressure field can equilibrate to
   *  the prescribed chi pattern before the pattern is allowed to answer back. */
  unsigned chi_freeze_steps = 0;

  /** Memory relaxation time. Should sit well above the LOCAL acoustic time sqrt(3)*L_P set
   *  by the pressure correlation length -- NOT by the box size, which overestimates it by
   *  roughly an order of magnitude. L_P is measurable only from a run, so this is checked by
   *  hand against a frozen-chi calibration, not at startup. */
  double tau_m = 1.;
  /** Threshold and smoothing width of the memory response g(P); width 0 is a sharp step */
  double pmem = 0., pmem_width = 0.;
  /** Initial memory */
  double m0 = 0.;

  /** Director initialisation: "uniform" (angle + noise) or "defect-pair" */
  std::string director_config = "uniform";
  /** Separation of the +/- 1/2 defect pair, lattice units (defect-pair mode) */
  double defect_sep = 0.;
  /** Initial amplitude of the nematic order */
  double init_order = 1.;

  /** THE VIDEO STREAM -- a second output channel, on its own clock.
   *
   * Every nvideo steps, four fields (|u|, P, m, chi) are block-averaged onto a
   * video_stride-times coarser lattice, quantised to uint8 against FIXED limits, and
   * APPENDED to four raw byte streams; a companion CSV carries the exact (double)
   * per-frame scalars. At L = 800, stride 2 this is 160 kB per field per frame against
   * ~170 MB for a full frame, i.e. 1000x cheaper, which is what makes a 0.2*tau_c video
   * affordable in the same run as 5*tau_c analysis frames.
   *
   * WHY THE LIMITS ARE FIXED AND EXPLICIT rather than auto-scaled per run. A video whose
   * colour scale follows its own field animates a quantity that grows by decades as if it
   * were constant, and two runs at different tau_m then cannot be compared by eye at all.
   * The campaign convention is one scale for every run, taken from the calibrated
   * sigma_P(zeta) and u_rms(zeta).
   *
   * WHY THE STORED RANGE IS NOT THE DISPLAYED RANGE. video_p_scale is the CLIPPING limit
   * of the stored byte, not the colour-bar limit; the renderer maps a narrower window
   * (the campaign convention is +/- 3 sigma_P) inside it. Storing wider than displaying
   * costs only quantisation resolution -- 128 levels instead of 255 across the displayed
   * window, still far finer than the eye or than any statistic computed from it -- and it
   * buys the one thing that cannot be recovered afterwards: the tails are not destroyed.
   * That matters because these streams are also the TIME SERIES the analysis uses (a
   * 5*tau_c frame cadence cannot resolve a lag of 3(tau_m + tau_chi) at small tau_m), and
   * a threshold statistic like f = fraction of time above pmem is meaningless if the
   * field was clipped near the threshold.
   *
   * Streams are opened on the first write and closed by the destructor; the per-frame
   * clipped fractions are recorded in the CSV so a badly chosen scale is visible in the
   * data rather than only in a washed-out picture. */
  unsigned nvideo = 0, video_stride = 2;
  double video_p_scale = 0., video_u_scale = 0.;
  /** Write only the seven fields the analysis reads, dropping the LB populations.
   *
   * A full frame is 27 field-equivalents, of which the nine ff components are by far the
   * largest item and the only reason to keep them is restartability. The analysis reads
   * exactly ux_mat, uy_mat, pressure, chi, m, QQxx, QQyx, so a light frame is ~4x smaller
   * -- 850 GB -> 220 GB over a 135-run campaign -- at the price of a frame one cannot
   * restart from. Default 0 (the full, restartable frame). */
  int frame_light = 0;
  std::ofstream vid_u, vid_p, vid_m, vid_chi, vid_meta;
  bool video_open = false;

  /** LAGRANGIAN TRACERS.
   *
   * The memory obeys D_t m = (g(P) - m)/tau_m along a MATERIAL path, so the clock tau_m
   * has to be compared against is the decorrelation time of P following a cell -- not the
   * Eulerian one at a fixed node, and not the instantaneous 1/lambda that cw_calib calls
   * t_eddy (a strain rate is not a correlation time). Tracers measure it directly.
   *
   * WHY IN THE MODEL AND NOT IN POST-PROCESSING. Integrating a trajectory needs the
   * velocity at every step. The full frames are 5 t_eddy apart, and the video stream is
   * block-averaged and quantised to uint8 -- the quantisation noise integrates into a
   * random walk and the block average removes exactly the small scales that decorrelate a
   * path. Storing u densely enough to do it afterwards costs more than the simulation
   * itself, writing already being ~92% of the runtime. Here it is one interpolation per
   * tracer per sub-stage against 10^6 nodes, i.e. free.
   *
   * ntracer is the sampling interval in steps; tracer-count is rounded to a near-square
   * grid. Tracers are NOT serialized: a restart re-seeds them. */
  unsigned ntracer = 0, tracer_count = 0;
  std::vector<double> tr_x, tr_y, tr_xN, tr_yN;
  unsigned tracer_nx = 0, tracer_ny = 0;
  std::ofstream trc_dat, trc_meta;
  bool tracer_open = false;

  /** Steps completed, i.e. calls to Step(). Drives chi_freeze_steps.
   *
   * Equals the runcard time only at nsubsteps = 1, which every runcard in this project
   * uses; with nsubsteps > 1 this counts sub-steps and chi-freeze-steps means sub-steps
   * too. Initialize() rejects nsubsteps > 1 together with a non-zero freeze rather than
   * leaving that ambiguous. */
  unsigned nstep_done = 0;

  /** Cached tau_chi > 0; keeps the branch out of the inner loop */
  bool use_switching = false;

  /** Free energy, molecular field, stress and pressures at one node */
  void UpdateNodeQuantities(unsigned);
  /** u_mat = [u_code + div(sigma)/(2n)]/[1 + friction/(2n)] at one node */
  void ComputeNodeMaterialVelocity(unsigned);
  /** Predictor-corrector update of Q, chi, m and the LB populations at one node */
  void UpdateNodeFields(unsigned, bool);
  /** Whole-lattice sweep for the material velocity */
  void ComputeMaterialVelocity();
  /** Set up the initial phenotype and memory */
  void ConfigurePhenotype();
  /** Diffuse a field until its correlation length is `length` (the noise initialisers) */
  void SmoothToLength(ScalarField&, double length);
  /** chi*(m), the target of the phenotype relaxation */
  double ChiStar(double) const;
  /** g(P), the target of the memory relaxation */
  double MemoryTarget(unsigned) const;
  /** Clamp a field back into [0,1] */
  static void ClampUnit(ScalarField&, unsigned);
  /** zeta_eff at one node: the constant zeta_open in open loop, the floored law otherwise */
  double Activity(unsigned) const;
  /** Block-average one field onto the video lattice and quantise it into `out` */
  void VideoPack(std::vector<unsigned char>& out, const ScalarField&,
                 double lo, double hi) const;
  /** Open the video streams (truncating) on the first write */
  void VideoOpen(const std::string& dir);
  /** Seed the tracers on a stratified near-square grid (deterministic, no seed) */
  void ConfigureTracers();
  /** Advance the tracers with the SAME stage structure as chi and m */
  void AdvanceTracers(bool first);
  /** Bilinear interpolation of a field at an off-lattice point (periodic wrap) */
  double InterpolateField(const ScalarField&, double x, double y) const;
  /** Open the tracer streams (truncating) on the first write */
  void TracerOpen(const std::string& dir);
  /** Sample P and m along the tracers, on the ntracer clock */
  void WriteTracers(const std::string& dir, unsigned t);

  virtual void UpdateQuantities();
  virtual void UpdateFields(bool);
  virtual void BoundaryConditionsFields();
  virtual void BoundaryConditionsFields2();

public:
  ConfluentWet(unsigned, unsigned, unsigned);
  virtual ~ConfluentWet();

  // functions from base class Model
  virtual void Initialize();
  virtual void Configure();
  virtual void ConfigureAtNode(unsigned);
  virtual void Step();
  virtual void RuntimeChecks();
  virtual void WriteAuxiliary(const std::string&, unsigned);
  virtual option_list GetOptions();

  /** Serialization of parameters (do not change; append only) */
  template<class Archive>
  void serialize_params(Archive& ar)
  {
    Nematic::serialize_params(ar);

    ar & auto_name(Dbio)
       & auto_name(tau_chi)
       & auto_name(chi_width)
       & auto_name(mc)
       & auto_name(switch_sign)
       & auto_name(chi_config)
       & auto_name(chi0)
       & auto_name(chi_noise)
       & auto_name(chi_length)
       & auto_name(tau_m)
       & auto_name(pmem)
       & auto_name(pmem_width)
       & auto_name(m0)
       & auto_name(director_config)
       & auto_name(defect_sep)
       & auto_name(init_order)
       // appended 2026-09-01 (memory campaign); never reorder the entries above
       & auto_name(zeta0_frac)
       & auto_name(open_loop)
       & auto_name(zeta_open)
       & auto_name(chi_lo)
       & auto_name(chi_hi)
       & auto_name(m_lo)
       & auto_name(m_hi)
       & auto_name(chi_block)
       & auto_name(chi_freeze_steps)
       & auto_name(nvideo)
       & auto_name(video_stride)
       & auto_name(video_p_scale)
       & auto_name(video_u_scale)
       & auto_name(frame_light)
       // appended 2026-09-01 (step-switch campaign); never reorder the entries above
       & auto_name(chi_seed);
  }

  /** Serialization of the current frame (time snapshot)
   *
   * ff is kept because it is the only complete state of the LB fluid and hence the only
   * honest restart point. n, ux/uy (the bare LB moment) and ux_mat/uy_mat (the physical
   * velocity) are all stored on top of it: reconstructing them in python is slow, and
   * ux_mat in particular cannot be reconstructed from ff alone -- it needs div(sigma). */
  template<class Archive>
  void serialize_frame(Archive& ar)
  {
    // The seven fields every analysis reads. Always written.
    ar & auto_name(QQxx)
       & auto_name(QQyx)
       & auto_name(chi)
       & auto_name(m)
       & auto_name(ux_mat)
       & auto_name(uy_mat)
       & auto_name(pressure);

    // Everything else: the restart state and the stress decomposition. Skipped by
    // frame_light, which is safe because the frame archive is write-only (the iarchive
    // side of serialize_frame is declared but unimplemented, and nothing in the C++ or the
    // python ever reads a frame back into a model).
    if(!frame_light)
      ar & auto_name(ff)
         & auto_name(n)
         & auto_name(ux)
         & auto_name(uy)
         & auto_name(sigmaXX)
         & auto_name(sigmaYY)
         & auto_name(sigmaXY)
         & auto_name(sigmaYX)
         & auto_name(sigma_bulk)
         & auto_name(sigma_active_xx)
         & auto_name(sigma_active_yx)
         & auto_name(pressure_lb);
  }
};

#endif//MODELS_CONFLUENT_WET_HPP_
