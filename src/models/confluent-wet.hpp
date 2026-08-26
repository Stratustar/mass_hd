#ifndef MODELS_CONFLUENT_WET_HPP_
#define MODELS_CONFLUENT_WET_HPP_

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
 *   active stress = -zeta (1 - chi) Q
 *   LB     : f_v <- f_v + (f^eq_v - f_v)/tau + w_v (F.c_v) , F = div(sigma) - friction u
 *
 *   D_t chi = (chi*(m) - chi)/tau_chi + Dbio lap chi , chi*(m)= .5(1+s tanh((m-mc)/chi_width))
 *   D_t m   = (g(P)    - m)/tau_m     + Dbio lap m   , g(P)   = .5(1+  tanh((P-pmem)/pmem_width))
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

  /** Biomass (cell-motility) diffusivity, applied IDENTICALLY to chi and m: a cell
   *  carries its phenotype and its memory together. Defaults to 0. */
  double Dbio = 0.;
  /** Phenotype relaxation time. tau_chi <= 0 disables switching (frozen-chi control). */
  double tau_chi = 0.;
  /** Width and centre of the tanh in chi*(m) */
  double chi_width = 0.15, mc = 0.5;
  /** +1: sustained compression drives chi UP, towards passive/grow (contact inhibition,
   *  a self-limiting negative feedback since the activity is zeta*(1-chi)). -1 is the
   *  opposite branch. */
  int switch_sign = 1;
  /** Phenotype initialisation: mode, mean, std and correlation length of the noise */
  std::string chi_config = "uniform";
  double chi0 = 0.5, chi_noise = 0., chi_length = 0.;

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
  /** chi*(m), the target of the phenotype relaxation */
  double ChiStar(double) const;
  /** g(P), the target of the memory relaxation */
  double MemoryTarget(unsigned) const;
  /** Clamp a field back into [0,1] */
  static void ClampUnit(ScalarField&, unsigned);

  virtual void UpdateQuantities();
  virtual void UpdateFields(bool);
  virtual void BoundaryConditionsFields();
  virtual void BoundaryConditionsFields2();

public:
  ConfluentWet(unsigned, unsigned, unsigned);

  // functions from base class Model
  virtual void Initialize();
  virtual void Configure();
  virtual void ConfigureAtNode(unsigned);
  virtual void Step();
  virtual void RuntimeChecks();
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
       & auto_name(init_order);
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
    ar & auto_name(ff)
       & auto_name(QQxx)
       & auto_name(QQyx)
       & auto_name(chi)
       & auto_name(m)
       & auto_name(n)
       & auto_name(ux)
       & auto_name(uy)
       & auto_name(ux_mat)
       & auto_name(uy_mat)
       & auto_name(sigmaXX)
       & auto_name(sigmaYY)
       & auto_name(sigmaXY)
       & auto_name(sigmaYX)
       & auto_name(sigma_bulk)
       & auto_name(sigma_active_xx)
       & auto_name(sigma_active_yx)
       & auto_name(pressure)
       & auto_name(pressure_lb);
  }
};

#endif//MODELS_CONFLUENT_WET_HPP_
