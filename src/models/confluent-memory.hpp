#ifndef MODELS_CONFLUENT_MEMORY_HPP_
#define MODELS_CONFLUENT_MEMORY_HPP_

#include "models/lyotropic-with-division.hpp"

/** Confluent proliferating active nematic with mechanical memory.
 *
 * A standalone model for the "cell confluent layer with memory" project. It is
 * deliberately INDEPENDENT of the go-or-grow line (GoOrGrow / DryGoOrGrow): it
 * derives only from LyotropicWithDivision, i.e. from the grid / field / boundary
 * / serialization / option infrastructure and the division-death masks, and
 * writes its own free energy, stress, velocity and field updates. Changes to the
 * go-or-grow models therefore cannot perturb this one, and vice versa.
 *
 * CONFLUENT means: no free boundary, no vacuum, no phase field. The double well
 * AA*phi^2*(1-phi)^2 and the interface term KK*|grad phi|^2 are dropped (both are
 * forced to zero and a non-zero value is rejected), because their minima at phi=0
 * and phi=1 are exactly what tears a dense layer into holes. The density is held
 * together by the one-sided crowding modulus B alone; rarefied patches are filled
 * back in by the pressure gradient. The free-energy density is
 *
 *   f = 1/2 B [(phi - phi_critical)_+]^2        crowding (one-sided)
 *     + 1/2 CC (Snem*phi - q^2)^2               nematic bulk, density coupled
 *     + 1/2 LL |grad Q|^2                       Frank elasticity
 *
 * DRY: the flow is overdamped, friction*u = div(sigma); there is no lattice
 * Boltzmann step at all. The inherited ff/fn arrays are allocated by the base
 * Initialize() but never read or written here.
 *
 * THREE TRANSPORTED SCALARS. phi is the cell density; the phenotype chi and the
 * mechanical memory m are intensive fields carried by it. Both are integrated in
 * their EXTENSIVE form so that transport is exactly conservative:
 *
 *   d_t phi     + div J           = R                                  , J = phi*u - GammaP*grad(mu)
 *   d_t(phi*chi) + div(chi*J)     = S_switch + Dchi*(...) + chi_growth*R
 *   d_t(phi*m)   + div(m*J)       = phi*(g(P) - m)/tau_m + m*R
 *
 * with P = -sigma_bulk the mechanical pressure (with AA=KK=0 there is no surface
 * contribution, so sigma_bulk is already the surface-free half-trace) and
 *
 *   g(P) = 1/2 [1 + tanh((P - pmem)/pmem_width)]      (pmem_width=0: sharp step)
 *
 * of unit amplitude, so m lies in [0,1] and is the exponentially weighted fraction
 * of recent time the material element spent above the pressure threshold. The
 * memory feeds back on the phenotype through
 *
 *   tau_chi * D_t chi = chi*(m) - chi,   chi*(m) = 1/2 [1 + s*tanh((m-mc)/chi_width)]
 *
 * i.e. the SAME first-order low-pass form as the memory itself, so the loop
 * P -> m -> chi -> activity -> flow -> P is a cascade of two lags whose times
 * tau_m and tau_chi are exact and directly comparable to each other and to the
 * flow time.
 *
 * Because the layer is confluent there is no vacuum, so none of the material
 * masking / neighbour-fallback machinery of GoOrGrow is needed: chi and m are
 * simply phichi/phi and phim/phi everywhere. If phi ever falls below phi-floor
 * the confluent assumption has broken and the run stops with an error rather than
 * silently producing garbage.
 */
class ConfluentMemory : public LyotropicWithDivision
{
protected:
  /** Phenotype: transported density phi*chi, its predictor and buffer, and the
   *  derived grow fraction chi = phichi/phi. */
  ScalarField phichi, phichiN, phichi_tmp, chi;
  /** Mechanical memory: transported density phi*m, its predictor and buffer, and
   *  the derived memory m = phim/phi. */
  ScalarField phim, phimN, phim_tmp, m;

  /** Crowding (one-sided) modulus; the only thing holding the layer together */
  double B = 1.;
  /** Preferred value of 1/2 Tr(Q^2) per unit density in the nematic bulk energy */
  double Snem = 1.;
  /** Phenotype diffusion coefficient */
  double Dchi = 0.;
  /** Phenotype switching, written as a RELAXATION so that its timescale is exact.
   *
   * V(chi,m) = (1/2)(chi - chi*(m))^2 / tau_chi gives
   *     tau_chi * D_t chi = chi*(m) - chi,
   *     chi*(m) = .5*(1 + switch_sign*tanh((m - mc)/chi_width)),
   * the exact mirror of the memory equation tau_m * D_t m = g(P) - m, so the two lags
   * of the loop P -> m -> chi -> activity -> flow -> P are the same kind of object and
   * can be compared directly. The earlier linear form V = Ochi*(m-mc)*chi made dV/dchi
   * independent of chi, so chi drifted at a constant rate and only stopped at the
   * [0,phi] clamp; its effective timescale then depended on |m-mc| and was measured to
   * be 2.2-9.0x longer than 1/|Ochi|, with the discrepancy itself correlated with
   * tau_m -- fatal for a study whose whole point is the ordering of the timescales.
   * switch_sign = +1: sustained compression drives chi UP, towards the passive/grow
   * phenotype (contact inhibition); since the activity is zeta*phi*(1-chi) that is a
   * negative, self-limiting feedback. switch_sign = -1 is the opposite branch.
   * tau_chi <= 0 disables switching entirely (the frozen-chi control). */
  double tau_chi = 0., chi_width = 0.15, mc = 0.5;
  int switch_sign = 1;
  /** Whether the phi*chi growth source is chi*R (chi-preserving) or R */
  int growTogether = 1;
  /** Phenotype initialization: mean, std and correlation length of the noise */
  double chi0 = 0.5, chi_noise = 0., chi_length = 0.;
  /** Phenotype initialization mode ("noise" or "uniform") */
  std::string chi_config = "noise";
  /** Memory relaxation time */
  double tau_m = 1.;
  /** Pressure threshold and smoothing width of the memory response g(P) */
  double pmem = 0., pmem_width = 0.;
  /** Initial memory */
  double m0 = 0.;
  /** Density below which the confluent assumption is declared broken */
  double phi_floor = 1e-3;
  /** Optional snapshot to seed Q, phi, phichi and phim from */
  std::string init_frame = "";
  /** With init_frame, reset chi uniform (=chi0) instead of loading it */
  int init_frame_uniform_chi = 0;

  /** Cached tau_chi > 0; keeps the string-free branch out of the inner loop */
  bool use_switching = false;

  /** Smallest phi seen in the last UpdateQuantities (confluence diagnostic) */
  double phi_min = 0.;

  /** Free energy, chemical potential, molecular field and stress at one node */
  void UpdateNodeQuantities(unsigned);
  /** Overdamped velocity friction*u = div(sigma) at one node */
  void ComputeNodeVelocity(unsigned);
  /** Predictor-corrector update of Q, phi, phichi and phim at one node */
  void UpdateNodeFields(unsigned, bool);
  /** Whole-lattice sweeps */
  void ComputeVelocity();
  /** Recompute chi and m from the transported densities */
  void UpdateDerivedFields();
  /** Set up the correlated-random initial phenotype and the initial memory */
  void ConfigurePhenotype();
  /** Seed Q, phi, phichi and phim from a saved frame */
  void ConfigureFromFrame();
  /** phi-weighted memory source phi*(g(P) - m)/tau_m */
  double MemorySource(unsigned, double) const;
  /** phi-weighted switching source phi*(chi*(m) - chi)/tau_chi */
  double SwitchingSource(double, double) const;
  /** Upwind value of an intensive field on the face (k, neighbour) */
  double UpwindFace(const ScalarField&, unsigned, unsigned, double) const;
  /** Clamp a phi-carried density back into [0, phi] */
  void ProjectDensity(ScalarField&, const ScalarField&);

  virtual void UpdateQuantities();
  virtual void UpdateFields(bool);
  virtual void BoundaryConditionsFields();

public:
  ConfluentMemory(unsigned, unsigned, unsigned);

  // functions from base class Model
  virtual void Initialize();
  virtual void Configure();
  virtual void Step();
  virtual option_list GetOptions();

  /** Serialization of parameters (do not change; append only) */
  template<class Archive>
  void serialize_params(Archive& ar)
  {
    LyotropicWithDivision::serialize_params(ar);

    ar & auto_name(B)
       & auto_name(Snem)
       & auto_name(Dchi)
       & auto_name(tau_chi)
       & auto_name(chi_width)
       & auto_name(switch_sign)
       & auto_name(mc)
       & auto_name(growTogether)
       & auto_name(chi0)
       & auto_name(chi_noise)
       & auto_name(chi_length)
       & auto_name(chi_config)
       & auto_name(tau_m)
       & auto_name(pmem)
       & auto_name(pmem_width)
       & auto_name(m0)
       & auto_name(phi_floor)
       & auto_name(init_frame)
       & auto_name(init_frame_uniform_chi);
  }

  /** Serialization of the current frame (time snapshot)
   *
   * Written explicitly rather than by chaining to the base class: this model has
   * no lattice Boltzmann step, so the nine-component ff array -- by far the
   * largest item in a Lyotropic frame -- carries no information and is omitted. */
  template<class Archive>
  void serialize_frame(Archive& ar)
  {
    ar & auto_name(QQxx)
       & auto_name(QQyx)
       & auto_name(phi)
       & auto_name(phichi)
       & auto_name(chi)
       & auto_name(phim)
       & auto_name(m)
       & auto_name(ux)
       & auto_name(uy)
       & auto_name(sigmaXX)
       & auto_name(sigmaYY)
       & auto_name(sigmaXY)
       & auto_name(sigmaYX)
       & auto_name(sigma_bulk)
       & auto_name(sigma_active_xx)
       & auto_name(sigma_active_yx)
       & auto_name(division_mask)
       & auto_name(death_mask);
  }
};

#endif//MODELS_CONFLUENT_MEMORY_HPP_
