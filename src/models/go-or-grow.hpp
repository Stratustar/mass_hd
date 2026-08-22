#ifndef MODELS_GO_OR_GROW_HPP_
#define MODELS_GO_OR_GROW_HPP_

#include "models/lyotropic-with-division.hpp"

class GoOrGrow : public LyotropicWithDivision
{
protected:
  /** Crowding/compressibility penalty strength */
  double B = 0.;
  /** Preferred value of 1/2 Tr(Q^2) in the nematic bulk free energy */
  double Snem = 1.;
  /** Phenotype diffusion coefficient */
  double Dchi = 0.;
  /** Phenotype switching potential parameters */
  double Achi = 0., Ochi = 0., pswitch = 0.;
  /** Whether phichi growth follows the local grow fraction */
  int growTogether = 0;
  /** Whether phase-field surface terms contribute to stress and velocity */
  int surface_stress = 1;
  /** Grow-type density and derived phenotype fraction */
  ScalarField phichi, phichiN, phichi_tmp, chi;
  /** Mechanical-memory density (phi*m) and the derived memory field m
   *
   * m is an intensive, biomass-carried memory of mechanical loading obeying
   *   tau_m * D_t m = g(P) - m,   g(P) = .5*(1+tanh((P-pmem)/pmem_width)),
   * so with the g amplitude fixed to 1 the memory m lies in [0,1] and measures
   * the recent fraction of time this material element spent above the pressure
   * threshold. D_t is the BIOMASS material derivative: the transported quantity
   * is the extensive phim = phi*m, carried by the full phi flux
   * J = phi*v - GammaP*grad(mu), exactly as phichi = phi*chi is. */
  ScalarField phim, phimN, phim_tmp, m;
  /** Isotropic stress (half-trace) WITHOUT the phase-field surface contribution
   *
   * sigma_bulk minus (f_surface - mu_surface*phi). Since mu_surface carries the
   * only non-local piece (-KK*lap(phi)), this is a purely local function of
   * (Qxx, Qyx, phi) at the node. P_mem = -sigma_bulk_nosurf is the pressure that
   * drives the mechanical-memory source g(P). */
  ScalarField sigma_bulk_nosurf;
  /** Which variable drives phenotype switching: "pressure" or "memory" */
  std::string switch_drive = "pressure";
  /** Mechanical-memory relaxation time */
  double tau_m = 1.;
  /** Pressure threshold and smoothing width of the memory source g(P) */
  double pmem = 0., pmem_width = 0.;
  /** Memory threshold in the switching potential V = Ochi*(m-mc)*chi */
  double mc = 0.5;
  /** Initial memory */
  double m0 = 0.;
  /** Cached switch_drive=="memory"; a string compare per node would be far too slow */
  bool use_memory = false;
  /** Phenotype initial mean, variance and correlation length */
  double chi0 = 0., chi_noise = 0., chi_length = 0.;
  /** Correlated-noise phi init: std of multiplicative modulation and its length */
  double phi_noise = 0., phi_length = 0.;
  /** Phenotype initialization mode */
  std::string chi_config = "noise";
  /** Optional snapshot to seed Q, phi (and chi) from, replacing the procedural
   *  random initialization. Empty string = normal initialization. */
  std::string init_frame = "";
  /** When seeding from init_frame, set chi uniform (=chi0) instead of loading it
   *  from the snapshot (used by the homogeneous control). */
  int init_frame_uniform_chi = 0;
  /** Dry free-energy relaxation before official dynamics */
  unsigned relax_steps = 0;
  double relax_dt = 1.;
  int relax_phi = 0, relax_Q = 0;

  /** Compute chemical potential, stress and derivatives */
  virtual void UpdateQuantities();
  /** UpdateQuantities() implementation */
  void UpdateQuantitiesAtNode(unsigned);
  /** Update fields using predictor-corrector method */
  virtual void UpdateFields(bool);
  /** Boundary Conditions for the fields */
  virtual void BoundaryConditionsFields();
  /** Setup the spatially correlated phenotype field */
  void ConfigurePhenotype();
  /** Seed Q, phi and (unless init_frame_uniform_chi) chi from a saved frame */
  void ConfigureFromFrame();
  /** Overlay a multiplicative correlated-noise modulation on the initial phi */
  void ConfigurePhiNoise();
  /** Update derived phenotype fraction */
  void UpdatePhenotypeQuantities();
  /** Density threshold below which chi is only an unphysical placeholder */
  double PhenotypePhiEpsilon() const;
  /** Whether this lattice node contains phenotype-carrying material */
  bool HasPhenotypeMaterial(unsigned) const;
  /** Whether a face can carry phi/phichi transport because either side has material */
  bool HasMaterialFace(unsigned, unsigned) const;
  /** Remove numerical phi transport across vacuum-vacuum faces */
  double MaskedPhiFaceIncrement(unsigned, unsigned, double) const;
  /** Material value of an intensive field nearby, if any node carries material */
  double LocalMaterialValue(const ScalarField&, unsigned, bool&) const;
  /** Use the preferred node only if it carries material; otherwise fall back */
  double MaterialValue(const ScalarField&, unsigned, unsigned) const;
  /** Intensive value carried by a phi transport increment into/out of a face */
  double TransportFaceValue(const ScalarField&, unsigned, unsigned, double) const;
  /** Wrappers of the above for the two intensive fields carried by phi */
  double MaterialChi(unsigned a, unsigned b) const
  { return MaterialValue(chi, a, b); }
  double MaterialM(unsigned a, unsigned b) const
  { return MaterialValue(m, a, b); }
  double TransportFaceChi(unsigned k, unsigned nb, double dphi) const
  { return TransportFaceValue(chi, k, nb, dphi); }
  double TransportFaceM(unsigned k, unsigned nb, double dphi) const
  { return TransportFaceValue(m, k, nb, dphi); }
  /** phi-weighted memory source phi*(g(P)-m)/tau_m, with P = -sigma_bulk_nosurf */
  double MemorySource(unsigned, double) const;
  /** Project a phi-carried density back into [0, phi] */
  void ProjectDensity(ScalarField&);
  /** Dry free-energy relaxation of phi and Q without hydrodynamics or growth */
  void RelaxFreeEnergy();
  /** Reset LB and velocity fields after dry relaxation */
  void ResetHydrodynamics();

public:
  GoOrGrow(unsigned, unsigned, unsigned);

  // functions from base class Model
  virtual void Initialize();
  virtual void Configure();
  virtual option_list GetOptions();

  /** Serialization of parameters (do not change) */
  template<class Archive>
  void serialize_params(Archive& ar)
  {
    // serialize from base class
    LyotropicWithDivision::serialize_params(ar);

    // serialize new parameters
    ar & auto_name(B)
       & auto_name(Snem)
       & auto_name(Dchi)
       & auto_name(Achi)
       & auto_name(Ochi)
       & auto_name(pswitch)
       & auto_name(growTogether)
       & auto_name(surface_stress)
       & auto_name(chi0)
       & auto_name(chi_noise)
       & auto_name(chi_length)
       & auto_name(chi_config)
       & auto_name(relax_steps)
       & auto_name(relax_dt)
       & auto_name(relax_phi)
       & auto_name(relax_Q)
       & auto_name(phi_noise)
       & auto_name(phi_length)
       & auto_name(init_frame)
       & auto_name(init_frame_uniform_chi)
       & auto_name(switch_drive)
       & auto_name(tau_m)
       & auto_name(pmem)
       & auto_name(pmem_width)
       & auto_name(mc)
       & auto_name(m0);
  }

  /** Serialization of the current frame (time snapshot) */
  template<class Archive>
  void serialize_frame(Archive& ar)
  {
    LyotropicWithDivision::serialize_frame(ar);
    ar & auto_name(phichi);
    ar & auto_name(chi);
    ar & auto_name(phim);
    ar & auto_name(m);
    ar & auto_name(sigma_bulk_nosurf);
  }
};

#endif//MODELS_GO_OR_GROW_HPP_
