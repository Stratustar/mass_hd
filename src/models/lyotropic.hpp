#ifndef MODELS_LYOTROPIC_HPP_
#define MODELS_LYOTROPIC_HPP_

#include "models.hpp"

class Lyotropic : public Model
{
protected:
  /** Lattice Boltzmann distribution
   *
   * Written as ff[k][v], where k is node and v is the direction.
   * */
  LBField ff, fn, ff_tmp, fn_tmp;
  /** Q-Tensor */
  ScalarField QQxx, QNxx, QQyx, QNyx;
  /** Binary order */
  ScalarField phi, phn, phi_tmp;

  /** Velocity */
  ScalarField ux, uy, ux_phi, uy_phi;
  /** Density */
  ScalarField n;
  /** Molecular field */
  ScalarField HHxx, HHyx, MU;
  /** Derivatives */
  ScalarField dxQQxx, dyQQxx, dxQQyx, dyQQyx;

  /** Nematic region level, concentration, director inclination, inital noise,
   * intial radius of the circle */
  double level, conc, angle_deg, angle, noise, radius;
  /** Total binary phase value*/
  double totalphi=0., countphi=0.;
  /** Fluid density */
  double rho = 40.;
  /** Fluid parameters */
  double GammaP, GammaQ, xi, tauNem, tauIso, friction, LL, KK, AA, CC, zeta, zetaI, phiJ = 1.;
  /** Intial configuration */
  std::string init_config;
  /** Number of correction steps in the predictor/corrector scheme */
  unsigned npc = 1;
  /** Sum of f (for checking purposes) */
  double ftot = 0;
  /** Total phi (for checking purposes) */
  double ptot = 0;
  /** Flag indicating if we need to conserve the phi field */
  bool conserve_phi = false;

  /** Update fields using predictor-corrector method
   *
   * Because of the way the predictor-corrector is implemented this function
   * can be called many times in a row in order to import the numerical
   * stability of the algorithm. Only the first call needs to have the parameter
   * set to true.
   * */
  virtual void UpdateFields(bool);
  /** Compute chemical potential, stress and derivatives */
  virtual void UpdateQuantities();
  /** UpdateFields() implementation */
  void UpdateFieldsAtNode(unsigned, bool);
  /** UpdateQuantities() implementation */
  void UpdateQuantitiesAtNode(unsigned);

  /** Boundary Conditions for the flow */
  virtual void BoundaryConditionsLB();
  /** Boundary Conditions for the fields */
  virtual void BoundaryConditionsFields();
  /** Boundary Conditions for the secondary fields */
  virtual void BoundaryConditionsFields2();
  /** Move the LB particles */
  void Move();
  /** Isotropic crowding pressure added on top of the bulk stress */
  double GetCrowdingPressure(double) const;

public:
  Lyotropic() = default;
  Lyotropic(unsigned, unsigned, unsigned);
  Lyotropic(unsigned, unsigned, unsigned, GridType);
  /** Stress tensor */
  ScalarField sigmaXX, sigmaYY, sigmaYX, sigmaXY;
  ScalarField sigma_bulk, sigma_elastic_xx, sigma_elastic_yx, sigma_phase_field_xx, sigma_phase_field_yx, sigma_active_xx, sigma_active_yx;

  /** Configure a single node
   *
   * This allows to change the way the arrays are configured in derived
   * classes, see for example LyotropicFreeBoundary.
   * */
  virtual void ConfigureAtNode(unsigned);

  // functions from base class Model
  virtual void Initialize();
  virtual void Step();
  virtual void Configure();
  virtual void RuntimeChecks();
  virtual option_list GetOptions();

  /** Serialization of parameters (do not change) */
  template<class Archive>
  void serialize_params(Archive& ar)
  {
    ar & auto_name(level)
       & auto_name(conc)
       & auto_name(angle)
       & auto_name(noise)
       & auto_name(totalphi)
       & auto_name(rho)
       & auto_name(GammaP)
       & auto_name(GammaQ)
       & auto_name(xi)
       & auto_name(zeta)
       & auto_name(zetaI)
       & auto_name(phiJ)
       & auto_name(tauNem)
       & auto_name(tauIso)
       & auto_name(friction)
       & auto_name(LL)
       & auto_name(KK)
       & auto_name(init_config)
       & auto_name(AA)
       & auto_name(CC);
  }

  /** Serialization of the current frame (time snapshot) */
  template<class Archive>
  void serialize_frame(Archive& ar)
  {
    ar & auto_name(ff)
       & auto_name(QQxx)
       & auto_name(QQyx)
       & auto_name(phi)
       & auto_name(sigmaXX)
       & auto_name(sigmaXY)
       & auto_name(sigmaYX)
       & auto_name(sigmaYY)
       & auto_name(sigma_bulk)
       & auto_name(sigma_elastic_xx)
       & auto_name(sigma_elastic_yx)
       & auto_name(sigma_phase_field_xx)
       & auto_name(sigma_phase_field_yx)
       & auto_name(sigma_active_xx)
       & auto_name(sigma_active_yx);
  }
};

#endif//MODELS_LYOTROPIC_HPP_
