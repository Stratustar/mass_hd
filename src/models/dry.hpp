#ifndef MODELS_DRY_HPP_
#define MODELS_DRY_HPP_

#include "models.hpp"

class Dry : public Model
{
protected:
  /** Lattice Boltzmann distribution
   *
   * Written as ff[k][v], where k is node and v is the direction.
   * */
  LBField ff;
  /** Q-Tensor */
  ScalarField QQxx, QNxx, QQyx, QNyx;

  /** Velocity */
  ScalarField ux, uy;
  /** Density */
  ScalarField n;
  /** Molecular field */
  ScalarField HHxx, HHyx;
  /** Derivatives */
  ScalarField dxQQxx, dyQQxx, dxQQyx, dyQQyx, del2QQxx, del2QQyx;
  /** Stress tensor */
  ScalarField sigmaXX, sigmaYY, sigmaYX, sigmaXY;

  /** Initial angle and noise */
  double angle_deg, angle, noise;
  /** Fluid density */
  double rho = 40.;
  /** Fluid parameters */
  double Gamma, xi, tau, friction, LL, CC, zeta, eta, gam;
  /** Intial configuration */
  std::string init_config;
  /** Number of correction steps in the predictor/corrector scheme */
  unsigned npc = 1;
  /** Sum of f (for checking purposes) */
  double ftot = 0;

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
  /** Compute velocity from force balance */
  void ComputeVelocity();
  /** UpdateFields() implementation */
  void UpdateFieldsAtNode(unsigned, bool);
  /** UpdateQuantities() implementation */
  void UpdateQuantitiesAtNode(unsigned);
  /** ComputeVelocity() implementation */
  void ComputeVelocityAtNode(unsigned);

public:
  Dry(unsigned, unsigned, unsigned);

  /** Configure a single node
   *
   * This allows to change the way the arrays are configured in derived
   * classes, see for example DryFreeBoundary.
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
    ar
       & auto_name(angle)
       & auto_name(noise)
       & auto_name(rho)
       & auto_name(Gamma)
       & auto_name(eta)
       & auto_name(xi)
       & auto_name(zeta)
       & auto_name(tau)
       & auto_name(friction)
       & auto_name(LL);
  }

  /** Serialization of the current frame (time snapshot) */
  template<class Archive>
  void serialize_frame(Archive& ar)
  {
    ar & auto_name(ff)
       & auto_name(QQxx)
       & auto_name(QQyx);
  }
};

#endif//MODELS_DRY_HPP_
