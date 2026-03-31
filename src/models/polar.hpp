#ifndef MODELS_POLAR_HPP_
#define MODELS_POLAR_HPP_

#include "models.hpp"

class Polar : public Model
{
protected:
  /** Lattice Boltzmann distribution
   *
   * Written as ff[k][v], where k is node and v is the direction.
   * */
  LBField ff, fn, ff_tmp, fn_tmp;
  /** P-Tensor */
  ScalarField Px, Py, PNx, PNy, Px_nem, Py_nem;

  /** Velocity */
  ScalarField ux, uy;
  /** Density */
  ScalarField n;
  /** Molecular field */
  ScalarField Hx, Hy;
  /** Derivatives */
  ScalarField dxPx, dyPx, dxPy, dyPy;
  /** Stress tensor */
  ScalarField SigmaXX, SigmaXY, SigmaYY;

  /** Nematic region level, concentration, director inclination, inital noise,
   * intial radius of the circle */
  double angle_deg, angle, noise, init_order = 1;
  /** Fluid properties */
  double rho = 40., tau = 1.;
  /** Polarisation dynamics */
  double xi, eta, nu, gamma, K, J, A, zeta, alpha;
  /** Number of correction steps in the predictor/corrector scheme */
  unsigned npc = 1;
  /** Sum of f (for checking purposes) */
  double ftot = 0, fcheck;

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
  virtual void BoundaryConditionsPol();
  /** Boundary Conditions for the secondary fields */
  virtual void BoundaryConditionsFields();
  /** Move the LB particles */
  void Move();

  /** Project polar to nematic at node */
  void Project(unsigned);

public:
  Polar(unsigned, unsigned, unsigned);

  /** Configure a single node
   *
   * This allows to change the way the arrays are configured in derived
   * classes, see for example PolarFreeBoundary.
   * */
  virtual void ConfigureAtNode(unsigned);

  // functions from base class Model
  virtual void Initialize();
  virtual void Step();
  virtual void Configure();
  virtual void RuntimeChecks();
  virtual void RuntimeStats();
  virtual option_list GetOptions();

  /** Serialization of parameters (do not change) */
  template<class Archive>
  void serialize_params(Archive& ar)
  {
    ar
       & auto_name(angle)
       & auto_name(noise)
       & auto_name(rho)
       & auto_name(gamma)
       & auto_name(xi)
       & auto_name(zeta)
       & auto_name(alpha)
       & auto_name(tau)
       & auto_name(K)
       & auto_name(init_order);
  }

  /** Serialization of the current frame (time snapshot) */
  template<class Archive>
  void serialize_frame(Archive& ar)
  {
    ar & auto_name(ff)
       & auto_name(Px)
       & auto_name(Py);
  }
};

#endif//MODELS_POLAR_HPP_
