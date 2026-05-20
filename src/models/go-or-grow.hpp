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
  double Achi = 0., Ochi = 0., phiswitch = 0.;
  /** Whether m growth follows the local grow fraction */
  int growTogether = 0;
  /** Grow-type density and derived phenotype fraction */
  ScalarField m, mN, m_tmp, chi;
  /** Phenotype initial mean, variance and correlation length */
  double chi0 = 0., chi_noise = 0., chi_length = 0.;
  /** Phenotype initialization mode */
  std::string chi_config = "noise";
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
  /** Update derived phenotype fraction */
  void UpdatePhenotypeQuantities();
  /** Project m back to a density consistent with phi */
  void ProjectM();
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
       & auto_name(phiswitch)
       & auto_name(growTogether)
       & auto_name(chi0)
       & auto_name(chi_noise)
       & auto_name(chi_length)
       & auto_name(chi_config)
       & auto_name(relax_steps)
       & auto_name(relax_dt)
       & auto_name(relax_phi)
       & auto_name(relax_Q);
  }

  /** Serialization of the current frame (time snapshot) */
  template<class Archive>
  void serialize_frame(Archive& ar)
  {
    LyotropicWithDivision::serialize_frame(ar);
    ar & auto_name(m);
    ar & auto_name(chi);
  }
};

#endif//MODELS_GO_OR_GROW_HPP_
