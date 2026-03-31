#ifndef NEMATIC_FREE_BOUNDARY_HPP
#define NEMATIC_FREE_BOUNDARY_HPP

#include "models/nematic.hpp"

class NematicFreeBoundary : public Nematic
{
protected:
  /** Mask for arbitrary boundary conditions
   *
   * This vector is basically telling whether a certain point lies inside or
   * outside the active region. Only important for BC==2.
   * */
  std::vector<bool> outside;

  /** For testing, sum over flux and laplacian in phi*/
#ifdef DEBUG
  double countphi_flux=0., countphi_laplace=0;
#endif//DEBUG

public:
  NematicFreeBoundary(unsigned, unsigned, unsigned);

  /** Initialize outside/inside region */
  void ConfigureBoundaries();
  /** Initialize outside and calls Lyotropic::Configure() */
  virtual void Configure();
  /** Same as in Lyotropic, just skip outside nodes */
  virtual void UpdateQuantities();
  /** Same as in Lyotropic, just skip outside nodes */
  void UpdateFields(bool first);
  /** Same as in Lyotropic, just skip outside nodes */
  void RuntimeChecks();
  void BoundaryConditionsLB();
  void BoundaryConditionsFields();
  void BoundaryConditionsFields2();

    /** Serialization of parameters (do not change) */
  template<class Archive>
  void serialize_params(Archive& ar)
  {
    // serialize from base class
    Nematic::serialize_params(ar);

    // serialize new parameters
    ar & auto_name(outside);
  }

  /** Serialization of the current frame (time snapshot) */
  template<class Archive>
  void serialize_frame(Archive& ar)
  {
    //To have correct boundary conditions in written data.
    BoundaryConditionsLB();
    BoundaryConditionsFields();
    BoundaryConditionsFields2();
    // serialize from base class
    Nematic::serialize_frame(ar);
  }
};

#endif//NEMATIC_FREE_BOUNDARY_HPP
