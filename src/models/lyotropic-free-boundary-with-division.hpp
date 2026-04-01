#ifndef MODELS_LYOTROPIC_FREE_BOUNDARY_WITH_DIVISION_HPP_
#define MODELS_LYOTROPIC_FREE_BOUNDARY_WITH_DIVISION_HPP_

#include "models/lyotropic-free-boundary.hpp"

class LyotropicFreeBoundaryWithDivision : public LyotropicFreeBoundary
{
protected:
  /** Mask for growth
   *
   * This is the term that enters the phi equation as growth_term * phi. This
   * term depends on the lattice site and is not homogeneous across the domain.
   * The masks define what regions are growing/shrinking.
   * */
  std::vector<bool> division_mask;
  /** Weirdest variable name ever */
  std::vector<bool> death_mask;

  /** Strength of division and death */
  double alpha, beta;
  /** Division rate */
  double division_rate;
  /** Life-time of a division events */
  double division_time;
  /** Radius of a patch created by a division event */
  double division_radius;
  /** Death rate */
  double death_rate;
  /** Life-time of a death events */
  double death_time;
  /** Radius of a patch created by a death event */
  double death_radius;
  /** Critical crowding pressure at which division stops */
  double P_critical;

  /** Time counters for shuffeling the patches */
  unsigned division_count, death_count;

public:
  LyotropicFreeBoundaryWithDivision(unsigned, unsigned, unsigned);

  // functions from base class Model
  virtual void Initialize();
  virtual void Step();
  virtual option_list GetOptions();
  virtual void UpdateFields(bool);

  /** Update the masks for division and death events
   *
   * The masks are created by choosing randomly N*density sites and creating
   * a circular patch around it. Update occurs every division_time or death_time
   * steps.
   * */
  void SetMasks();

  /** Serialization of parameters (do not change) */
  template<class Archive>
  void serialize_params(Archive& ar)
  {
    // serialize from base class
    LyotropicFreeBoundary::serialize_params(ar);

    // serialize new parameters
    ar & auto_name(alpha)
       & auto_name(beta)
       & auto_name(division_rate)
       & auto_name(division_time)
       & auto_name(division_radius)
       & auto_name(death_rate)
       & auto_name(death_time)
       & auto_name(death_radius)
       & auto_name(P_critical);
  }

  /** Serialization of the current frame (time snapshot) */
  template<class Archive>
  void serialize_frame(Archive& ar)
  {
    // serialize from base class
    LyotropicFreeBoundary::serialize_frame(ar);
    // serialize masks as well, for testing purposes
    ar & auto_name(division_mask)
       & auto_name(death_mask);
  }
};

#endif//MODELS_LYOTROPIC_FREE_BOUNDARY_WITH_DIVISION_HPP_
