#ifndef MODELS_NEMATIC_COMPRESSIBLE_DIVISION_HPP_
#define MODELS_NEMATIC_COMPRESSIBLE_DIVISION_HPP_

#include "models/nematic.hpp"

class NematicCompressibleDivision : public Nematic
{
protected:
  /** Mask that marks sites belonging to active growth patches. */
  std::vector<bool> division_mask;

  /** Growth strength at zero pressure. */
  double alpha = 0.;
  /** Probability to seed a new growth center when masks are refreshed. */
  double division_rate = 0.;
  /** Number of timesteps during which the current mask is kept. */
  unsigned division_time = 1u;
  /** Radius of a circular growth patch. */
  double division_radius = 0.;
  /** Pressure threshold above which growth switches off. */
  double P_critical = 1.;
  /** Amplitude of the additional compressive pressure term. */
  double pressure_A = 0.;
  /** Density scale entering the exponential pressure term. */
  double rho_critical = 1.;
  /** Background density outside the seeded active region. */
  double rho_background = 1e-3;
  /** Width of the diffuse interface used by config=circle. */
  double interface_width = 2.;

  /** Counter used to decide when a new set of growth patches is created. */
  unsigned division_count = 0u;

  /** Create or refresh the current set of growth patches. */
  void SetMasks();
  /** Extra isotropic pressure beyond the default LB equation of state. */
  virtual double GetExtraPressure(double) const;
  /** Total local pressure used to suppress growth. */
  double GetPressure(double) const;
  /** Local effective growth rate inside active patches. */
  double GetGrowthRate(unsigned) const;
  /** Update one node, adding the extra isotropic pressure contribution. */
  void UpdateQuantitiesAtNodeDivision(unsigned);
  /** Update one node and inject newly created mass at the local velocity. */
  void UpdateFieldsAtNodeDivision(unsigned, bool);

  virtual void UpdateQuantities() override;
  virtual void UpdateFields(bool) override;

public:
  NematicCompressibleDivision(unsigned, unsigned, unsigned);

  virtual void ConfigureAtNode(unsigned) override;
  virtual void Initialize() override;
  virtual void Step() override;
  virtual void RuntimeChecks() override;
  virtual option_list GetOptions() override;

  template<class Archive>
  void serialize_params(Archive& ar)
  {
    Nematic::serialize_params(ar);

    ar & auto_name(alpha)
       & auto_name(division_rate)
       & auto_name(division_time)
       & auto_name(division_radius)
       & auto_name(P_critical)
       & auto_name(pressure_A)
       & auto_name(rho_critical)
       & auto_name(rho_background)
       & auto_name(interface_width);
  }

  template<class Archive>
  void serialize_frame(Archive& ar)
  {
    Nematic::serialize_frame(ar);
    ar & auto_name(division_mask);
  }
};

#endif//MODELS_NEMATIC_COMPRESSIBLE_DIVISION_HPP_
