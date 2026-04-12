#ifndef MODELS_NEMATIC_COMPRESSIBLE_DIVISION_SOFTPRESSURE_HPP_
#define MODELS_NEMATIC_COMPRESSIBLE_DIVISION_SOFTPRESSURE_HPP_

#include "models/nematic-compressible-division.hpp"

class NematicCompressibleDivisionSoftPressure
  : public NematicCompressibleDivision
{
protected:
  virtual double GetExtraPressure(double) const override;

public:
  NematicCompressibleDivisionSoftPressure(unsigned, unsigned, unsigned);
};

#endif//MODELS_NEMATIC_COMPRESSIBLE_DIVISION_SOFTPRESSURE_HPP_
