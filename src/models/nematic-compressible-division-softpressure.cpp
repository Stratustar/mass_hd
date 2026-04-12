#include "header.hpp"
#include "models/nematic-compressible-division-softpressure.hpp"
#include "error_msg.hpp"

using namespace std;

NematicCompressibleDivisionSoftPressure::NematicCompressibleDivisionSoftPressure(
    unsigned LX, unsigned LY, unsigned BC)
  : NematicCompressibleDivision(LX, LY, BC)
{}

double NematicCompressibleDivisionSoftPressure::GetExtraPressure(double nn) const
{
  const double denom = 1. - nn/rho_critical;
  if(denom <= 0.)
    throw error_msg("soft pressure singularity reached: rho >= rho_critical (",
                    nn, " >= ", rho_critical, ")");

  return pressure_A/denom;
}
