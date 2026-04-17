#include "header.hpp"
#include "models/nematic-compressible-division.hpp"
#include "error_msg.hpp"
#include "random.hpp"
#include "lb.hpp"
#include "tools.hpp"

using namespace std;
namespace opt = boost::program_options;

// from main.cpp:
extern unsigned nthreads;
extern double time_step;

NematicCompressibleDivision::NematicCompressibleDivision(unsigned LX,
                                                         unsigned LY,
                                                         unsigned BC)
  : Nematic(LX, LY, BC)
{}

void NematicCompressibleDivision::ConfigureAtNode(unsigned k)
{
  double nematicOrder = 1.;
  double densityOrder = 1.;
  double theta = angle + noise*M_PI*(random_real()-.5);
  const unsigned x = GetXPosition(k);
  const unsigned y = GetYPosition(k);

  if(init_config.empty() || init_config=="uniform")
  {
    nematicOrder = 1.;
    densityOrder = 1.;
  }
  else if(init_config=="circle")
  {
    const double dx = static_cast<double>(diff(x, LX/2));
    const double dy = static_cast<double>(diff(y, LY/2));
    const double rr = sqrt(dx*dx + dy*dy);

    if(interface_width > 0.)
    {
      const double profile = 0.5*(1. - tanh((rr - radius)/interface_width));
      nematicOrder = profile;
      densityOrder = profile;
    }
    else
    {
      const bool inside = (rr <= radius);
      nematicOrder = inside ? 1. : 0.;
      densityOrder = inside ? 1. : 0.;
    }
    theta = M_PI*random_real();
  }
  else if(init_config=="half")
  {
    const bool inside = (y < level);
    nematicOrder = inside ? 1. : 0.;
    densityOrder = inside ? 1. : 0.;
  }
  else
    throw error_msg("error: initial configuration '", init_config, "' unknown.");

  QQxx[k] = nematicOrder*cos(2*theta);
  QQyx[k] = nematicOrder*sin(2*theta);
  ux[k] = uy[k] = 0.;
  n[k] = densityOrder*rho + (1. - densityOrder)*rho_background;
  ff[k] = GetEquilibriumDistribution(ux[k], uy[k], n[k]);
  ftot = accumulate(begin(ff[k]), end(ff[k]), ftot);
}

void NematicCompressibleDivision::Initialize()
{
  Nematic::Initialize();

  division_mask.assign(DomainSize, false);
  division_count = 0u;

  if(P_critical <= 0.)
    throw error_msg("P_critical must be > 0.");
  if(rho_critical <= 0.)
    throw error_msg("rho_critical must be > 0.");
  if(rho_background < 0.)
    throw error_msg("rho_background must be >= 0.");
  if(interface_width < 0.)
    throw error_msg("interface_width must be >= 0.");
  if(division_time == 0u)
    throw error_msg("division_time must be >= 1.");
}

double NematicCompressibleDivision::GetExtraPressure(double nn) const
{
  return pressure_A*exp(nn/rho_critical);
}

double NematicCompressibleDivision::GetPressure(double nn) const
{
  return nn/3. + GetExtraPressure(nn);
}

double NematicCompressibleDivision::GetGrowthRate(unsigned k) const
{
  if(!division_mask[k])
    return 0.;

  const double pressure = GetPressure(n[k]);
  if(pressure >= P_critical)
    return 0.;

  return alpha*(1. - pressure/P_critical);
}

void NematicCompressibleDivision::SetMasks()
{
  if(!division_count--)
  {
    division_mask.assign(DomainSize, false);
    division_count = division_time;

    vector<array<unsigned, 2>> centers;
    for(unsigned k=0; k<DomainSize; ++k)
    {
      const double order = QQxx[k]*QQxx[k] + QQyx[k]*QQyx[k];

      if(order > 1e-6
         && GetPressure(n[k]) < P_critical
         && random_real() < division_rate)
        centers.push_back({ GetXPosition(k), GetYPosition(k) });
    }

    const int radius = static_cast<int>(ceil(division_radius));
    const int radius_sq = static_cast<int>(ceil(division_radius*division_radius));

    for(const auto& c : centers)
    {
      for(int dx=-radius; dx<=radius; ++dx)
      {
        for(int dy=-radius; dy<=radius; ++dy)
        {
          if(dx*dx + dy*dy > radius_sq)
            continue;

          int x = static_cast<int>(c[0]) + dx;
          int y = static_cast<int>(c[1]) + dy;

          if(BC == 0)
          {
            x = modu(x, static_cast<int>(LX));
            y = modu(y, static_cast<int>(LY));
          }
          else if(BC == 1)
          {
            x = modu(x, static_cast<int>(LX));
            if(y < 0 || y >= static_cast<int>(LY))
              continue;
          }
          else if(x < 0 || x >= static_cast<int>(LX)
               || y < 0 || y >= static_cast<int>(LY))
          {
            continue;
          }

          division_mask[GetDomainIndex(static_cast<unsigned>(x),
                                       static_cast<unsigned>(y))] = true;
        }
      }
    }
  }
}

void NematicCompressibleDivision::UpdateQuantitiesAtNodeDivision(unsigned k)
{
  Nematic::UpdateQuantitiesAtNode(k);

  const double pressure_extra = GetExtraPressure(n[k]);
  sigmaXX[k] -= pressure_extra;
  sigmaYY[k] -= pressure_extra;
}

void NematicCompressibleDivision::UpdateFieldsAtNodeDivision(unsigned k,
                                                             bool first)
{
  const double nn = n[k];
  const double vx = ux[k];
  const double vy = uy[k];
  const double growth_rate = GetGrowthRate(k);

  Nematic::UpdateFieldsAtNode(k, first);

  if(growth_rate == 0.)
    return;

  const double delta_rho = time_step*growth_rate*nn;
  if(delta_rho == 0.)
    return;

  const auto growth_distribution =
    GetEquilibriumDistribution(vx, vy, delta_rho);

  if(first)
  {
    for(unsigned v=0; v<lbq; ++v)
    {
      fn[k][v] += growth_distribution[v];
      ff[k][v] += growth_distribution[v];
    }
  }
  else
  {
    for(unsigned v=0; v<lbq; ++v)
      ff[k][v] += growth_distribution[v];
  }
}

void NematicCompressibleDivision::UpdateQuantities()
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    UpdateQuantitiesAtNodeDivision(k);
}

void NematicCompressibleDivision::UpdateFields(bool first)
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    UpdateFieldsAtNodeDivision(k, first);
}

void NematicCompressibleDivision::Step()
{
  SetMasks();
  Nematic::Step();
}

void NematicCompressibleDivision::RuntimeChecks()
{
}

option_list NematicCompressibleDivision::GetOptions()
{
  auto options = Nematic::GetOptions();

  options[0].add_options()
    ("alpha", opt::value<double>(&alpha),
     "zero-pressure proliferation strength")
    ("division-rate", opt::value<double>(&division_rate),
     "probability to seed a proliferation patch when masks are refreshed")
    ("division-time", opt::value<unsigned>(&division_time),
     "number of steps for which a proliferation mask is kept")
    ("division-radius", opt::value<double>(&division_radius),
     "radius of a proliferation patch")
    ("P_critical", opt::value<double>(&P_critical)->default_value(1.0),
     "critical pressure above which proliferation stops")
    ("pressure-A", opt::value<double>(&pressure_A),
     "amplitude of the additional pressure term A exp(rho/rho_critical)")
    ("rho_critical", opt::value<double>(&rho_critical)->default_value(1.0),
     "density scale entering the additional pressure term")
    ("rho_background", opt::value<double>(&rho_background)->default_value(1e-3),
     "background density outside the seeded active region")
    ("interface_width", opt::value<double>(&interface_width)->default_value(2.0),
     "width of the diffuse interface used for config=circle");

  return options;
}
