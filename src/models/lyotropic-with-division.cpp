#include "models/lyotropic-with-division.hpp"
#include "error_msg.hpp"
#include "random.hpp"
#include "lb.hpp"
#include "tools.hpp"

using namespace std;
namespace opt = boost::program_options;

// from main.cpp:
extern unsigned nthreads, nsubsteps;
extern double time_step;

LyotropicWithDivision::LyotropicWithDivision(unsigned LX_, unsigned LY_, unsigned BC)
  : Lyotropic(LX_, LY_, BC)
{}

void LyotropicWithDivision::Initialize()
{
  // base class init
  Lyotropic::Initialize();

  division_mask.assign(DomainSize, false);
  death_mask.assign(DomainSize, false);

  // initialize the counters
  division_count = death_count = 0u;

  // check if phi should be conserved: yes if either strength/rate=0 for division and death, else no
  // (time=0 applies masks for 1 step, radius=0 applies masks for 1 lattice point, so phi is not conserved)
  conserve_phi = ( (alpha==0 || division_rate==0) && ( beta==0 || death_rate==0 ) );

  if(P_critical<=0)
    throw error_msg("P_critical must be > 0.");
}

void LyotropicWithDivision::SetMasks()
{
  if(not division_count--)
  {
    // reinit masks and counters
    division_mask.assign(DomainSize, false);
    division_count = division_time;
    // create a temporary list of centers
    vector<array<unsigned, 2>> centers;
    for(unsigned k=0; k<DomainSize; ++k)
    {
      if(phi[k]>=.5)
      {
        const double pressure = GetCrowdingPressure(phi[k]);
        const double adjusted_division_rate = pressure < P_critical ? division_rate * phi[k] : 0.;
        if(random_real()<adjusted_division_rate)
          centers.emplace_back(array<unsigned,2>{GetXPosition(k), GetYPosition(k)});
      }
    }
    // draw the circles
    for(unsigned k=0; k<DomainSize; ++k)
    {
      // corrdinates of the point
      const unsigned yk = GetYPosition(k);
      const unsigned xk = GetXPosition(k);

      // set the corresponding mask to 1 if the point falls is close enough from
      // a center
      for(const auto& c : centers)
      {
        if(pow(wrap(diff(yk, c[1]), LY), 2) + pow(wrap(diff(xk, c[0]), LX), 2)
            <= ceil(division_radius*division_radius))
        {
          division_mask[k] = true;
          break;
        }
      }
    }
  }

  // yes this is shamelessly copy-pasted
  if(not death_count--)
  {
    // reinit masks and counters
    death_mask.assign(DomainSize, false);
    death_count = death_time;
    // create a temporary list of centers
    vector<array<unsigned, 2>> centers;
    for(unsigned k=0; k<DomainSize; ++k)
    {
      if(phi[k]>=.5)
      {
        const double adjusted_death_rate = death_rate * phi[k];
        if(random_real()<adjusted_death_rate)
          centers.emplace_back(array<unsigned,2>{GetXPosition(k), GetYPosition(k)});
      }
    }
    // draw the circles
    for(unsigned k=0; k<DomainSize; ++k)
    {
      // corrdinates of the point
      const unsigned yk = GetYPosition(k);
      const unsigned xk = GetXPosition(k);
      // set the corresponding mask to 1 if the point falls is close enough from
      // a center
      for(const auto& c : centers)
      {
        if(pow(wrap(diff(yk, c[1]), LY), 2) + pow(wrap(diff(xk, c[0]), LX), 2)
            <= ceil(death_radius*death_radius))
        {
          death_mask[k] = true;
          break;
        }
      }
    }
  }
}

void LyotropicWithDivision::UpdateFields(bool first)
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
  {
    const double pressure = GetCrowdingPressure(phi[k]);
    const double division_m = pressure < P_critical ? ( 1. ? division_mask[k] : 0. ) : 0.;
    const double death_m = ( 1. ? death_mask[k] : 0. );

    // normal update
    Lyotropic::UpdateFieldsAtNode(k, first);

    // we then correct for division and death
    phi_tmp[k] += phi[k] * ( alpha*division_m - beta*death_m );
  }

  swap(phi, phi_tmp);
}

option_list LyotropicWithDivision::GetOptions()
{
  // get options from base model
  auto options = Lyotropic::GetOptions();

  // add new model options
  options[0].add_options()
    ("alpha", opt::value<double>(&alpha), "division strength")
    ("beta", opt::value<double>(&beta), "death strength")
    ("division-rate", opt::value<double>(&division_rate), "division rate")
    ("division-time", opt::value<double>(&division_time), "division life time")
    ("division-radius", opt::value<double>(&division_radius), "division radius")
    ("death-rate", opt::value<double>(&death_rate), "death rate")
    ("death-time", opt::value<double>(&death_time), "death life time")
    ("death-radius", opt::value<double>(&death_radius), "death radius")
    ("P_critical", opt::value<double>(&P_critical)->default_value(1.0), "critical crowding pressure above which division stops");

  return options;
}

void LyotropicWithDivision::Step()
{
  // set the growth rate
  SetMasks();

  // normal step
  Lyotropic::Step();
}
