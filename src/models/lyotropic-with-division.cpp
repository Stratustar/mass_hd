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
        const double crowding_factor = division_phi_critical_factor
          ? (1 - phi[k]/phi_critical) : 1.;
        const double adjusted_division_rate = division_rate * phi[k] * crowding_factor;
        if(random_real()<adjusted_division_rate)
          centers.emplace_back(array<unsigned,2>{GetXPosition(k), GetYPosition(k)});
      }
    }
    // draw the circles: stamp a disc of radius division_radius around each
    // center locally, instead of testing every grid node against every center.
    // Identical result to the O(DomainSize * n_centers) scan, but O(n_centers *
    // radius^2) -- the previous version dominated runtime once the colony grew.
    const double thr = ceil(division_radius*division_radius);
    const int Rbox = static_cast<int>(ceil(sqrt(thr)));
    for(const auto& c : centers)
    {
      const long int cx = c[0], cy = c[1];
      for(int dx=-Rbox; dx<=Rbox; ++dx)
        for(int dy=-Rbox; dy<=Rbox; ++dy)
          if(dx*dx + dy*dy <= thr)
          {
            const unsigned x = modu(cx+dx, static_cast<long int>(LX));
            const unsigned y = modu(cy+dy, static_cast<long int>(LY));
            division_mask[GetDomainIndex(x, y)] = true;
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
    // draw the circles: stamp a disc of radius death_radius around each center
    // locally (see the division block above for the rationale).
    const double thr = ceil(death_radius*death_radius);
    const int Rbox = static_cast<int>(ceil(sqrt(thr)));
    for(const auto& c : centers)
    {
      const long int cx = c[0], cy = c[1];
      for(int dx=-Rbox; dx<=Rbox; ++dx)
        for(int dy=-Rbox; dy<=Rbox; ++dy)
          if(dx*dx + dy*dy <= thr)
          {
            const unsigned x = modu(cx+dx, static_cast<long int>(LX));
            const unsigned y = modu(cy+dy, static_cast<long int>(LY));
            death_mask[GetDomainIndex(x, y)] = true;
          }
    }
  }
}

void LyotropicWithDivision::UpdateFields(bool first)
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
  {
    const double division_m = division_mask[k] ? 1. : 0.;
    const double death_m = death_mask[k] ? 1. : 0.;

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
    ("phi-critical", opt::value<double>(&phi_critical), "phi critical")
    ("division-phi-critical-factor", opt::value<int>(&division_phi_critical_factor),
     "whether division sampling uses the (1 - phi/phi-critical) factor");

  return options;
}

void LyotropicWithDivision::Step()
{
  // set the growth rate
  SetMasks();

  // normal step
  Lyotropic::Step();
}
