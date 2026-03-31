#include "models/lyotropic-free-boundary-with-division.hpp"
#include "error_msg.hpp"
#include "random.hpp"
#include "lb.hpp"
#include "tools.hpp"

using namespace std;
namespace opt = boost::program_options;

// from main.cpp:
extern unsigned nthreads, nsubsteps;
extern double time_step;

LyotropicFreeBoundaryWithDivision::LyotropicFreeBoundaryWithDivision(unsigned LX_, unsigned LY_, unsigned BC_)
  : LyotropicFreeBoundary(LX_, LY_, BC_)
{}

void LyotropicFreeBoundaryWithDivision::Initialize()
{
  // base class init
  LyotropicFreeBoundary::Initialize();

  division_mask.assign(DomainSize, false);
  death_mask.assign(DomainSize, false);

  // initialize the counters
  division_count = death_count = 0u;

  // check if phi should be conserved: yes if either strength/rate=0 for division and death, else no
  // (time=0 applies masks for 1 step, radius=0 applies masks for 1 lattice point, so phi is not conserved)
  conserve_phi = false; //( (alpha==0 || division_rate==0) && ( beta==0 || death_rate==0 ) );
}

void LyotropicFreeBoundaryWithDivision::SetMasks()
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
      if(phi[k]>1e-3 and random_real()<division_rate*phi[k] and !outside[k] and GetYPosition(k)<level)
      {
        centers.emplace_back(array<unsigned,2>{GetXPosition(k), GetYPosition(k)});
      }
    }

    // draw the circles
    for(unsigned k=0; k<DomainSize; ++k)
    {
      if(outside[k]) continue;
      // coordinates of the point
      const unsigned yk = GetYPosition(k);
      const unsigned xk = GetXPosition(k);

      // set the corresponding mask to 1 if the point falls is close enough from
      // a center
      for(const auto& c : centers)
      {
        if (BC==0 or BC==3 or BC==6)
        {
          if(pow(wrap(diff(yk, c[1]), LY), 2) + pow(wrap(diff(xk, c[0]), LX), 2)
              <= ceil(division_radius*division_radius))
          {
            division_mask[k] = true;
            break;
          }
        }
        else if (BC==2 or BC==5 or BC==11 or BC==12)
        {
          if(pow(diff(yk, c[1]), 2) + pow(wrap(diff(xk, c[0]), LX), 2)
              <= ceil(division_radius*division_radius))
          {
            division_mask[k] = true;
            break;
          }
        }
        else
        {
          if(pow(diff(yk, c[1]), 2) + pow(diff(xk, c[0]), 2)
              <= ceil(division_radius*division_radius))
          {
            division_mask[k] = true;
            break;
          }
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
      if(phi[k]>=.5 and !outside[k])
      {
        const double adjusted_death_rate = death_rate * phi[k];
        if(random_real()<adjusted_death_rate)
          centers.emplace_back(array<unsigned,2>{GetXPosition(k), GetYPosition(k)});
      }
    }
    // draw the circles
    for(unsigned k=0; k<DomainSize; ++k)
    {
      if(outside[k]) continue;
      // corrdinates of the point
      const unsigned yk = GetYPosition(k);
      const unsigned xk = GetXPosition(k);
      // set the corresponding mask to 1 if the point falls is close enough from
      // a center
      for(const auto& c : centers)
      {
        if (BC==0 or BC==3 or BC==6)
        {
          if(pow(wrap(diff(yk, c[1]), LY), 2) + pow(wrap(diff(xk, c[0]), LX), 2)
              <= ceil(death_radius*death_radius))
          {
            death_mask[k] = true;
            break;
          }
        }
        else if (BC==2 or BC==5 or BC==11 or BC==12)
        {
          if(pow(diff(yk, c[1]), 2) + pow(wrap(diff(xk, c[0]), LX), 2)
              <= ceil(death_radius*death_radius))
          {
            death_mask[k] = true;
            break;
          }
        }
        else
        {
          if(pow(diff(yk, c[1]), 2) + pow(diff(xk, c[0]), 2)
              <= ceil(death_radius*death_radius))
          {
            death_mask[k] = true;
            break;
          }
        }
      }
    }
  }
}

void LyotropicFreeBoundaryWithDivision::UpdateFields(bool first)
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
  {
    if(outside[k]) continue;

    const double division_m = ( 1. ? division_mask[k] : 0. );
    const double death_m = ( 1. ? death_mask[k] : 0. );

    // normal update
    Lyotropic::UpdateFieldsAtNode(k, first);

    // we then correct for division and death
    phi_tmp[k] += phi[k]*alpha*division_m*(1.-phi[k]/phi_critical) - phi[k]*beta*death_m;
  }

  swap(phi.get_data(), phi_tmp.get_data());
}

option_list LyotropicFreeBoundaryWithDivision::GetOptions()
{
  // get options from base model
  auto options = LyotropicFreeBoundary::GetOptions();

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
    ("phi-critical", opt::value<double>(&phi_critical), "phi critical");

  return options;
}

void LyotropicFreeBoundaryWithDivision::Step()
{
  // set the growth rate
  SetMasks();
  // the rest is just a normal Lyotropic::Step()
  LyotropicFreeBoundary::Step();
}
