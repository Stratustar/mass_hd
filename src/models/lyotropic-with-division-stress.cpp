#include "models/lyotropic-with-division-stress.hpp"
#include "error_msg.hpp"
#include "random.hpp"
#include "lb.hpp"
#include "tools.hpp"

using namespace std;
namespace opt = boost::program_options;

// from main.cpp:
extern unsigned nthreads, nsubsteps;
extern double time_step;

LyotropicWithDivisionStress::LyotropicWithDivisionStress(unsigned LX_, unsigned LY_, unsigned BC)
  : Lyotropic(LX_, LY_, BC)
{}

void LyotropicWithDivisionStress::Initialize()
{
  // base class init
  Lyotropic::Initialize();

  division_mask.assign(DomainSize, false);
  death_mask.assign(DomainSize, false);

  // initialize the counters
  death_count = division_rate_count = 0u;

  // check if phi should be conserved: yes if either strength/rate=0 for division and death, else no
  // (time=0 applies masks for 1 step, radius=0 applies masks for 1 lattice point, so phi is not conserved)
  conserve_phi = ( (alpha==0 || division_rate==0) && ( beta==0 || death_rate==0 ) );

  if(P_critical<=0)
    throw error_msg("P_critical must be > 0.");
}
  // Here we find the maximum position of J2, and use that as our cell division point. 
void LyotropicWithDivisionStress::SetDivision()
{
  if(not division_rate_count--)
  {
    const double pressure_threshold = P_critical;
    //Here we reset the counter
    division_rate_count = division_rate;

    const unsigned margin = 10;
    const unsigned min_x = margin;
    const unsigned max_x = LX - margin;
    const unsigned min_y = margin;
    const unsigned max_y = LY - margin;

    //We calcuate the stress
    double I1, J2, J2t;
    unsigned kt=0;
    J2=0.0;
    for(unsigned k=0; k<DomainSize; ++k)
    {
      I1=sigmaXX[k]+sigmaYY[k];
      J2t=0.5 * (sigmaXX[k]-I1/2.0) * (sigmaXX[k]-I1/2.0)+0.5 * (sigmaYY[k]-I1/2.0) * (sigmaYY[k]-I1/2.0)+sigmaXY[k] * sigmaYX[k];
      
      const unsigned xk = GetXPosition(k);
      const unsigned yk = GetYPosition(k);
      
      const double pressure = GetCrowdingPressure(phi[k]);
      if (pressure < pressure_threshold && J2t > J2 && xk >= min_x && xk <= max_x && yk >= min_y && yk <= max_y)
      {
        J2=J2t;
        kt=k;
      }
    }
    //Division happens at point of largest stress if the local crowding pressure permits it.
    if(GetCrowdingPressure(phi[kt]) < pressure_threshold)
    {
      division_count.push_back(division_time);
      centers.push_back(array<unsigned,2>{GetXPosition(kt), GetYPosition(kt)});
    }
    
  }

}

void LyotropicWithDivisionStress::SetMasks()
{
  //if(not division_count--)
  {
    // reinit masks and counters
    division_mask.assign(DomainSize, false);
    // draw the circles
    for(unsigned k=0; k<DomainSize; ++k)
    {
      // corrdinates of the point
      const unsigned yk = GetYPosition(k);
      const unsigned xk = GetXPosition(k);

      // set the corresponding mask to 1 if the point falls is close enough from
      // a center
      //for(const auto& c : centers)
          
      for (unsigned t=0; t<division_count.size(); ++t)
      {
        if (division_count[t]>0){
            if(pow(wrap(diff(yk, centers[t][1]), LY), 2) + pow(wrap(diff(xk, centers[t][0]), LX), 2)
                <= ceil(division_radius*division_radius))
            {     
            division_mask[k] = true;
            break;
            }
        }
      }
    }
    //Then count back to make sure that division only takes so long. 
    for (unsigned t=0; t<division_count.size(); ++t)
    {
      if (division_count[t]>0){
        division_count[t]=division_count[t]-1;
      }
    }
  }

  // yes this is shamelessly copy-pasted. In this version, we ignore death.
  if(not death_count--)
  {
  /*  // reinit masks and counters
    death_mask.assign(DomainSize, false);
    death_count = death_time;
    // create a temporary list of centers
    //vector<array<unsigned, 2>> centers;
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
    }*/
  }
}

void LyotropicWithDivisionStress::UpdateFields(bool first)
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

option_list LyotropicWithDivisionStress::GetOptions()
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

void LyotropicWithDivisionStress::Step()
{
  // set the growth rate
  SetDivision();
  SetMasks();

  // normal step
  Lyotropic::Step();
}
