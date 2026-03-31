#include "header.hpp"
#include "models/dry.hpp"
#include "error_msg.hpp"
#include "random.hpp"
#include "lb.hpp"
#include "tools.hpp"

using namespace std;
namespace opt = boost::program_options;

// from main.cpp:
extern unsigned nthreads, nsubsteps;
extern double time_step;

Dry::Dry(unsigned LX, unsigned LY, unsigned BC)
  : Model(LX, LY, 0, GridType::Periodic)
{
  if(BC!=0)
    throw error_msg("model requires PBC.");
}

void Dry::Initialize()
{
  // initialize variables
  angle = angle_deg*M_PI/180.;

  // allocate memory
  ff.SetSize(LX, LY, Type);
  QQxx.SetSize(LX, LY, Type);
  QQyx.SetSize(LX, LY, Type);
  QNxx.SetSize(LX, LY, Type);
  QNyx.SetSize(LX, LY, Type);
  n.SetSize(LX, LY, Type);
  ux.SetSize(LX, LY, Type);
  uy.SetSize(LX, LY, Type);
  HHxx.SetSize(LX, LY, Type);
  HHyx.SetSize(LX, LY, Type);
  dxQQxx.SetSize(LX, LY, Type);
  dyQQxx.SetSize(LX, LY, Type);
  dxQQyx.SetSize(LX, LY, Type);
  dyQQyx.SetSize(LX, LY, Type);
  del2QQxx.SetSize(LX, LY, Type);
  del2QQyx.SetSize(LX, LY, Type);

  if(nsubsteps>1)
    throw error_msg("time stepping not implemented for this model"
                    ", please set nsubsteps=1.");
}

void Dry::ConfigureAtNode(unsigned k)
{
  // add noise (for meta-stable configs)
  // theta is the angle of the director
  const double theta = angle + noise*M_PI*(random_real()-.5);
  QQxx[k] = cos(2*theta);
  QQyx[k] = sin(2*theta);
  // equilibrium dist
  ux[k] = uy[k] = 0;
  n[k]  = rho;
  ff[k] = GetEquilibriumDistribution(ux[k], uy[k], n[k]);
  // compute totals for later checks
  ftot  = accumulate(begin(ff[k]), end(ff[k]), ftot);
}

void Dry::Configure()
{
  for(unsigned k=0; k<DomainSize; ++k)
    ConfigureAtNode(k);
}

void Dry::UpdateQuantitiesAtNode(unsigned k)
{
  // array placeholders for current node
  const auto& d = get_neighbours(k);
  //const auto& f = ff[k];
  // Q-tensor
  const double Qxx = QQxx[k];
  const double Qyx = QQyx[k];

  // compute derivatives etc.
  const double del2Qxx  = laplacian(QQxx,  d, sD);
  const double dxQxx    = derivX   (QQxx,  d, sB);
  const double dyQxx    = derivY   (QQxx,  d, sB);
  const double del2Qyx  = laplacian(QQyx,  d, sD);
  const double dxQyx    = derivX   (QQyx,  d, sB);
  const double dyQyx    = derivY   (QQyx,  d, sB);

  // computation of the chemical potential and molecular field...
  // ...term that couples the binary phase to the degree of nematic order
  const double term = 1. - Qxx*Qxx - Qyx*Qyx;
  // ...molecular field
  const double Hxx = CC*term*Qxx + LL*del2Qxx;
  const double Hyx = CC*term*Qyx + LL*del2Qyx;

  // transfer to arrays
  n[k]        = 1;
  HHxx[k]     = Hxx;
  HHyx[k]     = Hyx;
  dxQQxx[k]   = dxQxx;
  dxQQyx[k]   = dxQyx;
  dyQQxx[k]   = dyQxx;
  dyQQyx[k]   = dyQyx;
  del2QQxx[k] = del2Qxx;
  del2QQyx[k] = del2Qyx;
}

void Dry::ComputeVelocityAtNode(unsigned k)
{
  const auto& d = get_neighbours(k);
  const double dxQxx = dxQQxx[k];
  const double dyQxx = dyQQxx[k];
  const double dxQyx = dxQQyx[k];
  const double dyQyx = dyQQyx[k];
  const double deldxQxx = laplacian(dxQQxx, d, sD);
  const double deldyQxx = laplacian(dyQQxx, d, sD);
  const double deldxQyx = laplacian(dxQQyx, d, sD);
  const double deldyQyx = laplacian(dyQQyx, d, sD);

  // velocity is set from simple force balance
  const double vx = -zeta/friction*(dxQxx + dyQyx
                      + eta/friction*(deldxQxx + deldyQyx));
  const double vy = -zeta/friction*(dxQyx - dyQxx
                      + eta/friction*(deldxQyx - deldyQxx));
  // we use the f-field for the output
  ff[k] = GetEquilibriumDistribution(vx, vy, 1);

  ux[k]      =  vx;
  uy[k]      =  vy;
}

void Dry::UpdateFieldsAtNode(unsigned k, bool first)
{
  // pointer to neighbours
  const auto& d = get_neighbours(k);

  // store data in arrays in local varibles
  const double vx = ux[k];
  const double vy = uy[k];
  const double Qxx = QQxx[k];
  const double Qyx = QQyx[k];
  const double Hxx = HHxx[k];
  const double Hyx = HHyx[k];
  const double dxQxx = dxQQxx[k];
  const double dyQxx = dyQQxx[k];
  const double dxQyx = dxQQyx[k];
  const double dyQyx = dyQQyx[k];
  const double dxdxQxx = derivX(dxQQxx, d, sB);
  const double dxdyQxx = derivX(dyQQxx, d, sB);
  const double dydyQxx = derivY(dyQQxx, d, sB);
  const double dxdxQyx = derivX(dxQQyx, d, sB);
  const double dxdyQyx = derivX(dyQQyx, d, sB);
  const double dydyQyx = derivY(dyQQyx, d, sB);

  // laplacian of the chem pot
  // derivatives of velocity
  const double dxux = derivX(ux, d, sB);
  const double dyux = derivY(ux, d, sB);
  const double dxuy = derivX(uy, d, sB);
  const double dyuy = derivY(uy, d, sB);

  // expansion/compression term
  const double expansion = dxux + dyuy;
  // shear
  const double shear     = .5*(dxuy + dyux);
  // vorticity
  const double vorticity = .5*(dxuy - dyux);
  // trace
  const double traceQL   = Qxx*(dxux - dyuy) + 2*Qyx*shear;
  // that new term
  const double Gxx = -zeta*(Qxx*(dxdxQxx-dydyQxx)+2*Qyx*dxdyQxx);
  const double Gyx = -zeta*(Qxx*(dxdxQyx-dydyQyx)+2*Qyx*dxdyQyx);

  // xx component of the Beris-Edwards equation
  const double Dxx = Gamma*(Hxx - gam*laplacian(del2QQxx, d, sD)) + Gxx - 0*vx*dxQxx - 0*vy*dyQxx - 2*vorticity*Qyx
    + xi*((Qxx+1)*(2*dxux-traceQL) +2*Qyx*shear -expansion);
  // zz component of the Beris-Edwards equation
  const double Dyx = Gamma*(Hyx - gam*laplacian(del2QQyx, d, sD)) + Gyx - 0*vx*dxQyx - 0*vy*dyQyx + 2*vorticity*Qxx
    + xi*( Qyx*(expansion-traceQL) + 2*shear);

  if(first)
  {
    QNxx[k] = QQxx[k] + .5*Dxx;
    QNyx[k] = QQyx[k] + .5*Dyx;
    QQxx[k] = QQxx[k] +    Dxx;
    QQyx[k] = QQyx[k] +    Dyx;
  }
  else
  {
    QQxx[k] = QNxx[k] + .5*Dxx;
    QQyx[k] = QNyx[k] + .5*Dyx;
  }
}

void Dry::UpdateQuantities()
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    UpdateQuantitiesAtNode(k);
}

void Dry::ComputeVelocity()
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    ComputeVelocityAtNode(k);
}

void Dry::UpdateFields(bool first)
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    UpdateFieldsAtNode(k, first);
}

void Dry::Step()
{
  // predictor step
  UpdateQuantities();
  ComputeVelocity();
  // update all the fields
  this->UpdateFields(true);
  // corrector steps
  for(unsigned n=1; n<=npc; ++n)
  {
    // same thing (no advection)
    UpdateQuantities();
    ComputeVelocity();
    this->UpdateFields(false);
  }
}

void Dry::RuntimeChecks()
{
}

option_list Dry::GetOptions()
{
  // model specific options
  opt::options_description model_options("Model options");
  model_options.add_options()
    ("Gamma", opt::value<double>(&Gamma),
     "Mobility")
    ("xi", opt::value<double>(&xi),
     "tumbling/aligning parameter")
    ("tau", opt::value<double>(&tau),
     "viscosity")
    ("rho", opt::value<double>(&rho),
     "fluid density")
    ("eta", opt::value<double>(&eta),
     "viscosity parameter")
    ("gamma", opt::value<double>(&gam),
     "regularisation parameter (Delta^2 Q)")
    ("friction", opt::value<double>(&friction),
     "friction")
    ("CC", opt::value<double>(&CC),
     "coupling constant")
    ("LL", opt::value<double>(&LL),
     "elastic constant")
    ("zeta", opt::value<double>(&zeta),
     "activity parameter")
    ("npc", opt::value<unsigned>(&npc),
     "number of correction steps for the predictor-corrector method");

  opt::options_description config_options("Initial configuration options");
  config_options.add_options()
    ("angle", opt::value<double>(&angle_deg),
     "initial angle to x direction (in degrees)")
    ("noise", opt::value<double>(&noise),
     "size of initial variations");

  return { model_options, config_options };
}
