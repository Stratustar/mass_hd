#include "header.hpp"
#include "models/nematic.hpp"
#include "error_msg.hpp"
#include "random.hpp"
#include "lb.hpp"
#include "tools.hpp"

using namespace std;
namespace opt = boost::program_options;

// from main.cpp:
extern unsigned nthreads, nsubsteps;
extern double time_step;

Nematic::Nematic(unsigned LX, unsigned LY, unsigned BC)
  : Model(LX, LY, BC, BC==0 ? GridType::Periodic : GridType::Layer)
{}

void Nematic::Initialize()
{
  // initialize variables
  angle = angle_deg*M_PI/180.;

  // allocate memory
  ff.SetSize(LX, LY, Type);
  fn.SetSize(LX, LY, Type);
  ff_tmp.SetSize(LX, LY, Type);
  fn_tmp.SetSize(LX, LY, Type);
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
  sigmaXX.SetSize(LX, LY, Type);
  sigmaYY.SetSize(LX, LY, Type);
  sigmaYX.SetSize(LX, LY, Type);
  sigmaXY.SetSize(LX, LY, Type);

  if(nsubsteps>1)
    throw error_msg("time stepping not implemented for this model"
                    ", please set nsubsteps=1.");

  Q_fluct = (Q_kBT!=0);
  u_fluct = (u_kBT!=0);
}

void Nematic::ConfigureAtNode(unsigned k)
{
  double nematicOrder = 1.;
  const unsigned x = GetXPosition(k);
  const unsigned y = GetYPosition(k);

  if(init_config.empty() || init_config=="uniform")
  {
    nematicOrder = 1.;
  }
  else if(init_config=="circle")
  {
    nematicOrder = (pow(diff(x, LX/2), 2) + pow(diff(y, LY/2), 2)
                   <= radius*radius) ? 1. : 0.;
  }
  else if(init_config=="half")
  {
    nematicOrder = (y < level) ? 1. : 0.;
  }
  else
    throw error_msg("error: initial configuration '", init_config, "' unknown.");

  // theta is the angle of the director inside the nematic region
  const double theta = angle + noise*M_PI*(random_real()-.5);
  QQxx[k] = nematicOrder*cos(2*theta);
  QQyx[k] = nematicOrder*sin(2*theta);
  // equilibrium dist
  ux[k] = uy[k] = 0;
  n[k]  = rho;
  ff[k] = GetEquilibriumDistribution(ux[k], uy[k], n[k]);
  // compute totals for later checks
  ftot  = accumulate(begin(ff[k]), end(ff[k]), ftot);
}

void Nematic::Configure()
{
  for(unsigned k=0; k<DomainSize; ++k)
    ConfigureAtNode(k);
}

void Nematic::UpdateQuantitiesAtNode(unsigned k)
{
  // array placeholders for current node
  const auto& d = get_neighbours(k);
  const auto& f = ff[k];
  // Q-tensor
  const double Qxx = QQxx[k];
  const double Qyx = QQyx[k];

  // compute velocities
  const double nn = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];
  const double vx = nn>1e-12
    ? (f[1] - f[2] + f[5] - f[6] - f[7] + f[8])/nn
    : 0.;
  const double vy = nn>1e-12
    ? (f[3] - f[4] + f[5] - f[6] + f[7] - f[8])/nn
    : 0.;

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

  // computation of sigma...
  // ... on-diagonal stress components
  const double sigmaB = .5*CC*term*term;
  const double sigmaF = 2*xi*( (Qxx*Qxx-1.)*Hxx + Qxx*Qyx*Hyx )
                        - zeta*Qxx
                        + LL*(dyQxx*dyQxx+dyQyx*dyQyx-dxQxx*dxQxx-dxQyx*dxQyx);
  // .. off-diagonal stress components
  const double sigmaS = 2*xi*(Qyx*Qxx*Hxx + (Qyx*Qyx-1)*Hyx) - zeta*Qyx
                        - 2*LL*(dxQxx*dyQxx+dxQyx*dyQyx);
  const double sigmaA = 2*(Qxx*Hyx - Qyx*Hxx);

  // transfer to arrays
  n[k]       =  nn;
  ux[k]      =  vx;
  uy[k]      =  vy;
  HHxx[k]    =  Hxx;
  HHyx[k]    =  Hyx;
  dxQQxx[k]  =  dxQxx;
  dxQQyx[k]  =  dxQyx;
  dyQQxx[k]  =  dyQxx;
  dyQQyx[k]  =  dyQyx;
  sigmaXX[k] =  sigmaF + sigmaB;
  sigmaYY[k] = -sigmaF + sigmaB;
  sigmaXY[k] =  sigmaS + sigmaA;
  sigmaYX[k] =  sigmaS - sigmaA;
}

void Nematic::UpdateFieldsAtNode(unsigned k, bool first)
{
  // pointer to neighbours
  const auto& d = get_neighbours(k);

  // store data in arrays in local varibles
  const double nn = n[k];
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

  // laplacian of the chem pot
  // derivatives of velocity
  const double dxux = derivX(ux, d, sB);
  const double dyux = derivY(ux, d, sB);
  const double dxuy = derivX(uy, d, sB);
  const double dyuy = derivY(uy, d, sB);
  // derivatives of stress field
  const double dxSxx = derivX(sigmaXX, d, sB);
  const double dySxy = derivY(sigmaXY, d, sB);
  const double dxSyx = derivX(sigmaYX, d, sB);
  const double dySyy = derivY(sigmaYY, d, sB);

  // expansion/compression term
  const double expansion = dxux + dyuy;
  // shear
  const double shear     = .5*(dxuy + dyux);
  // vorticity
  const double vorticity = .5*(dxuy - dyux);
  // trace
  const double traceQL   = Qxx*(dxux - dyuy) + 2*Qyx*shear;

  // xx component of the Beris-Edwards equation
  const double Dxx = Gamma*Hxx - vx*dxQxx - vy*dyQxx - 2*vorticity*Qyx
    + xi*((Qxx+1)*(2*dxux-traceQL) +2*Qyx*shear -expansion);
  // zz component of the Beris-Edwards equation
  const double Dyx = Gamma*Hyx - vx*dxQyx - vy*dyQyx + 2*vorticity*Qxx
    + xi*( Qyx*(expansion-traceQL) + 2*shear);

  // forcing term
  const double Fx = dxSxx + dySxy - friction*vx;
  const double Fy = dxSyx + dySyy - friction*vy;

  // calculate the equilibrium distribution fe
  const auto fe = GetEquilibriumDistribution(vx, vy, nn);

  if(first)
  {
    double Qxx_noise = 0., Qxy_noise = 0.;
    LBNode ff_noise = {0.};

    if(Q_fluct)
    {
      static const double Q_stren = sqrt(Gamma*Q_kBT);
      Qxx_noise = Q_stren*random_real();
      Qxy_noise = Q_stren*random_real();
    }

    if(u_fluct)
    {
      const auto ff_stren = sqrt(3.*nn*u_kBT*(2.*tau-1.)/tau/tau);
      ff_noise = GenerateNoiseDistribution(ff_stren);
    }

    QNxx[k] = QQxx[k] + .5*Dxx + Qxx_noise;
    QNyx[k] = QQyx[k] + .5*Dyx + Qxy_noise;
    QQxx[k] = QNxx[k] + .5*Dxx;
    QQyx[k] = QNyx[k] + .5*Dyx;

    for(unsigned v=0; v<lbq; ++v)
    {
      fn[k][v] = ff[k][v] + .5*((fe[v]-ff[k][v])/tau
          + w[v]*(Fx*xdir(v) + Fy*ydir(v)))
          + ff_noise[v];
      ff[k][v] = fn[k][v] + .5*((fe[v]-ff[k][v])/tau
          + w[v]*(Fx*xdir(v) + Fy*ydir(v)));
    }
  }
  else
  {
    QQxx[k] = QNxx[k] + .5*Dxx;
    QQyx[k] = QNyx[k] + .5*Dyx;

    for(unsigned v=0; v<lbq; ++v)
      ff[k][v] = fn[k][v] + .5*((fe[v]-ff[k][v])/tau
          + w[v]*(Fx*xdir(v) + Fy*ydir(v)));
  }
}

void Nematic::UpdateQuantities()
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    UpdateQuantitiesAtNode(k);
}

void Nematic::UpdateFields(bool first)
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    UpdateFieldsAtNode(k, first);
}

void Nematic::Move()
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<TotalSize; ++k)
  {
    for(unsigned v=0; v<lbq; ++v)
    {
      // advect particles
      ff_tmp[next(k, v)][v] = ff[k][v];
      fn_tmp[next(k, v)][v] = fn[k][v];
    }
  }
  // swap temp variables
  swap(ff, ff_tmp);
  swap(fn, fn_tmp);
}

void Nematic::BoundaryConditionsLB()
{
  switch(BC)
  {
    // pbc without bdry layer (nothing to do)
    case 0:
      break;
    // free-slip channel
    case 1:
    {
      auto apply_bc = [](LBField& field) {
        // pbc on the left and right walls
        field.ApplyPBC(PBCWall::LeftRight);
        // Free-slip on the front and back walls
        field.ApplyFreeSlip(Wall::Front);
        field.ApplyFreeSlip(Wall::Back);
        // corners
        field.ApplyFreeSlip(Corner::RightBack, Wall::Back);
        field.ApplyFreeSlip(Corner::RightFront, Wall::Front);
        field.ApplyFreeSlip(Corner::LeftBack, Wall::Back);
        field.ApplyFreeSlip(Corner::LeftFront, Wall::Front);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    // pbc with boundary layer
    default:
      ff.ApplyPBC();
      fn.ApplyPBC();
  }
}

void Nematic::BoundaryConditionsFields()
{
  switch(BC)
  {
    // pbc without bdry layer (nothing to do)
    case 0:
      break;
    // free-slip channel
    case 1:
    {
      auto apply_bc = [](ScalarField& field) {
        // pbc on the left and right walls
        field.ApplyPBC(PBCWall::LeftRight);
        // Neumann on the fron and back walls
        field.ApplyNeumann(Wall::Front);
        field.ApplyNeumann(Wall::Back);
        // corners
        field.ApplyNeumann(Corner::RightBack, Wall::Back);
        field.ApplyNeumann(Corner::RightFront, Wall::Front);
        field.ApplyNeumann(Corner::LeftBack, Wall::Back);
        field.ApplyNeumann(Corner::LeftFront, Wall::Front);
      };

      apply_bc(QQxx);
      apply_bc(QQyx);
      break;
    }
    // pbc with bdry layer
    default:
      QQxx.ApplyPBC();
      QQyx.ApplyPBC();
  }
}

void Nematic::BoundaryConditionsFields2()
{
  switch(BC)
  {
    // pbc without bdry layer (nothing to do)
    case 0:
      break;
    // free-slip channel
    case 1:
    {
      auto apply_bc = [](ScalarField& field) {
        // pbc on the left and right walls
        field.ApplyPBC(PBCWall::LeftRight);
        // Neumann on the fron and back walls
        field.ApplyNeumann(Wall::Front);
        field.ApplyNeumann(Wall::Back);
        // corners
        field.ApplyNeumann(Corner::RightBack, Wall::Back);
        field.ApplyNeumann(Corner::RightFront, Wall::Front);
        field.ApplyNeumann(Corner::LeftBack, Wall::Back);
        field.ApplyNeumann(Corner::LeftFront, Wall::Front);
      };

      apply_bc(ux);
      apply_bc(uy);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }
    // pbc with bdry layer
    default:
      ux.ApplyPBC();
      uy.ApplyPBC();
      sigmaXX.ApplyPBC();
      sigmaYY.ApplyPBC();
      sigmaYX.ApplyPBC();
      sigmaXY.ApplyPBC();
  }
}

void Nematic::Step()
{
  // boundary conditions for primary fields
  BoundaryConditionsFields();
  // predictor step
  UpdateQuantities();
  // boundary conditions for the secondary
  // fields MU, u, sigma, and phi
  BoundaryConditionsFields2();
  // update all the fields
  this->UpdateFields(true);
  // boundary conditions for
  // the flow before advection
  BoundaryConditionsLB();
  // move LB particles
  Move();
  // corrector steps
  for(unsigned n=1; n<=npc; ++n)
  {
    // same thing (no advection)
    BoundaryConditionsFields();
    UpdateQuantities();
    BoundaryConditionsFields2();
    this->UpdateFields(false);
  }
}

void Nematic::RuntimeChecks()
{
  // check that the sum of f is constant
  {
    double fcheck = 0;
    for(unsigned k=0; k<DomainSize; ++k)
        fcheck = accumulate(begin(ff[k]), end(ff[k]), fcheck);
    if(abs(ftot-fcheck)>1)
      throw error_msg("f is not conserved (", ftot, "/", fcheck, ")");
  }
}

option_list Nematic::GetOptions()
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
    ("Q_kBT", opt::value<double>(&Q_kBT),
     "hydrodynamic fluctations strength")
    ("u_kBT", opt::value<double>(&u_kBT),
     "nematic fluctuations strength")
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
    ("config", opt::value<string>(&init_config)->default_value("uniform"),
     "initial configuration: uniform, circle, half")
    ("level", opt::value<double>(&level),
     "height of the nematic region for config=half")
    ("angle", opt::value<double>(&angle_deg),
     "initial angle to x direction (in degrees)")
    ("radius", opt::value<double>(&radius),
     "radius of the initial circle for config=circle")
    ("noise", opt::value<double>(&noise),
     "size of initial variations");

  return { model_options, config_options };
}
