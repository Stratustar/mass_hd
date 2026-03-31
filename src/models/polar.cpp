#include "header.hpp"
#include "models/polar.hpp"
#include "error_msg.hpp"
#include "random.hpp"
#include "lb.hpp"
#include "tools.hpp"

using namespace std;
namespace opt = boost::program_options;

// from main.cpp:
extern unsigned nthreads, nsubsteps;
extern double time_step;

Polar::Polar(unsigned LX, unsigned LY, unsigned BC)
  : Model(LX, LY, BC, BC==1 ? GridType::Layer : GridType::Periodic)
{}

void Polar::Initialize()
{
  // initialize variables
  angle = angle_deg*M_PI/180.;

  // allocate memory
  ff.SetSize(LX, LY, Type);
  fn.SetSize(LX, LY, Type);
  ff_tmp.SetSize(LX, LY, Type);
  fn_tmp.SetSize(LX, LY, Type);
  Px.SetSize(LX, LY, Type);
  Py.SetSize(LX, LY, Type);
  //Px_nem.SetSize(LX, LY, Type);
  //Py_nem.SetSize(LX, LY, Type);
  PNx.SetSize(LX, LY, Type);
  PNy.SetSize(LX, LY, Type);
  n.SetSize(LX, LY, Type);
  ux.SetSize(LX, LY, Type);
  uy.SetSize(LX, LY, Type);
  Hx.SetSize(LX, LY, Type);
  Hy.SetSize(LX, LY, Type);
  dxPx.SetSize(LX, LY, Type);
  dyPx.SetSize(LX, LY, Type);
  dxPy.SetSize(LX, LY, Type);
  dyPy.SetSize(LX, LY, Type);
  SigmaXX.SetSize(LX, LY, Type);
  SigmaXY.SetSize(LX, LY, Type);
  SigmaYY.SetSize(LX, LY, Type);

  if(nsubsteps>1)
    throw error_msg("time stepping not implemented for this model"
                    ", please set nsubsteps=1.");
}

void Polar::ConfigureAtNode(unsigned k)
{
  // add noise (for meta-stable configs)
  // theta is the angle of the director
  const double theta = angle + noise*M_PI*2*(random_real() - .5);
  Px[k] = init_order*cos(theta);
  Py[k] = init_order*sin(theta);
  //Project(k);
  // equilibrium dist
  ux[k] = uy[k] = 0;
  n[k]  = rho;
  ff[k] = GetEquilibriumDistribution(ux[k], uy[k], n[k]);
  // compute totals for later checks
  fcheck = ftot  = accumulate(begin(ff[k]), end(ff[k]), ftot);
}

void Polar::Configure()
{
  for(unsigned k=0; k<DomainSize; ++k)
    ConfigureAtNode(k);
}

void Polar::UpdateQuantitiesAtNode(unsigned k)
{
  // array placeholders for current node
  const auto& d = get_neighbours(k);
  const auto& f = ff[k];
  // Polarisation
  const double px = Px[k];
  const double py = Py[k];
  const double ps = px*px+py*py;

  // compute velocities
  const double nn = f[0]+f[1]+f[2]+f[3]+f[4]+f[5]+f[6]+f[7]+f[8];
  const double vx = (f[1]-f[2]+f[5]-f[6]-f[7]+f[8])/nn;
  const double vy = (f[3]-f[4]+f[5]-f[6]+f[7]-f[8])/nn;

  // compute derivatives etc.
  const double dxpx = derivX(Px, d, sB);
  const double dypx = derivY(Px, d, sB);
  const double dxpy = derivX(Py, d, sB);
  const double dypy = derivY(Py, d, sB);

  // computation of the molecular field...
  const double hx = A*(1.-ps)*px + K*laplacian(Px, d, sD)
    + J*( 2.*ps*laplacian(Px, d, sD)
         +2.*px*(dxpx*dxpx+dypx*dypx-dxpy*dxpy-dypy*dypy)
         +4.*py*(dxpx*dxpy+dypx*dypy));
  const double hy = A*(1.-ps)*py + K*laplacian(Py, d, sD)
    + J*( 2.*ps*laplacian(Py, d, sD)
         +2.*py*(dxpy*dxpy+dypy*dypy-dxpx*dxpx-dypx*dypx)
         +4.*px*(dxpx*dxpy+dypx*dypy));

  // computation of sigma...
  const double sigmaXX = -zeta*.5*(px*px-py*py) + nu*.5*(px*hx-py*hy);
  const double sigmaYY = -zeta*.5*(py*py-px*px) + nu*.5*(py*hy-px*hx);
  const double sigmaXY = -zeta*px*py + .5*(nu+1)*px*hy + .5*(nu-1)*py*hx;

  // transfer to arrays
  n[k]  = nn;
  ux[k] = vx;
  uy[k] = vy;
  Hx[k] = hx;
  Hy[k] = hy;
  dxPx[k] = dxpx;
  dxPy[k] = dxpy;
  dyPx[k] = dypx;
  dyPy[k] = dypy;
  SigmaXX[k] = sigmaXX;
  SigmaXY[k] = sigmaXY;
  SigmaYY[k] = sigmaYY;
}

void Polar::UpdateFieldsAtNode(unsigned k, bool first)
{
  // pointer to neighbours
  const auto& d = get_neighbours(k);

  // store data in arrays in local varibles
  const double nn = n[k];
  const double vx = ux[k];
  const double vy = uy[k];
  const double px = Px[k];
  const double py = Py[k];
  const double hx = Hx[k];
  const double hy = Hy[k];
  const double dxpx = dxPx[k];
  const double dypx = dyPx[k];
  const double dxpy = dxPy[k];
  const double dypy = dyPy[k];

  // derivatives of velocity
  const double dxux = derivX(ux, d, sB);
  const double dyux = derivY(ux, d, sB);
  const double dxuy = derivX(uy, d, sB);
  const double dyuy = derivY(uy, d, sB);
  // derivatives of stress field
  const double dxSxx = derivX(SigmaXX, d, sB);
  const double dySyy = derivY(SigmaYY, d, sB);
  const double dySxy = derivY(SigmaXY, d, sB);
  const double dxSyx = derivX(SigmaXY, d, sB);

  // corrections to the polarisation
  const double Dx = hx/gamma - vx*dxpx - vy*dypx
                    -nu*(dxux*px + .5*(dyux+dxuy)*py) - .5*(dyux-dxuy)*py;
  const double Dy = hy/gamma - vx*dxpy - vy*dypy
                    -nu*(dyuy*py + .5*(dxuy+dyux)*px) - .5*(dxuy-dyux)*px;

  // forcing term
  const double Fx = dxSxx + dySxy - xi*vx + alpha*px;
  const double Fy = dxSyx + dySyy - xi*vy + alpha*py;

  // calculate the equilibrium distribution fe
  const auto fe = GetEquilibriumDistribution(vx, vy, nn);

  if(first)
  {
    PNx[k] = Px[k] + .5*Dx;
    PNy[k] = Py[k] + .5*Dy;

    Px[k] = Px[k] + Dx;
    Py[k] = Py[k] + Dy;

    for(unsigned v=0; v<lbq; ++v)
    {
      fn[k][v] = ff[k][v] + .5*((fe[v]-ff[k][v])/tau
          + w[v]*(Fx*xdir(v) + Fy*ydir(v)));
      ff[k][v] = ff[k][v] +    ((fe[v]-ff[k][v])/tau
          + w[v]*(Fx*xdir(v) + Fy*ydir(v)));
    }
  }
  else
  {
    Px[k] = PNx[k] + .5*Dx;
    Py[k] = PNy[k] + .5*Dy;

    for(unsigned v=0; v<lbq; ++v)
      ff[k][v] = fn[k][v] + .5*((fe[v]-ff[k][v])/tau
          + w[v]*(Fx*xdir(v) + Fy*ydir(v)));
  }
}

void Polar::UpdateQuantities()
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    UpdateQuantitiesAtNode(k);
}

void Polar::UpdateFields(bool first)
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    UpdateFieldsAtNode(k, first);
}

void Polar::Move()
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

void Polar::BoundaryConditionsLB()
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

void Polar::BoundaryConditionsPol()
{
  switch(BC)
  {
    // pbc without bdry layer (nothing to do)
    case 0:
      break;
    // free-slip channel
    case 1:
    {
      const auto apply_bc = [](ScalarField& field) {
        // pbc on the left and right walls
        field.ApplyPBC(PBCWall::LeftRight);
        field.ApplyNeumann(Wall::Front);
        field.ApplyNeumann(Wall::Back);
        // corners
        field.ApplyNeumann(Corner::RightBack,  Wall::Back);
        field.ApplyNeumann(Corner::RightFront, Wall::Front);
        field.ApplyNeumann(Corner::LeftBack,   Wall::Back);
        field.ApplyNeumann(Corner::LeftFront,  Wall::Front);
        /*
        // Dirichlet on the front and back walls
        field.ApplyDirichlet(Wall::Front, value);
        field.ApplyDirichlet(Wall::Back, value);
        // corners
        field.ApplyPBC(Corner::RightBack,  Wall::Right);
        field.ApplyPBC(Corner::RightFront, Wall::Right);
        field.ApplyPBC(Corner::LeftBack,   Wall::Left);
        field.ApplyPBC(Corner::LeftFront,  Wall::Left);
        */
      };

      //const double theta_anch = Pi/2.;
      apply_bc(Px/*, cos(theta_anch)*/);
      apply_bc(Py/*, sin(theta_anch)*/);
      break;
    }
    // pbc
    default:
      Px.ApplyPBC();
      Py.ApplyPBC();
  }
}

void Polar::BoundaryConditionsFields()
{
  switch(BC)
  {
    // pbc without bdry layer (nothing to do)
    case 0:
      break;
    // free-slip channel
    case 1:
    {
      const auto apply_bc = [](ScalarField& field) {
        // pbc on the left and right walls
        field.ApplyPBC(PBCWall::LeftRight);
        // Neumann on the front and back walls
        field.ApplyNeumann(Wall::Front);
        field.ApplyNeumann(Wall::Back);
        // corners
        field.ApplyNeumann(Corner::RightBack,  Wall::Back);
        field.ApplyNeumann(Corner::RightFront, Wall::Front);
        field.ApplyNeumann(Corner::LeftBack,   Wall::Back);
        field.ApplyNeumann(Corner::LeftFront,  Wall::Front);
      };

      apply_bc(ux);
      apply_bc(uy);
      apply_bc(SigmaXX);
      apply_bc(SigmaYY);
      apply_bc(SigmaXY);
      break;
    }
    // pbc
    default:
      ux.ApplyPBC();
      uy.ApplyPBC();
      SigmaXX.ApplyPBC();
      SigmaYY.ApplyPBC();
      SigmaXY.ApplyPBC();
  }
}

void Polar::Step()
{
  // boundary conditions for the polarisation
  BoundaryConditionsPol();
  // predictor step
  UpdateQuantities();
  // boundary conditions for the other fields
  BoundaryConditionsFields();
  // update all the fields
  this->UpdateFields(true);
  // boundary conditions for the flow
  BoundaryConditionsLB();
  // move LB particles
  Move();
  // corrector steps
  for(unsigned n=1; n<=npc; ++n)
  {
    // same thing (no advection)
    BoundaryConditionsPol();
    UpdateQuantities();
    BoundaryConditionsFields();
    this->UpdateFields(false);
  }
}

void Polar::RuntimeStats()
{
  // compute sum of LB particles
  fcheck = 0;
  for(unsigned k=0; k<DomainSize; ++k)
    fcheck = accumulate(begin(ff[k]), end(ff[k]), fcheck);

  std::cout << "Total LB particles: " << fcheck << '\n';
}

void Polar::RuntimeChecks()
{
  // check that the sum of f is constant
  if(abs(ftot-fcheck)>1)
    throw error_msg("f is not conserved (", ftot, "/", fcheck, ")");
}

option_list Polar::GetOptions()
{
  // model specific options
  opt::options_description model_options("Model options");
  model_options.add_options()
    ("gamma", opt::value<double>(&gamma),
     "polarisation mobility")
    ("nu", opt::value<double>(&nu),
     "tumbling/aligning parameter")
    ("tau", opt::value<double>(&tau),
     "fluid viscosity")
    ("rho", opt::value<double>(&rho),
     "fluid density")
    ("xi", opt::value<double>(&xi),
     "substrate friction")
    ("K", opt::value<double>(&K),
     "polarisation elastic constant")
    ("J", opt::value<double>(&J),
     "polarisation nematic elastic constant")
    ("A", opt::value<double>(&A),
     "Free energy strength")
    ("zeta", opt::value<double>(&zeta),
     "activity parameter")
    ("alpha", opt::value<double>(&alpha),
     "motility parameter")
    ("npc", opt::value<unsigned>(&npc),
     "number of correction steps for the predictor-corrector method");

  // init config options
  opt::options_description config_options("Initial configuration options");
  config_options.add_options()
    ("angle", opt::value<double>(&angle_deg),
     "initial angle to x direction (in degrees)")
    ("noise", opt::value<double>(&noise),
     "size of initial variations")
    ("initial-order", opt::value<double>(&init_order),
     "initial order of the polarisation field");

  return { model_options, config_options };
}
