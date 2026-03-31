#include "header.hpp"
#include "models/phases.hpp"
#include "error_msg.hpp"
#include "random.hpp"
#include "lb.hpp"
#include "tools.hpp"

using namespace std;
namespace opt = boost::program_options;

// defined in main.cpp:
extern double time_step;
extern unsigned nsubsteps;
extern bool no_warning, stop_at_warning;
extern unsigned verbose;

Phases::Phases(unsigned LX, unsigned LY, unsigned BC)
  : Model(LX, LY, BC, GridType::Periodic)
{}

void Phases::InitializeFields()
{
  phi.resize(nphases, {LX, LY, GridType::Periodic});
  phi_old.resize(nphases, {LX, LY, GridType::Periodic});
  V.resize(nphases, {LX, LY, GridType::Periodic});
  potential.resize(nphases, {LX, LY, GridType::Periodic});
  potential_old.resize(nphases, {LX, LY, GridType::Periodic});
  area.resize(nphases, 0.);
  area_cnt.resize(nphases, 0.);
  domain_min.resize(nphases, {0, 0});
  domain_max.resize(nphases, {LX, LY});
  com.resize(nphases, {0., 0.});
  com_prev.resize(nphases, {0., 0.});
  pol.resize(nphases, {0., 0.});
  velp.resize(nphases, {0., 0.});
  velc.resize(nphases, {0., 0.});
  velf.resize(nphases, {0., 0.});
  vel.resize(nphases, {0., 0.});
  com_x.resize(nphases, 0.);
  com_y.resize(nphases, 0.);
  c.resize(nphases, 0);
  S00.resize(nphases, 0.);
  S01.resize(nphases, 0.);
  S_order.resize(nphases, 0.);
  S_angle.resize(nphases, 0.);
  theta.resize(nphases, 0.);
  dtheta.resize(nphases, 0.);
  //torque.resize(nphases, 0.);
}

void Phases::Initialize()
{
  // initialize memory
  walls.SetSize(LX, LY, GridType::Periodic);
  walls_dx.SetSize(LX, LY, GridType::Periodic);
  walls_dy.SetSize(LX, LY, GridType::Periodic);
  walls_laplace.SetSize(LX, LY, GridType::Periodic);
  sum.SetSize(LX, LY, GridType::Periodic);
  sum_cnt.SetSize(LX, LY, GridType::Periodic);
  square.SetSize(LX, LY, GridType::Periodic);
  square_cnt.SetSize(LX, LY, GridType::Periodic);
  Px.SetSize(LX, LY, GridType::Periodic);
  Py.SetSize(LX, LY, GridType::Periodic);
  Theta.SetSize(LX, LY, GridType::Periodic);
  Q00.SetSize(LX, LY, GridType::Periodic);
  Q01.SetSize(LX, LY, GridType::Periodic);
  Px_cnt.SetSize(LX, LY, GridType::Periodic);
  Py_cnt.SetSize(LX, LY, GridType::Periodic);
  Theta_cnt.resize(DomainSize, 0.);
  Q00_cnt.resize(DomainSize, 0.);
  Q01_cnt.resize(DomainSize, 0.);
  P.resize(DomainSize, 0.);
  P_cnt.resize(DomainSize, 0.);
  InitializeFields();

  // pre-compute coefficients
  C1 = 60./lambda/lambda;
  C2 = Pi*R*R;
  C3 = C1/lambda/lambda;

  // extend the parameters with the last given value
  gamma.resize(nphases, gamma.back());
  mu.resize(nphases, mu.back());

  // compute tables
  for(unsigned i=0; i<LX; ++i)
    com_x_table.push_back(exp(complex<double>(0., 1.)*(-Pi+2.*Pi*i/LX)));
  for(unsigned i=0; i<LY; ++i)
    com_y_table.push_back(exp(complex<double>(0., 1.)*(-Pi+2.*Pi*i/LY)));

  // check birth boundaries
  if(birth_bdries.size()==0)
    birth_bdries = {0, LX, 0, LY};
  else if(birth_bdries.size()!=4)
    throw error_msg("Birth boundaries have wrong format, see help.");

  // set flags
  division = (division_rate!=0.);
}

void Phases::AddCell(unsigned n, const array<unsigned, 2>& c)
{
  // we create smaller cells that will then relax
  // this improves greatly the stability at the first steps
  const auto radius = max(R/2., 4.);

  // create the cells at the centers we just computed
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
  {
    const unsigned xk = GetXPosition(k);
    const unsigned yk = GetYPosition(k);

    // round shape (do not wrap if no PBC)
    if(
        (BC==0 and pow(wrap(diff(yk, c[1]), LY), 2)
         + pow(wrap(diff(xk, c[0]), LX), 2)<=ceil(radius*radius))
        or
        (BC>=1 and pow(diff(yk, c[1]), 2)
         + pow(diff(xk, c[0]), 2)<=ceil(radius*radius))
      )
    {
      phi[n][k]     = 1.;
      phi_old[n][k] = 1.;
      area[n]      += 1.;
      square[k]    += 1.;
      sum[k]       += 1.;
    }
    else
    {
      phi[n][k]     = 0.;
      phi_old[n][k] = 0.;
    }
  }

  this->c[n] = c0;
  theta[n]   = 2.*Pi*random_real();
  dtheta[n]  = 0;

  if(tracking)
  {
    domain_min[n][0] = (c[0]+LX-margin)%LX;
    domain_max[n][0] = (c[0]+margin)%LX;
    domain_min[n][1] = (c[1]+LY-margin)%LY;
    domain_max[n][1] = (c[1]+margin)%LY;
  }
  else
  {
    domain_min[n][0] = 0u;
    domain_max[n][0] = LX;
    domain_min[n][1] = 0u;
    domain_max[n][1] = LY;
  }
}

void Phases::Configure()
{
  // ===========================================================================
  // adding cells at random while trying to keep their center non overlapping
  if(init_config=="random" and BC==0)
  {
    // target radius for spacing between cells
    unsigned radius = sqrt(double(LX*LY/nphases)/Pi);
    // list of all centers
    vector<array<unsigned, 2>> centers;

    for(unsigned n=0; n<nphases; ++n)
    {
      // generate new center while trying to keep safe distance
      while(true)
      {
        const array<unsigned, 2> center = {
          static_cast<unsigned>(birth_bdries[0]
              +random_real()*(birth_bdries[1]-birth_bdries[0])),
          static_cast<unsigned>(birth_bdries[2]
              +random_real()*(birth_bdries[3]-birth_bdries[2]))
        };

        bool is_overlapping = false;
        for(const auto& c : centers)
        {
          if(pow(wrap(diff(center[0], c[0]), LX), 2)
              + pow(wrap(diff(center[1], c[1]), LY), 2) < 0.9*radius*radius)
          {
            is_overlapping = true;
            break;
          }
        }

        if(!is_overlapping)
        {
          centers.emplace_back(center);
          break;
        }
      }

      // add cell
      AddCell(n, centers.back());
    }
  }
  // ===========================================================================
  // same but with walls: we need to be careful not to create cells on the wall
  else if(init_config=="random" and BC>=1)
  {
    // target radius for spacing between cells
    unsigned radius = sqrt(double(LX*LY/nphases)/Pi);
    // list of all centers
    vector<array<unsigned, 2>> centers;

    for(unsigned n=0; n<nphases; ++n)
    {
      // generate new center while trying to keep safe distance
      while(true)
      {
        const array<unsigned, 2> center = {
          static_cast<unsigned>(birth_bdries[0]
              +random_real()*(birth_bdries[1]-birth_bdries[0])),
          static_cast<unsigned>(birth_bdries[2]
              +random_real()*(birth_bdries[3]-birth_bdries[2]))
        };

        // detect walls
        // ... box only
        if(BC==1)
          if(center[0]<0.9*R or LX-center[0]<0.9*R) continue;
        // ... box and channel
        if(BC==1 or BC==2)
          if(center[1]<0.9*R or LY-center[1]<0.9*R) continue;
        // ... ellipse
        if(BC==3)
        {
          // compute distance from the elliptic wall
          // ... angle of the current point (from center of the domain)
          const auto angle = atan2(LY/2.-center[1], LX/2.-center[0]);
          // ... small helper function to compute radius
          auto rad = [](double x, double y) { return sqrt(x*x + y*y); };
          // ... distance is the difference between wall and current point
          const auto d = rad(LX/2.*cos(angle), LY/2.*sin(angle))
                        -rad(LX/2.-center[0], LY/2.-center[1]);

          if(d<0.9*R) continue;
        }

        // overlapp between cells
        bool is_overlapping = false;
        for(const auto& c : centers)
        {
          if(pow(wrap(diff(center[0], c[0]), LX), 2)
              + pow(wrap(diff(center[1], c[1]), LY), 2) < 0.9*radius*radius)
          {
            is_overlapping = true;
            break;
          }
        }

        if(!is_overlapping)
        {
          centers.emplace_back(center);
          break;
        }
      }

      // add cell
      AddCell(n, centers.back());
    }
  }
  // ===========================================================================
  // cluster of close cells in the center
  else if(init_config=="cluster")
  {
    const double theta  = 2*Pi/nphases;
    const double radius = R + nphases - 2;

    for(unsigned n=0; n<nphases; ++n)
      AddCell(n, {unsigned(LX/2+radius*(cos(n*theta)+noise*random_real())),
                   unsigned(LY/2+radius*(sin(n*theta)+noise*random_real())) });
  }
  // ===========================================================================
  // single cell in the middle
  else if(init_config=="single")
  {
    if(nphases!=1)
      throw error_msg("error: initial conditions require "
                      "nphases=1.");

    AddCell(0, {LX/2, LY/2});
  }
  else throw error_msg("error: initial configuration '",
      init_config, "' unknown.");

  // ===========================================================================
  // Initialize the walls
  ConfigureWalls();
}

void Phases::ConfigureWalls()
{
  switch(BC)
  {
  case 0:
    // no walls (pbc)
    for(unsigned k=0; k<DomainSize; ++k)
      walls[k] = 0;
    break;
  case 1:
    // Exponentially falling phase-field:
    for(unsigned k=0; k<DomainSize; ++k)
    {
      const double x = GetXPosition(k);
      const double y = GetYPosition(k);

      // this is the easies way: each wall contributes as an exponentially
      // falling potential and we do not care about overalps
      walls[k] =   exp(-y/wall_thickness)
                 + exp(-x/wall_thickness)
                 + exp(-(LX-1-x)/wall_thickness)
                 + exp(-(LY-1-y)/wall_thickness);
    }
    break;
  // Same as above but channel.
  case 2:
    for(unsigned k=0; k<DomainSize; ++k)
    {
      const auto y = GetYPosition(k);

      // exponentially falling on both sides
      walls[k] = exp(-double(y)/wall_thickness)
        + exp(-double(LY-y-1)/wall_thickness);
    }
    break;
  // ellipse!
  case 3:
    for(unsigned k=0; k<DomainSize; ++k)
    {
      const auto x = GetXPosition(k);
      const auto y = GetYPosition(k);

      // compute distance from the elliptic wall
      // ... angle of the current point (from center of the domain)
      const auto theta = atan2(LY/2.-y, LX/2.-x);
      // ... small helper function to compute radius
      const auto rad = [](double x, double y) { return sqrt(x*x + y*y); };
      // ... distance is the difference between wall and current point
      const auto d = rad(LX/2.*cos(theta), LY/2.*sin(theta))
                    -rad(LX/2.-x, LY/2.-y);
      // set the wall
      if(d<0)
        walls[k] = 1.;
      else
        walls[k] = exp(-d/wall_thickness);
    }
    break;
  default:
    throw error_msg("boundary condition unknown.");
  }

  // pre-compute derivatives
  for(unsigned k=0; k<DomainSize; ++k)
  {
    const auto& d  = get_neighbours(k);

    walls_dx[k] = derivX(walls, d, sB);
    walls_dy[k] = derivY(walls, d, sB);
    walls_laplace[k] = laplacian(walls, d, sD);
  }
}

void Phases::Pre()
{
  // we make the system relax (without activity)
  if(relax_time>0)
  {
    double save_c0 = 0.; swap(c0, save_c0);
    double save_zeta = 0.; swap(zeta, save_zeta);
    double save_beta = 0.; swap(beta, save_beta);
    double save_alpha = 0.; swap(alpha, save_alpha);
    double save_division_rate = 0.; swap(save_division_rate, division_rate);

    if(relax_nsubsteps) swap(nsubsteps, relax_nsubsteps);

    for(unsigned i=0; i<relax_time*nsubsteps; ++i) Step();

    if(relax_nsubsteps) swap(nsubsteps, relax_nsubsteps);

    swap(c0, save_c0);
    swap(zeta, save_zeta);
    swap(beta, save_beta);
    swap(alpha, save_alpha);
    swap(save_division_rate, division_rate);
  }
}

void Phases::PreRunStats()
{
  // packing fraction
  {
    double packing = nphases*Pi*R*R;
    if(BC<=2)
      packing/= (birth_bdries[1]-birth_bdries[0])
               *(birth_bdries[3]-birth_bdries[2]);
    if(BC==3) packing /= Pi*LX*LY/4.;

    cout << "Packing fraction = " << packing << endl;
  }
}

void Phases::RuntimeStats()
{
  // TBD
}

void Phases::RuntimeChecks()
{
  // check that the area is more or less conserved (15%)
  for(const auto a : area)
    if(abs(1.-a/C2)>.15)
      throw warning_msg("area is not conserved.");

  for(unsigned n=0; n<nphases; ++n)
  {
    // check that the cells are not leaking, i.e. that at least 95% of the
    // contributions to the area comes from inside the cell (>1/2).
    double a = 0.;
    // compute area of points outside the cell (<1/2)
    for(const auto v : phi[n]) if(v<.5) a += v*v;
    // check that it is less than 5%
    if(a/area[n]>.05)
      throw warning_msg("your cells are leaking!");

    // check that the phase fields stay between 0 and 1
    for(const auto& p : phi[n])
      if(p<-0.5 or p>1.5)
        throw warning_msg("phase-field is not in [0,1]!");
  }
}

inline void Phases::UpdateFieldsAtNode(unsigned n, unsigned k)
{
  const auto& d  = get_neighbours(k);
  // cell properties
  const auto& p  = phi[n][k];
  const auto  l  = laplacian(phi[n], d, sD);
  const auto  dx = derivX(phi[n], d, sB);
  const auto  dy = derivY(phi[n], d, sB);
  const auto  a  = area[n];
  // all-cells properties
  const auto  ls = laplacian(sum, d, sD);
  const auto  dxs  = derivX(sum, d, sB);
  const auto  dys  = derivY(sum, d, sB);
  const auto  dxp0 = derivX(Px, d, sB);
  const auto  dyp0 = derivY(Px, d, sB);
  const auto  dxp1 = derivX(Py, d, sB);
  const auto  dyp1 = derivY(Py, d, sB);

  const auto  lw = walls_laplace[k];
  const auto dxw = walls_dx[k];
  const auto dyw = walls_dy[k];

  // delta F / delta phi
  const double force = (
      + C1*(
        // repulsion term
        + kappa*p*(square[k]-p*p)
        + wall_kappa*p*walls[k]*walls[k]
        )
      - C3*(
        // adhesion term
        + omega*(ls-l)
        + wall_omega*lw
        )
      + (
        // elasticity properties
        + C1*gamma[n]*p*(1.-p)*(1.-2.*p)
        - 2.*gamma[n]*l
        - 2.*mu[n]*(1.-a/C2)*2.*p
        )
      );

  // potential
  V[n][k]     = force;
  // passive force
  velp[n][0] += dx*force;
  velp[n][1] += dy*force;
  // contractility force
  velc[n][0] += ( (P[k]+zeta*Q00[k])*dx + zeta*Q01[k]*dy );
  velc[n][1] += ( zeta*Q01[k]*dx + (P[k]-zeta*Q00[k])*dy );
  // friction force
  velf[n][0] += + f*alpha/xi*pol[n][0]*(dx*(dxs-dx)+dy*(dys-dy))
                - f*alpha/xi*(dx*(dxp0-pol[n][0]*dx)+dy*(dyp0-pol[n][0]*dy))
                //+ f_walls*alpha/xi*pol[n][0]*(dx*dxw+dy*dyw);
                - dyw*f_walls/xi*(pol[n][1]*dxw-pol[n][0]*dyw)*(dx*dxw+dy*dyw);
  velf[n][1] += + f*alpha/xi*pol[n][1]*(dx*(dxs-dx)+dy*(dys-dy))
                - f*alpha/xi*(dx*(dxp1-pol[n][1]*dx)+dy*(dyp1-pol[n][1]*dy))
                //+ f_walls*alpha/xi*pol[n][1]*(dx*dxw+dy*dyw);
                + dxw*f_walls/xi*(pol[n][1]*dxw-pol[n][0]*dyw)*(dx*dxw+dy*dyw);

  // these are different alignment torques and must be cleaned up once we decide
  // which one is the best
  //
  // alignment torque (orientational)
  //const auto dxQ00 = derivX(Q00, d, sB);
  //const auto dxQ01 = derivX(Q01, d, sB);
  //const auto dyQ00 = derivY(Q00, d, sB);
  //const auto dyQ01 = derivY(Q01, d, sB);
  //torque[n] += - 4*J1*(Q00[n]*(dx*dxQ01+dy*dyQ01)-Q01[n]*(dx*dxQ00+dy*dyQ00));
  //
  // alignment torque (directional)
  //const auto  dxt  = derivX(Theta, d, sB);
  //const auto  dyt  = derivY(Theta, d, sB);
  //torque[n] -= J1*(dx*dxt+dy*dyt - theta[n]*(dx*dxs+dy*dys));
  //
  // alignment torque (shear stress)
  //torque[n] -= J1*zeta*(2*Q00[n]*dx*dy + Q01[n]*(dx*dx-dy*dy));
}

inline void Phases::UpdateAtNode(unsigned n, unsigned k)
{
  // compute potential
  {
    const auto& d  = get_neighbours(k);
    const auto  dx = derivX(phi[n], d, sB);
    const auto  dz = derivY(phi[n], d, sB);

    potential[n][k] = (
      // free energy term
      -.5*V[n][k]
      // advection term
      - (alpha*pol[n][0]+velp[n][0]/*+velc[n][0]+velf[n][0]*/)*dx/xi
      - (alpha*pol[n][1]+velp[n][1]/*+velc[n][1]+velf[n][1]*/)*dz/xi
      );
  }

  // store values
  if(store)
  {
    potential_old[n][k] = potential[n][k];
    phi_old[n][k]       = phi[n][k];
  }

  // predictor-corrector
  {
    const double p = phi_old[n][k]
                     + time_step*.5*(potential[n][k] + potential_old[n][k]);

    // update for next call
    phi[n][k]    = p;
    com_x[n]    += com_x_table[GetXPosition(k)]*p;
    com_y[n]    += com_y_table[GetYPosition(k)]*p;
    area_cnt[n] += p*p;
  }
}

void Phases::UpdatePolarization(unsigned n)
{
  // align to the sum of passive and viscous forces
  array<double, 2> v = {
    velp[n][0] + velf[n][0] /*+ velc[n][0]*/,
    velp[n][1] + velf[n][1] /*+ velc[n][1]*/
  };

  // ...the norm of the passive velocity
  const double ni = sqrt(abs(v[0]*v[0]+v[1]*v[1]));

  // alignement torque
  double torque = -ni*atan2(v[0]*pol[n][1]-v[1]*pol[n][0],
                            v[0]*pol[n][0]+v[1]*pol[n][1]);

  // ...euler-marijuana update
  dtheta[n] += time_step*(torque-dtheta[n])/tauc;
  theta[n]  += time_step*J*dtheta[n] + sqrt(time_step)*D*random_normal();

  // update polarisation and contractility
  pol[n] = {cos(theta[n]), sin(theta[n])};

  // dynamics of the contractility: needs some cleaning up once settled
  //
  //c[n]  += -time_step*(c[n]-c0)/tauc + beta*c0*(area_cnt[n]-area[n]);
  //c[n]  += time_step*( -(c[n]-c0) - beta*(1.-area_cnt[n]/C2/.88171))/tauc;
  //c[n]  -= time_step*(c[n]-beta*c0*area_cnt[n]/C2)/tauc;
}

void Phases::ComputeCoM(unsigned n)
{
  // the strategy to deal with the periodic boundary conditions is to compute
  // all the integrals in Fourier space and come back at the end. This way the
  // periodicity of the domain is automatically taken into account.
  const auto mx = arg(com_x[n]/static_cast<double>(DomainSize)) + Pi;
  const auto my = arg(com_y[n]/static_cast<double>(DomainSize)) + Pi;
  com[n] = { mx/2./Pi*LX, my/2./Pi*LY };
}

void Phases::UpdateWindow(unsigned n)
{
  // update walls
  domain_min[n][0] = (static_cast<unsigned>(round(com[n][0])) + LX - margin)%LX;
  domain_min[n][1] = (static_cast<unsigned>(round(com[n][1])) + LY - margin)%LY;
  domain_max[n][0] = (static_cast<unsigned>(round(com[n][0])) + margin)%LX;
  domain_max[n][1] = (static_cast<unsigned>(round(com[n][1])) + margin)%LY;
}

inline void Phases::UpdateStructureTensorAtNode(unsigned n, unsigned k)
{
  const auto& d  = get_neighbours(k);
  const auto  p  = phi[n][k];
  const auto  dx = 2*p*derivX(phi[n], d, sB);
  const auto  dy = 2*p*derivY(phi[n], d, sB);
  S00[n] += 0.5*(dx*dx-dy*dy);
  S01[n] += dx*dy;
}

void Phases::ComputeShape(unsigned n)
{
  // shape: we remap the 2x2 traceless symmetric matrix to polar coord for ease
  // of manipulation
  S_order[n] = sqrt(S00[n]*S00[n]+S01[n]*S01[n]);
  S_angle[n] = atan2(S01[n], S00[n])/2.;
}

inline void Phases::SquareAndSumAtNode(unsigned n, unsigned k)
{
  const auto p = phi[n][k];

  // we swap counters and values afterwards
  sum_cnt[k]    += p;
  square_cnt[k] += p*p;
  P_cnt[k]      -= p*c[n];
  Q00_cnt[k]    -= p*.5*(pol[n][0]*pol[n][0]-pol[n][1]*pol[n][1]);
  Q01_cnt[k]    -= p*pol[n][0]*pol[n][1];
  Px_cnt[k]     += p*pol[n][0];
  Py_cnt[k]     += p*pol[n][1];
  Theta_cnt[k]  += p*theta[n];

  // coupling to shape
  //Q00_cnt[k] += p*p*(-c[n]+zeta*S_order[n]*cos(2*S_angle[n]));
  //Q01_cnt[k] += p*p*(     +zeta*S_order[n]*sin(2*S_angle[n]));
  //Q00_cnt[k] += p*p*(-c[n]-zeta*S_order[n]*cos(2*S_angle[n]));
}

inline void Phases::ReinitSquareAndSumAtNode(unsigned k)
{
  // reinit values
  sum[k]    = 0;
  square[k] = 0;
  P[k]      = 0;
  Theta[k]  = 0;
  Q00[k]    = 0;
  Q01[k]    = 0;
  Px[k]     = 0;
  Py[k]     = 0;
}

void Phases::Divide(unsigned n)
{
  // work in progress

  const double r = random_real();

  if(r<nsubsteps*division_rate)
  //if(not pine--)
  {
    cout << "DIVIDE" << endl;
    // add a new cell and extend memory
    const auto m = nphases++;
    InitializeFields();
    gamma.resize(nphases, gamma[n]);
    mu.resize(nphases, mu[n]);
    theta[m] = theta[n];
    // pick random axis
    //const double theta0 = 2*Pi*r/nsubsteps/division_rate;
    // set division axis to be perp to the polarization
    const double theta0 = theta[n] + Pi;

    #pragma omp parallel for num_threads(nthreads) if(nthreads)
    for(unsigned k=0; k<DomainSize; ++k)
    {
      const auto x = wrap(diff(unsigned(com[n][0]), GetXPosition(k)), LX);
      const auto y = wrap(diff(unsigned(com[n][1]), GetXPosition(k)), LY);
      const auto theta = atan2(y, x);

      // distance from cutting line
      const auto d = sqrt(x*x+y*y)*sin(theta-theta0);
      // current size of the cell (more or less)
      const auto l = 30;//5*sqrt(area[n]/Pi);

      const auto p = phi[n][k];
      phi[n][k] = .5*(1.+tanh( d/l))*p;
      phi[m][k] = .5*(1.-tanh(-d/l))*p;
    }
  }
}

void Phases::Update()
{
  // 1) Compute induced force and passive velocity
  //
  // We need to loop over all nodes once before updating the phase fields
  // because the passive component of the velocity requires a integral
  // involving a derivative. This means that the field phi must be fully
  // updated first.

  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned n=0; n<nphases; ++n)
  {
    velp[n] = {0., 0.};
    velc[n] = {0., 0.};
    velf[n] = {0., 0.};
    //torque[n] = 0.;

    // update in restricted domain only
    UpdateDomain(&Phases::UpdateFieldsAtNode, n);
  }

  // 2) Predictor-corrector function for updating the phase fields
  //
  // The predictor corrector is such that it can be used many time in a row in
  // order to give you better precision, effectively giving higher order
  // approximations.

  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned n=0; n<nphases; ++n)
  {
    // only update fields in the restricted domain of field n
    UpdateDomain(&Phases::UpdateAtNode, n);
    // because the polarisation dynamics is first
    // order (euler-maruyama) we need to update only
    // once per predictor-corrector step
    if(store) UpdatePolarization(n);
    // update center of mass
    ComputeCoM(n);
    // update domain walls
    if(tracking) UpdateWindow(n);
    // update structure tensor
    UpdateDomain(&Phases::UpdateStructureTensorAtNode, n);
    // and get shape
    ComputeShape(n);

    // reinit for next round
    com_x[n] = com_y[n] = area[n] = 0.;
    S00[n] = S01[n] = 0.;
  }

  // 3) Compute square and sum
  //
  // We need yet another loop here for parallelization, because we can not send
  // each phi to a different core when computing the square and the sum of all
  // phase fields. This is much faster than using an atomic portion in the
  // previous loop.

  for(unsigned n=0; n<nphases; ++n)
    // update only domain (in parallel, each node to a different core)
    UpdateDomainP(&Phases::SquareAndSumAtNode, n);

  // 4) Reinit counters and swap to get correct values
  //
  // We use this construction because it is much faster with OpenMP: the single
  // threaded portion of the code consists only of these swaps!

  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k) ReinitSquareAndSumAtNode(k);

  swap(sum.get_data(), sum_cnt.get_data());
  swap(square.get_data(), square_cnt.get_data());
  swap(P, P_cnt);
  swap(Theta.get_data(), Theta_cnt);
  swap(Q00.get_data(), Q00_cnt);
  swap(Q01.get_data(), Q01_cnt);
  swap(Px.get_data(), Px_cnt.get_data());
  swap(Py.get_data(), Py_cnt.get_data());
  swap(area, area_cnt);
}

void Phases::Step()
{
  // first sweeps produces estimate of values
  store = true; // ok this is really bad :-(
  Update();
  store = false;

  // subsequent sweeps produce corrected values
  for(unsigned i=0; i<npc; ++i)
    Update();

  // division
  const auto m = nphases; // this is needed because we might increment nphases
  for(unsigned n=0; n<m; ++n) Divide(n);

  // compute center-of-mass velocity
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned n=0; n<nphases; ++n)
  {
    vel[n]      = { (com[n][0]-com_prev[n][0])/time_step,
                    (com[n][1]-com_prev[n][1])/time_step };
    com_prev[n] = com[n];
  }
}

option_list Phases::GetOptions()
{
  // model specific options
  opt::options_description model_options("Model options");
  model_options.add_options()
    ("nphases", opt::value<unsigned>(&nphases),
      "Number of phases")
    ("gamma", opt::value<vector<double>>(&gamma),
      "Elastic constant of each phase (array)")
    ("mu", opt::value<vector<double>>(&mu),
      "Energy penalty for area of each phase (array)")
    ("lambda", opt::value<double>(&lambda),
      "Interface thickness parameter")
    ("kappa", opt::value<double>(&kappa),
      "Interaction strength")
    ("alpha", opt::value<double>(&alpha),
      "Migration alpha (passive)")
    ("beta", opt::value<double>(&beta),
      "..........")
    ("npc", opt::value<unsigned>(&npc)->default_value(1u),
      "Number of predictor-corrector steps")
    ("margin", opt::value<unsigned>(&margin)->default_value(0u),
      "Margin for the definition of restricted domains (if 0: update full box)")
    ("friction", opt::value<double>(&f),
      "Cell-cell friction parameter")
    ("friction-walls", opt::value<double>(&f_walls),
      "Cell-wall friction parameter")
    ("xi", opt::value<double>(&xi),
      "Substrate friction parameter")
    ("c0", opt::value<double>(&c0),
     "Base level for cellular contractility")
    ("tauc", opt::value<double>(&tauc),
     "Relaxation time for the contractility")
    ("zeta", opt::value<double>(&zeta),
     "Friction parameter")
    ("omega", opt::value<double>(&omega),
      "Adhesion parameter")
    ("J", opt::value<double>(&J),
      "Strength of alignment torque")
    ("D", opt::value<double>(&D),
      "Rotational diffusion constant")
    ("tracking", opt::bool_switch(&tracking),
     "enable tracking of cells in subregions only")
    ("wall-thickness", opt::value<double>(&wall_thickness),
      "Wall thickness (typical decay length)")
    ("wall-kappa", opt::value<double>(&wall_kappa)->default_value(kappa),
      "Wall repulsion")
    ("wall-omega", opt::value<double>(&wall_omega)->default_value(0.),
      "Wall adhesion")
    ("division-rate",  opt::value<double>(&division_rate),
      "Rate of division")
    ("R", opt::value<double>(&R),
      "Preferred radius (defines area Pi*R*R)");

  // init config options
  opt::options_description config_options("Initial configuration options");
  config_options.add_options()
    ("config", opt::value<string>(&init_config),
      "Initial configuration")
    ("relax-time", opt::value<unsigned>(&relax_time)->default_value(0u),
      "Relaxation time steps at initialization.")
    ("noise", opt::value<double>(&noise),
      "Noise level")
    ("birth-boundaries", opt::value<vector<unsigned>>(&birth_bdries)->multitoken(),
     "Boundaries in which the cells are created "
     "when the initial configuration 'random' is choosed. "
     "Format: {min x, max, x, min y, max y}")
    ("relax-nsubsteps", opt::value<unsigned>(&relax_nsubsteps)->default_value(0u),
      "Value of nsubsteps to use at initial relaxation (0 means use nsubsteps).");

  return { model_options, config_options };
}
