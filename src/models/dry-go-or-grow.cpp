#include "header.hpp"
#include "models/dry-go-or-grow.hpp"
#include "error_msg.hpp"
#include "lb.hpp"

using namespace std;

// from main.cpp:
extern unsigned nthreads;

DryGoOrGrow::DryGoOrGrow(unsigned LX_, unsigned LY_, unsigned BC)
  : GoOrGrow(LX_, LY_, BC)
{}

void DryGoOrGrow::Initialize()
{
  GoOrGrow::Initialize();

  if(friction <= 0)
    throw error_msg("dry-go-or-grow requires friction > 0.");
}

void DryGoOrGrow::ComputeDryVelocityAtNode(unsigned k)
{
  const auto& d = get_neighbours(k);

  const double dxSxx = derivX(sigmaXX, d, sB);
  const double dySxy = derivY(sigmaXY, d, sB);
  const double dxSyx = derivX(sigmaYX, d, sB);
  const double dySyy = derivY(sigmaYY, d, sB);

  const double vx = (dxSxx + dySxy)/friction;
  const double vy = (dxSyx + dySyy)/friction;
  const double p = phi[k];
  const auto fe = GetEquilibriumDistribution(vx, vy, rho);

  ux[k] = vx;
  uy[k] = vy;
  ux_phi[k] = vx*p;
  uy_phi[k] = vy*p;
  n[k] = rho;
  ff[k] = fe;
  fn[k] = fe;
  ff_tmp[k] = fe;
  fn_tmp[k] = fe;
}

void DryGoOrGrow::ComputeDryVelocity()
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    ComputeDryVelocityAtNode(k);
}

void DryGoOrGrow::UpdateDryFieldsAtNode(unsigned k, bool first)
{
  const auto& d = get_neighbours(k);

  const double vx = ux[k];
  const double vy = uy[k];
  const double Qxx = QQxx[k];
  const double Qyx = QQyx[k];
  const double p = phi[k];
  const double Hxx = HHxx[k];
  const double Hyx = HHyx[k];
  const double dxQxx = dxQQxx[k];
  const double dyQxx = dyQQxx[k];
  const double dxQyx = dxQQyx[k];
  const double dyQyx = dyQQyx[k];

  const double del2mu = laplacian(MU, d, sD);
  const double dxux = derivX(ux, d, sB);
  const double dyux = derivY(ux, d, sB);
  const double dxuy = derivX(uy, d, sB);
  const double dyuy = derivY(uy, d, sB);
  const double pFlux = derivX(ux_phi, d, sB) + derivY(uy_phi, d, sB);

  const double expansion = dxux + dyuy;
  const double shear = .5*(dxuy + dyux);
  const double vorticity = .5*(dxuy - dyux);
  const double traceQL = Qxx*(dxux - dyuy) + 2*Qyx*shear;

  const double Dxx = GammaQ*Hxx - vx*dxQxx - vy*dyQxx - 2*vorticity*Qyx
    + xi*((Qxx+1)*(2*dxux-traceQL) + 2*Qyx*shear - expansion);
  const double Dyx = GammaQ*Hyx - vx*dxQyx - vy*dyQyx + 2*vorticity*Qxx
    + xi*(Qyx*(expansion-traceQL) + 2*shear);
  const double Dp = GammaP*del2mu - pFlux
    - (conserve_phi ? (countphi-totalphi)/DomainSize : 0);

  const double division_m = division_mask[k] ? 1. : 0.;
  const double death_m = death_mask[k] ? 1. : 0.;
  const double growth_rate = alpha*division_m - beta*death_m;
  const double chi_eff = chi[k];
  const double R = chi_eff*p*growth_rate;

  const double mFlux =
    .5*(ux[d[1]]*m[d[1]] - ux[d[2]]*m[d[2]])
  + .5*(uy[d[3]]*m[d[3]] - uy[d[4]]*m[d[4]]);

  const double chi_px = chi[d[1]];
  const double chi_mx = chi[d[2]];
  const double chi_py = chi[d[3]];
  const double chi_my = chi[d[4]];
  const double diffusiveMFlux =
    .5*(chi_px + chi_eff)*(MU[d[1]] - MU[k])
  - .5*(chi_eff + chi_mx)*(MU[k] - MU[d[2]])
  + .5*(chi_py + chi_eff)*(MU[d[3]] - MU[k])
  - .5*(chi_eff + chi_my)*(MU[k] - MU[d[4]]);
  const double phenotypeDiffusion =
    .5*(phi[d[1]] + p)*(chi_px - chi_eff)
  - .5*(p + phi[d[2]])*(chi_eff - chi_mx)
  + .5*(phi[d[3]] + p)*(chi_py - chi_eff)
  - .5*(p + phi[d[4]])*(chi_eff - chi_my);
  const double dVdChi =
    2*Achi*chi_eff*(1-chi_eff)*(1-2*chi_eff)
  + Ochi*(p-phiswitch);
  const double Sswitch = -p*dVdChi;
  const double mGrowth = growTogether ? chi_eff*R : R;
  const double Dm = GammaP*diffusiveMFlux + Dchi*phenotypeDiffusion
    - mFlux + Sswitch + mGrowth;

  if(first)
  {
    QNxx[k] = QQxx[k] + .5*Dxx;
    QNyx[k] = QQyx[k] + .5*Dyx;
    phn[k] = phi[k] + .5*Dp;

    QQxx[k] = QQxx[k] + Dxx;
    QQyx[k] = QQyx[k] + Dyx;
    phi_tmp[k] = phi[k] + Dp;

    mN[k] = m[k] + .5*Dm;
    m_tmp[k] = m[k] + Dm;
  }
  else
  {
    QQxx[k] = QNxx[k] + .5*Dxx;
    QQyx[k] = QNyx[k] + .5*Dyx;
    phi_tmp[k] = phn[k] + .5*Dp;
    m_tmp[k] = mN[k] + .5*Dm;
  }

  phi_tmp[k] += R;

  const double upper = phi_tmp[k] > 0 ? phi_tmp[k] : 0.;
  if(m_tmp[k] < 0)
    m_tmp[k] = 0;
  else if(m_tmp[k] > upper)
    m_tmp[k] = upper;
}

void DryGoOrGrow::UpdateFields(bool first)
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    UpdateDryFieldsAtNode(k, first);

  swap(phi.get_data(), phi_tmp.get_data());
  swap(m.get_data(), m_tmp.get_data());
  UpdatePhenotypeQuantities();
}

void DryGoOrGrow::Step()
{
  SetMasks();

  BoundaryConditionsFields();
  UpdateQuantities();
  Lyotropic::BoundaryConditionsFields2();
  ComputeDryVelocity();
  Lyotropic::BoundaryConditionsFields2();
  this->UpdateFields(true);

  for(unsigned n=1; n<=npc; ++n)
  {
    BoundaryConditionsFields();
    UpdateQuantities();
    Lyotropic::BoundaryConditionsFields2();
    ComputeDryVelocity();
    Lyotropic::BoundaryConditionsFields2();
    this->UpdateFields(false);
  }
}
