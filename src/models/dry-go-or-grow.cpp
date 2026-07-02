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

  const double dxux = derivX(ux, d, sB);
  const double dyux = derivY(ux, d, sB);
  const double dxuy = derivX(uy, d, sB);
  const double dyuy = derivY(uy, d, sB);

  const double expansion = dxux + dyuy;
  const double shear = .5*(dxuy + dyux);
  const double vorticity = .5*(dxuy - dyux);
  const double traceQL = Qxx*(dxux - dyuy) + 2*Qyx*shear;

  const double Dxx = GammaQ*Hxx - vx*dxQxx - vy*dyQxx - 2*vorticity*Qyx
    + xi*((Qxx+1)*(2*dxux-traceQL) + 2*Qyx*shear - expansion);
  const double Dyx = GammaQ*Hyx - vx*dxQyx - vy*dyQyx + 2*vorticity*Qxx
    + xi*(Qyx*(expansion-traceQL) + 2*shear);
  const double dphi_px_raw = GammaP*(MU[d[1]] - MU[k]) - .5*ux_phi[d[1]];
  const double dphi_mx_raw = GammaP*(MU[d[2]] - MU[k]) + .5*ux_phi[d[2]];
  const double dphi_py_raw = GammaP*(MU[d[3]] - MU[k]) - .5*uy_phi[d[3]];
  const double dphi_my_raw = GammaP*(MU[d[4]] - MU[k]) + .5*uy_phi[d[4]];
  const double dphi_px = MaskedPhiFaceIncrement(k, d[1], dphi_px_raw);
  const double dphi_mx = MaskedPhiFaceIncrement(k, d[2], dphi_mx_raw);
  const double dphi_py = MaskedPhiFaceIncrement(k, d[3], dphi_py_raw);
  const double dphi_my = MaskedPhiFaceIncrement(k, d[4], dphi_my_raw);
  const double Dp = dphi_px + dphi_mx + dphi_py + dphi_my
    - (conserve_phi ? (countphi-totalphi)/DomainSize : 0);

  const double division_m = division_mask[k] ? 1. : 0.;
  const double death_m = death_mask[k] ? 1. : 0.;
  const double growth_rate = alpha*division_m - beta*death_m;
  const double chi_eff = MaterialChi(k, k);
  const double R = chi_eff*p*growth_rate;
  const double Dphi = Dp + R;

  const double mTransport =
    TransportFaceChi(k, d[1], dphi_px)*dphi_px
  + TransportFaceChi(k, d[2], dphi_mx)*dphi_mx
  + TransportFaceChi(k, d[3], dphi_py)*dphi_py
  + TransportFaceChi(k, d[4], dphi_my)*dphi_my
  - chi_eff*(conserve_phi ? (countphi-totalphi)/DomainSize : 0.);

  const double chi_px = MaterialChi(d[1], k);
  const double chi_mx = MaterialChi(d[2], k);
  const double chi_py = MaterialChi(d[3], k);
  const double chi_my = MaterialChi(d[4], k);
  const double phenotypeDiffusion =
    .5*(phi[d[1]] + p)*(chi_px - chi_eff)
  - .5*(p + phi[d[2]])*(chi_eff - chi_mx)
  + .5*(phi[d[3]] + p)*(chi_py - chi_eff)
  - .5*(p + phi[d[4]])*(chi_eff - chi_my);
  const double press = -sigma_bulk[k];
  // Dead-band phenotype switching: convert go<->grow only when the pressure lies
  // more than switch_cost from pswitch. The rate grows with that excess (further
  // from pswitch -> faster) and with the available reservoir (mass-action: grow->go
  // slows as chi->0, go->grow slows as chi->1). Inside the dead band Sswitch=0 so
  // chi is carried by advection alone -- hysteresis without a bistable attractor.
  const double excess_low  = press < pswitch - switch_cost ? (pswitch - switch_cost) - press : 0.;
  const double excess_high = press > pswitch + switch_cost ? press - (pswitch + switch_cost) : 0.;
  const double Sswitch = p*( switch_gogrow*(1-chi_eff)*excess_low
                           - switch_growgo*chi_eff*excess_high );
  const double mGrowth = growTogether ? chi_eff*R : R;
  const double Dm = mTransport + Dchi*phenotypeDiffusion + Sswitch + mGrowth;

  if(first)
  {
    QNxx[k] = QQxx[k] + .5*Dxx;
    QNyx[k] = QQyx[k] + .5*Dyx;
    phn[k] = phi[k] + .5*Dphi;

    QQxx[k] = QQxx[k] + Dxx;
    QQyx[k] = QQyx[k] + Dyx;
    phi_tmp[k] = phi[k] + Dphi;

    mN[k] = m[k] + .5*Dm;
    m_tmp[k] = m[k] + Dm;
  }
  else
  {
    QQxx[k] = QNxx[k] + .5*Dxx;
    QQyx[k] = QNyx[k] + .5*Dyx;
    phi_tmp[k] = phn[k] + .5*Dphi;
    m_tmp[k] = mN[k] + .5*Dm;
  }

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
