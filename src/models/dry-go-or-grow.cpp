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
  // Advective face flux uses the shared (pair-symmetric) face momentum .5*(u_phi[k]+u_phi[n])
  // so each per-face increment is antisymmetric (dphi_px@k = -dphi_mx@n). The extra .5*u_phi[k]
  // terms cancel within Dp, so phi evolution is unchanged to round-off; but weighting these now
  // antisymmetric increments by the upwind chi (phichiTransport below) conserves
  // Sum(phichi) instead of leaking.
  const double dphi_px_raw = GammaP*(MU[d[1]] - MU[k]) - .5*(ux_phi[k] + ux_phi[d[1]]);
  const double dphi_mx_raw = GammaP*(MU[d[2]] - MU[k]) + .5*(ux_phi[k] + ux_phi[d[2]]);
  const double dphi_py_raw = GammaP*(MU[d[3]] - MU[k]) - .5*(uy_phi[k] + uy_phi[d[3]]);
  const double dphi_my_raw = GammaP*(MU[d[4]] - MU[k]) + .5*(uy_phi[k] + uy_phi[d[4]]);
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

  const double phichiTransport =
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
  // Mechanical memory: phim = phi*m is carried by the SAME face increments as phichi
  // (i.e. by the full phi flux J = phi*v - GammaP*grad(mu), not by phi*v alone), so
  //   d_t(phi*m) + div(m*J) = phi*(g(P)-m)/tau_m + m*R,
  // which is the conservative form of tau_m * D_t m = g(P) - m along the biomass.
  // The m*R term makes new biomass inherit the memory of the material it grew from.
  const double m_eff = use_memory ? MaterialM(k, k) : 0.;
  double Dphim = 0.;
  if(use_memory)
  {
    const double phimTransport =
      TransportFaceM(k, d[1], dphi_px)*dphi_px
    + TransportFaceM(k, d[2], dphi_mx)*dphi_mx
    + TransportFaceM(k, d[3], dphi_py)*dphi_py
    + TransportFaceM(k, d[4], dphi_my)*dphi_my
    - m_eff*(conserve_phi ? (countphi-totalphi)/DomainSize : 0.);
    Dphim = phimTransport + MemorySource(k, m_eff) + m_eff*R;
  }

  // Phenotype switching is driven either by the instantaneous local pressure or by
  // the accumulated mechanical memory: V(chi,m) = Ochi*(m-mc)*chi.
  const double press = -sigma_bulk[k];
  const double dVdChi =
    2*Achi*chi_eff*(1-chi_eff)*(1-2*chi_eff)
  + Ochi*(use_memory ? (m_eff-mc) : (press-pswitch));
  const double Sswitch = -p*dVdChi;
  const double phichiGrowth = growTogether ? chi_eff*R : R;
  const double Dphichi = phichiTransport + Dchi*phenotypeDiffusion + Sswitch + phichiGrowth;

  if(first)
  {
    QNxx[k] = QQxx[k] + .5*Dxx;
    QNyx[k] = QQyx[k] + .5*Dyx;
    phn[k] = phi[k] + .5*Dphi;

    QQxx[k] = QQxx[k] + Dxx;
    QQyx[k] = QQyx[k] + Dyx;
    phi_tmp[k] = phi[k] + Dphi;

    phichiN[k] = phichi[k] + .5*Dphichi;
    phichi_tmp[k] = phichi[k] + Dphichi;

    phimN[k] = phim[k] + .5*Dphim;
    phim_tmp[k] = phim[k] + Dphim;
  }
  else
  {
    QQxx[k] = QNxx[k] + .5*Dxx;
    QQyx[k] = QNyx[k] + .5*Dyx;
    phi_tmp[k] = phn[k] + .5*Dphi;
    phichi_tmp[k] = phichiN[k] + .5*Dphichi;
    phim_tmp[k] = phimN[k] + .5*Dphim;
  }

  // Both carried densities obey 0 <= . <= phi, since chi and m each lie in [0,1].
  const double upper = phi_tmp[k] > 0 ? phi_tmp[k] : 0.;
  if(phichi_tmp[k] < 0)
    phichi_tmp[k] = 0;
  else if(phichi_tmp[k] > upper)
    phichi_tmp[k] = upper;
  if(phim_tmp[k] < 0)
    phim_tmp[k] = 0;
  else if(phim_tmp[k] > upper)
    phim_tmp[k] = upper;
}

void DryGoOrGrow::UpdateFields(bool first)
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    UpdateDryFieldsAtNode(k, first);

  swap(phi.get_data(), phi_tmp.get_data());
  swap(phichi.get_data(), phichi_tmp.get_data());
  swap(phim.get_data(), phim_tmp.get_data());
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
