#include "header.hpp"
#include "models/go-or-grow.hpp"
#include "error_msg.hpp"
#include "lb.hpp"
#include "random.hpp"

using namespace std;
namespace opt = boost::program_options;

// from main.cpp:
extern unsigned nthreads;

GoOrGrow::GoOrGrow(unsigned LX_, unsigned LY_, unsigned BC)
  : LyotropicWithDivision(LX_, LY_, BC)
{}

void GoOrGrow::Initialize()
{
  LyotropicWithDivision::Initialize();

  if(chi_noise < 0)
    throw error_msg("chi-noise must be non-negative.");
  if(chi_length < 0)
    throw error_msg("chi-length must be non-negative.");

  chi.SetSize(LX, LY, Type);
  m.SetSize(LX, LY, Type);
  mN.SetSize(LX, LY, Type);
  m_tmp.SetSize(LX, LY, Type);
}

void GoOrGrow::Configure()
{
  LyotropicWithDivision::Configure();
  ConfigurePhenotype();
}

void GoOrGrow::ConfigurePhenotype()
{
  if(chi_config=="noise")
  {
    for(unsigned k=0; k<DomainSize; ++k)
      chi[k] = random_normal();

    BoundaryConditionsFields();

    constexpr double smooth_rate = 0.2;
    const unsigned smooth_steps =
      chi_length > 0 ? static_cast<unsigned>(ceil(chi_length*chi_length/(2*smooth_rate))) : 0u;

    for(unsigned s=0; s<smooth_steps; ++s)
    {
      for(unsigned k=0; k<DomainSize; ++k)
      {
        const auto& d = get_neighbours(k);
        m_tmp[k] = chi[k] + smooth_rate*laplacian(chi, d, sD);
      }

      swap(chi.get_data(), m_tmp.get_data());
      BoundaryConditionsFields();
    }

    double mean = 0.;
    for(unsigned k=0; k<DomainSize; ++k)
      mean += chi[k];
    mean /= DomainSize;

    double variance = 0.;
    for(unsigned k=0; k<DomainSize; ++k)
    {
      chi[k] -= mean;
      variance += chi[k]*chi[k];
    }
    variance /= DomainSize;

    const double amplitude = sqrt(chi_noise);
    const double norm = variance > 0 ? amplitude/sqrt(variance) : 0.;
    for(unsigned k=0; k<DomainSize; ++k)
      chi[k] = chi0 + norm*chi[k];
  }
  else if(chi_config=="front")
  {
    const double length = chi_length > 0 ? chi_length : 1.;
    const double center = .5*(LX-1);

    for(unsigned k=0; k<DomainSize; ++k)
    {
      const double x = GetXPosition(k);
      chi[k] = .5*(1. - tanh((x-center)/length));
    }
  }
  else
    throw error_msg("error: chi configuration '", chi_config, "' unknown.");

  for(unsigned k=0; k<DomainSize; ++k)
    m[k] = phi[k]*chi[k];

  ProjectM();
  UpdatePhenotypeQuantities();
  BoundaryConditionsFields();
}

void GoOrGrow::UpdatePhenotypeQuantities()
{
  constexpr double phi_epsilon = 1e-12;

  for(unsigned k=0; k<DomainSize; ++k)
    chi[k] = phi[k] > phi_epsilon ? m[k]/phi[k] : 0.;

  BoundaryConditionsFields();
}

void GoOrGrow::ProjectM()
{
  for(unsigned k=0; k<DomainSize; ++k)
  {
    const double upper = phi[k] > 0 ? phi[k] : 0.;
    if(m[k] < 0)
      m[k] = 0;
    else if(m[k] > upper)
      m[k] = upper;
  }
}

void GoOrGrow::UpdateQuantitiesAtNode(unsigned k)
{
  // array placeholders for current node
  const auto& d = get_neighbours(k);
  const auto& f = ff[k];
  // Q-tensor and binary phase order
  const double Qxx = QQxx[k];
  const double Qyx = QQyx[k];
  const double p   = phi[k];

  // compute velocities
  const double nn = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];
  const double vx = (f[1] - f[2] + f[5] - f[6] - f[7] + f[8])/nn;
  const double vy = (f[3] - f[4] + f[5] - f[6] + f[7] - f[8])/nn;

  // compute derivatives etc.
  const double del2Qxx  = laplacian(QQxx,  d, sD);
  const double dxQxx    = derivX   (QQxx,  d, sB);
  const double dyQxx    = derivY   (QQxx,  d, sB);
  const double del2Qyx  = laplacian(QQyx,  d, sD);
  const double dxQyx    = derivX   (QQyx,  d, sB);
  const double dyQyx    = derivY   (QQyx,  d, sB);
  const double del2p    = laplacian(phi, d, sD);
  const double dxPhi    = derivX   (phi, d, sB);
  const double dyPhi    = derivY   (phi, d, sB);

  // computation of the chemical potential and molecular field...
  // ...term controlling the preferred nematic order magnitude
  const double q2 = Qxx*Qxx + Qyx*Qyx;
  const double term = Snem - q2;
  // ...molecular field
  const double Hxx = 2*CC*term*Qxx + LL*del2Qxx;
  const double Hyx = 2*CC*term*Qyx + LL*del2Qyx;
  // ...compressibility/crowding contribution
  const double dp_critical = p - phi_critical;
  const double f_compress = dp_critical > 0 ? .5*B*dp_critical*dp_critical : 0.;
  const double mu_compress = dp_critical > 0 ? B*dp_critical : 0.;
  // ...chemical potential
  const double mu = AA*p*(1-p)*(1-2*p) + mu_compress - KK*del2p;

  // computation of sigma...
  // ... on-diagonal stress components
  const double sigmaB = .5*AA*p*p*(1-p)*(1-p)
    + f_compress + .5*CC*term*term - mu*p;
  const double active_prefactor = zeta*p*(1-chi[k]);
  const double sigmaF = 2*xi*( (Qxx*Qxx-1)*Hxx + Qxx*Qyx*Hyx )
    - active_prefactor*Qxx + .5*KK*(dyPhi*dyPhi-dxPhi*dxPhi)
    + LL*(dyQxx*dyQxx+dyQyx*dyQyx-dxQxx*dxQxx-dxQyx*dxQyx);
  // .. off-diagonal stress components
  const double sigmaS = 2*xi*(Qyx*Qxx*Hxx + (Qyx*Qyx-1)*Hyx)
    - active_prefactor*Qyx  - KK*dxPhi*dyPhi
    - 2*LL*(dxQxx*dyQxx+dxQyx*dyQyx);
  const double sigmaA = 2*(Qxx*Hyx - Qyx*Hxx);

  // transfer to arrays
  ux[k]      =  vx;
  uy[k]      =  vy;
  ux_phi[k]  =  vx*p;
  uy_phi[k]  =  vy*p;
  n[k]       =  nn;
  HHxx[k]    =  Hxx;
  HHyx[k]    =  Hyx;
  MU[k]      =  mu;
  dxQQxx[k]  =  dxQxx;
  dxQQyx[k]  =  dxQyx;
  dyQQxx[k]  =  dyQxx;
  dyQQyx[k]  =  dyQyx;
  sigmaXX[k] =  sigmaF + sigmaB + zetaI * (conc-p);
  sigmaYY[k] = -sigmaF + sigmaB + zetaI * (conc-p);
  sigmaXY[k] =  sigmaS + sigmaA;
  sigmaYX[k] =  sigmaS - sigmaA;
  sigma_bulk[k] = sigmaB + zetaI * (conc-p);
  sigma_elastic_xx[k] = sigmaF - .5*KK*(dyPhi*dyPhi-dxPhi*dxPhi) + active_prefactor*Qxx;
  sigma_elastic_yx[k] = sigmaS + KK*dxPhi*dyPhi + active_prefactor*Qyx;
  sigma_phase_field_xx[k] = .5*KK*(dyPhi*dyPhi-dxPhi*dxPhi);
  sigma_phase_field_yx[k] = -KK*dxPhi*dyPhi;
  sigma_active_xx[k] = -active_prefactor*Qxx;
  sigma_active_yx[k] = -active_prefactor*Qyx;
}

void GoOrGrow::UpdateQuantities()
{
  // sum -> countphi
  double sum = 0;

  #pragma omp parallel for reduction (+:sum) num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
  {
    // sums the binary order to check for discrepancy
    sum = sum + phi[k];

    // do the job
    UpdateQuantitiesAtNode(k);
  }

  countphi = sum;
  UpdatePhenotypeQuantities();
}

void GoOrGrow::UpdateFields(bool first)
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
  {
    const auto& d = get_neighbours(k);

    const double division_m = ( 1. ? division_mask[k] : 0. );
    const double death_m = ( 1. ? death_mask[k] : 0. );
    const double growth_rate = alpha*division_m - beta*death_m;
    const double chi_eff = chi[k];
    const double R = chi_eff*phi[k]*growth_rate;

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
    const double Dm = GammaP*diffusiveMFlux - mFlux + chi_eff*R;

    // normal lyotropic update
    Lyotropic::UpdateFieldsAtNode(k, first);

    // same division/death correction as LyotropicWithDivision
    phi_tmp[k] += R;

    if(first)
    {
      mN[k] = m[k] + .5*Dm;
      m_tmp[k] = m[k] + Dm;
    }
    else
      m_tmp[k] = mN[k] + .5*Dm;

    const double upper = phi_tmp[k] > 0 ? phi_tmp[k] : 0.;
    if(m_tmp[k] < 0)
      m_tmp[k] = 0;
    else if(m_tmp[k] > upper)
      m_tmp[k] = upper;
  }

  swap(phi.get_data(), phi_tmp.get_data());
  swap(m.get_data(), m_tmp.get_data());
  UpdatePhenotypeQuantities();
}

void GoOrGrow::BoundaryConditionsFields()
{
  LyotropicWithDivision::BoundaryConditionsFields();

  switch(BC)
  {
    // pbc without bdry layer (nothing to do)
    case 0:
      break;
    // channel
    case 1:
    case 2:
      m.ApplyNeumannChannel();
      chi.ApplyNeumannChannel();
      break;
    // box
    case 3:
    case 4:
      m.ApplyNeumann();
      chi.ApplyNeumann();
      break;
    // pbc with bdry layer
    default:
      m.ApplyPBC();
      chi.ApplyPBC();
  }
}

option_list GoOrGrow::GetOptions()
{
  // get options from base model
  auto options = LyotropicWithDivision::GetOptions();

  // add new model options
  options[0].add_options()
    ("B", opt::value<double>(&B),
     "crowding/compressibility penalty strength")
    ("Snem", opt::value<double>(&Snem),
     "preferred value of one half Tr(Q^2)")
    ("chi-config", opt::value<string>(&chi_config),
     "phenotype initialization mode: noise or front")
    ("chi0", opt::value<double>(&chi0),
     "initial phenotype mean")
    ("chi-noise", opt::value<double>(&chi_noise),
     "initial phenotype variance")
    ("chi-length", opt::value<double>(&chi_length),
     "initial phenotype correlation length");

  return options;
}
