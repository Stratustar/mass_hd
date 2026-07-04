#include "header.hpp"
#include "models/go-or-grow.hpp"
#include "error_msg.hpp"
#include "serialization.hpp"
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
  if(relax_dt < 0)
    throw error_msg("relax-dt must be non-negative.");

  chi.SetSize(LX, LY, Type);
  m.SetSize(LX, LY, Type);
  mN.SetSize(LX, LY, Type);
  m_tmp.SetSize(LX, LY, Type);
}

void GoOrGrow::Configure()
{
  LyotropicWithDivision::Configure();

  if(init_frame.empty())
  {
    // Standard procedural/random initialization.
    RelaxFreeEnergy();
    ResetHydrodynamics();
    ConfigurePhiNoise();
    ConfigurePhenotype();
  }
  else
  {
    // Seed the dynamical fields from a saved snapshot instead.
    ConfigureFromFrame();
  }
}

void GoOrGrow::ConfigureFromFrame()
{
  // Read the snapshot once, then overwrite the fields set by the base Configure.
  const string content = read_file_to_string(init_frame);

  auto load_into = [&](const string& nm, ScalarField& fld)
  {
    const auto vals = read_frame_field(content, nm);
    if(vals.size() != DomainSize)
      throw error_msg("init-frame field '", nm, "' has ", vals.size(),
                      " values but this run has DomainSize=", DomainSize,
                      " (grid-size mismatch between snapshot and run).");
    auto& data = fld.get_data();
    for(unsigned k=0; k<DomainSize; ++k)
      data[k] = vals[k];
  };

  // Director and interface come from the snapshot in every case.
  load_into("QQxx", QQxx);
  load_into("QQyx", QQyx);
  load_into("phi",  phi);

  // Phenotype: carry over the snapshot's chi (via m), or reset it uniform (=chi0)
  // for the homogeneous control.
  if(init_frame_uniform_chi)
    for(unsigned k=0; k<DomainSize; ++k)
      m[k] = phi[k]*chi0;
  else
    load_into("m", m);

  BoundaryConditionsFields();
  ProjectM();
  UpdatePhenotypeQuantities();

  // Start hydrodynamics from rest, as after dry relaxation on the normal path.
  ResetHydrodynamics();

  // Retarget the mass reference used by conserve_phi / diagnostics.
  double total = 0.;
  for(unsigned k=0; k<DomainSize; ++k)
    total += phi[k];
  totalphi = countphi = ptot = total;
}

void GoOrGrow::RelaxFreeEnergy()
{
  if(relax_steps == 0 || (!relax_phi && !relax_Q))
    return;

  for(unsigned step=0; step<relax_steps; ++step)
  {
    Lyotropic::BoundaryConditionsFields();

    #pragma omp parallel for num_threads(nthreads) if(nthreads)
    for(unsigned k=0; k<DomainSize; ++k)
    {
      const auto& d = get_neighbours(k);
      const double Qxx = QQxx[k];
      const double Qyx = QQyx[k];
      const double p = phi[k];
      const double q2 = Qxx*Qxx + Qyx*Qyx;
      const double term = Snem*p - q2;
      const double dp_critical = p - phi_critical;
      const double mu_compress = dp_critical > 0 ? B*dp_critical : 0.;

      HHxx[k] = 2*CC*term*Qxx + LL*laplacian(QQxx, d, sD);
      HHyx[k] = 2*CC*term*Qyx + LL*laplacian(QQyx, d, sD);
      MU[k] = AA*p*(1-p)*(1-2*p) + mu_compress
        + CC*Snem*term - KK*laplacian(phi, d, sD);
    }

    switch(BC)
    {
      case 0:
        break;
      case 1:
      case 2:
        MU.ApplyNeumannChannel();
        break;
      case 3:
      case 4:
        MU.ApplyNeumann();
        break;
      default:
        MU.ApplyPBC();
    }

    #pragma omp parallel for num_threads(nthreads) if(nthreads)
    for(unsigned k=0; k<DomainSize; ++k)
    {
      const auto& d = get_neighbours(k);

      if(relax_phi)
        phi_tmp[k] = phi[k] + relax_dt*GammaP*laplacian(MU, d, sD);
      else
        phi_tmp[k] = phi[k];

      if(relax_Q)
      {
        QNxx[k] = QQxx[k] + relax_dt*GammaQ*HHxx[k];
        QNyx[k] = QQyx[k] + relax_dt*GammaQ*HHyx[k];
      }
      else
      {
        QNxx[k] = QQxx[k];
        QNyx[k] = QQyx[k];
      }
    }

    swap(phi.get_data(), phi_tmp.get_data());
    swap(QQxx.get_data(), QNxx.get_data());
    swap(QQyx.get_data(), QNyx.get_data());
  }

  Lyotropic::BoundaryConditionsFields();

  double relaxed_phi = 0.;
  for(unsigned k=0; k<DomainSize; ++k)
    relaxed_phi += phi[k];
  totalphi = relaxed_phi;
  countphi = relaxed_phi;
  ptot = relaxed_phi;
}

void GoOrGrow::ResetHydrodynamics()
{
  double ftotal = 0.;

  #pragma omp parallel for reduction(+:ftotal) num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
  {
    ux[k] = uy[k] = ux_phi[k] = uy_phi[k] = 0.;
    n[k] = rho;

    const auto fe = GetEquilibriumDistribution(0., 0., rho);
    ff[k] = fe;
    fn[k] = fe;
    ff_tmp[k] = fe;
    fn_tmp[k] = fe;

    ftotal = accumulate(begin(fe), end(fe), ftotal);
  }

  ftot = ftotal;
}

void GoOrGrow::ConfigurePhiNoise()
{
  // Overlay a multiplicative, spatially correlated modulation phi *= max(1+xi, 0)
  // on the relaxed initial phi, where xi is a zero-mean Gaussian field of standard
  // deviation phi_noise smoothed to correlation length phi_length. Reuses the chi
  // noise recipe; chi is used as workspace (overwritten by ConfigurePhenotype next).
  if(phi_noise <= 0.)
    return;

  for(unsigned k=0; k<DomainSize; ++k)
    chi[k] = random_normal();

  BoundaryConditionsFields();

  constexpr double smooth_rate = 0.2;
  const unsigned smooth_steps =
    phi_length > 0 ? static_cast<unsigned>(ceil(phi_length*phi_length/(2*smooth_rate))) : 0u;

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

  const double norm = variance > 0 ? phi_noise/sqrt(variance) : 0.;

  for(unsigned k=0; k<DomainSize; ++k)
  {
    const double factor = 1. + norm*chi[k];
    phi[k] *= factor > 0 ? factor : 0.;
  }

  // The modulation changes the total mass; retarget the conserve_phi reference.
  totalphi = 0.;
  for(unsigned k=0; k<DomainSize; ++k)
    totalphi += phi[k];
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
  for(unsigned k=0; k<DomainSize; ++k)
  {
    if(!HasPhenotypeMaterial(k))
    {
      phi[k] = 0.;
      m[k] = 0.;
      chi[k] = 0.;
    }
    else
      chi[k] = m[k]/phi[k];
  }

  BoundaryConditionsFields();
}

double GoOrGrow::PhenotypePhiEpsilon() const
{
  return 1e-12;
}

bool GoOrGrow::HasPhenotypeMaterial(unsigned k) const
{
  return phi[k] > PhenotypePhiEpsilon();
}

bool GoOrGrow::HasMaterialFace(unsigned k, unsigned neighbour) const
{
  return HasPhenotypeMaterial(k) || HasPhenotypeMaterial(neighbour);
}

double GoOrGrow::MaskedPhiFaceIncrement(unsigned k, unsigned neighbour, double increment) const
{
  return HasMaterialFace(k, neighbour) ? increment : 0.;
}

double GoOrGrow::LocalMaterialChi(unsigned k, bool& found) const
{
  if(HasPhenotypeMaterial(k))
  {
    found = true;
    return chi[k];
  }

  const auto& d = get_neighbours(k);
  double sum = 0.;
  unsigned count = 0;

  for(unsigned q=1; q<9; ++q)
    if(HasPhenotypeMaterial(d[q]))
    {
      sum += chi[d[q]];
      ++count;
    }

  found = count > 0;
  return found ? sum/count : 0.;
}

double GoOrGrow::MaterialChi(unsigned preferred, unsigned fallback) const
{
  if(HasPhenotypeMaterial(preferred))
    return chi[preferred];
  if(HasPhenotypeMaterial(fallback))
    return chi[fallback];

  bool found = false;
  const double preferred_chi = LocalMaterialChi(preferred, found);
  if(found)
    return preferred_chi;

  const double fallback_chi = LocalMaterialChi(fallback, found);
  if(found)
    return fallback_chi;

  return 0.;
}

double GoOrGrow::TransportFaceChi(unsigned k, unsigned neighbour, double dphi_k) const
{
  return dphi_k >= 0 ? MaterialChi(neighbour, k) : MaterialChi(k, neighbour);
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
  const double term = Snem*p - q2;
  // ...molecular field
  const double Hxx = 2*CC*term*Qxx + LL*del2Qxx;
  const double Hyx = 2*CC*term*Qyx + LL*del2Qyx;
  // ...compressibility/crowding contribution
  const double dp_critical = p - phi_critical;
  const double f_surface = .5*AA*p*p*(1-p)*(1-p);
  const double mu_surface = AA*p*(1-p)*(1-2*p) - KK*del2p;
  const double f_compress = dp_critical > 0 ? .5*B*dp_critical*dp_critical : 0.;
  const double mu_compress = dp_critical > 0 ? B*dp_critical : 0.;
  // ...chemical potential
  const double mu = mu_surface + mu_compress + CC*Snem*term;

  // computation of sigma...
  // ... on-diagonal stress components
  const double sigma_surface_bulk = f_surface - mu_surface*p;
  const double sigmaB = f_compress + .5*CC*term*term
    - (mu_compress + CC*Snem*term)*p
    + (surface_stress ? sigma_surface_bulk : 0.);
  const double active_prefactor = zeta*p*(1-MaterialChi(k, k));
  const double phase_xx = surface_stress ? .5*KK*(dyPhi*dyPhi-dxPhi*dxPhi) : 0.;
  const double phase_yx = surface_stress ? -KK*dxPhi*dyPhi : 0.;
  const double sigmaF = 2*xi*( (Qxx*Qxx-1)*Hxx + Qxx*Qyx*Hyx )
    - active_prefactor*Qxx + phase_xx
    + LL*(dyQxx*dyQxx+dyQyx*dyQyx-dxQxx*dxQxx-dxQyx*dxQyx);
  // .. off-diagonal stress components
  const double sigmaS = 2*xi*(Qyx*Qxx*Hxx + (Qyx*Qyx-1)*Hyx)
    - active_prefactor*Qyx + phase_yx
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
  sigma_elastic_xx[k] = sigmaF - phase_xx + active_prefactor*Qxx;
  sigma_elastic_yx[k] = sigmaS - phase_yx + active_prefactor*Qyx;
  sigma_phase_field_xx[k] = phase_xx;
  sigma_phase_field_yx[k] = phase_yx;
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

    const double division_m = division_mask[k] ? 1. : 0.;
    const double death_m = death_mask[k] ? 1. : 0.;
    const double growth_rate = alpha*division_m - beta*death_m;
    const double chi_eff = MaterialChi(k, k);
    const double R = chi_eff*phi[k]*growth_rate;

    const double dphi_px_raw = GammaP*(MU[d[1]] - MU[k]) - .5*ux_phi[d[1]];
    const double dphi_mx_raw = GammaP*(MU[d[2]] - MU[k]) + .5*ux_phi[d[2]];
    const double dphi_py_raw = GammaP*(MU[d[3]] - MU[k]) - .5*uy_phi[d[3]];
    const double dphi_my_raw = GammaP*(MU[d[4]] - MU[k]) + .5*uy_phi[d[4]];
    const double dphi_px = MaskedPhiFaceIncrement(k, d[1], dphi_px_raw);
    const double dphi_mx = MaskedPhiFaceIncrement(k, d[2], dphi_mx_raw);
    const double dphi_py = MaskedPhiFaceIncrement(k, d[3], dphi_py_raw);
    const double dphi_my = MaskedPhiFaceIncrement(k, d[4], dphi_my_raw);
    const double phiRawTransport = dphi_px_raw + dphi_mx_raw + dphi_py_raw + dphi_my_raw
      - (conserve_phi ? (countphi-totalphi)/DomainSize : 0.);
    const double phiTransport = dphi_px + dphi_mx + dphi_py + dphi_my
      - (conserve_phi ? (countphi-totalphi)/DomainSize : 0.);
    const double phiTransportCorrection = phiTransport - phiRawTransport;
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
      .5*(phi[d[1]] + phi[k])*(chi_px - chi_eff)
    - .5*(phi[k] + phi[d[2]])*(chi_eff - chi_mx)
    + .5*(phi[d[3]] + phi[k])*(chi_py - chi_eff)
    - .5*(phi[k] + phi[d[4]])*(chi_eff - chi_my);
    const double press = -sigma_bulk[k];
    const double dVdChi =
      2*Achi*chi_eff*(1-chi_eff)*(1-2*chi_eff)
    + Ochi*(press-pswitch);
    const double Sswitch = -phi[k]*dVdChi;
    const double mGrowth = growTogether ? chi_eff*R : R;
    const double Dm = mTransport + Dchi*phenotypeDiffusion + Sswitch + mGrowth;

    // normal lyotropic update
    Lyotropic::UpdateFieldsAtNode(k, first);

    // Apply growth with the same predictor-corrector timing used for m.
    if(first)
    {
      phn[k] += .5*(phiTransportCorrection + R);
      phi_tmp[k] += phiTransportCorrection + R;
    }
    else
      phi_tmp[k] += .5*(phiTransportCorrection + R);

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
    ("Dchi", opt::value<double>(&Dchi),
     "phenotype diffusion coefficient")
    ("Achi", opt::value<double>(&Achi),
     "phenotype switching double-well barrier")
    ("Ochi", opt::value<double>(&Ochi),
     "phenotype switching bias strength")
    ("pswitch", opt::value<double>(&pswitch),
     "pressure threshold for phenotype switching bias")
    ("growTogether", opt::value<int>(&growTogether),
     "m growth source: 0 uses R, 1 uses chi*R")
    ("surface-stress", opt::value<int>(&surface_stress),
     "whether phase-field surface terms contribute to stress and velocity")
    ("chi-config", opt::value<string>(&chi_config),
     "phenotype initialization mode: noise or front")
    ("chi0", opt::value<double>(&chi0),
     "initial phenotype mean")
    ("chi-noise", opt::value<double>(&chi_noise),
     "initial phenotype variance")
    ("chi-length", opt::value<double>(&chi_length),
     "initial phenotype correlation length")
    ("init-frame", opt::value<string>(&init_frame),
     "seed Q, phi (and chi) from this frame*.json snapshot instead of random init")
    ("init-frame-uniform-chi", opt::value<int>(&init_frame_uniform_chi),
     "with init-frame, set chi uniform (=chi0) instead of loading it from the snapshot")
    ("phi-noise", opt::value<double>(&phi_noise),
     "std of multiplicative correlated-noise phi modulation")
    ("phi-length", opt::value<double>(&phi_length),
     "correlation length of phi noise modulation")
    ("relax-steps", opt::value<unsigned>(&relax_steps),
     "dry free-energy relaxation steps before official dynamics")
    ("relax-dt", opt::value<double>(&relax_dt),
     "time step multiplier for dry free-energy relaxation")
    ("relax-phi", opt::value<int>(&relax_phi),
     "relax phi by Cahn-Hilliard free-energy descent before official dynamics")
    ("relax-Q", opt::value<int>(&relax_Q),
     "relax Q by molecular-field free-energy descent before official dynamics");

  return options;
}
