#ifndef MODELS_PHASES_HPP_
#define MODELS_PHASES_HPP_

#include "models.hpp"

/** Phase field simulation of cellular monolayers
 *
 * We follow reference doi:10.1038/srep11745. Each cell is represented by a
 * unique phase field and the dynamics is given by
 *
 *    \partial_t \phi + v \partial_x \phi = -\frac 1 2 V
 *
 * where V = \delta F / \delta \phi is computed from a free energy. The
 * velocity v contains both a passive and an acitve contribution. As it stands
 * this is a dry model, but hydrodynamics can be easily added later. Please
 * refer to the original paper for the definition of the coefficients and more
 * detailed discussion.
 *
 * Solution of the evolution equation is performed using a predictor corrector
 * method (with an arbitrary high number of correction steps).
 *
 * -----------------------------------------------------------------------------
 *
 * In order to reduce the amount of computations associated with a high number
 * of different phases we use the fact that the phases are usually well
 * located (blobs that represent cells) and that it is unnecessary to take a
 * phase field into account far outside this domain. Hence we track each cell
 * and update only the region around it.
 *
 * The domains are obtained by computing the center of mass (com) of each cell
 * and by adding a 'margin' of nodes on each side.
 *
 * The walls are implemented as if they were another cell... more to come here.
 */
class Phases : public Model
{
protected:
  /** Phases */
  std::vector<ScalarField> phi;
  /** Predicted phi in a PC^n step */
  std::vector<ScalarField> phi_old;
  /** V = delta F / delta phi */
  std::vector<ScalarField> V;
  /** Potential: -V/2 - adv. term */
  std::vector<ScalarField> potential;
  /** Predicted potential in a P(C)^n step */
  std::vector<ScalarField> potential_old;
  /** Active velocity angle and angular acceleration */
  std::vector<double> theta, dtheta;
  /** Area associated with a phase field */
  std::vector<double> area;
  /** Counter for computing the area */
  std::vector<double> area_cnt;
  /** Sum of phi at each node */
  ScalarField sum, sum_cnt;
  /** Sum of square phi at each node */
  ScalarField square, square_cnt;
  /** Phase-field for the walls */
  ScalarField walls, walls_dx, walls_dy, walls_laplace;
  /** Cell polarity */
  std::vector<std::array<double, 2>> pol;
  /** Passive velocities */
  std::vector<std::array<double, 2>> velp;
  /** Contraction velocities */
  std::vector<std::array<double, 2>> velc;
  /** Friction velocities */
  std::vector<std::array<double, 2>> velf;
  /** Total com velocity */
  std::vector<std::array<double, 2>> vel;
  /** Alignment torque */
  // std::vector<double> torque;
  /** Center-of-mass */
  std::vector<std::array<double, 2>> com, com_prev;
  /** Overall polarization of the tissue */
  ScalarField Px, Py, Px_cnt, Py_cnt;
  /** Parameters for the polarisation dynamics */
  double D=1., J=1.;
  /** Contractility parameters */
  double c0, tauc = 1;
  /** Contractility */
  std::vector<double> c;
  /** Shape parameter: order */
  std::vector<double> S_order;
  /** Shape parameter: angle */
  std::vector<double> S_angle;
  /** Structure tensor */
  std::vector<double> S00, S01;
  /** Polarity tensor */
  ScalarField Q00, Q01;
  ScalarField Theta;
  std::vector<double> Theta_cnt;
  /** Counters for polarity tensor */
  std::vector<double> Q00_cnt, Q01_cnt;
  /** Internal pressure */
  std::vector<double> P, P_cnt;

  /** Number of predictor-corrector steps */
  unsigned npc = 1;
  /** Relaxation time at initialization */
  unsigned relax_time;
  /** Value of nsubstep to use for initialization */
  unsigned relax_nsubsteps;

  /** Enable tracking? */
  bool tracking = false;
  /** Min of the boundaries of the domains and center of mass */
  std::vector<std::array<unsigned, 2>> domain_min;
  /** Max of the boundaries of the domains and center of mass */
  std::vector<std::array<unsigned, 2>> domain_max;
  /** Counter to compute com in Fourier space */
  std::vector<std::complex<double>> com_x;
  /** Counter to compute com in Fourier space */
  std::vector<std::complex<double>> com_y;
  /** Margin for the definition of domains */
  unsigned margin;
  /** Precomputed tables for sin and cos (as complex numbers) used in the
   * computation of the com.
   * */
  std::vector<std::complex<double>> com_x_table, com_y_table;

  /** Number of phases */
  unsigned nphases;
  /** Elasticity */
  std::vector<double> gamma;
  /** Energy penalty for area */
  std::vector<double> mu;
  /** Interface thickness */
  double lambda;
  /**  Interaction stength */
  double kappa;
  /** Friction parameter */
  double zeta = 0.;
  /** Cell-cell friction parameter */
  double f;
  /** Cell-wall friction parameter */
  double f_walls;
  /** Substrate friction parameter */
  double xi;
  /** Adhesion parameter */
  double omega;
  /** Prefered radius (area = pi*R*R) */
  double R;
  /** Migration speed */
  double alpha;
  /** Coupling between area and contractility */
  double beta;
  /** Intial configuration */
  std::string init_config;
  /** Needa store? (yes this is dirty) */
  bool store;
  /** Noise level for initial configurations */
  double noise = 0;

  /** Division flag */
  bool division = false;
  /** Division rate */
  double division_rate = 0.;
  /** Relaxation time for division
   *
   * This is the time we give the just-divided cells to relax with all the other
   * cells fixed. See Divide().
   * */
  unsigned division_relax_time = 100;

  /** Wall thickness */
  double wall_thickness = 1.;
  /** Repuslion by the wall */
  double wall_kappa;
  /** Adhesion on the wall */
  double wall_omega;

  /** Boudaries for cell generation
   *
   * These are the boundaries (min and max x and y components) of the domain in
   * which the cells are created when the initial config 'random' is choosen.
   * */
  std::vector<unsigned> birth_bdries;

  /** Pre-computed coefficients */
  double C1, C2, C3;

  /** Helper function
   *
   * Update the fields in a square domain that is entirely contained inside the
   * domain, i.e. that is not wrapping around the borders.
   * */
  template<typename Ret, typename ...Args>
  void UpdateSubDomain(Ret (Phases::*)(unsigned, unsigned, Args...),
                       unsigned, unsigned, unsigned, unsigned,
                       unsigned, Args&&... args);
  /** Parallel version */
  template<typename Ret, typename ...Args>
  void UpdateSubDomainP(Ret (Phases::*)(unsigned, unsigned, Args...),
                        unsigned, unsigned, unsigned, unsigned,
                        unsigned, Args&&... args);
  /** Helper function
   *
   * This function is used to updated the fields only in a restricted domain
   * around the cell center. One needs to be careful because of the periodic
   * boundary conditions. The template argument is the function used to update
   * the fields at each node (called ***AtNode()).
   * */
  template<typename R, typename ...Args>
  void UpdateDomain(R (Phases::*)(unsigned, unsigned, Args...),
                    unsigned, Args&&... args);
  /** Parallel version */
  template<typename R, typename ...Args>
  void UpdateDomainP(R (Phases::*)(unsigned, unsigned, Args...),
                     unsigned, Args&&... args);

  /** Small helper function to add a cell at a certain point */
  void AddCell(unsigned, const std::array<unsigned, 2>&);

  /** Compute center of mass of a given phase field */
  inline void ComputeCoM(unsigned);

  /** Initialize all qties that depend on the number of phases
   *
   * We put this in a separate function because it is reused when a cell
   * divides.
   * */
  void InitializeFields();

public:
  Phases(unsigned, unsigned, unsigned);

  // functions from base class Model
  virtual void Initialize();
  virtual void Step();
  virtual void Configure();
  virtual void Pre();
  virtual void PreRunStats();
  virtual void RuntimeChecks();
  virtual void RuntimeStats();
  virtual option_list GetOptions();

  /** Configure boundary conditions, i.e. the walls */
  void ConfigureWalls();

  /** Predictor-corrector function for updating the system */
  void Update();

  /** Subfunction for update */
  inline void UpdateAtNode(unsigned, unsigned);
  /** Subfunction for update */
  inline void UpdateFieldsAtNode(unsigned, unsigned);
  /** Subfunction for update */
  inline void UpdateStructureTensorAtNode(unsigned, unsigned);
  /** Subfunction for update */
  inline void UpdateFrictionForceAtNode(unsigned, unsigned);
  /** Subfunction for update */
  inline void SquareAndSumAtNode(unsigned, unsigned);
  /** Subfunction for update */
  inline void ReinitSquareAndSumAtNode(unsigned);

  /** Make a cell divide
   *
   * The strategy for implementing division is to chose a division axis randomly
   * then divide the given cell in two while fixing all the other cells. We then
   * let the two new cells relax while fixing the other cells such that they are
   * able to create a common interface.
   * */
  void Divide(unsigned i);

  /** Update polarisation of a given field
   *
   * This function updates the polarisation of the cell which give the direction
   * of the active velocity of the cell. We follow reference 10.1101/095133 and
   * define the dynamics as
   *
   *    d theta / dt = J_r torque + 2 D_r eta
   *
   * where eta is gaussian white noise with zero mean and unit variance, see
   * paper for more details. Note that we use the euler-maruyama method, instead
   * of a predictor-corrector method.
   * */
  void UpdatePolarization(unsigned);

  /** Update friction force
   *
   * This function updates the friction force and needs to be in a separate
   * loop because it must be computed after the passive part of the velocity has
   * been computed fully. See paper for more details.
   * */
  void UpdateFrictionForce();

  /** Compute shape parameters
   *
   * This function effectively computes the second moment of area, which ca n be used to
   * fit the shape of a cell to an ellipse.
   * */
  void ComputeShape(unsigned);

  /** Update the window for tracking */
  void UpdateWindow(unsigned);

  /** Serialization of parameters */
  template<class Archive>
  void serialize_params(Archive& ar)
  {
    ar & auto_name(gamma)
       & auto_name(mu)
       & auto_name(nphases)
       & auto_name(lambda)
       & auto_name(kappa)
       & auto_name(alpha)
       & auto_name(R)
       & auto_name(xi)
       & auto_name(omega)
       & auto_name(init_config)
       & auto_name(zeta)
       & auto_name(D)
       & auto_name(J)
       & auto_name(f)
       & auto_name(f_walls)
       & auto_name(wall_thickness)
       & auto_name(wall_kappa)
       & auto_name(wall_omega)
       & auto_name(walls);

    ar & auto_name(tracking)
       & auto_name(margin);
  }

  /** Serialization of the current frame */
  template<class Archive>
  void serialize_frame(Archive& ar)
  {
    ar & auto_name(phi)
       & auto_name(area)
       & auto_name(com)
       & auto_name(S_order)
       & auto_name(S_angle)
       & auto_name(pol)
       & auto_name(velp)
       & auto_name(velf)
       & auto_name(velc)
       & auto_name(vel);
    if(tracking) ar
       & auto_name(domain_min)
       & auto_name(domain_max);
  }
};

// =============================================================================
// Implementation

extern unsigned nthreads;

template<typename Ret, typename ...Args>
void Phases::UpdateSubDomain(Ret (Phases::*fun)(unsigned, unsigned, Args...),
                             unsigned n,
                             unsigned m0, unsigned m1,
                             unsigned M0, unsigned M1,
                             Args&&... args)
{
  // only update on the subregion
  for(unsigned i=m0; i<M0; ++i)
    for(unsigned j=m1; j<M1; ++j)
      // if you want to look it up, this is called a pointer
      // to member function and is an obscure C++ feature...
      (this->*fun)(n, GetDomainIndex(i, j),
                   std::forward<Args>(args)...);
}

template<typename Ret, typename ...Args>
void Phases::UpdateSubDomainP(Ret (Phases::*fun)(unsigned, unsigned, Args...),
                              unsigned n,
                              unsigned m0, unsigned m1,
                              unsigned M0, unsigned M1,
                              Args&&... args)
{
  // same but with openmp
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<(M0-m0)*(M1-m1); ++k)
      // if you want to look it up, this is called a pointer
      // to member function and is an obscure C++ feature...
      (this->*fun)(n, GetDomainIndex(m0+k%(M0-m0), m1+k/(M0-m0)),
                   std::forward<Args>(args)...);
}

template<typename Ret, typename ...Args>
void Phases::UpdateDomainP(Ret (Phases::*fun)(unsigned, unsigned, Args...),
                           unsigned n, Args&&... args)
{
  if(domain_min[n][0]>=domain_max[n][0] and
     domain_min[n][1]>=domain_max[n][1])
  {
    // domain is across the corners
    UpdateSubDomainP(fun, n, domain_min[n][0], domain_min[n][1], LX, LY);
    UpdateSubDomainP(fun, n, 0u, 0u, domain_max[n][0], domain_max[n][1]);
    UpdateSubDomainP(fun, n, domain_min[n][0], 0u, LX, domain_max[n][1]);
    UpdateSubDomainP(fun, n, 0u, domain_min[n][1], domain_max[n][0], LY);
  }
  else if(domain_min[n][0]>=domain_max[n][0])
  {
    // domain is across the left/right border
    UpdateSubDomainP(fun, n, domain_min[n][0], domain_min[n][1], LX, domain_max[n][1]);
    UpdateSubDomainP(fun, n, 0u, domain_min[n][1], domain_max[n][0], domain_max[n][1]);
  }
  else if(domain_min[n][1]>=domain_max[n][1])
  {
    // domain is across the up/down border
    UpdateSubDomainP(fun, n, domain_min[n][0], domain_min[n][1], domain_max[n][0], LY);
    UpdateSubDomainP(fun, n, domain_min[n][0], 0u, domain_max[n][0], domain_max[n][1]);
  }
  else
    // domain is in the middle
    UpdateSubDomainP(fun, n, domain_min[n][0], domain_min[n][1], domain_max[n][0], domain_max[n][1]);
}

template<typename Ret, typename ...Args>
void Phases::UpdateDomain(Ret (Phases::*fun)(unsigned, unsigned, Args...),
                          unsigned n, Args&&... args)
{
  if(domain_min[n][0]>=domain_max[n][0] and
     domain_min[n][1]>=domain_max[n][1])
  {
    // domain is across the corners
    UpdateSubDomain(fun, n, domain_min[n][0], domain_min[n][1], LX, LY);
    UpdateSubDomain(fun, n, 0u, 0u, domain_max[n][0], domain_max[n][1]);
    UpdateSubDomain(fun, n, domain_min[n][0], 0u, LX, domain_max[n][1]);
    UpdateSubDomain(fun, n, 0u, domain_min[n][1], domain_max[n][0], LY);
  }
  else if(domain_min[n][0]>=domain_max[n][0])
  {
    // domain is across the left/right border
    UpdateSubDomain(fun, n, domain_min[n][0], domain_min[n][1], LX, domain_max[n][1]);
    UpdateSubDomain(fun, n, 0u, domain_min[n][1], domain_max[n][0], domain_max[n][1]);
  }
  else if(domain_min[n][1]>=domain_max[n][1])
  {
    // domain is across the up/down border
    UpdateSubDomain(fun, n, domain_min[n][0], domain_min[n][1], domain_max[n][0], LY);
    UpdateSubDomain(fun, n, domain_min[n][0], 0u, domain_max[n][0], domain_max[n][1]);
  }
  else
    // domain is in the middle
    UpdateSubDomain(fun, n, domain_min[n][0], domain_min[n][1], domain_max[n][0], domain_max[n][1]);
}

#endif//MODELS_PHASES_HPP_
