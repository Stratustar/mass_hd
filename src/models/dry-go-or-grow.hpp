#ifndef MODELS_DRY_GO_OR_GROW_HPP_
#define MODELS_DRY_GO_OR_GROW_HPP_

#include "models/go-or-grow.hpp"

class DryGoOrGrow : public GoOrGrow
{
protected:
  /** Compute velocity from overdamped force balance */
  void ComputeDryVelocity();
  /** Compute velocity at one node from div(sigma) = friction*u */
  void ComputeDryVelocityAtNode(unsigned);
  /** Update Q, phi, and m without LB collision/streaming */
  void UpdateDryFieldsAtNode(unsigned, bool);

  /** Update fields using predictor-corrector method */
  virtual void UpdateFields(bool);

public:
  DryGoOrGrow(unsigned, unsigned, unsigned);

  // functions from base class Model
  virtual void Initialize();
  virtual void Step();
};

#endif//MODELS_DRY_GO_OR_GROW_HPP_
