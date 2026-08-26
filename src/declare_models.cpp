// ======================================================================
// Model declaration (gets compiled in models.cpp)
// model headers and declare_model must be consistent!


// model headers
#include "models/minimal.hpp"
#include "models/lyotropic.hpp"
#include "models/lyotropic-free-boundary.hpp"
#include "models/lyotropic-with-division.hpp"
#include "models/lyotropic-with-division-stress.hpp"
#include "models/go-or-grow.hpp"
#include "models/dry-go-or-grow.hpp"
#include "models/confluent-memory.hpp"
#include "models/confluent-wet.hpp"
#include "models/nematic.hpp"
#include "models/dry.hpp"
#include "models/polar.hpp"
#include "models/phases.hpp"
#include "models/lyotropic-free-boundary-with-division.hpp"
#include "models/nematic-free-boundary.hpp"

void DeclareModels()
{
  declare_model<Minimal>(
     "minimal",
      "This is just an example model showing a minimal implementation. "
     "This does exactly nothing (and not particularly fast)."
    );

  declare_model<Lyotropic>(
      "lyotropic",
      "Biphasic, lyotropic, nematic model as presented as described in "
      "10.1103/PhysRevLett.113.248303. We refer the user to this reference "
      "for further information."
      );

  declare_model<LyotropicFreeBoundary>(
      "lyotropic-free-boundary",
      "The lyotropic model with more complicated boundary conditions. In "
      "principle the implementation allows for arbitrary boundary conditions. "
      "In principle."
      );

  declare_model<LyotropicWithDivision>(
      "lyotropic-with-division",
      "As the name suggests, this model is exactly the same as the lyotropic "
      "model with cell division (and death) added. Division and death are "
      "modeled separately by introducing random patches where the growth rate "
      "is locally changed."
      );
      
   declare_model<LyotropicWithDivisionStress>(
        "lyotropic-with-division-stress",
        "As the name suggests, this model is exactly the same as the lyotropic "
        "division model with cell division (and death) added. "
        "However, the division regions are now based on the local stress"
        );      

  declare_model<GoOrGrow>(
      "go-or-grow",
      "Lyotropic model with cell division and a crowding/compressibility "
      "free-energy penalty for phi above phi-critical."
      );

  declare_model<DryGoOrGrow>(
      "dry-go-or-grow",
      "Go-or-grow model with the same Q, phi, and phenotype dynamics, but "
      "with velocity solved from overdamped force balance instead of LB "
      "hydrodynamics."
      );

  declare_model<ConfluentMemory>(
      "confluent-memory",
      "Confluent proliferating active nematic with mechanical memory. Dry "
      "(overdamped) and phase-field-free: no free boundary, no vacuum, the "
      "layer is held together by the one-sided crowding modulus B alone. "
      "Carries a phenotype chi and a mechanical memory m, both transported "
      "with the cell density; m obeys tau_m*D_t m = g(P) - m with g a smoothed "
      "step in the pressure, and drives phenotype switching through "
      "V(chi,m) = Ochi*(m-mc)*chi."
      );

  declare_model<ConfluentWet>(
      "confluent-wet",
      "Confluent, incompressible active nematic with phenotype and mechanical "
      "memory on a wet (lattice Boltzmann) foundation. Rebuilt WITHOUT a phase "
      "field: derived from the traditional wet nematic model, so there is no phi, "
      "no double well, no interface term, no Cahn-Hilliard mobility, no crowding "
      "modulus and no proliferation. The layer is the whole domain. Because the "
      "flow is incompressible the pressure is the Lagrange multiplier of the "
      "constraint rather than an equation of state: the physical field is "
      "P = -1/2 Tr(Pi) = n/3 - sigma_bulk, which obeys a Poisson equation sourced "
      "by the deviatoric stress alone. The phenotype chi and the memory m are "
      "INTENSIVE scalars carried by the material derivative and advected with the "
      "true material velocity u + F/(2n), not the bare lattice-Boltzmann moment; "
      "chi modulates the activity as zeta*(1-chi) and is driven by m, which is "
      "driven by P."
      );

  declare_model<Nematic>(
      "nematic",
      "Pure nematic model with LdG free energy."
      );

  declare_model<Polar>(
      "polar",
      "Polar model with both motility and contractility. This is version with "
      "full hydrodynamic of the model described in arXiv:1705.00501."
      );

  declare_model<Phases>(
      "phases",
      "Phase-field model for the simulation of active cellular monolayers as "
      "described in an upcoming high-impact publication."
      );

  declare_model<LyotropicFreeBoundaryWithDivision>(
      "lyotropic-free-boundary-with-division",
      "As the name suggests, this model is exactly the same as the lyotropic "
      "model adapted to more complicated boundary conditions and with cell division"
      "(and death) added. Division and death are modeled separately by introducing"
      "random patches where the growth rate is locally changed."
      );

  declare_model<Dry>(
      "dry",
      "Dry version of the nematic model, where v is solved from the force "
      "balance at every time step."
      );

   declare_model<NematicFreeBoundary>(
     "nematic-free-boundary",
     "The nematic model with more complicated boundary conditions. In "
     "principle the implementation allows for arbitrary boundary conditions. "
     "In principle."
     );
   // add your models here....
}
