// ======================================================================
// Model declaration (gets compiled in models.cpp)
// model headers and declare_model must be consistent!


// model headers
#include "models/minimal.hpp"
#include "models/lyotropic.hpp"
#include "models/lyotropic-free-boundary.hpp"
#include "models/lyotropic-with-division.hpp"
#include "models/lyotropic-with-division-stress.hpp"
#include "models/nematic.hpp"
#include "models/dry.hpp"
#include "models/polar.hpp"
#include "models/phases.hpp"
#include "models/lyotropic-free-boundary-with-division.hpp"
#include "models/nematic-free-boundary.hpp"
#include "models/nematic-compressible-division.hpp"
#include "models/nematic-compressible-division-softpressure.hpp"

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

   declare_model<NematicCompressibleDivision>(
     "nematic_compressible_division",
     "Compressible wet nematic model with patch-based proliferation and an "
     "additional density-dependent pressure term."
     );

   declare_model<NematicCompressibleDivisionSoftPressure>(
     "nematic_compressible_division_softpressure",
     "Compressible wet nematic model with patch-based proliferation and a "
     "softer rational density-dependent pressure term."
     );
   // add your models here....
}
