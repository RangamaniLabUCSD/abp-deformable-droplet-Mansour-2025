# Modified Cytosim Source Code

## Origin
The base Cytosim version can be found at https://gitlab.com/f-nedelec/cytosim. The modifications introduced by the authors to the dynamic_ellipse fork of Cytosim are documented in the commit history from 2024 onward and build upon the extensive preexisting development that produced the underlying Cytosim version.

## Capping Protein
Capping protein was implemented implicitly as a special type of growing fiber introduced in this work. The following files are the major source code modifications:
- Cytosim/src/sim/fibers/capping_fiber.cc
- Cytosim/src/sim/fibers/capping_fiber.h
- Cytosim/src/sim/fibers/capping_fiber_prop.cc
- Cytosim/src/sim/fibers/capping_fiber_prop.h

## Dynamic Ellipsoidal Deformation
The implementation of dynamic ellipsoidal deformation was built as part of the dynamic_ellipse fork of the base Cytosim version and presented in Dmitrieff et al. 2017. The following files are the major modifications that introduced the simple deformable boundary to Cytosim:
- Cytosim/src/sim/spaces/space_dynamic_ellipse.cc
- Cytosim/src/sim/spaces/space_dynamic_ellipse.h
- Cytosim/src/sim/spaces/space_dynamic_prop.cc
- Cytosim/src/sim/spaces/space_dynamic_prop.h

## Dynamic Multimerization
Dynamically multimerizing crosslinkers was implemented as an extension of the dynamic dimers introduced in Walker et al. 2025. The following files are the major source code modifications:
- Cytosim/src/sim/hands/dimerizer.cc
- Cytosim/src/sim/hands/dimerizer.h
- Cytosim/src/sim/hands/dimerizer_prop.cc
- Cytosim/src/sim/hands/dimerizer_prop.h