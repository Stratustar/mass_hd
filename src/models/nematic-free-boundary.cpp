#include "header.hpp"
#include "models/nematic-free-boundary.hpp"
#include "error_msg.hpp"
#include "random.hpp"
#include "lb.hpp"
#include "tools.hpp"

using namespace std;
namespace opt = boost::program_options;

// from main.cpp:
extern unsigned nthreads, nsubsteps, verbose;

NematicFreeBoundary::NematicFreeBoundary(unsigned LX, unsigned LY, unsigned BC_):Nematic( LX,  LY,  BC_)
{
  BC=BC_;

  //In certain cases (when we need additional nodes), the grid has to be reset.
  //This also redefines the neighbour list.

  if (BC==101 or BC==102 or BC==103)
  {
    SetSize(LX, LY, GridType::SquareAn);
  }
  else
  {
    SetSize(LX, LY, GridType::Periodic);
  }

  outside.resize(DomainSize);
}

void NematicFreeBoundary::ConfigureBoundaries()
{
  // define the inside region
  for(unsigned k=0; k<DomainSize; ++k)
  {
    const unsigned x = GetXPosition(k);
    const unsigned y = GetYPosition(k);

    // pbc
    if(BC==0)
      outside[k] = false;
    // normal boxes with free- or no-slip and von Neumann or more complex bc
    else if(BC==1)
      outside[k] = ( x<=0 or x>=LX-1 or y<=0 or y>=LY-1 );

//BC 101 and 102 are meant to be FreeSlip. 103 is No Slip.
//BC 101 is not meant to work in theory
    else if (BC==101 or BC==102 or BC==103)
    {
//        outside[k] = ( x<=0 or x>=LX-1 or y<=0 or y>=LY-1 or (x>=SqX_left and x<= SqX_right and y<= SqY_back and y>= SqY_front));
      if (x==0 or x== LX-1 )
        {outside[k] = true;}
      else if ( y==0 or y==LY-1)
        {outside[k] = true;}
      else if (x>=(SqX_left ) and x<= (SqX_right) and y<= (SqY_back) and y>= (SqY_front))
        {outside[k] = true;}
      else
        {outside[k]=false;}
    }

    else
    {
      throw error_msg("Unknown boundary condition bc=", BC);
    }
  }
  //for (unsigned k = 0; k < GetSize(BoundaryLayer_enum::BoundaryLayer); ++k)
  //{
  //  unsigned index= GetIndex(ExtraNode_enum::ExtraNode, k);
  //  outside[index]=true;
  //}
}

void NematicFreeBoundary::Configure()
{
  // intialise the outside array
  ConfigureBoundaries();

  for(unsigned k=0; k<DomainSize; ++k)
  {
    // skip outside cells
    if(outside[k]) continue;
    // see Lyotropic::Configure()
    ConfigureAtNode(k);
  }
}

void NematicFreeBoundary::UpdateQuantities()
{

  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
  {
    // skip outside cells
    if(outside[k]) continue;
    // see Nematic::UpdateQuantities()
    UpdateQuantitiesAtNode(k);
  }

}

void NematicFreeBoundary::UpdateFields(bool first)
{
#ifdef DEBUG
  double sum3=0;
  double sum2=0;
  #pragma omp parallel for reduction (+:sum2,sum3) num_threads(nthreads) if(nthreads)
#else
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
#endif

  for(unsigned k=0; k<DomainSize; ++k)
  {
    // skip outside cells
    if(outside[k]) continue;
    // see Lyotropic::UpdateFields()
    UpdateFieldsAtNode(k, first);
  }
#ifdef DEBUG
  countphi_laplace=sum3;
  countphi_flux=sum2;
#endif
}

void NematicFreeBoundary::BoundaryConditionsLB()
{
  switch(BC)
  {
    // pbc without bdry layer (nothing to do)
    case 0:
      break;
    // free-slip box
    case 1:
    {
      auto apply_bc = [this](LBField& field) {
        // walls
        field.ApplyFreeSlipFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyFreeSlipFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyFreeSlipFreeWall(1, 0, LX-2, Wall::Front);
        field.ApplyFreeSlipFreeWall(1, LY-1, LX-2, Wall::Back);
        // corners
        field.ApplyFreeSlipFreeConcaveCorner(LX-1,LY-1,Corner::RightBack);
        field.ApplyFreeSlipFreeConcaveCorner(LX-1,0,Corner::RightFront);
        field.ApplyFreeSlipFreeConcaveCorner(0,LY-1,Corner::LeftBack);
        field.ApplyFreeSlipFreeConcaveCorner(0,0,Corner::LeftFront);

      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }

        case 101:
        case 102:
    {
      auto apply_bc = [this](LBField& field) {
        // walls
        field.ApplyFreeSlipFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyFreeSlipFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyFreeSlipFreeWall(1, 0, LX-2, Wall::Front);
        field.ApplyFreeSlipFreeWall(1, LY-1, LX-2, Wall::Back);
        // corners
        field.ApplyFreeSlipFreeConcaveCorner(LX-1,LY-1,Corner::RightBack);
        field.ApplyFreeSlipFreeConcaveCorner(LX-1,0,Corner::RightFront);
        field.ApplyFreeSlipFreeConcaveCorner(0,LY-1,Corner::LeftBack);
        field.ApplyFreeSlipFreeConcaveCorner(0,0,Corner::LeftFront);

        // walls
        field.ApplyFreeSlipFreeWall(SqX_left, SqY_front+1, SqWallHeight, Wall::Right);
        field.ApplyFreeSlipFreeWall(SqX_right, SqY_front+1 , SqWallHeight, Wall::Left);
        field.ApplyFreeSlipFreeWall(SqX_left+1, SqY_back , SqWallLength, Wall::Front);
        field.ApplyFreeSlipFreeWall(SqX_left+1,  SqY_front, SqWallLength, Wall::Back);
        // corners
        field.ApplyFreeSlipFreeConvexCorner(SqX_right , SqY_back ,Corner::RightBack);
        field.ApplyFreeSlipFreeConvexCorner(SqX_right ,SqY_front ,Corner::RightFront);
        field.ApplyFreeSlipFreeConvexCorner(SqX_left ,SqY_back ,Corner::LeftBack);
        field.ApplyFreeSlipFreeConvexCorner(SqX_left ,SqY_front ,Corner::LeftFront);


      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }

        case 103:
    {
      auto apply_bc = [this](LBField& field) {
        // walls
        field.ApplyNoSlipFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNoSlipFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyNoSlipFreeWall(1, 0, LX-2, Wall::Front);
        field.ApplyNoSlipFreeWall(1, LY-1, LX-2, Wall::Back);
        // corners
        field.ApplyNoSlipFreeConcaveCorner(LX-1,LY-1,Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-1,0,Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(0,LY-1,Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(0,0,Corner::LeftFront);

        // walls
        field.ApplyNoSlipFreeWall(SqX_left, SqY_front+1, SqWallHeight, Wall::Right);
        field.ApplyNoSlipFreeWall(SqX_right, SqY_front+1 , SqWallHeight, Wall::Left);
        field.ApplyNoSlipFreeWall(SqX_left+1, SqY_back , SqWallLength, Wall::Front);
        field.ApplyNoSlipFreeWall(SqX_left+1,  SqY_front, SqWallLength, Wall::Back);
        // corners
        field.ApplyNoSlipFreeConvexCorner(SqX_right , SqY_back ,Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(SqX_right ,SqY_front ,Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(SqX_left ,SqY_back ,Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(SqX_left ,SqY_front ,Corner::LeftFront);


      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }

}}

void NematicFreeBoundary::BoundaryConditionsFields()
{
  switch(BC)
  {
    // pbc without bdry layer (nothing to do)
    case 101:
          {
      auto apply_bc = [this](ScalarField& field) {
        // walls
        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyNeumannFreeWall(1, 0, LX-2, Wall::Front);
        field.ApplyNeumannFreeWall(1, LY-1, LX-2, Wall::Back);
        // corners
        field.ApplyNeumannFreeConcaveCorner(LX-1,LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1,0, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(0,LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(0,0, Corner::LeftFront);

// walls
        field.ApplyNeumannFreeWall(SqX_left, SqY_front+1, SqWallHeight, Wall::Right);
        field.ApplyNeumannFreeWall(SqX_right, SqY_front+1 , SqWallHeight, Wall::Left);
        field.ApplyNeumannFreeWall(SqX_left+1, SqY_back , SqWallLength, Wall::Front);
        field.ApplyNeumannFreeWall(SqX_left+1,  SqY_front, SqWallLength, Wall::Back);
        // corners
        field.ApplyNeumannFreeConvexCorner(SqX_right , SqY_back ,Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(SqX_right ,SqY_front ,Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(SqX_left ,SqY_back ,Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(SqX_left ,SqY_front ,Corner::LeftFront);


      };
      apply_bc(QQxx);
      apply_bc(QQyx);
      break;
    }

        case 102:
          {
      auto apply_bc = [this](ScalarField& field) {
        // walls
        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyNeumannFreeWall(1, 0, LX-2, Wall::Front);
        field.ApplyNeumannFreeWall(1, LY-1, LX-2, Wall::Back);
        // corners
        field.ApplyNeumannFreeConcaveCorner(LX-1,LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1,0, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(0,LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(0,0, Corner::LeftFront);

// walls
        field.ApplyNeumannFreeWall(SqX_left, SqY_front+1, SqWallHeight, Wall::Right);
        field.ApplyNeumannFreeWall(SqX_right, SqY_front+1 , SqWallHeight, Wall::Left);
        field.ApplyNeumannFreeWall(SqX_left+1, SqY_back , SqWallLength, Wall::Front);
        field.ApplyNeumannFreeWall(SqX_left+1,  SqY_front, SqWallLength, Wall::Back);
        // corners
        field.ApplyNeumannFreeConvexCorner(SqX_right , SqY_back ,Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(SqX_right ,SqY_front ,Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(SqX_left ,SqY_back ,Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(SqX_left ,SqY_front ,Corner::LeftFront);


      };
      apply_bc(QQxx);
      apply_bc(QQyx);
      break;
    }

            case 103:
          {
      auto apply_bc = [this](ScalarField& field) {
        // walls
        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyNeumannFreeWall(1, 0, LX-2, Wall::Front);
        field.ApplyNeumannFreeWall(1, LY-1, LX-2, Wall::Back);
        // corners
        field.ApplyNeumannFreeConcaveCorner(LX-1,LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1,0, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(0,LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(0,0, Corner::LeftFront);

// walls
        field.ApplyNeumannFreeWall(SqX_left, SqY_front+1, SqWallHeight, Wall::Right);
        field.ApplyNeumannFreeWall(SqX_right, SqY_front+1 , SqWallHeight, Wall::Left);
        field.ApplyNeumannFreeWall(SqX_left+1, SqY_back , SqWallLength, Wall::Front);
        field.ApplyNeumannFreeWall(SqX_left+1,  SqY_front, SqWallLength, Wall::Back);
        // corners
        field.ApplyNeumannFreeConvexCorner(SqX_right , SqY_back ,Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(SqX_right ,SqY_front ,Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(SqX_left ,SqY_back ,Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(SqX_left ,SqY_front ,Corner::LeftFront);

      };
      apply_bc(QQxx);
      apply_bc(QQyx);

      break;
    }


    // von neumann box
    case 1:
    {
      auto apply_bc = [this](ScalarField& field) {
        // walls
        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyNeumannFreeWall(1, 0, LX-2, Wall::Front);
        field.ApplyNeumannFreeWall(1, LY-1, LX-2, Wall::Back);
        // corners
        field.ApplyNeumannFreeConcaveCorner(LX-1,LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1,0, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(0,LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(0,0, Corner::LeftFront);
      };
      apply_bc(QQxx);
      apply_bc(QQyx);

      break;
    }
  }
}

void NematicFreeBoundary::BoundaryConditionsFields2()
{
  switch(BC)
  {
    // pbc without bdry layer (nothing to do)
    case 0:
      break;

    case 101:
      {
      auto apply_bc = [this](ScalarField& field) {
        // walls
        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyNeumannFreeWall(1, 0, LX-2, Wall::Front);
        field.ApplyNeumannFreeWall(1, LY-1, LX-2, Wall::Back);
        // corners (don't matter, but this depends on implementation of derivative)
        field.ApplyNeumannFreeConcaveCorner(LX-1,LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1,0, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(0,LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(0,0, Corner::LeftFront);

        // walls
        field.ApplyNeumannFreeWall(SqX_left, SqY_front+1, SqWallHeight, Wall::Right);
        field.ApplyNeumannFreeWall(SqX_right, SqY_front+1 , SqWallHeight, Wall::Left);
        field.ApplyNeumannFreeWall(SqX_left+1, SqY_front , SqWallLength, Wall::Back);
        field.ApplyNeumannFreeWall(SqX_left+1,  SqY_back, SqWallLength, Wall::Front);
        // corners
        field.ApplyNeumannFreeConvexCorner(SqX_right , SqY_back ,Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(SqX_right ,SqY_front ,Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(SqX_left ,SqY_back ,Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(SqX_left ,SqY_front ,Corner::LeftFront);

      };
      auto apply_u_bc = [this](ScalarField& field) {
        // walls
        field.CopyDerivativeFreeWall(0, 1, LY-2, Wall::Left);
        field.CopyDerivativeFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.CopyDerivativeFreeWall(1, 0, LX-2, Wall::Front);
        field.CopyDerivativeFreeWall(1, LY-1, LX-2, Wall::Back);
        // corners (don't matter, but this depends on implementation of derivative)
        field.CopyDerivativeFreeConcaveCorner(LX-1,LY-1, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-1,0, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(0,LY-1, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(0,0, Corner::LeftFront);
        // walls
        field.CopyDerivativeFreeWall(SqX_left, SqY_front+1, SqWallHeight, Wall::Right);
        field.CopyDerivativeFreeWall(SqX_right, SqY_front+1 , SqWallHeight, Wall::Left);
        field.CopyDerivativeFreeWall(SqX_left+1, SqY_front , SqWallLength, Wall::Back);
        field.CopyDerivativeFreeWall(SqX_left+1,  SqY_back, SqWallLength, Wall::Front);
        // corners
        field.CopyDerivativeFreeConvexCorner(SqX_right , SqY_back ,Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(SqX_right ,SqY_front ,Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(SqX_left ,SqY_back ,Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(SqX_left ,SqY_front ,Corner::LeftFront);

      };

      apply_u_bc(ux);
      apply_u_bc(uy);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }

     case 102:
    {
      auto apply_bc = [this](ScalarField& field) {
               // walls
        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyNeumannFreeWall(1, 0, LX-2, Wall::Front);
        field.ApplyNeumannFreeWall(1, LY-1, LX-2, Wall::Back);
        // corners (don't matter, but this depends on implementation of derivative)
        field.ApplyNeumannFreeConcaveCorner(LX-1,LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1,0, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(0,LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(0,0, Corner::LeftFront);

        // walls
        field.ApplyNeumannFreeWall(SqX_left, SqY_front+1, SqWallHeight, Wall::Right);
        field.ApplyNeumannFreeWall(SqX_right, SqY_front+1 , SqWallHeight, Wall::Left);
        field.ApplyNeumannFreeWall(SqX_left+1, SqY_front , SqWallLength, Wall::Back);
        field.ApplyNeumannFreeWall(SqX_left+1,  SqY_back, SqWallLength, Wall::Front);
        // corners
        field.ApplyNeumannFreeConvexCorner(SqX_right , SqY_back ,Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(SqX_right ,SqY_front ,Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(SqX_left ,SqY_back ,Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(SqX_left ,SqY_front ,Corner::LeftFront);


      };

//NEED NEW FUNCTIONS FOR UX AND UY INDIVIDUALLY
      auto apply_u_bc = [this](ScalarField& field) {
        // walls
        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.CopyDerivativeFreeWall(1, 0, LX-2, Wall::Front);
        field.CopyDerivativeFreeWall(1, LY-1, LX-2, Wall::Back);

        field.CopyDerivativeFreeConcaveCorner(LX-1,LY-1, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-1,0, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(0,LY-1, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(0,0, Corner::LeftFront);
        // walls
        field.ApplyNeumannFreeWall(SqX_left, SqY_front+1, SqWallHeight, Wall::Right);
        field.ApplyNeumannFreeWall(SqX_right, SqY_front+1 , SqWallHeight, Wall::Left);
        field.CopyDerivativeFreeWall(SqX_left+1, SqY_front , SqWallLength, Wall::Back);
        field.CopyDerivativeFreeWall(SqX_left+1,  SqY_back, SqWallLength, Wall::Front);
        // corners
        field.CopyDerivativeFreeConvexCorner(SqX_right , SqY_back ,Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(SqX_right ,SqY_front ,Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(SqX_left ,SqY_back ,Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(SqX_left ,SqY_front ,Corner::LeftFront);
      };

      apply_u_bc(ux);
      apply_u_bc(uy);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }


     case 103:
         {
      auto apply_bc = [this](ScalarField& field) {
        // walls
        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyNeumannFreeWall(1, 0, LX-2, Wall::Front);
        field.ApplyNeumannFreeWall(1, LY-1, LX-2, Wall::Back);
        // corners (don't matter, but this depends on implementation of derivative)
        field.ApplyNeumannFreeConcaveCorner(LX-1,LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1,0, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(0,LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(0,0, Corner::LeftFront);

        // walls
        field.ApplyNeumannFreeWall(SqX_left, SqY_front+1, SqWallHeight, Wall::Right);
        field.ApplyNeumannFreeWall(SqX_right, SqY_front+1 , SqWallHeight, Wall::Left);
        field.ApplyNeumannFreeWall(SqX_left+1, SqY_front , SqWallLength, Wall::Back);
        field.ApplyNeumannFreeWall(SqX_left+1,  SqY_back, SqWallLength, Wall::Front);
        // corners
        field.ApplyNeumannFreeConvexCorner(SqX_right , SqY_back ,Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(SqX_right ,SqY_front ,Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(SqX_left ,SqY_back ,Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(SqX_left ,SqY_front ,Corner::LeftFront);

      };


      auto apply_u_bc = [this](ScalarField& field) {
        // walls
        field.CopyDerivativeFreeWall(0, 1, LY-2, Wall::Left);
        field.CopyDerivativeFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.CopyDerivativeFreeWall(1, 0, LX-2, Wall::Front);
        field.CopyDerivativeFreeWall(1, LY-1, LX-2, Wall::Back);
        // corners (don't matter, but this depends on implementation of derivative)
        field.CopyDerivativeFreeConcaveCorner(LX-1,LY-1, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-1,0, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(0,LY-1, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(0,0, Corner::LeftFront);
        // walls
        field.CopyDerivativeFreeWall(SqX_left, SqY_front+1, SqWallHeight, Wall::Right);
        field.CopyDerivativeFreeWall(SqX_right, SqY_front+1 , SqWallHeight, Wall::Left);
        field.CopyDerivativeFreeWall(SqX_left+1, SqY_front , SqWallLength, Wall::Back);
        field.CopyDerivativeFreeWall(SqX_left+1,  SqY_back, SqWallLength, Wall::Front);
        // corners
        field.CopyDerivativeFreeConvexCorner(SqX_right , SqY_back ,Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(SqX_right ,SqY_front ,Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(SqX_left ,SqY_back ,Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(SqX_left ,SqY_front ,Corner::LeftFront);

      };

      apply_u_bc(ux);
      apply_u_bc(uy);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }
    // box with Neumann
    case 1:
    {
      auto apply_bc = [this](ScalarField& field) {
        // walls
        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyNeumannFreeWall(1, 0, LX-2, Wall::Front);
        field.ApplyNeumannFreeWall(1, LY-1, LX-2, Wall::Back);
        // corners (don't matter, but this depends on implementation of derivative)
        field.ApplyNeumannFreeConcaveCorner(LX-1,LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1,0, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(0,LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(0,0, Corner::LeftFront);
      };

      auto apply_u_bc = [this](ScalarField& field) {
        // walls
        field.CopyDerivativeFreeWall(0, 1, LY-2, Wall::Left);
        field.CopyDerivativeFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.CopyDerivativeFreeWall(1, 0, LX-2, Wall::Front);
        field.CopyDerivativeFreeWall(1, LY-1, LX-2, Wall::Back);
        // corners (don't matter, but this depends on implementation of derivative)
        field.CopyDerivativeFreeConcaveCorner(LX-1,LY-1, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-1,0, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(0,LY-1, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(0,0, Corner::LeftFront);
      };

      apply_u_bc(ux);
      apply_u_bc(uy);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }

  }
}

void NematicFreeBoundary::RuntimeChecks()
{
  //check derivatives of the fields at the boundaries
#ifdef DEBUG
  {
    BoundaryConditionsFields();
    BoundaryConditionsFields2();

    array<double,4> test = {0., 0., 0., 0.};

    //walls
    switch(BC)
    {
      case 1:
      {
        test[0]=QQxx.CheckNeumannFreeWall(0, 1, LY-2, Wall::Left);
        test[1]=QQxx.CheckNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        test[2]=QQxx.CheckNeumannFreeWall(1, 0, LX-2, Wall::Front);
        test[3]=QQxx.CheckNeumannFreeWall(1, LY-1, LX-2, Wall::Back);

        std::cout << "\nCheck of bc for field (incl corners):\nLeft: "
                  << test[0] << ", Right: " << test[1] << ", Front: "
                  << test[2] << ", Back: " << test[3];
      break;
      }    }
  }
#endif //DEBUG

  // check that the sum of f is constant
  {
    double fcheck = 0;
    for(unsigned k=0; k<DomainSize; ++k)
      if(!outside[k]) fcheck = accumulate(begin(ff[k]), end(ff[k]), fcheck);
    std::cout << "\nftot=" << ftot << ", difference: " << ftot-fcheck;
    if(abs(ftot-fcheck)>1)
      throw error_msg("f is not conserved (", ftot, "/", fcheck, ")");
  }

  cout << "\n";
}
