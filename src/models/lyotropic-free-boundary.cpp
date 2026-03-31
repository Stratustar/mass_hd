#include "header.hpp"
#include "models/lyotropic-free-boundary.hpp"
#include "error_msg.hpp"
#include "random.hpp"
#include "lb.hpp"
#include "tools.hpp"

using namespace std;
namespace opt = boost::program_options;

// from main.cpp:
extern unsigned nthreads, nsubsteps, verbose;

LyotropicFreeBoundary::LyotropicFreeBoundary(unsigned LX, unsigned LY, unsigned BC_)
{
  BC=BC_;

  //In certain cases (when we need additional nodes), the grid has to be reset.
  //This also redefines the neighbour list.
  if (BC==7 or BC==8 or BC==11 or BC==12)
  {
    SetSize(LX, LY, GridType::Custom1);
  }
  else if (BC==9 or BC==10)
  {
    SetSize(LX, LY, GridType::Custom2);
  }
  else if (BC==13 or BC==14)
  {
    SetSize(LX, LY, GridType::Custom3);
  }
  else if (BC==15 or BC==16 or BC==17 or BC==18 or BC==19)
  {
    SetSize(LX, LY, GridType::Custom4);
  }
  else if (BC==20)
  {
    SetSize(LX, LY, GridType::Chimney);
  }
  else if (BC==26)
  {
    SetSize(LX, LY, GridType::Chimney4);
  }
  else if (BC==27)
  {
    SetSize(LX, LY, GridType::Chimney7);
  }
  else if (BC==28)
  {
    SetSize(LX, LY, GridType::Chimney10);
  }
  else if (BC==29)
  {
    SetSize(LX, LY, GridType::Chimney13);
  }
  else if (BC==30)
  {
    SetSize(LX, LY, GridType::Chimney16);
  }
  else if (BC==31)
  {
    SetSize(LX, LY, GridType::ChimneyExit);
  }
  else if (BC==32)
  {
    SetSize(LX, LY, GridType::ChimneyExit4);
  }
  else if (BC==33)
  {
    SetSize(LX, LY, GridType::ChimneyExit7);
  }
  else if (BC==34)
  {
    SetSize(LX, LY, GridType::ChimneyExit10);
  }
  else if (BC==35)
  {
    SetSize(LX, LY, GridType::ChimneyExit13);
  }
  else if (BC==36)
  {
    SetSize(LX, LY, GridType::ChimneyExit16);
  }
  else if (BC==101 or BC==102 or BC==103)
  {
    SetSize(LX, LY, GridType::SquareAn);        //Should be Layer?
  }
  else
  {
    SetSize(LX, LY, GridType::Periodic);
  }

  outside.resize(DomainSize);
}

void LyotropicFreeBoundary::ConfigureBoundaries()
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
    else if(BC==1 or BC==4 or BC==21 or BC==22 or BC==23 or BC==24 or BC==25)
      outside[k] = ( x<=0 or x>=LX-1 or y<=0 or y>=LY-1 );
    //BC 101 and 102 are meant to be FreeSlip. 103 is No Slip.
    //BC 101 is not meant to work in theory
    else if (BC==101 or BC==102 or BC==103)
    {
      //outside[k] = ( x<=0 or x>=LX-1 or y<=0 or y>=LY-1 or (x>=SqX_left and x<= SqX_right and y<= SqY_back and y>= SqY_front));
      if (x==0 or x== LX-1 )
        {outside[k] = true;}
      else if ( y==0 or y==LY-1)
        {outside[k] = true;}
      else if (x>=(SqX_left ) and x<= (SqX_right) and y<= (SqY_back) and y>= (SqY_front))
        {outside[k] = true;}
      else
        {outside[k]=false;}

    }
    // channel with free- or no-slip and von Neumann
    else if (BC==2 or BC==5)
      outside[k] = ( y==0 or y==LY-1 );
    // channel in y-direction with free- or no-slip and von Neumann
    else if (BC==201 or BC==501)
      outside[k] = ( x==0 or x==LX-1 );
    // 'inverted' box with free- or no-slip and von Neumann, merely for testing
    else if (BC==3 or BC==6)
      outside[k] = not ( x<10 or x>40 or y<10 or y>30 );
    // channel that tightens and widens again with free- or no-slip and von Neumann
    else if (BC==7 or BC==8)
    {
      outside[k]= x==0 or x==LX-1 or y==0 or y==LY-1 or
                  (o<=y and y<=o2-1 and (x<=d or x>=LX-d-1));
    }
    // channel that widens with free- or no-slip and von Neumann.
    // In principle as BC==7/BC==8 but without the lower reservoir.
    else if (BC==9 or BC==10)
    {
      outside[k]= x==0 or x==LX-1 or y==0 or y==LY-1 or
                  (y<=l and (x<=d or x>=LX-d-1));
    }
    // channel that tightens and widens again with free- or no-slip and von Neumann
    //, but periodic at the reservoirs
    else if (BC==11 or BC==12)
    {
      outside[k]= y==0 or y==LY-1 or
                  (o<=y and y<=o2-1 and (x<=d or x>=LX-d-1));
    }
    // channel that tightens and widens again with free- or no-slip and von Neumann
    // and with rounded corners
    else if (BC==13 or BC==14)
    {
      if (x==0 or x==LX-1 or y==0 or y==LY-1)
      {
        outside[k]=true;
      }
      else if ( (x<=d-5 or x>=LX-d+4) and (y==o or y==o2-1) )
      {
        outside[k]=true;
      }
      else if ( (x<=d-4 or x>=LX-d+3) and (y==o+1 or y==o2-2) )
      {
        outside[k]=true;
      }
      else if ( (x<=d-2 or x>=LX-d+1) and (y==o+2 or y==o+3 or y==o2-3 or y==o2-4) )
      {
        outside[k]=true;
      }
      else if ( (x<=d-1 or x>=LX-d  ) and (y==o+4 or y==o2-5) )
      {
        outside[k]=true;
      }
      else if ( (x<=d or x>=LX-d-1  ) and y>=o+5 and y<=o2-6 )
      {
        outside[k]=true;
      }
    }
    else if (BC==15 or BC==16)
    {
      //automatically generated code
      if (x<=0 and y>=1 and y<=1+o-1)
      {
        outside[k]=true;
      }
      if (x<=d and y>=o+16 and y<=o+16+l-32)
      {
        outside[k]=true;
      }
      if (x<=0 and y>=o2 and y<=o2+r-1)
      {
        outside[k]=true;
      }
      if (x>=LX-1 and y>=1 and y<=1+o-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1 and y>=o+16 and y<=o+16+l-32)
      {
        outside[k]=true;
      }
      if (x>=LX-1 and y>=o2 and y<=o2+r-1)
      {
        outside[k]=true;
      }
      if (y==0 and x>=1 and x<=1+LX-2)
      {
        outside[k]=true;
      }
      if (y==o2-1 and x>=1 and x<=1+d-14)
      {
        outside[k]=true;
      }
      if (y==o2-1 and x>=LX-d+15 and x<=LX-d+15+d-14)
      {
        outside[k]=true;
      }
      if (y==LY-1 and x>=1 and x<=1+LX-2)
      {
        outside[k]=true;
      }
      if (y==o and x>=1 and x<=1+d-14)
      {
        outside[k]=true;
      }
      if (y==o and x>=LX-d+15 and x<=LX-d+15+d-14)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o+8)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+10)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+11)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+13)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+14)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+15)
      {
        outside[k]=true;
      }
      if (x<=d-12 and y==o+0)
      {
        outside[k]=true;
      }
      if (x<=d-9 and y==o+1)
      {
        outside[k]=true;
      }
      if (x<=d-7 and y==o+2)
      {
        outside[k]=true;
      }
      if (x<=d-6 and y==o+3)
      {
        outside[k]=true;
      }
      if (x<=d-5 and y==o+4)
      {
        outside[k]=true;
      }
      if (x<=d-4 and y==o+5)
      {
        outside[k]=true;
      }
      if (x<=d-3 and y==o+6)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o+7)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+9)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+12)
      {
        outside[k]=true;
      }
      if (x<=d-12 and y==o+1)
      {
        outside[k]=true;
      }
      if (x<=d-9 and y==o+2)
      {
        outside[k]=true;
      }
      if (x<=d-7 and y==o+3)
      {
        outside[k]=true;
      }
      if (x<=d-6 and y==o+4)
      {
        outside[k]=true;
      }
      if (x<=d-5 and y==o+5)
      {
        outside[k]=true;
      }
      if (x<=d-4 and y==o+6)
      {
        outside[k]=true;
      }
      if (x<=d-3 and y==o+7)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o+9)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+12)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+15)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+14)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+13)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+11)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+10)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o+8)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+12)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+9)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o+7)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+3 and y==o+6)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+4 and y==o+5)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+5 and y==o+4)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+6 and y==o+3)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+7 and y==o+2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+9 and y==o+1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+12 and y==o+0)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+12)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o+9)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+3 and y==o+7)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+4 and y==o+6)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+5 and y==o+5)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+6 and y==o+4)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+7 and y==o+3)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+9 and y==o+2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+12 and y==o+1)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o2-1-15)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o2-1-14)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o2-1-13)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o2-1-11)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o2-1-10)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o2-1-8)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o2-1-12)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o2-1-9)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o2-1-7)
      {
        outside[k]=true;
      }
      if (x<=d-3 and y==o2-1-6)
      {
        outside[k]=true;
      }
      if (x<=d-4 and y==o2-1-5)
      {
        outside[k]=true;
      }
      if (x<=d-5 and y==o2-1-4)
      {
        outside[k]=true;
      }
      if (x<=d-6 and y==o2-1-3)
      {
        outside[k]=true;
      }
      if (x<=d-7 and y==o2-1-2)
      {
        outside[k]=true;
      }
      if (x<=d-9 and y==o2-1-1)
      {
        outside[k]=true;
      }
      if (x<=d-12 and y==o2-1+0)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o2-1-12)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o2-1-9)
      {
        outside[k]=true;
      }
      if (x<=d-3 and y==o2-1-7)
      {
        outside[k]=true;
      }
      if (x<=d-4 and y==o2-1-6)
      {
        outside[k]=true;
      }
      if (x<=d-5 and y==o2-1-5)
      {
        outside[k]=true;
      }
      if (x<=d-6 and y==o2-1-4)
      {
        outside[k]=true;
      }
      if (x<=d-7 and y==o2-1-3)
      {
        outside[k]=true;
      }
      if (x<=d-9 and y==o2-1-2)
      {
        outside[k]=true;
      }
      if (x<=d-12 and y==o2-1-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o2-1-8)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o2-1-10)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o2-1-11)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o2-1-13)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o2-1-14)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o2-1-15)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+12 and y==o2-1+0)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+9 and y==o2-1-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+7 and y==o2-1-2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+6 and y==o2-1-3)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+5 and y==o2-1-4)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+4 and y==o2-1-5)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+3 and y==o2-1-6)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o2-1-7)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o2-1-9)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o2-1-12)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+12 and y==o2-1-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+9 and y==o2-1-2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+7 and y==o2-1-3)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+6 and y==o2-1-4)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+5 and y==o2-1-5)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+4 and y==o2-1-6)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+3 and y==o2-1-7)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o2-1-9)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o2-1-12)
      {
        outside[k]=true;
      }
      if (x<=0 and y==0)
      {
        outside[k]=true;
      }
      if (x<=0 and y==o2-1)
      {
        outside[k]=true;
      }
      if (x<=0 and y==o)
      {
        outside[k]=true;
      }
      if (x<=0 and y==LY-1)
      {
        outside[k]=true;
      }
      if (x>=LX-1 and y==LY-1)
      {
        outside[k]=true;
      }
      if (x>=LX-1 and y==o)
      {
        outside[k]=true;
      }
      if (x>=LX-1 and y==0)
      {
        outside[k]=true;
      }
      if (x>=LX-1 and y==o2-1)
      {
        outside[k]=true;
      }
    }
    else if (BC==17 or BC==18 or BC==19)
    {
      //automatically generated code
      if (x<=d and y>=o+16 and y<=o+16+l-32)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1 and y>=o+16 and y<=o+16+l-32)
      {
        outside[k]=true;
      }
      if (y==0 and x<=0+LX)
      {
        outside[k]=true;
      }
      if (y==o2-1 and x<=0+d-15)
      {
        outside[k]=true;
      }
      if (y==o2-1 and x>=LX-d+15 and x<=LX-d+15+d-15)
      {
        outside[k]=true;
      }
      if (y==LY-1 and x<=0+LX)
      {
        outside[k]=true;
      }
      if (y==o and x<=0+d-15)
      {
        outside[k]=true;
      }
      if (y==o and x>=LX-d+15 and x<=LX-d+15+d-15)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o+8)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+10)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+11)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+13)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+14)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+15)
      {
        outside[k]=true;
      }
      if (x<=d-12 and y==o+0)
      {
        outside[k]=true;
      }
      if (x<=d-9 and y==o+1)
      {
        outside[k]=true;
      }
      if (x<=d-7 and y==o+2)
      {
        outside[k]=true;
      }
      if (x<=d-6 and y==o+3)
      {
        outside[k]=true;
      }
      if (x<=d-5 and y==o+4)
      {
        outside[k]=true;
      }
      if (x<=d-4 and y==o+5)
      {
        outside[k]=true;
      }
      if (x<=d-3 and y==o+6)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o+7)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+9)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+12)
      {
        outside[k]=true;
      }
      if (x<=d-12 and y==o+1)
      {
        outside[k]=true;
      }
      if (x<=d-9 and y==o+2)
      {
        outside[k]=true;
      }
      if (x<=d-7 and y==o+3)
      {
        outside[k]=true;
      }
      if (x<=d-6 and y==o+4)
      {
        outside[k]=true;
      }
      if (x<=d-5 and y==o+5)
      {
        outside[k]=true;
      }
      if (x<=d-4 and y==o+6)
      {
        outside[k]=true;
      }
      if (x<=d-3 and y==o+7)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o+9)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+12)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+15)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+14)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+13)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+11)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+10)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o+8)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+12)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+9)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o+7)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+3 and y==o+6)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+4 and y==o+5)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+5 and y==o+4)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+6 and y==o+3)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+7 and y==o+2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+9 and y==o+1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+12 and y==o+0)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+12)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o+9)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+3 and y==o+7)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+4 and y==o+6)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+5 and y==o+5)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+6 and y==o+4)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+7 and y==o+3)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+9 and y==o+2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+12 and y==o+1)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o2-1-15)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o2-1-14)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o2-1-13)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o2-1-11)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o2-1-10)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o2-1-8)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o2-1-12)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o2-1-9)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o2-1-7)
      {
        outside[k]=true;
      }
      if (x<=d-3 and y==o2-1-6)
      {
        outside[k]=true;
      }
      if (x<=d-4 and y==o2-1-5)
      {
        outside[k]=true;
      }
      if (x<=d-5 and y==o2-1-4)
      {
        outside[k]=true;
      }
      if (x<=d-6 and y==o2-1-3)
      {
        outside[k]=true;
      }
      if (x<=d-7 and y==o2-1-2)
      {
        outside[k]=true;
      }
      if (x<=d-9 and y==o2-1-1)
      {
        outside[k]=true;
      }
      if (x<=d-12 and y==o2-1+0)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o2-1-12)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o2-1-9)
      {
        outside[k]=true;
      }
      if (x<=d-3 and y==o2-1-7)
      {
        outside[k]=true;
      }
      if (x<=d-4 and y==o2-1-6)
      {
        outside[k]=true;
      }
      if (x<=d-5 and y==o2-1-5)
      {
        outside[k]=true;
      }
      if (x<=d-6 and y==o2-1-4)
      {
        outside[k]=true;
      }
      if (x<=d-7 and y==o2-1-3)
      {
        outside[k]=true;
      }
      if (x<=d-9 and y==o2-1-2)
      {
        outside[k]=true;
      }
      if (x<=d-12 and y==o2-1-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o2-1-8)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o2-1-10)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o2-1-11)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o2-1-13)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o2-1-14)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o2-1-15)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+12 and y==o2-1+0)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+9 and y==o2-1-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+7 and y==o2-1-2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+6 and y==o2-1-3)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+5 and y==o2-1-4)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+4 and y==o2-1-5)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+3 and y==o2-1-6)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o2-1-7)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o2-1-9)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o2-1-12)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+12 and y==o2-1-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+9 and y==o2-1-2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+7 and y==o2-1-3)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+6 and y==o2-1-4)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+5 and y==o2-1-5)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+4 and y==o2-1-6)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+3 and y==o2-1-7)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o2-1-9)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o2-1-12)
      {
        outside[k]=true;
      }
    }
    // reservoir with a single chimney
    else if (BC==20)
    {
      //automatically generated code
      if (x<=d and y>=o+1 and y<=o+1+LY-o-2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1 and y>=o+1 and y<=o+1+LY-o-2)
      {
        outside[k]=true;
      }
      if (y==0 and x<=0+LX)
      {
        outside[k]=true;
      }
      if (y==LY-1 and x>=d+1 and x<=d+1+LX-2*d-2)
      {
        outside[k]=true;
      }
      if (y==o and x>=LX-d and x<=LX-d+d)
      {
        outside[k]=true;
      }
      if (y==o and x<=0+d)
      {
        outside[k]=true;
      }
      if (x<=d and y==LY-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1 and y==LY-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1 and y==o)
      {
        outside[k]=true;
      }
      if (x<=d and y==o)
      {
        outside[k]=true;
      }
    }
    //reservoir with chimney, radius corners=4
    else if (BC==26)
    {
      //automatically generated code
      if (x<=d and y>=o+5 and y<=o+5+LY-o-6)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1 and y>=o+5 and y<=o+5+LY-o-6)
      {
        outside[k]=true;
      }
      if (y==0 and x<=0+LX)
      {
        outside[k]=true;
      }
      if (y==LY-1 and x>=d+1 and x<=d+1+LX-2*d-2)
      {
        outside[k]=true;
      }
      if (x<=d and y==LY-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1 and y==LY-1)
      {
        outside[k]=true;
      }
      if (y==o and x>=LX-d+4 and x<=LX-d+4+d-4)
      {
        outside[k]=true;
      }
      if (y==o and x<=0+d-4)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+4)
      {
        outside[k]=true;
      }
      if (x<=d-3 and y==o+0)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o+1)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+2)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+3)
      {
        outside[k]=true;
      }
      if (x<=d-3 and y==o+1)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o+2)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+3)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+4)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+3)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o+1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+3 and y==o+0)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+3)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o+2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+3 and y==o+1)
      {
        outside[k]=true;
      }
    }
    //reservoir with chimney, radius corners=7
    else if (BC==27)
    {
      //automatically generated code
      if (x<=d and y>=o+8 and y<=o+8+LY-o-9)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1 and y>=o+8 and y<=o+8+LY-o-9)
      {
        outside[k]=true;
      }
      if (y==0 and x<=0+LX)
      {
        outside[k]=true;
      }
      if (y==LY-1 and x>=d+1 and x<=d+1+LX-2*d-2)
      {
        outside[k]=true;
      }
      if (x<=d and y==LY-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1 and y==LY-1)
      {
        outside[k]=true;
      }
      if (y==o and x>=LX-d+7 and x<=LX-d+7+d-7)
      {
        outside[k]=true;
      }
      if (y==o and x<=0+d-7)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+4)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+6)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+7)
      {
        outside[k]=true;
      }
      if (x<=d-5 and y==o+0)
      {
        outside[k]=true;
      }
      if (x<=d-3 and y==o+1)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o+2)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+3)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+5)
      {
        outside[k]=true;
      }
      if (x<=d-5 and y==o+1)
      {
        outside[k]=true;
      }
      if (x<=d-3 and y==o+2)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o+3)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+5)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+7)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+6)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+4)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+5)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+3)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o+2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+3 and y==o+1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+5 and y==o+0)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+5)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o+3)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+3 and y==o+2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+5 and y==o+1)
      {
        outside[k]=true;
      }
    }
    //reservoir with chimney, radius corners=10
    else if (BC==28)
    {
      //automatically generated code
      if (x<=d and y>=o+11 and y<=o+11+LY-o-12)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1 and y>=o+11 and y<=o+11+LY-o-12)
      {
        outside[k]=true;
      }
      if (y==0 and x<=0+LX)
      {
        outside[k]=true;
      }
      if (y==LY-1 and x>=d+1 and x<=d+1+LX-2*d-2)
      {
        outside[k]=true;
      }
      if (x<=d and y==LY-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1 and y==LY-1)
      {
        outside[k]=true;
      }
      if (y==o and x>=LX-d+10 and x<=LX-d+10+d-10)
      {
        outside[k]=true;
      }
      if (y==o and x<=0+d-10)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+6)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+8)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+9)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+10)
      {
        outside[k]=true;
      }
      if (x<=d-7 and y==o+0)
      {
        outside[k]=true;
      }
      if (x<=d-5 and y==o+1)
      {
        outside[k]=true;
      }
      if (x<=d-4 and y==o+2)
      {
        outside[k]=true;
      }
      if (x<=d-3 and y==o+3)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o+4)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+5)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+7)
      {
        outside[k]=true;
      }
      if (x<=d-7 and y==o+1)
      {
        outside[k]=true;
      }
      if (x<=d-5 and y==o+2)
      {
        outside[k]=true;
      }
      if (x<=d-4 and y==o+3)
      {
        outside[k]=true;
      }
      if (x<=d-3 and y==o+4)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o+5)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+7)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+10)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+9)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+8)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+6)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+7)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+5)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o+4)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+3 and y==o+3)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+4 and y==o+2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+5 and y==o+1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+7 and y==o+0)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+7)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o+5)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+3 and y==o+4)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+4 and y==o+3)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+5 and y==o+2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+7 and y==o+1)
      {
        outside[k]=true;
      }
    }
    //reservoir with chimney, radius corners=13
    else if (BC==29)
    {
      //automatically generated code
      if (x<=d and y>=o+14 and y<=o+14+LY-o-15)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1 and y>=o+14 and y<=o+14+LY-o-15)
      {
        outside[k]=true;
      }
      if (y==0 and x<=0+LX)
      {
        outside[k]=true;
      }
      if (y==LY-1 and x>=d+1 and x<=d+1+LX-2*d-2)
      {
        outside[k]=true;
      }
      if (x<=d and y==LY-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1 and y==LY-1)
      {
        outside[k]=true;
      }
      if (y==o and x>=LX-d+13 and x<=LX-d+13+d-13)
      {
        outside[k]=true;
      }
      if (y==o and x<=0+d-13)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+8)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+9)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+11)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+12)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+13)
      {
        outside[k]=true;
      }
      if (x<=d-10 and y==o+0)
      {
        outside[k]=true;
      }
      if (x<=d-7 and y==o+1)
      {
        outside[k]=true;
      }
      if (x<=d-6 and y==o+2)
      {
        outside[k]=true;
      }
      if (x<=d-5 and y==o+3)
      {
        outside[k]=true;
      }
      if (x<=d-4 and y==o+4)
      {
        outside[k]=true;
      }
      if (x<=d-3 and y==o+5)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o+6)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+7)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+10)
      {
        outside[k]=true;
      }
      if (x<=d-10 and y==o+1)
      {
        outside[k]=true;
      }
      if (x<=d-7 and y==o+2)
      {
        outside[k]=true;
      }
      if (x<=d-6 and y==o+3)
      {
        outside[k]=true;
      }
      if (x<=d-5 and y==o+4)
      {
        outside[k]=true;
      }
      if (x<=d-4 and y==o+5)
      {
        outside[k]=true;
      }
      if (x<=d-3 and y==o+6)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o+7)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+10)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+13)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+12)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+11)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+9)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+8)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+10)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+7)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o+6)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+3 and y==o+5)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+4 and y==o+4)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+5 and y==o+3)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+6 and y==o+2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+7 and y==o+1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+10 and y==o+0)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+10)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o+7)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+3 and y==o+6)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+4 and y==o+5)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+5 and y==o+4)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+6 and y==o+3)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+7 and y==o+2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+10 and y==o+1)
      {
        outside[k]=true;
      }

    }
    //reservoir with chimney, radius corners=16
    else if (BC==30)
    {
      //automatically generated code
      if (x<=d and y>=o+17 and y<=o+17+LY-o-18)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1 and y>=o+17 and y<=o+17+LY-o-18)
      {
        outside[k]=true;
      }
      if (y==0 and x<=0+LX)
      {
        outside[k]=true;
      }
      if (y==LY-1 and x>=d+1 and x<=d+1+LX-2*d-2)
      {
        outside[k]=true;
      }
      if (x<=d and y==LY-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1 and y==LY-1)
      {
        outside[k]=true;
      }
      if (y==o and x>=LX-d+16 and x<=LX-d+16+d-16)
      {
        outside[k]=true;
      }
      if (y==o and x<=0+d-16)
      {
        outside[k]=true;
      }
      if (x<=d-4 and y==o+6)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o+9)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+11)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+12)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+14)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+15)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+16)
      {
        outside[k]=true;
      }
      if (x<=d-13 and y==o+0)
      {
        outside[k]=true;
      }
      if (x<=d-10 and y==o+1)
      {
        outside[k]=true;
      }
      if (x<=d-8 and y==o+2)
      {
        outside[k]=true;
      }
      if (x<=d-7 and y==o+3)
      {
        outside[k]=true;
      }
      if (x<=d-5 and y==o+4)
      {
        outside[k]=true;
      }
      if (x<=d-4 and y==o+5)
      {
        outside[k]=true;
      }
      if (x<=d-3 and y==o+7)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o+8)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+10)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+13)
      {
        outside[k]=true;
      }
      if (x<=d-13 and y==o+1)
      {
        outside[k]=true;
      }
      if (x<=d-10 and y==o+2)
      {
        outside[k]=true;
      }
      if (x<=d-8 and y==o+3)
      {
        outside[k]=true;
      }
      if (x<=d-7 and y==o+4)
      {
        outside[k]=true;
      }
      if (x<=d-5 and y==o+5)
      {
        outside[k]=true;
      }
      if (x<=d-4 and y==o+7)
      {
        outside[k]=true;
      }
      if (x<=d-3 and y==o+8)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o+10)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+13)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+16)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+15)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+14)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+12)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+11)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o+9)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+4 and y==o+6)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+13)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+10)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o+8)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+3 and y==o+7)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+4 and y==o+5)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+5 and y==o+4)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+7 and y==o+3)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+8 and y==o+2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+10 and y==o+1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+13 and y==o+0)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+13)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o+10)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+3 and y==o+8)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+4 and y==o+7)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+5 and y==o+5)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+7 and y==o+4)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+8 and y==o+3)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+10 and y==o+2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+13 and y==o+1)
      {
        outside[k]=true;
      }
    }
    //reservoir with chimney exit
    else if (BC==31)
    {
      //automatically generated code
      if (x>=LX-d-1 and y==0)
      {
        outside[k]=true;
      }
      if (x<=d and y==0)
      {
        outside[k]=true;
      }
      if (y==0 and x>=d+1 and x<=d+1+LX-2*d-2)
      {
        outside[k]=true;
      }
      if (x<=d and y>=1 and y<=1+o)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1 and y>=1 and y<=1+o)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1 and y==o+1)
      {
        outside[k]=true;
      }
      if (x<=d and y==o+1)
      {
        outside[k]=true;
      }
      if (y==o+1 and x>=LX-d and x<=LX-d+d)
      {
        outside[k]=true;
      }
      if (y==o+1 and x<=0+d)
      {
        outside[k]=true;
      }
      if (y==LY-1 and x<=0+LX)
      {
        outside[k]=true;
      }
    }
    //reservoir with chimney exit, radius corners=4
    else if (BC==32)
    {
      //automatically generated code
      if (x<=d+0 and y==o+1-4)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+1-3)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+1-2)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o+1-1)
      {
        outside[k]=true;
      }
      if (x<=d-3 and y==o+1+0)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+1-3)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o+1-2)
      {
        outside[k]=true;
      }
      if (x<=d-3 and y==o+1-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+1-4)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+3 and y==o+1+0)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o+1-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+1-2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+1-3)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+3 and y==o+1-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o+1-2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+1-3)
      {
        outside[k]=true;
      }
      if (y==0 and x>=d+1 and x<=d+1+LX-2*d-2)
      {
        outside[k]=true;
      }
      if (x<=d and y>=1 and y<=1+o-4)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1 and y>=1 and y<=1+o-4)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1 and y==0)
      {
        outside[k]=true;
      }
      if (x<=d and y==0)
      {
        outside[k]=true;
      }
      if (y==o+1 and x>=LX-d+4 and x<=LX-d+4+d-4)
      {
        outside[k]=true;
      }
      if (y==o+1 and x<=0+d-4)
      {
        outside[k]=true;
      }
      if (y==LY-1 and x<=0+LX)
      {
        outside[k]=true;
      }
    }
    //reservoir with chimney exit, radius corners=7
    else if (BC==33)
    {
      //automatically generated code
      if (x<=d+0 and y==o+1-7)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+1-6)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+1-4)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+1-5)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+1-3)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o+1-2)
      {
        outside[k]=true;
      }
      if (x<=d-3 and y==o+1-1)
      {
        outside[k]=true;
      }
      if (x<=d-5 and y==o+1+0)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+1-5)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o+1-3)
      {
        outside[k]=true;
      }
      if (x<=d-3 and y==o+1-2)
      {
        outside[k]=true;
      }
      if (x<=d-5 and y==o+1-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+1-4)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+1-6)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+1-7)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+5 and y==o+1+0)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+3 and y==o+1-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o+1-2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+1-3)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+1-5)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+5 and y==o+1-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+3 and y==o+1-2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o+1-3)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+1-5)
      {
        outside[k]=true;
      }
      if (y==0 and x>=d+1 and x<=d+1+LX-2*d-2)
      {
        outside[k]=true;
      }
      if (x<=d and y>=1 and y<=1+o-7)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1 and y>=1 and y<=1+o-7)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1 and y==0)
      {
        outside[k]=true;
      }
      if (x<=d and y==0)
      {
        outside[k]=true;
      }
      if (y==o+1 and x>=LX-d+7 and x<=LX-d+7+d-7)
      {
        outside[k]=true;
      }
      if (y==o+1 and x<=0+d-7)
      {
        outside[k]=true;
      }
      if (y==LY-1 and x<=0+LX)
      {
        outside[k]=true;
      }
    }
    //reservoir with chimney exit, radius corners=10
    else if (BC==34)
    {
      //automatically generated code
      if (x<=d+0 and y==o+1-10)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+1-9)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+1-8)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+1-6)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+1-7)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+1-5)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o+1-4)
      {
        outside[k]=true;
      }
      if (x<=d-3 and y==o+1-3)
      {
        outside[k]=true;
      }
      if (x<=d-4 and y==o+1-2)
      {
        outside[k]=true;
      }
      if (x<=d-5 and y==o+1-1)
      {
        outside[k]=true;
      }
      if (x<=d-7 and y==o+1+0)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+1-7)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o+1-5)
      {
        outside[k]=true;
      }
      if (x<=d-3 and y==o+1-4)
      {
        outside[k]=true;
      }
      if (x<=d-4 and y==o+1-3)
      {
        outside[k]=true;
      }
      if (x<=d-5 and y==o+1-2)
      {
        outside[k]=true;
      }
      if (x<=d-7 and y==o+1-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+1-6)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+1-8)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+1-9)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+1-10)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+7 and y==o+1+0)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+5 and y==o+1-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+4 and y==o+1-2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+3 and y==o+1-3)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o+1-4)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+1-5)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+1-7)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+7 and y==o+1-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+5 and y==o+1-2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+4 and y==o+1-3)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+3 and y==o+1-4)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o+1-5)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+1-7)
      {
        outside[k]=true;
      }
      if (y==0 and x>=d+1 and x<=d+1+LX-2*d-2)
      {
        outside[k]=true;
      }
      if (x<=d and y>=1 and y<=1+o-10)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1 and y>=1 and y<=1+o-10)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1 and y==0)
      {
        outside[k]=true;
      }
      if (x<=d and y==0)
      {
        outside[k]=true;
      }
      if (y==o+1 and x>=LX-d+10 and x<=LX-d+10+d-10)
      {
        outside[k]=true;
      }
      if (y==o+1 and x<=0+d-10)
      {
        outside[k]=true;
      }
      if (y==LY-1 and x<=0+LX)
      {
        outside[k]=true;
      }
    }
    //reservoir with chimney exit, radius corners=13
    else if (BC==35)
    {
      //automatically generated code
      if (x<=d+0 and y==o+1-13)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+1-12)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+1-11)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+1-9)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+1-8)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+1-10)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+1-7)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o+1-6)
      {
        outside[k]=true;
      }
      if (x<=d-3 and y==o+1-5)
      {
        outside[k]=true;
      }
      if (x<=d-4 and y==o+1-4)
      {
        outside[k]=true;
      }
      if (x<=d-5 and y==o+1-3)
      {
        outside[k]=true;
      }
      if (x<=d-6 and y==o+1-2)
      {
        outside[k]=true;
      }
      if (x<=d-7 and y==o+1-1)
      {
        outside[k]=true;
      }
      if (x<=d-10 and y==o+1+0)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+1-10)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o+1-7)
      {
        outside[k]=true;
      }
      if (x<=d-3 and y==o+1-6)
      {
        outside[k]=true;
      }
      if (x<=d-4 and y==o+1-5)
      {
        outside[k]=true;
      }
      if (x<=d-5 and y==o+1-4)
      {
        outside[k]=true;
      }
      if (x<=d-6 and y==o+1-3)
      {
        outside[k]=true;
      }
      if (x<=d-7 and y==o+1-2)
      {
        outside[k]=true;
      }
      if (x<=d-10 and y==o+1-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+1-8)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+1-9)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+1-11)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+1-12)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+1-13)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+10 and y==o+1+0)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+7 and y==o+1-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+6 and y==o+1-2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+5 and y==o+1-3)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+4 and y==o+1-4)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+3 and y==o+1-5)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o+1-6)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+1-7)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+1-10)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+10 and y==o+1-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+7 and y==o+1-2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+6 and y==o+1-3)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+5 and y==o+1-4)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+4 and y==o+1-5)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+3 and y==o+1-6)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o+1-7)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+1-10)
      {
        outside[k]=true;
      }
      if (y==0 and x>=d+1 and x<=d+1+LX-2*d-2)
      {
        outside[k]=true;
      }
      if (x<=d and y>=1 and y<=1+o-13)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1 and y>=1 and y<=1+o-13)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1 and y==0)
      {
        outside[k]=true;
      }
      if (x<=d and y==0)
      {
        outside[k]=true;
      }
      if (y==o+1 and x>=LX-d+13 and x<=LX-d+13+d-13)
      {
        outside[k]=true;
      }
      if (y==o+1 and x<=0+d-13)
      {
        outside[k]=true;
      }
      if (y==LY-1 and x<=0+LX)
      {
        outside[k]=true;
      }
    }
    //reservoir with chimney exit, radius corners=16
    else if (BC==36)
    {
      //automatically generated code
      if (x<=d+0 and y==o+1-16)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+1-15)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+1-14)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+1-12)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+1-11)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o+1-9)
      {
        outside[k]=true;
      }
      if (x<=d-4 and y==o+1-6)
      {
        outside[k]=true;
      }
      if (x<=d+0 and y==o+1-13)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+1-10)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o+1-8)
      {
        outside[k]=true;
      }
      if (x<=d-3 and y==o+1-7)
      {
        outside[k]=true;
      }
      if (x<=d-4 and y==o+1-5)
      {
        outside[k]=true;
      }
      if (x<=d-5 and y==o+1-4)
      {
        outside[k]=true;
      }
      if (x<=d-7 and y==o+1-3)
      {
        outside[k]=true;
      }
      if (x<=d-8 and y==o+1-2)
      {
        outside[k]=true;
      }
      if (x<=d-10 and y==o+1-1)
      {
        outside[k]=true;
      }
      if (x<=d-13 and y==o+1+0)
      {
        outside[k]=true;
      }
      if (x<=d-1 and y==o+1-13)
      {
        outside[k]=true;
      }
      if (x<=d-2 and y==o+1-10)
      {
        outside[k]=true;
      }
      if (x<=d-3 and y==o+1-8)
      {
        outside[k]=true;
      }
      if (x<=d-4 and y==o+1-7)
      {
        outside[k]=true;
      }
      if (x<=d-5 and y==o+1-5)
      {
        outside[k]=true;
      }
      if (x<=d-7 and y==o+1-4)
      {
        outside[k]=true;
      }
      if (x<=d-8 and y==o+1-3)
      {
        outside[k]=true;
      }
      if (x<=d-10 and y==o+1-2)
      {
        outside[k]=true;
      }
      if (x<=d-13 and y==o+1-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+4 and y==o+1-6)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o+1-9)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+1-11)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+1-12)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+1-14)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+1-15)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+1-16)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+13 and y==o+1+0)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+10 and y==o+1-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+8 and y==o+1-2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+7 and y==o+1-3)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+5 and y==o+1-4)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+4 and y==o+1-5)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+3 and y==o+1-7)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o+1-8)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+1-10)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+0 and y==o+1-13)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+13 and y==o+1-1)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+10 and y==o+1-2)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+8 and y==o+1-3)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+7 and y==o+1-4)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+5 and y==o+1-5)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+4 and y==o+1-7)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+3 and y==o+1-8)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+2 and y==o+1-10)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1+1 and y==o+1-13)
      {
        outside[k]=true;
      }
      if (y==0 and x>=d+1 and x<=d+1+LX-2*d-2)
      {
        outside[k]=true;
      }
      if (x<=d and y>=1 and y<=1+o-16)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1 and y>=1 and y<=1+o-16)
      {
        outside[k]=true;
      }
      if (x>=LX-d-1 and y==0)
      {
        outside[k]=true;
      }
      if (x<=d and y==0)
      {
        outside[k]=true;
      }
      if (y==o+1 and x>=LX-d+16 and x<=LX-d+16+d-16)
      {
        outside[k]=true;
      }
      if (y==o+1 and x<=0+d-16)
      {
        outside[k]=true;
      }
      if (y==LY-1 and x<=0+LX)
      {
        outside[k]=true;
      }
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

void LyotropicFreeBoundary::Configure()
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

void LyotropicFreeBoundary::UpdateQuantities()
{
  // sum -> countphi
  double sum = 0;

  #pragma omp parallel for reduction (+:sum) num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
  {
    // skip outside cells
    if(outside[k]) continue;
    // see Lyotropic::UpdateQuantities()
    sum = sum + phi[k];
    UpdateQuantitiesAtNode(k);
  }

  countphi = sum;
}

void LyotropicFreeBoundary::UpdateFields(bool first)
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
#ifdef DEBUG
    const auto& d = get_neighbours(k);
    sum3 = sum3 + laplacian(MU, d, sD);
    sum2= sum2+flux  (phi, ux, uy, d, 0);
#endif
    // see Lyotropic::UpdateFields()
    UpdateFieldsAtNode(k, first);
  }
#ifdef DEBUG
  countphi_laplace=sum3;
  countphi_flux=sum2;
#endif
  swap(phi.get_data(), phi_tmp.get_data());
}

void LyotropicFreeBoundary::BoundaryConditionsLB()
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
    // free slip and von neumann in a channel
    case 2:
    {
      auto apply_bc = [this](LBField& field) {
        // pbc on the left and right walls
        field.ApplyFreeSlipFreeWall(0, 0, LX, Wall::Front);
        field.ApplyFreeSlipFreeWall(0, LY-1, LX, Wall::Back);
        // no corners
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    // free slip and von neumann in a channel in y-direction
    case 201:
    {
      auto apply_bc = [this](LBField& field) {
        // pbc on the left and right walls
        field.ApplyFreeSlipFreeWall(0, 0, LY, Wall::Left);
        field.ApplyFreeSlipFreeWall(LX-1, 0, LY, Wall::Right);
        // no corners
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    // 'inverted' free-slip box (for testing purposes)
    case 3:
    {
      auto apply_bc = [this](LBField& field) {
        // walls
        field.ApplyFreeSlipFreeWall(40, 11, 19, Wall::Left);
        field.ApplyFreeSlipFreeWall(10, 11, 19, Wall::Right);
        field.ApplyFreeSlipFreeWall(11, 30, 29, Wall::Front);
        field.ApplyFreeSlipFreeWall(11, 10, 29, Wall::Back);
        // corners
        field.ApplyFreeSlipFreeConvexCorner(10,10,Corner::LeftFront);
        field.ApplyFreeSlipFreeConvexCorner(10,30,Corner::LeftBack);
        field.ApplyFreeSlipFreeConvexCorner(40,10,Corner::RightFront);
        field.ApplyFreeSlipFreeConvexCorner(40,30,Corner::RightBack);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    // no-slip box
    case 4:
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
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    //no slip and von neumann in a channel
    case 5:
    {
      auto apply_bc = [this](LBField& field) {
        // pbc on the left and right walls
        field.ApplyNoSlipFreeWall(0, 0, LX, Wall::Front);
        field.ApplyNoSlipFreeWall(0, LY-1, LX, Wall::Back);
        // no corners
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    //no slip and von neumann in a channel in y-direction
    case 501:
    {
      auto apply_bc = [this](LBField& field) {
        // pbc on the left and right walls
        field.ApplyNoSlipFreeWall(0, 0, LY, Wall::Left);
        field.ApplyNoSlipFreeWall(LX-1, 0, LY, Wall::Right);
        // no corners
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    // 'inverted' no-slip box (for testing purposes)
    case 6:
    {
      auto apply_bc = [this](LBField& field) {
        // walls
        field.ApplyNoSlipFreeWall(40, 11, 19, Wall::Left);
        field.ApplyNoSlipFreeWall(10, 11, 19, Wall::Right);
        field.ApplyNoSlipFreeWall(11, 30, 29, Wall::Front);
        field.ApplyNoSlipFreeWall(11, 10, 29, Wall::Back);
        // corners
        field.ApplyNoSlipFreeConvexCorner(10,10,Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(10,30,Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(40,10,Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(40,30,Corner::RightBack);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    // channel that tightens and widens again with no-slip and von Neumann
    case 7:
    {
      auto apply_bc = [this](LBField& field) {
        // walls
        field.ApplyNoSlipFreeWall(0, 1,   o-1, Wall::Left);
        field.ApplyNoSlipFreeWall(d, o+1, l-2, Wall::Left);
        field.ApplyNoSlipFreeWall(0, o2,  r-1, Wall::Left);

        field.ApplyNoSlipFreeWall(LX-1,   1,   o-1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1, o+1, l-2, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-1,   o2,  r-1, Wall::Right);

        field.ApplyNoSlipFreeWall(1,    0,    LX-2, Wall::Front);
        field.ApplyNoSlipFreeWall(1,    o2-1, d-1,  Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d, o2-1, d-1,  Wall::Front);

        field.ApplyNoSlipFreeWall(1,    LY-1, LX-2, Wall::Back);
        field.ApplyNoSlipFreeWall(1,    o,    d-1,  Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d, o,    d-1,  Wall::Back);

        // concave corners
        field.ApplyNoSlipFreeConcaveCorner(0, 0,    Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(0, o2-1, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(0, o,    Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(0, LY-1, Corner::LeftBack);

        field.ApplyNoSlipFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-1, o,    Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-1, 0,    Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-1, o2-1, Corner::RightFront);

        //convex corners
        field.ApplyNoSlipFreeConvexCorner(LX-d-1, o,   Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(d,      o,   Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1, o2-1, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(d,      o2-1, Corner::RightBack);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    // channel that tightens and widens again with free-slip and von Neumann
    case 8:
    {
      auto apply_bc = [this](LBField& field) {
        // walls
        field.ApplyFreeSlipFreeWall(0, 1,   o-1, Wall::Left);
        field.ApplyFreeSlipFreeWall(d, o+1, l-2, Wall::Left);
        field.ApplyFreeSlipFreeWall(0, o2,  r-1, Wall::Left);

        field.ApplyFreeSlipFreeWall(LX-1,   1,   o-1, Wall::Right);
        field.ApplyFreeSlipFreeWall(LX-d-1, o+1, l-2, Wall::Right);
        field.ApplyFreeSlipFreeWall(LX-1,   o2,  r-1, Wall::Right);

        field.ApplyFreeSlipFreeWall(1,    0,    LX-2, Wall::Front);
        field.ApplyFreeSlipFreeWall(1,    o2-1, d-1,  Wall::Front);
        field.ApplyFreeSlipFreeWall(LX-d, o2-1, d-1,  Wall::Front);

        field.ApplyFreeSlipFreeWall(1,    LY-1, LX-2, Wall::Back);
        field.ApplyFreeSlipFreeWall(1,    o,    d-1,  Wall::Back);
        field.ApplyFreeSlipFreeWall(LX-d, o,    d-1,  Wall::Back);

        // concave corners
        field.ApplyFreeSlipFreeConcaveCorner(0, 0,    Corner::LeftFront);
        field.ApplyFreeSlipFreeConcaveCorner(0, o2-1, Corner::LeftFront);
        field.ApplyFreeSlipFreeConcaveCorner(0, o,    Corner::LeftBack);
        field.ApplyFreeSlipFreeConcaveCorner(0, LY-1, Corner::LeftBack);

        field.ApplyFreeSlipFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.ApplyFreeSlipFreeConcaveCorner(LX-1, o,    Corner::RightBack);
        field.ApplyFreeSlipFreeConcaveCorner(LX-1, 0,    Corner::RightFront);
        field.ApplyFreeSlipFreeConcaveCorner(LX-1, o2-1, Corner::RightFront);

        //convex corners
        field.ApplyFreeSlipFreeConvexCorner(LX-d-1, o,    Corner::LeftFront);
        field.ApplyFreeSlipFreeConvexCorner(d,      o,    Corner::RightFront);
        field.ApplyFreeSlipFreeConvexCorner(LX-d-1, o2-1, Corner::LeftBack);
        field.ApplyFreeSlipFreeConvexCorner(d,      o2-1, Corner::RightBack);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    // channel that widens with no-slip and von Neumann.
    // In principle as BC==7/BC==8 but without the lower reservoir.
    case 9:
    {
      auto apply_bc = [this](LBField& field) {
        // walls
        field.ApplyNoSlipFreeWall(d, 1,   l-1, Wall::Left);
        field.ApplyNoSlipFreeWall(0, l+1, r-1, Wall::Left);

        field.ApplyNoSlipFreeWall(LX-d-1, 1,   l-1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-1,   l+1, r-1, Wall::Right);

        field.ApplyNoSlipFreeWall(d+1,  0, LX-2-2*d, Wall::Front);
        field.ApplyNoSlipFreeWall(1,    l, d-1,      Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d, l, d-1,      Wall::Front);

        field.ApplyNoSlipFreeWall(1,    LY-1, LX-2, Wall::Back);

        // concave corners
        field.ApplyNoSlipFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(0, l, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(0, LY-1, Corner::LeftBack);

        field.ApplyNoSlipFreeConcaveCorner(LX-1,   LY-1, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1, 0,    Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-1,   l,    Corner::RightFront);

        //convex corners
        field.ApplyNoSlipFreeConvexCorner(LX-d-1, l, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(d,      l, Corner::RightBack);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    // channel that widens with free-slip and von Neumann.
    // In principle as BC==7/BC==8 but without the lower reservoir.
    case 10:
    {
      auto apply_bc = [this](LBField& field) {
        // walls
        field.ApplyFreeSlipFreeWall(d, 1,   l-1, Wall::Left);
        field.ApplyFreeSlipFreeWall(0, l+1, r-1, Wall::Left);

        field.ApplyFreeSlipFreeWall(LX-d-1, 1,   l-1, Wall::Right);
        field.ApplyFreeSlipFreeWall(LX-1,   l+1, r-1, Wall::Right);

        field.ApplyFreeSlipFreeWall(d+1,  0, LX-2-2*d, Wall::Front);
        field.ApplyFreeSlipFreeWall(1,    l, d-1,      Wall::Front);
        field.ApplyFreeSlipFreeWall(LX-d, l, d-1,      Wall::Front);

        field.ApplyFreeSlipFreeWall(1,    LY-1, LX-2, Wall::Back);

        // concave corners
        field.ApplyFreeSlipFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.ApplyFreeSlipFreeConcaveCorner(0, l, Corner::LeftFront);
        field.ApplyFreeSlipFreeConcaveCorner(0, LY-1, Corner::LeftBack);

        field.ApplyFreeSlipFreeConcaveCorner(LX-1,   LY-1, Corner::RightBack);
        field.ApplyFreeSlipFreeConcaveCorner(LX-d-1, 0,    Corner::RightFront);
        field.ApplyFreeSlipFreeConcaveCorner(LX-1,   l,    Corner::RightFront);

        //convex corners
        field.ApplyFreeSlipFreeConvexCorner(LX-d-1, l, Corner::LeftBack);
        field.ApplyFreeSlipFreeConvexCorner(d,      l, Corner::RightBack);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    // channel that tightens and widens again with no-slip and von Neumann
    // Similar to C==7/BC==8, but periodic at the reservoirs
    case 11:
    {
      auto apply_bc = [this](LBField& field) {
        // walls
        field.ApplyNoSlipFreeWall(d, o+1, l-2, Wall::Left);
        field.ApplyNoSlipFreeWall(LX-d-1, o+1, l-2, Wall::Right);

        field.ApplyNoSlipFreeWall(0,    0,    LX , Wall::Front);
        field.ApplyNoSlipFreeWall(0,    o2-1, d  ,  Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d, o2-1, d  ,  Wall::Front);

        field.ApplyNoSlipFreeWall(0,    LY-1, LX , Wall::Back);
        field.ApplyNoSlipFreeWall(0,    o,    d  ,  Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d, o,    d  ,  Wall::Back);

        //convex corners
        field.ApplyNoSlipFreeConvexCorner(LX-d-1, o,   Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(d,      o,   Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1, o2-1, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(d,      o2-1, Corner::RightBack);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    // channel that tightens and widens again with free-slip and von Neumann
    // Similar to C==7/BC==8,but periodic at the reservoirs
    case 12:
    {
      auto apply_bc = [this](LBField& field) {
        // walls
        field.ApplyFreeSlipFreeWall(d, o+1, l-2, Wall::Left);
        field.ApplyFreeSlipFreeWall(LX-d-1, o+1, l-2, Wall::Right);

        field.ApplyFreeSlipFreeWall(0,    0,    LX, Wall::Front);
        field.ApplyFreeSlipFreeWall(0,    o2-1, d ,  Wall::Front);
        field.ApplyFreeSlipFreeWall(LX-d, o2-1, d ,  Wall::Front);

        field.ApplyFreeSlipFreeWall(0,    LY-1, LX, Wall::Back);
        field.ApplyFreeSlipFreeWall(0,    o,    d ,  Wall::Back);
        field.ApplyFreeSlipFreeWall(LX-d, o,    d ,  Wall::Back);

        //convex corners
        field.ApplyFreeSlipFreeConvexCorner(LX-d-1, o,    Corner::LeftFront);
        field.ApplyFreeSlipFreeConvexCorner(d,      o,    Corner::RightFront);
        field.ApplyFreeSlipFreeConvexCorner(LX-d-1, o2-1, Corner::LeftBack);
        field.ApplyFreeSlipFreeConvexCorner(d,      o2-1, Corner::RightBack);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    // channel that tightens and widens again with no-slip and von Neumann
    // and with rounded corners
    case 13:
    {
      auto apply_bc = [this](LBField& field) {
        // walls
        field.ApplyNoSlipFreeWall(0, 1,   o-1, Wall::Left);
        field.ApplyNoSlipFreeWall(d, o+6, l-12, Wall::Left);
        field.ApplyNoSlipFreeWall(0, o2,  r-1, Wall::Left);

        field.ApplyNoSlipFreeWall(LX-1,   1,   o-1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1, o+6, l-12, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-1,   o2,  r-1, Wall::Right);

        field.ApplyNoSlipFreeWall(1,      0,    LX-2, Wall::Front);
        field.ApplyNoSlipFreeWall(1,      o2-1, d-6,  Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d+5, o2-1, d-6,  Wall::Front);

        field.ApplyNoSlipFreeWall(1,      LY-1, LX-2, Wall::Back);
        field.ApplyNoSlipFreeWall(1,      o,    d-6,  Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d+5, o,    d-6,  Wall::Back);

        // concave corners
        field.ApplyNoSlipFreeConcaveCorner(0, 0,    Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(0, o2-1, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(0, o,    Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(0, LY-1, Corner::LeftBack);

        field.ApplyNoSlipFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-1, o,    Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-1, 0,    Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-1, o2-1, Corner::RightFront);

        //rounded convex corners
        //Lower left corner
        field.ApplyNoSlipFreeConvexCorner(d-5,      o,     Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-4,      o+1,   Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-2,      o+2,   Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-1,      o+4,   Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d  ,      o+5,   Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(d-5,      o+1,   Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-4,      o+2,   Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-2,      o+4,   Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-1,      o+5,   Corner::LeftBack);
        field.ApplyNoSlipFreeWall(d-3, o+2, 1,  Wall::Back);
        field.ApplyNoSlipFreeWall(d-2, o+3, 1,  Wall::Left);


        //Lower right corner
        field.ApplyNoSlipFreeConvexCorner(LX-d+4, o,     Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d+3, o+1,   Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d+1, o+2,   Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d  , o+4,   Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1, o+5,   Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d+4, o+1,   Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d+3, o+2,   Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d+1, o+4,   Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d  , o+5,   Corner::RightBack);
        field.ApplyNoSlipFreeWall(LX-d+2, o+2, 1,  Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d+1, o+3, 1,  Wall::Right);

        //Upper left corner
        field.ApplyNoSlipFreeConvexCorner(d-5, o2-1, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-4, o2-2, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-2, o2-3, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-1, o2-5, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d  , o2-6, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(d-5, o2-2, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-4, o2-3, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-2, o2-5, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-1, o2-6, Corner::LeftFront);
        field.ApplyNoSlipFreeWall(d-3, o2-3, 1,  Wall::Front);
        field.ApplyNoSlipFreeWall(d-2, o2-4, 1,  Wall::Left);

        //Upper right corner
        field.ApplyNoSlipFreeConvexCorner(LX-d+4, o2-1, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d+3, o2-2, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d+1, o2-3, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d  , o2-5, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1, o2-6, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d+4, o2-2, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d+3, o2-3, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d+1, o2-5, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d  , o2-6, Corner::RightFront);
        field.ApplyNoSlipFreeWall(LX-d+2, o2-3, 1,  Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d+1, o2-4, 1,  Wall::Right);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    //other type of rounded corners
    case 15:
    {
      //automatically generated code
      auto apply_bc = [this](LBField& field) {

        field.ApplyNoSlipFreeWall(0, 1, o-1, Wall::Left);
        field.ApplyNoSlipFreeWall(d, o+16, l-32, Wall::Left);
        field.ApplyNoSlipFreeWall(0, o2, r-1, Wall::Left);
        field.ApplyNoSlipFreeWall(LX-1, 1, o-1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1, o+16, l-32, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-1, o2, r-1, Wall::Right);
        field.ApplyNoSlipFreeWall(1, 0, LX-2, Wall::Front);
        field.ApplyNoSlipFreeWall(1, o2-1, d-14, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d+15, o2-1, d-14, Wall::Front);
        field.ApplyNoSlipFreeWall(1, LY-1, LX-2, Wall::Back);
        field.ApplyNoSlipFreeWall(1, o, d-14, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d+15, o, d-14, Wall::Back);
        field.ApplyNoSlipFreeWall(d-2, o+8, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-1, o+10, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-1, o+11, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+13, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+14, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+15, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-15, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-14, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-13, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-11, o+1, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-10, o+1, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-8, o+2, 1, Wall::Back);
        field.ApplyNoSlipFreeConvexCorner(d-12, o+0, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-9, o+1, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-7, o+2, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-6, o+3, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-5, o+4, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-4, o+5, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-3, o+6, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-2, o+7, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-1, o+9, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d+0, o+12, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(d-12, o+1, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-9, o+2, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-7, o+3, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-6, o+4, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-5, o+5, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-4, o+6, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-3, o+7, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-2, o+9, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-1, o+12, Corner::LeftBack);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+15, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+14, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+13, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o+11, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o+10, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+2, o+8, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+8, o+2, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+10, o+1, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+11, o+1, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+13, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+14, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+15, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+0, o+12, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+1, o+9, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+2, o+7, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+3, o+6, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+4, o+5, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+5, o+4, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+6, o+3, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+7, o+2, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+9, o+1, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+12, o+0, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+1, o+12, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+2, o+9, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+3, o+7, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+4, o+6, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+5, o+5, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+6, o+4, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+7, o+3, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+9, o+2, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+12, o+1, Corner::RightBack);
        field.ApplyNoSlipFreeWall(d+0, o2-1-15, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o2-1-14, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o2-1-13, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-1, o2-1-11, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-1, o2-1-10, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-2, o2-1-8, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-8, o2-1-2, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-10, o2-1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-11, o2-1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-13, o2-1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-14, o2-1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-15, o2-1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeConvexCorner(d+0, o2-1-12, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-1, o2-1-9, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-2, o2-1-7, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-3, o2-1-6, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-4, o2-1-5, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-5, o2-1-4, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-6, o2-1-3, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-7, o2-1-2, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-9, o2-1-1, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-12, o2-1+0, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(d-1, o2-1-12, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-2, o2-1-9, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-3, o2-1-7, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-4, o2-1-6, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-5, o2-1-5, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-6, o2-1-4, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-7, o2-1-3, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-9, o2-1-2, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-12, o2-1-1, Corner::LeftFront);
        field.ApplyNoSlipFreeWall(LX-d-1+2, o2-1-8, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o2-1-10, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o2-1-11, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o2-1-13, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o2-1-14, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o2-1-15, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+15, o2-1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+14, o2-1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+13, o2-1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+11, o2-1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+10, o2-1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+8, o2-1-2, 1, Wall::Front);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+12, o2-1+0, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+9, o2-1-1, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+7, o2-1-2, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+6, o2-1-3, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+5, o2-1-4, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+4, o2-1-5, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+3, o2-1-6, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+2, o2-1-7, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+1, o2-1-9, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+0, o2-1-12, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+12, o2-1-1, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+9, o2-1-2, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+7, o2-1-3, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+6, o2-1-4, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+5, o2-1-5, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+4, o2-1-6, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+3, o2-1-7, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+2, o2-1-9, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+1, o2-1-12, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(0, 0, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(0, o2-1, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(0, o, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(0, LY-1, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-1, o, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-1, 0, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-1, o2-1, Corner::RightFront);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    //other type of rounded corners with periodic reservoir
    case 17:
    {
      //automatically generated code
      auto apply_bc = [this](LBField& field) {

        field.ApplyNoSlipFreeWall(d, o+16, l-32, Wall::Left);
        field.ApplyNoSlipFreeWall(LX-d-1, o+16, l-32, Wall::Right);
        field.ApplyNoSlipFreeWall(0, 0, LX, Wall::Front);
        field.ApplyNoSlipFreeWall(0, o2-1, d-15, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d+15, o2-1, d-15, Wall::Front);
        field.ApplyNoSlipFreeWall(0, LY-1, LX, Wall::Back);
        field.ApplyNoSlipFreeWall(0, o, d-15, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d+15, o, d-15, Wall::Back);
        field.ApplyNoSlipFreeWall(d-2, o+8, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-1, o+10, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-1, o+11, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+13, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+14, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+15, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-15, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-14, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-13, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-11, o+1, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-10, o+1, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-8, o+2, 1, Wall::Back);
        field.ApplyNoSlipFreeConvexCorner(d-12, o+0, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-9, o+1, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-7, o+2, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-6, o+3, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-5, o+4, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-4, o+5, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-3, o+6, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-2, o+7, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-1, o+9, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d+0, o+12, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(d-12, o+1, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-9, o+2, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-7, o+3, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-6, o+4, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-5, o+5, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-4, o+6, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-3, o+7, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-2, o+9, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-1, o+12, Corner::LeftBack);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+15, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+14, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+13, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o+11, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o+10, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+2, o+8, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+8, o+2, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+10, o+1, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+11, o+1, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+13, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+14, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+15, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+0, o+12, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+1, o+9, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+2, o+7, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+3, o+6, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+4, o+5, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+5, o+4, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+6, o+3, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+7, o+2, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+9, o+1, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+12, o+0, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+1, o+12, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+2, o+9, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+3, o+7, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+4, o+6, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+5, o+5, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+6, o+4, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+7, o+3, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+9, o+2, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+12, o+1, Corner::RightBack);
        field.ApplyNoSlipFreeWall(d+0, o2-1-15, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o2-1-14, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o2-1-13, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-1, o2-1-11, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-1, o2-1-10, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-2, o2-1-8, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-8, o2-1-2, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-10, o2-1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-11, o2-1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-13, o2-1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-14, o2-1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-15, o2-1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeConvexCorner(d+0, o2-1-12, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-1, o2-1-9, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-2, o2-1-7, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-3, o2-1-6, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-4, o2-1-5, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-5, o2-1-4, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-6, o2-1-3, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-7, o2-1-2, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-9, o2-1-1, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-12, o2-1+0, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(d-1, o2-1-12, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-2, o2-1-9, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-3, o2-1-7, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-4, o2-1-6, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-5, o2-1-5, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-6, o2-1-4, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-7, o2-1-3, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-9, o2-1-2, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-12, o2-1-1, Corner::LeftFront);
        field.ApplyNoSlipFreeWall(LX-d-1+2, o2-1-8, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o2-1-10, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o2-1-11, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o2-1-13, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o2-1-14, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o2-1-15, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+15, o2-1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+14, o2-1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+13, o2-1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+11, o2-1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+10, o2-1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+8, o2-1-2, 1, Wall::Front);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+12, o2-1+0, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+9, o2-1-1, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+7, o2-1-2, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+6, o2-1-3, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+5, o2-1-4, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+4, o2-1-5, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+3, o2-1-6, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+2, o2-1-7, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+1, o2-1-9, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+0, o2-1-12, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+12, o2-1-1, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+9, o2-1-2, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+7, o2-1-3, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+6, o2-1-4, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+5, o2-1-5, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+4, o2-1-6, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+3, o2-1-7, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+2, o2-1-9, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+1, o2-1-12, Corner::RightFront);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    //other type of rounded corners with periodic reservoir with von Neumann at upper end
    case 19:
    {
      //automatically generated code
      auto apply_bc = [this](LBField& field) {

        field.ApplyNoSlipFreeWall(d, o+16, l-32, Wall::Left);
        field.ApplyNoSlipFreeWall(LX-d-1, o+16, l-32, Wall::Right);
        field.ApplyNoSlipFreeWall(0, 0, LX, Wall::Front);
        field.ApplyNoSlipFreeWall(0, o2-1, d-15, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d+15, o2-1, d-15, Wall::Front);
        field.ApplyNeumannFreeWall(0, LY-1, LX, Wall::Back);
        field.ApplyNoSlipFreeWall(0, o, d-15, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d+15, o, d-15, Wall::Back);
        field.ApplyNoSlipFreeWall(d-2, o+8, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-1, o+10, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-1, o+11, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+13, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+14, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+15, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-15, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-14, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-13, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-11, o+1, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-10, o+1, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-8, o+2, 1, Wall::Back);
        field.ApplyNoSlipFreeConvexCorner(d-12, o+0, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-9, o+1, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-7, o+2, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-6, o+3, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-5, o+4, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-4, o+5, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-3, o+6, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-2, o+7, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-1, o+9, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d+0, o+12, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(d-12, o+1, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-9, o+2, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-7, o+3, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-6, o+4, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-5, o+5, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-4, o+6, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-3, o+7, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-2, o+9, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-1, o+12, Corner::LeftBack);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+15, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+14, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+13, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o+11, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o+10, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+2, o+8, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+8, o+2, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+10, o+1, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+11, o+1, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+13, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+14, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+15, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+0, o+12, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+1, o+9, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+2, o+7, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+3, o+6, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+4, o+5, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+5, o+4, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+6, o+3, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+7, o+2, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+9, o+1, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+12, o+0, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+1, o+12, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+2, o+9, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+3, o+7, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+4, o+6, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+5, o+5, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+6, o+4, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+7, o+3, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+9, o+2, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+12, o+1, Corner::RightBack);
        field.ApplyNoSlipFreeWall(d+0, o2-1-15, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o2-1-14, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o2-1-13, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-1, o2-1-11, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-1, o2-1-10, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-2, o2-1-8, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-8, o2-1-2, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-10, o2-1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-11, o2-1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-13, o2-1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-14, o2-1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-15, o2-1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeConvexCorner(d+0, o2-1-12, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-1, o2-1-9, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-2, o2-1-7, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-3, o2-1-6, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-4, o2-1-5, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-5, o2-1-4, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-6, o2-1-3, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-7, o2-1-2, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-9, o2-1-1, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-12, o2-1+0, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(d-1, o2-1-12, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-2, o2-1-9, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-3, o2-1-7, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-4, o2-1-6, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-5, o2-1-5, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-6, o2-1-4, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-7, o2-1-3, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-9, o2-1-2, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-12, o2-1-1, Corner::LeftFront);
        field.ApplyNoSlipFreeWall(LX-d-1+2, o2-1-8, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o2-1-10, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o2-1-11, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o2-1-13, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o2-1-14, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o2-1-15, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+15, o2-1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+14, o2-1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+13, o2-1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+11, o2-1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+10, o2-1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+8, o2-1-2, 1, Wall::Front);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+12, o2-1+0, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+9, o2-1-1, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+7, o2-1-2, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+6, o2-1-3, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+5, o2-1-4, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+4, o2-1-5, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+3, o2-1-6, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+2, o2-1-7, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+1, o2-1-9, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+0, o2-1-12, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+12, o2-1-1, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+9, o2-1-2, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+7, o2-1-3, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+6, o2-1-4, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+5, o2-1-5, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+4, o2-1-6, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+3, o2-1-7, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+2, o2-1-9, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+1, o2-1-12, Corner::RightFront);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    //reservoir with chimney
    case 20:
    {
      //automatically generated code
      auto apply_bc = [this](LBField& field) {

        field.ApplyNoSlipFreeWall(d, o+1, LY-o-2, Wall::Left);
        field.ApplyNoSlipFreeWall(LX-d-1, o+1, LY-o-2, Wall::Right);
        field.ApplyNoSlipFreeWall(0, 0, LX, Wall::Front);
        field.ApplyNoSlipFreeWall(d+1, LY-1, LX-2-d*2, Wall::Back);
        //field.ApplyPressureOutletFreeWall(d+1, LY-1, LX-2-d*2, Wall::Back, n, ux, uy);
        field.ApplyNoSlipFreeWall(LX-d, o, d, Wall::Back);
        field.ApplyNoSlipFreeWall(0, o, d, Wall::Back);
        field.ApplyNoSlipFreeConcaveCorner(d, LY-1, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1, o, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(d, o, Corner::RightFront);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    //box with outfulx at top
    case 21:
    {
      //automatically generated code
      auto apply_bc = [this](LBField& field) {

        field.ApplyNoSlipFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNoSlipFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyNoSlipFreeWall(1, 0, LX-2, Wall::Front);
        field.ApplyNeumannFreeWall(1, LY-1, LX-2, Wall::Back);
        field.ApplyNoSlipFreeConcaveCorner(0, 0, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(0, LY-1, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-1, 0, Corner::RightFront);
      };
      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    //box with influx at bottom and outfulx at top
    case 22:
    {
      //automatically generated code
      auto apply_bc = [this](LBField& field) {

        field.ApplyNoSlipFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNoSlipFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyInletFreeWall(1, 0, LX-2, Wall::Front,0,uin,n);
        field.ApplyOutletFreeWall(1, LY-1, LX-2, Wall::Back,uin);
        field.ApplyNoSlipFreeConcaveCorner(0, 0, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(0, LY-1, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-1, 0, Corner::RightFront);
      };
      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    //box with pressure boundary at top
    case 23:
    //box with pressure boundary at top and diffusive inlet at bottom
    case 24:
    {
      //automatically generated code
      auto apply_bc = [this](LBField& field) {

        field.ApplyNoSlipFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNoSlipFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyNoSlipFreeWall(1, 0, LX-2, Wall::Front);
        field.ApplyPressureOutletFreeWall(1, LY-1, LX-2, Wall::Back, n, ux, uy);
        field.ApplyNoSlipFreeConcaveCorner(0, 0, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(0, LY-1, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-1, 0, Corner::RightFront);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    //box with pressure boundary and inlet at bottom
    case 25:
    {
      auto apply_bc = [this](LBField& field) {

        field.ApplyNoSlipFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNoSlipFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyNoSlipFreeConcaveCorner(0, 0, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(0, LY-1, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-1, 0, Corner::RightFront);
        field.ApplyPressureOutletFreeWall(0, LY-1, LX-1, Wall::Back, n, ux, uy);
        field.ApplyInletFreeWall(0, 0, LX-1, Wall::Front,0,uin,n);

      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    //reservoir with chimney, radius corners=4
    case 26:
    {
      //automatically generated code
      auto apply_bc = [this](LBField& field) {

        field.ApplyNoSlipFreeWall(d, o+5, LY-o-6, Wall::Left);
        field.ApplyNoSlipFreeWall(LX-d-1, o+5, LY-o-6, Wall::Right);
        field.ApplyNoSlipFreeWall(0, 0, LX, Wall::Front);
        field.ApplyPressureOutletFreeWall(d+1, LY-1, LX-2-d*2, Wall::Back, n, ux, uy);
        field.ApplyNoSlipFreeConcaveCorner(d, LY-1, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack);
        field.ApplyNoSlipFreeWall(LX-d+4, o, d-4, Wall::Back);
        field.ApplyNoSlipFreeWall(0, o, d-4, Wall::Back);
        field.ApplyNoSlipFreeWall(d+0, o+4, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-4, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeConvexCorner(d-3, o+0, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-2, o+1, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-1, o+2, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d+0, o+3, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(d-3, o+1, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-2, o+2, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-1, o+3, Corner::LeftBack);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+4, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+4, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+0, o+3, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+1, o+2, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+2, o+1, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+3, o+0, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+1, o+3, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+2, o+2, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+3, o+1, Corner::RightBack);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    //reservoir with chimney, radius corners=7
    case 27:
    {
      //automatically generated code
      auto apply_bc = [this](LBField& field) {

        field.ApplyNoSlipFreeWall(d, o+8, LY-o-9, Wall::Left);
        field.ApplyNoSlipFreeWall(LX-d-1, o+8, LY-o-9, Wall::Right);
        field.ApplyNoSlipFreeWall(0, 0, LX, Wall::Front);
        field.ApplyPressureOutletFreeWall(d+1, LY-1, LX-2-d*2, Wall::Back, n, ux, uy);
        field.ApplyNoSlipFreeConcaveCorner(d, LY-1, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack);
        field.ApplyNoSlipFreeWall(LX-d+7, o, d-7, Wall::Back);
        field.ApplyNoSlipFreeWall(0, o, d-7, Wall::Back);
        field.ApplyNoSlipFreeWall(d-1, o+4, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+6, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+7, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-7, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-6, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-4, o+1, 1, Wall::Back);
        field.ApplyNoSlipFreeConvexCorner(d-5, o+0, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-3, o+1, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-2, o+2, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-1, o+3, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d+0, o+5, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(d-5, o+1, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-3, o+2, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-2, o+3, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-1, o+5, Corner::LeftBack);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+7, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+6, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o+4, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+4, o+1, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+6, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+7, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+0, o+5, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+1, o+3, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+2, o+2, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+3, o+1, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+5, o+0, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+1, o+5, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+2, o+3, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+3, o+2, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+5, o+1, Corner::RightBack);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    //reservoir with chimney, radius corners=10
    case 28:
    {
      //automatically generated code
      auto apply_bc = [this](LBField& field) {

        field.ApplyNoSlipFreeWall(d, o+11, LY-o-12, Wall::Left);
        field.ApplyNoSlipFreeWall(LX-d-1, o+11, LY-o-12, Wall::Right);
        field.ApplyNoSlipFreeWall(0, 0, LX, Wall::Front);
        field.ApplyPressureOutletFreeWall(d+1, LY-1, LX-2-d*2, Wall::Back, n, ux, uy);
        field.ApplyNoSlipFreeConcaveCorner(d, LY-1, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack);
        field.ApplyNoSlipFreeWall(LX-d+10, o, d-10, Wall::Back);
        field.ApplyNoSlipFreeWall(0, o, d-10, Wall::Back);
        field.ApplyNoSlipFreeWall(d-1, o+6, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+8, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+9, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+10, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-10, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-9, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-8, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-6, o+1, 1, Wall::Back);
        field.ApplyNoSlipFreeConvexCorner(d-7, o+0, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-5, o+1, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-4, o+2, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-3, o+3, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-2, o+4, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-1, o+5, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d+0, o+7, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(d-7, o+1, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-5, o+2, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-4, o+3, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-3, o+4, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-2, o+5, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-1, o+7, Corner::LeftBack);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+10, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+9, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+8, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o+6, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+6, o+1, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+8, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+9, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+10, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+0, o+7, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+1, o+5, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+2, o+4, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+3, o+3, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+4, o+2, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+5, o+1, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+7, o+0, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+1, o+7, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+2, o+5, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+3, o+4, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+4, o+3, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+5, o+2, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+7, o+1, Corner::RightBack);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    //reservoir with chimney, radius corners=13
    case 29:
    {
      //automatically generated code
      auto apply_bc = [this](LBField& field) {

        field.ApplyNoSlipFreeWall(d, o+14, LY-o-15, Wall::Left);
        field.ApplyNoSlipFreeWall(LX-d-1, o+14, LY-o-15, Wall::Right);
        field.ApplyNoSlipFreeWall(0, 0, LX, Wall::Front);
        field.ApplyPressureOutletFreeWall(d+1, LY-1, LX-2-d*2, Wall::Back, n, ux, uy);
        field.ApplyNoSlipFreeConcaveCorner(d, LY-1, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack);
        field.ApplyNoSlipFreeWall(LX-d+13, o, d-13, Wall::Back);
        field.ApplyNoSlipFreeWall(0, o, d-13, Wall::Back);
        field.ApplyNoSlipFreeWall(d-1, o+8, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-1, o+9, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+11, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+12, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+13, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-13, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-12, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-11, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-9, o+1, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-8, o+1, 1, Wall::Back);
        field.ApplyNoSlipFreeConvexCorner(d-10, o+0, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-7, o+1, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-6, o+2, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-5, o+3, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-4, o+4, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-3, o+5, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-2, o+6, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-1, o+7, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d+0, o+10, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(d-10, o+1, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-7, o+2, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-6, o+3, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-5, o+4, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-4, o+5, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-3, o+6, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-2, o+7, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-1, o+10, Corner::LeftBack);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+13, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+12, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+11, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o+9, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o+8, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+8, o+1, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+9, o+1, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+11, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+12, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+13, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+0, o+10, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+1, o+7, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+2, o+6, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+3, o+5, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+4, o+4, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+5, o+3, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+6, o+2, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+7, o+1, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+10, o+0, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+1, o+10, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+2, o+7, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+3, o+6, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+4, o+5, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+5, o+4, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+6, o+3, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+7, o+2, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+10, o+1, Corner::RightBack);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    //reservoir with chimney, radius corners=16
    case 30:
    {
      //automatically generated code
      auto apply_bc = [this](LBField& field) {

        field.ApplyNoSlipFreeWall(d, o+17, LY-o-18, Wall::Left);
        field.ApplyNoSlipFreeWall(LX-d-1, o+17, LY-o-18, Wall::Right);
        field.ApplyNoSlipFreeWall(0, 0, LX, Wall::Front);
        field.ApplyPressureOutletFreeWall(d+1, LY-1, LX-2-d*2, Wall::Back, n, ux, uy);
        field.ApplyNoSlipFreeConcaveCorner(d, LY-1, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack);
        field.ApplyNoSlipFreeWall(LX-d+16, o, d-16, Wall::Back);
        field.ApplyNoSlipFreeWall(0, o, d-16, Wall::Back);
        field.ApplyNoSlipFreeWall(d-4, o+6, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-2, o+9, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-1, o+11, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-1, o+12, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+14, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+15, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+16, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-16, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-15, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-14, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-12, o+1, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-11, o+1, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-9, o+2, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(d-6, o+4, 1, Wall::Back);
        field.ApplyNoSlipFreeConvexCorner(d-13, o+0, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-10, o+1, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-8, o+2, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-7, o+3, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-5, o+4, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-4, o+5, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-3, o+7, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-2, o+8, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d-1, o+10, Corner::RightFront);
        field.ApplyNoSlipFreeConvexCorner(d+0, o+13, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(d-13, o+1, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-10, o+2, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-8, o+3, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-7, o+4, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-5, o+5, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-4, o+7, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-3, o+8, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-2, o+10, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(d-1, o+13, Corner::LeftBack);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+16, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+15, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+14, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o+12, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o+11, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+2, o+9, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+4, o+6, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+6, o+4, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+9, o+2, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+11, o+1, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+12, o+1, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+14, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+15, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeWall(LX-d-1+16, o+0, 1, Wall::Back);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+0, o+13, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+1, o+10, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+2, o+8, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+3, o+7, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+4, o+5, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+5, o+4, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+7, o+3, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+8, o+2, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+10, o+1, Corner::LeftFront);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+13, o+0, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+1, o+13, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+2, o+10, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+3, o+8, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+4, o+7, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+5, o+5, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+7, o+4, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+8, o+3, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+10, o+2, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+13, o+1, Corner::RightBack);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    //reservoir with chimney exit
    case 31:
    {
      //automatically generated code
      auto apply_bc = [this](LBField& field) {

        field.ApplyNoSlipFreeConcaveCorner(LX-d-1, 0, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.ApplyNoSlipFreeWall(d+1, 0, LX-2*d-2, Wall::Front);
        field.ApplyNoSlipFreeWall(d, 1, o, Wall::Left);
        field.ApplyNoSlipFreeWall(LX-d-1, 1, o, Wall::Right);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1, o+1, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(d, o+1, Corner::RightBack);
        field.ApplyNoSlipFreeWall(LX-d, o+1, d, Wall::Front);
        field.ApplyNoSlipFreeWall(0, o+1, d, Wall::Front);
        field.ApplyPressureOutletFreeWall(0, LY-1, LX, Wall::Back, n, ux, uy);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }

    //reservoir with chimney exit, radius corners=4
    case 32:
    {
      //automatically generated code
      auto apply_bc = [this](LBField& field) {

        field.ApplyNoSlipFreeWall(d+0, o+1-4, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-4, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeConvexCorner(d+0, o+1-3, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-1, o+1-2, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-2, o+1-1, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-3, o+1+0, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(d-1, o+1-3, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-2, o+1-2, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-3, o+1-1, Corner::LeftFront);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+1-4, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+4, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+3, o+1+0, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+2, o+1-1, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+1, o+1-2, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+0, o+1-3, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+3, o+1-1, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+2, o+1-2, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+1, o+1-3, Corner::RightFront);
        field.ApplyNoSlipFreeWall(d+1, 0, LX-2*d-2, Wall::Front);
        field.ApplyNoSlipFreeWall(d, 1, o-4, Wall::Left);
        field.ApplyNoSlipFreeWall(LX-d-1, 1, o-4, Wall::Right);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1, 0, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.ApplyNoSlipFreeWall(LX-d+4, o+1, d-4, Wall::Front);
        field.ApplyNoSlipFreeWall(0, o+1, d-4, Wall::Front);
        field.ApplyPressureOutletFreeWall(0, LY-1, LX, Wall::Back, n, ux, uy);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    //reservoir with chimney exit, radius corners=7
    case 33:
    {
      //automatically generated code
      auto apply_bc = [this](LBField& field) {

        field.ApplyNoSlipFreeWall(d+0, o+1-7, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+1-6, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-1, o+1-4, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-4, o+1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-6, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-7, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeConvexCorner(d+0, o+1-5, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-1, o+1-3, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-2, o+1-2, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-3, o+1-1, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-5, o+1+0, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(d-1, o+1-5, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-2, o+1-3, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-3, o+1-2, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-5, o+1-1, Corner::LeftFront);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o+1-4, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+1-6, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+1-7, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+7, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+6, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+4, o+1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+5, o+1+0, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+3, o+1-1, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+2, o+1-2, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+1, o+1-3, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+0, o+1-5, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+5, o+1-1, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+3, o+1-2, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+2, o+1-3, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+1, o+1-5, Corner::RightFront);
        field.ApplyNoSlipFreeWall(d+1, 0, LX-2*d-2, Wall::Front);
        field.ApplyNoSlipFreeWall(d, 1, o-7, Wall::Left);
        field.ApplyNoSlipFreeWall(LX-d-1, 1, o-7, Wall::Right);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1, 0, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.ApplyNoSlipFreeWall(LX-d+7, o+1, d-7, Wall::Front);
        field.ApplyNoSlipFreeWall(0, o+1, d-7, Wall::Front);
        field.ApplyPressureOutletFreeWall(0, LY-1, LX, Wall::Back, n, ux, uy);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    //reservoir with chimney exit, radius corners=10
    case 34:
    {
      //automatically generated code
      auto apply_bc = [this](LBField& field) {

        field.ApplyNoSlipFreeWall(d+0, o+1-10, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+1-9, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+1-8, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-1, o+1-6, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-6, o+1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-8, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-9, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-10, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeConvexCorner(d+0, o+1-7, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-1, o+1-5, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-2, o+1-4, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-3, o+1-3, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-4, o+1-2, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-5, o+1-1, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-7, o+1+0, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(d-1, o+1-7, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-2, o+1-5, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-3, o+1-4, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-4, o+1-3, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-5, o+1-2, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-7, o+1-1, Corner::LeftFront);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o+1-6, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+1-8, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+1-9, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+1-10, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+10, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+9, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+8, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+6, o+1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+7, o+1+0, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+5, o+1-1, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+4, o+1-2, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+3, o+1-3, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+2, o+1-4, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+1, o+1-5, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+0, o+1-7, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+7, o+1-1, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+5, o+1-2, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+4, o+1-3, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+3, o+1-4, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+2, o+1-5, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+1, o+1-7, Corner::RightFront);
        field.ApplyNoSlipFreeWall(d+1, 0, LX-2*d-2, Wall::Front);
        field.ApplyNoSlipFreeWall(d, 1, o-10, Wall::Left);
        field.ApplyNoSlipFreeWall(LX-d-1, 1, o-10, Wall::Right);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1, 0, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.ApplyNoSlipFreeWall(LX-d+10, o+1, d-10, Wall::Front);
        field.ApplyNoSlipFreeWall(0, o+1, d-10, Wall::Front);
        field.ApplyPressureOutletFreeWall(0, LY-1, LX, Wall::Back, n, ux, uy);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    //reservoir with chimney, radius corners=13
    case 35:
    {
      //automatically generated code
      auto apply_bc = [this](LBField& field) {

        field.ApplyNoSlipFreeWall(d+0, o+1-13, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+1-12, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+1-11, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-1, o+1-9, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-1, o+1-8, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-8, o+1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-9, o+1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-11, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-12, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-13, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeConvexCorner(d+0, o+1-10, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-1, o+1-7, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-2, o+1-6, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-3, o+1-5, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-4, o+1-4, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-5, o+1-3, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-6, o+1-2, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-7, o+1-1, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-10, o+1+0, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(d-1, o+1-10, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-2, o+1-7, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-3, o+1-6, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-4, o+1-5, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-5, o+1-4, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-6, o+1-3, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-7, o+1-2, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-10, o+1-1, Corner::LeftFront);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o+1-8, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o+1-9, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+1-11, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+1-12, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+1-13, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+13, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+12, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+11, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+9, o+1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+8, o+1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+10, o+1+0, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+7, o+1-1, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+6, o+1-2, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+5, o+1-3, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+4, o+1-4, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+3, o+1-5, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+2, o+1-6, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+1, o+1-7, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+0, o+1-10, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+10, o+1-1, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+7, o+1-2, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+6, o+1-3, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+5, o+1-4, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+4, o+1-5, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+3, o+1-6, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+2, o+1-7, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+1, o+1-10, Corner::RightFront);
        field.ApplyNoSlipFreeWall(d+1, 0, LX-2*d-2, Wall::Front);
        field.ApplyNoSlipFreeWall(d, 1, o-13, Wall::Left);
        field.ApplyNoSlipFreeWall(LX-d-1, 1, o-13, Wall::Right);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1, 0, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.ApplyNoSlipFreeWall(LX-d+13, o+1, d-13, Wall::Front);
        field.ApplyNoSlipFreeWall(0, o+1, d-13, Wall::Front);
        field.ApplyPressureOutletFreeWall(0, LY-1, LX, Wall::Back, n, ux, uy);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    //reservoir with chimney exit, radius corners=16
    case 36:
    {
      //automatically generated code
      auto apply_bc = [this](LBField& field) {

        field.ApplyNoSlipFreeWall(d+0, o+1-16, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+1-15, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+1-14, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-1, o+1-12, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-1, o+1-11, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-2, o+1-9, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-4, o+1-6, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-6, o+1-4, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-9, o+1-2, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-11, o+1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-12, o+1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-14, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-15, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-16, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeConvexCorner(d+0, o+1-13, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-1, o+1-10, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-2, o+1-8, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-3, o+1-7, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-4, o+1-5, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-5, o+1-4, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-7, o+1-3, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-8, o+1-2, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-10, o+1-1, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-13, o+1+0, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(d-1, o+1-13, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-2, o+1-10, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-3, o+1-8, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-4, o+1-7, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-5, o+1-5, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-7, o+1-4, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-8, o+1-3, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-10, o+1-2, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-13, o+1-1, Corner::LeftFront);
        field.ApplyNoSlipFreeWall(LX-d-1+4, o+1-6, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+2, o+1-9, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o+1-11, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o+1-12, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+1-14, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+1-15, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+1-16, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+16, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+15, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+14, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+12, o+1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+11, o+1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+9, o+1-2, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+6, o+1-4, 1, Wall::Front);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+13, o+1+0, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+10, o+1-1, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+8, o+1-2, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+7, o+1-3, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+5, o+1-4, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+4, o+1-5, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+3, o+1-7, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+2, o+1-8, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+1, o+1-10, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+0, o+1-13, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+13, o+1-1, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+10, o+1-2, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+8, o+1-3, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+7, o+1-4, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+5, o+1-5, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+4, o+1-7, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+3, o+1-8, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+2, o+1-10, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+1, o+1-13, Corner::RightFront);
        field.ApplyNoSlipFreeWall(d+1, 0, LX-2*d-2, Wall::Front);
        field.ApplyNoSlipFreeWall(d, 1, o-16, Wall::Left);
        field.ApplyNoSlipFreeWall(LX-d-1, 1, o-16, Wall::Right);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1, 0, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.ApplyNoSlipFreeWall(LX-d+16, o+1, d-16, Wall::Front);
        field.ApplyNoSlipFreeWall(0, o+1, d-16, Wall::Front);
        field.ApplyPressureOutletFreeWall(0, LY-1, LX, Wall::Back, n, ux, uy);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
  }
}

void LyotropicFreeBoundary::BoundaryConditionsFields()
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
      apply_bc(phi);
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
      apply_bc(phi);
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
      apply_bc(phi);
      break;
    }


    // von neumann box
    case 1:
    case 4:
    case 21:
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
      apply_bc(phi);
      break;
    }
    case 25:
    {
      auto apply_QQxx_bc = [this](ScalarField& field) {
        // walls
        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front, cos(2*angle));
        field.CopyDerivativeFreeWall(1, LY-1, LX-2, Wall::Back);
        // corners
        field.ApplyNeumannFreeConcaveCorner(LX-1,LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1,0, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(0,LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(0,0, Corner::LeftFront);
      };
      auto apply_QQyx_bc = [this](ScalarField& field) {
        // walls
        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front, sin(2*angle));
        field.CopyDerivativeFreeWall(1, LY-1, LX-2, Wall::Back);
        // corners
        field.ApplyNeumannFreeConcaveCorner(LX-1,LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1,0, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(0,LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(0,0, Corner::LeftFront);
      };
      auto apply_phi_bc = [this](ScalarField& field) {
        // walls
        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front, conc);
        field.CopyDerivativeFreeWall(1, LY-1, LX-2, Wall::Back);
        // corners
        field.ApplyNeumannFreeConcaveCorner(LX-1,LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1,0, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(0,LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(0,0, Corner::LeftFront);
      };
      apply_QQxx_bc(QQxx);
      apply_QQyx_bc(QQyx);
      apply_phi_bc(phi);
      break;
    }
    // von neumann channel
    case 2:
    case 5:
    {
      auto apply_bc = [this](ScalarField& field) {
        // pbc on the left and right walls
        field.ApplyNeumannFreeWall(0, 0, LX, Wall::Front);
        field.ApplyNeumannFreeWall(0, LY-1, LX, Wall::Back);
        // no corners
      };
      apply_bc(QQxx);
      apply_bc(QQyx);
      apply_bc(phi);
      break;
    }
    // von neumann channel in y-direction
    case 201:
    case 501:
    {
      auto apply_bc = [this](ScalarField& field) {
        // pbc on the left and right walls
        field.ApplyNeumannFreeWall(0, 0, LY, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 0, LY, Wall::Right);
        // no corners
      };
      apply_bc(QQxx);
      apply_bc(QQyx);
      apply_bc(phi);
      break;
    }
    //not implemented at the moment
    case 3:
    case 6:
    {
      break;
    }
    //channel that tightens and widens again with free-slip or no-slip and von Neumann
    case 7:
    case 8:
    {
      auto apply_bc = [this](ScalarField& field) {
        // walls
        field.ApplyNeumannFreeWall(0, 1,   o-1, Wall::Left);
        field.ApplyNeumannFreeWall(d, o+1, l-2, Wall::Left);
        field.ApplyNeumannFreeWall(0, o2,  r-1, Wall::Left);

        field.ApplyNeumannFreeWall(LX-1,   1,   o-1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1, o+1, l-2, Wall::Right);
        field.ApplyNeumannFreeWall(LX-1,   o2,  r-1, Wall::Right);

        field.ApplyNeumannFreeWall(1,    0,    LX-2, Wall::Front);
        field.ApplyNeumannFreeWall(1,    o2-1, d-1,  Wall::Front);
        field.ApplyNeumannFreeWall(LX-d, o2-1, d-1,  Wall::Front);

        field.ApplyNeumannFreeWall(1,    LY-1, LX-2, Wall::Back);
        field.ApplyNeumannFreeWall(1,    o,    d-1,  Wall::Back);
        field.ApplyNeumannFreeWall(LX-d, o,    d-1,  Wall::Back);

        // concave corners
        field.ApplyNeumannFreeConcaveCorner(0, 0,    Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, o2-1, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, o,    Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(0, LY-1, Corner::LeftBack);

        field.ApplyNeumannFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, o,    Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, 0,    Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-1, o2-1, Corner::RightFront);

        //convex corners
        field.ApplyNeumannFreeConvexCorner(LX-d-1, o,    Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(d,      o,    Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1, o2-1, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(d,      o2-1, Corner::RightBack);
      };
      apply_bc(QQxx);
      apply_bc(QQyx);
      apply_bc(phi);
      break;
    }
    // channel that widens with free-slip or no-slip and von Neumann.
    // In principle as BC==7/BC==8 but without the lower reservoir.
    case 9:
    case 10:
    {
      auto apply_bc = [this](ScalarField& field) {
        // walls
        field.ApplyNeumannFreeWall(d, 1,   l-1, Wall::Left);
        field.ApplyNeumannFreeWall(0, l+1, r-1, Wall::Left);

        field.ApplyNeumannFreeWall(LX-d-1, 1,   l-1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-1,   l+1, r-1, Wall::Right);

        field.ApplyNeumannFreeWall(d+1,  0, LX-2-2*d, Wall::Front);
        field.ApplyNeumannFreeWall(1,    l, d-1,      Wall::Front);
        field.ApplyNeumannFreeWall(LX-d, l, d-1,      Wall::Front);

        field.ApplyNeumannFreeWall(1,    LY-1, LX-2, Wall::Back);

        // concave corners
        field.ApplyNeumannFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, l, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, LY-1, Corner::LeftBack);

        field.ApplyNeumannFreeConcaveCorner(LX-1,   LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1, 0,    Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-1,   l,    Corner::RightFront);

        //convex corners
        field.ApplyNeumannFreeConvexCorner(LX-d-1, l, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(d,      l, Corner::RightBack);
      };
      apply_bc(QQxx);
      apply_bc(QQyx);
      apply_bc(phi);
      break;
    }
    //channel that tightens and widens again with free-slip or no-slip and von Neumann
    // Similar to C==7/BC==8 , but periodic at the reservoirs
    case 11:
    case 12:
    {
      auto apply_bc = [this](ScalarField& field) {
        // walls
        field.ApplyNeumannFreeWall(d, o+1, l-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-d-1, o+1, l-2, Wall::Right);

        field.ApplyNeumannFreeWall(0,    0,    LX, Wall::Front);
        field.ApplyNeumannFreeWall(0,    o2-1, d ,  Wall::Front);
        field.ApplyNeumannFreeWall(LX-d, o2-1, d ,  Wall::Front);

        field.ApplyNeumannFreeWall(0,    LY-1, LX, Wall::Back);
        field.ApplyNeumannFreeWall(0,    o,    d ,  Wall::Back);
        field.ApplyNeumannFreeWall(LX-d, o,    d ,  Wall::Back);

        //convex corners
        field.ApplyNeumannFreeConvexCorner(LX-d-1, o,    Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(d,      o,    Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1, o2-1, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(d,      o2-1, Corner::RightBack);
      };
      apply_bc(QQxx);
      apply_bc(QQyx);
      apply_bc(phi);
      break;
    }
    // channel that tightens and widens again with von Neumann
    // and with rounded corners
    case 13:
    case 14:
    {
      auto apply_bc = [this](ScalarField& field) {
        // walls
        field.ApplyNeumannFreeWall(0, 1,   o-1, Wall::Left);
        field.ApplyNeumannFreeWall(d, o+6, l-12, Wall::Left);
        field.ApplyNeumannFreeWall(0, o2,  r-1, Wall::Left);

        field.ApplyNeumannFreeWall(LX-1,   1,   o-1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1, o+6, l-12, Wall::Right);
        field.ApplyNeumannFreeWall(LX-1,   o2,  r-1, Wall::Right);

        field.ApplyNeumannFreeWall(1,      0,    LX-2, Wall::Front);
        field.ApplyNeumannFreeWall(1,      o2-1, d-6,  Wall::Front);
        field.ApplyNeumannFreeWall(LX-d+5, o2-1, d-6,  Wall::Front);

        field.ApplyNeumannFreeWall(1,      LY-1, LX-2, Wall::Back);
        field.ApplyNeumannFreeWall(1,      o,    d-6,  Wall::Back);
        field.ApplyNeumannFreeWall(LX-d+5, o,    d-6,  Wall::Back);

        // concave corners
        field.ApplyNeumannFreeConcaveCorner(0, 0,    Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, o2-1, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, o,    Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(0, LY-1, Corner::LeftBack);

        field.ApplyNeumannFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, o,    Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, 0,    Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-1, o2-1, Corner::RightFront);

        //rounded convex corners
        //Lower left corner
        field.ApplyNeumannFreeConvexCorner(d-5,      o,     Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-4,      o+1,   Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-2,      o+2,   Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-1,      o+4,   Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d  ,      o+5,   Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(d-5,      o+1,   Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-4,      o+2,   Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-2,      o+4,   Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-1,      o+5,   Corner::LeftBack);
        field.ApplyNeumannFreeWall(d-3, o+2, 1,  Wall::Back);
        field.ApplyNeumannFreeWall(d-2, o+3, 1,  Wall::Left);

        //Lower right corner
        field.ApplyNeumannFreeConvexCorner(LX-d+4, o,     Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d+3, o+1,   Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d+1, o+2,   Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d  , o+4,   Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1, o+5,   Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d+4, o+1,   Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d+3, o+2,   Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d+1, o+4,   Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d  , o+5,   Corner::RightBack);
        field.ApplyNeumannFreeWall(LX-d+2, o+2, 1,  Wall::Back);
        field.ApplyNeumannFreeWall(LX-d+1, o+3, 1,  Wall::Right);

        //Upper left corner
        field.ApplyNeumannFreeConvexCorner(d-5, o2-1, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-4, o2-2, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-2, o2-3, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-1, o2-5, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d  , o2-6, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(d-5, o2-2, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-4, o2-3, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-2, o2-5, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-1, o2-6, Corner::LeftFront);
        field.ApplyNeumannFreeWall(d-3, o2-3, 1,  Wall::Front);
        field.ApplyNeumannFreeWall(d-2, o2-4, 1,  Wall::Left);

        //Upper right corner
        field.ApplyNeumannFreeConvexCorner(LX-d+4, o2-1, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d+3, o2-2, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d+1, o2-3, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d  , o2-5, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1, o2-6, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d+4, o2-2, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d+3, o2-3, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d+1, o2-5, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d  , o2-6, Corner::RightFront);
        field.ApplyNeumannFreeWall(LX-d+2, o2-3, 1,  Wall::Front);
        field.ApplyNeumannFreeWall(LX-d+1, o2-4, 1,  Wall::Right);

      };
      apply_bc(QQxx);
      apply_bc(QQyx);
      apply_bc(phi);
      break;
    }
    //other type of rounded corners
    case 15:
    case 16:
    {
      //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(0, 1, o-1, Wall::Left);
        field.ApplyNeumannFreeWall(d, o+16, l-32, Wall::Left);
        field.ApplyNeumannFreeWall(0, o2, r-1, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, o-1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1, o+16, l-32, Wall::Right);
        field.ApplyNeumannFreeWall(LX-1, o2, r-1, Wall::Right);
        field.ApplyNeumannFreeWall(1, 0, LX-2, Wall::Front);
        field.ApplyNeumannFreeWall(1, o2-1, d-14, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d+15, o2-1, d-14, Wall::Front);
        field.ApplyNeumannFreeWall(1, LY-1, LX-2, Wall::Back);
        field.ApplyNeumannFreeWall(1, o, d-14, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d+15, o, d-14, Wall::Back);
        field.ApplyNeumannFreeWall(d-2, o+8, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o+10, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o+11, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+13, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+14, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+15, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-15, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-14, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-13, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-11, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-10, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-8, o+2, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(d-12, o+0, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-9, o+1, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-7, o+2, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-6, o+3, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-5, o+4, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-4, o+5, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-3, o+6, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-2, o+7, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-1, o+9, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d+0, o+12, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(d-12, o+1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-9, o+2, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-7, o+3, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-6, o+4, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-5, o+5, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-4, o+6, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-3, o+7, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-2, o+9, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-1, o+12, Corner::LeftBack);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+15, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+14, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+13, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o+11, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o+10, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+2, o+8, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+8, o+2, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+10, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+11, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+13, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+14, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+15, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+0, o+12, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+1, o+9, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+2, o+7, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+3, o+6, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+4, o+5, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+5, o+4, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+6, o+3, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+7, o+2, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+9, o+1, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+12, o+0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+1, o+12, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+2, o+9, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+3, o+7, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+4, o+6, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+5, o+5, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+6, o+4, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+7, o+3, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+9, o+2, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+12, o+1, Corner::RightBack);
        field.ApplyNeumannFreeWall(d+0, o2-1-15, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o2-1-14, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o2-1-13, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o2-1-11, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o2-1-10, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-2, o2-1-8, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-8, o2-1-2, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-10, o2-1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-11, o2-1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-13, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-14, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-15, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeConvexCorner(d+0, o2-1-12, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-1, o2-1-9, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-2, o2-1-7, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-3, o2-1-6, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-4, o2-1-5, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-5, o2-1-4, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-6, o2-1-3, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-7, o2-1-2, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-9, o2-1-1, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-12, o2-1+0, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(d-1, o2-1-12, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-2, o2-1-9, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-3, o2-1-7, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-4, o2-1-6, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-5, o2-1-5, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-6, o2-1-4, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-7, o2-1-3, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-9, o2-1-2, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-12, o2-1-1, Corner::LeftFront);
        field.ApplyNeumannFreeWall(LX-d-1+2, o2-1-8, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o2-1-10, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o2-1-11, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o2-1-13, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o2-1-14, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o2-1-15, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+15, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+14, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+13, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+11, o2-1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+10, o2-1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+8, o2-1-2, 1, Wall::Front);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+12, o2-1+0, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+9, o2-1-1, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+7, o2-1-2, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+6, o2-1-3, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+5, o2-1-4, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+4, o2-1-5, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+3, o2-1-6, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+2, o2-1-7, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+1, o2-1-9, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+0, o2-1-12, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+12, o2-1-1, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+9, o2-1-2, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+7, o2-1-3, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+6, o2-1-4, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+5, o2-1-5, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+4, o2-1-6, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+3, o2-1-7, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+2, o2-1-9, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+1, o2-1-12, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(0, 0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, o2-1, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, o, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(0, LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, o, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, 0, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-1, o2-1, Corner::RightFront);
      };

      apply_bc(QQxx);
      apply_bc(QQyx);
      apply_bc(phi);
      break;
    }
    //other type of rounded corners with periodic reservoir
    case 17:
    case 18:
    case 19:
    {
      //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(d, o+16, l-32, Wall::Left);
        field.ApplyNeumannFreeWall(LX-d-1, o+16, l-32, Wall::Right);
        field.ApplyNeumannFreeWall(0, 0, LX, Wall::Front);
        field.ApplyNeumannFreeWall(0, o2-1, d-15, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d+15, o2-1, d-15, Wall::Front);
        field.ApplyNeumannFreeWall(0, LY-1, LX, Wall::Back);
        field.ApplyNeumannFreeWall(0, o, d-15, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d+15, o, d-15, Wall::Back);
        field.ApplyNeumannFreeWall(d-2, o+8, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o+10, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o+11, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+13, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+14, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+15, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-15, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-14, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-13, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-11, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-10, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-8, o+2, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(d-12, o+0, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-9, o+1, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-7, o+2, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-6, o+3, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-5, o+4, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-4, o+5, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-3, o+6, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-2, o+7, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-1, o+9, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d+0, o+12, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(d-12, o+1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-9, o+2, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-7, o+3, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-6, o+4, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-5, o+5, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-4, o+6, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-3, o+7, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-2, o+9, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-1, o+12, Corner::LeftBack);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+15, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+14, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+13, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o+11, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o+10, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+2, o+8, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+8, o+2, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+10, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+11, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+13, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+14, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+15, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+0, o+12, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+1, o+9, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+2, o+7, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+3, o+6, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+4, o+5, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+5, o+4, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+6, o+3, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+7, o+2, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+9, o+1, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+12, o+0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+1, o+12, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+2, o+9, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+3, o+7, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+4, o+6, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+5, o+5, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+6, o+4, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+7, o+3, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+9, o+2, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+12, o+1, Corner::RightBack);
        field.ApplyNeumannFreeWall(d+0, o2-1-15, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o2-1-14, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o2-1-13, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o2-1-11, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o2-1-10, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-2, o2-1-8, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-8, o2-1-2, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-10, o2-1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-11, o2-1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-13, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-14, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-15, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeConvexCorner(d+0, o2-1-12, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-1, o2-1-9, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-2, o2-1-7, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-3, o2-1-6, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-4, o2-1-5, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-5, o2-1-4, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-6, o2-1-3, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-7, o2-1-2, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-9, o2-1-1, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-12, o2-1+0, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(d-1, o2-1-12, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-2, o2-1-9, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-3, o2-1-7, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-4, o2-1-6, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-5, o2-1-5, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-6, o2-1-4, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-7, o2-1-3, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-9, o2-1-2, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-12, o2-1-1, Corner::LeftFront);
        field.ApplyNeumannFreeWall(LX-d-1+2, o2-1-8, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o2-1-10, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o2-1-11, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o2-1-13, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o2-1-14, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o2-1-15, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+15, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+14, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+13, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+11, o2-1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+10, o2-1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+8, o2-1-2, 1, Wall::Front);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+12, o2-1+0, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+9, o2-1-1, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+7, o2-1-2, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+6, o2-1-3, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+5, o2-1-4, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+4, o2-1-5, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+3, o2-1-6, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+2, o2-1-7, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+1, o2-1-9, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+0, o2-1-12, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+12, o2-1-1, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+9, o2-1-2, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+7, o2-1-3, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+6, o2-1-4, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+5, o2-1-5, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+4, o2-1-6, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+3, o2-1-7, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+2, o2-1-9, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+1, o2-1-12, Corner::RightFront);
      };

      apply_bc(QQxx);
      apply_bc(QQyx);
      apply_bc(phi);
      break;
    }
    //reservoir with chimney
    case 20:
    {
      //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(d, o+1, LY-o-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-d-1, o+1, LY-o-2, Wall::Right);
        field.ApplyNeumannFreeWall(0, 0, LX, Wall::Front);
        field.ApplyNeumannFreeWall(d+1, LY-1, LX-2*d-2, Wall::Back);
        //field.CopyDerivativeFreeWall(d+1, LY-1, LX-2*d-2, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d, o, d, Wall::Back);
        field.ApplyNeumannFreeWall(0, o, d, Wall::Back);
        field.ApplyNeumannFreeConcaveCorner(d, LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1, o, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(d, o, Corner::RightFront);
      };

      apply_bc(QQxx);
      apply_bc(QQyx);
      apply_bc(phi);
      break;
    }
    //long channel with inlet at bottom and outlet at top
    case 22:
    {
      //automatically generated code
      auto apply_phi_bc = [this](ScalarField& field) {
        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front,(1+noise*random_real())*conc);
        field.ApplyNeumannFreeWall(1, LY-1, LX-2, Wall::Back);
        field.ApplyNeumannFreeConcaveCorner(0, 0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, 0, Corner::RightFront);
      };

      double theta = angle + noise*M_PI*(random_real() - .5);

      auto apply_QQxx_bc = [this](ScalarField& field, double theta) {
        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front,cos(2*theta));
        field.ApplyNeumannFreeWall(1, LY-1, LX-2, Wall::Back);
        field.ApplyNeumannFreeConcaveCorner(0, 0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, 0, Corner::RightFront);
      };

      auto apply_QQyx_bc = [this](ScalarField& field, double theta) {
        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front,sin(2*theta));
        field.ApplyNeumannFreeWall(1, LY-1, LX-2, Wall::Back);
        field.ApplyNeumannFreeConcaveCorner(0, 0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, 0, Corner::RightFront);
      };

      apply_QQxx_bc(QQxx,theta);
      apply_QQyx_bc(QQyx,theta);
      apply_phi_bc(phi);
      break;
    }
    // box with pressure boundary on top
    case 23:
    {
      //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyNeumannFreeWall(1, 0, LX-2, Wall::Front);
        field.CopyDerivativeFreeWall(1, LY-1, LX-2, Wall::Back);
        field.ApplyNeumannFreeConcaveCorner(0, 0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, 0, Corner::RightFront);
      };
      apply_bc(QQxx);
      apply_bc(QQyx);
      apply_bc(phi);
      break;
    }
    // box with pressure boundary on top and diffusive phi inlet at the bottom
    case 24:
    {
      //automatically generated code
      //automatically generated code
      auto apply_phi_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front, conc);
        field.CopyDerivativeFreeWall(1, LY-1, LX-2, Wall::Back);
        field.ApplyNeumannFreeConcaveCorner(0, 0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, 0, Corner::RightFront);
      };

      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyNeumannFreeWall(1, 0, LX-2, Wall::Front);
        field.CopyDerivativeFreeWall(1, LY-1, LX-2, Wall::Back);
        field.ApplyNeumannFreeConcaveCorner(0, 0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, 0, Corner::RightFront);
      };
      apply_bc(QQxx);
      apply_bc(QQyx);
      apply_phi_bc(phi);
      break;
    }
    //reservoir with chimney, radius corners=4
    case 26:
    {
      //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(d, o+5, LY-o-6, Wall::Left);
        field.ApplyNeumannFreeWall(LX-d-1, o+5, LY-o-6, Wall::Right);
        field.ApplyNeumannFreeWall(0, 0, LX, Wall::Front);
        field.CopyDerivativeFreeWall(d+1, LY-1, LX-2*d-2, Wall::Back);
        field.ApplyNeumannFreeConcaveCorner(d, LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeWall(LX-d+4, o, d-4, Wall::Back);
        field.ApplyNeumannFreeWall(0, o, d-4, Wall::Back);
        field.ApplyNeumannFreeWall(d+0, o+4, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-4, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(d-3, o+0, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-2, o+1, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-1, o+2, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d+0, o+3, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(d-3, o+1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-2, o+2, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-1, o+3, Corner::LeftBack);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+4, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+4, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+0, o+3, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+1, o+2, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+2, o+1, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+3, o+0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+1, o+3, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+2, o+2, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+3, o+1, Corner::RightBack);
      };

      apply_bc(QQxx);
      apply_bc(QQyx);
      apply_bc(phi);
      break;
    }
    //reservoir with chimney, radius corners=7
    case 27:
    {
      //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(d, o+8, LY-o-9, Wall::Left);
        field.ApplyNeumannFreeWall(LX-d-1, o+8, LY-o-9, Wall::Right);
        field.ApplyNeumannFreeWall(0, 0, LX, Wall::Front);
        field.CopyDerivativeFreeWall(d+1, LY-1, LX-2*d-2, Wall::Back);
        field.ApplyNeumannFreeConcaveCorner(d, LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeWall(LX-d+7, o, d-7, Wall::Back);
        field.ApplyNeumannFreeWall(0, o, d-7, Wall::Back);
        field.ApplyNeumannFreeWall(d-1, o+4, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+6, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+7, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-7, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-6, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-4, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(d-5, o+0, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-3, o+1, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-2, o+2, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-1, o+3, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d+0, o+5, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(d-5, o+1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-3, o+2, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-2, o+3, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-1, o+5, Corner::LeftBack);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+7, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+6, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o+4, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+4, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+6, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+7, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+0, o+5, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+1, o+3, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+2, o+2, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+3, o+1, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+5, o+0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+1, o+5, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+2, o+3, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+3, o+2, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+5, o+1, Corner::RightBack);
      };

      apply_bc(QQxx);
      apply_bc(QQyx);
      apply_bc(phi);
      break;
    }
        //reservoir with chimney, radius corners=10
    case 28:
    {
      //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(d, o+11, LY-o-12, Wall::Left);
        field.ApplyNeumannFreeWall(LX-d-1, o+11, LY-o-12, Wall::Right);
        field.ApplyNeumannFreeWall(0, 0, LX, Wall::Front);
        field.CopyDerivativeFreeWall(d+1, LY-1, LX-2*d-2, Wall::Back);
        field.ApplyNeumannFreeConcaveCorner(d, LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeWall(LX-d+10, o, d-10, Wall::Back);
        field.ApplyNeumannFreeWall(0, o, d-10, Wall::Back);
        field.ApplyNeumannFreeWall(d-1, o+6, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+8, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+9, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+10, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-10, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-9, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-8, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-6, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(d-7, o+0, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-5, o+1, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-4, o+2, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-3, o+3, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-2, o+4, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-1, o+5, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d+0, o+7, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(d-7, o+1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-5, o+2, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-4, o+3, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-3, o+4, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-2, o+5, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-1, o+7, Corner::LeftBack);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+10, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+9, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+8, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o+6, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+6, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+8, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+9, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+10, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+0, o+7, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+1, o+5, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+2, o+4, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+3, o+3, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+4, o+2, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+5, o+1, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+7, o+0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+1, o+7, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+2, o+5, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+3, o+4, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+4, o+3, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+5, o+2, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+7, o+1, Corner::RightBack);
      };

      apply_bc(QQxx);
      apply_bc(QQyx);
      apply_bc(phi);
      break;
    }
    //reservoir with chimney, radius corners=13
    case 29:
    {
      //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(d, o+14, LY-o-15, Wall::Left);
        field.ApplyNeumannFreeWall(LX-d-1, o+14, LY-o-15, Wall::Right);
        field.ApplyNeumannFreeWall(0, 0, LX, Wall::Front);
        field.CopyDerivativeFreeWall(d+1, LY-1, LX-2*d-2, Wall::Back);
        field.ApplyNeumannFreeConcaveCorner(d, LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeWall(LX-d+13, o, d-13, Wall::Back);
        field.ApplyNeumannFreeWall(0, o, d-13, Wall::Back);
        field.ApplyNeumannFreeWall(d-1, o+8, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o+9, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+11, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+12, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+13, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-13, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-12, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-11, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-9, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-8, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(d-10, o+0, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-7, o+1, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-6, o+2, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-5, o+3, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-4, o+4, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-3, o+5, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-2, o+6, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-1, o+7, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d+0, o+10, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(d-10, o+1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-7, o+2, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-6, o+3, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-5, o+4, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-4, o+5, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-3, o+6, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-2, o+7, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-1, o+10, Corner::LeftBack);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+13, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+12, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+11, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o+9, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o+8, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+8, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+9, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+11, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+12, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+13, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+0, o+10, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+1, o+7, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+2, o+6, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+3, o+5, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+4, o+4, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+5, o+3, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+6, o+2, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+7, o+1, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+10, o+0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+1, o+10, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+2, o+7, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+3, o+6, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+4, o+5, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+5, o+4, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+6, o+3, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+7, o+2, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+10, o+1, Corner::RightBack);
      };

      apply_bc(QQxx);
      apply_bc(QQyx);
      apply_bc(phi);
      break;
    }
    //reservoir with chimney, radius corners=16
    case 30:
    {
      //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(d, o+17, LY-o-18, Wall::Left);
        field.ApplyNeumannFreeWall(LX-d-1, o+17, LY-o-18, Wall::Right);
        field.ApplyNeumannFreeWall(0, 0, LX, Wall::Front);
        field.CopyDerivativeFreeWall(d+1, LY-1, LX-2*d-2, Wall::Back);
        field.ApplyNeumannFreeConcaveCorner(d, LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeWall(LX-d+16, o, d-16, Wall::Back);
        field.ApplyNeumannFreeWall(0, o, d-16, Wall::Back);
        field.ApplyNeumannFreeWall(d-4, o+6, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-2, o+9, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o+11, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o+12, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+14, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+15, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+16, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-16, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-15, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-14, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-12, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-11, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-9, o+2, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-6, o+4, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(d-13, o+0, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-10, o+1, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-8, o+2, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-7, o+3, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-5, o+4, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-4, o+5, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-3, o+7, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-2, o+8, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-1, o+10, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d+0, o+13, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(d-13, o+1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-10, o+2, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-8, o+3, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-7, o+4, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-5, o+5, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-4, o+7, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-3, o+8, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-2, o+10, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-1, o+13, Corner::LeftBack);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+16, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+15, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+14, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o+12, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o+11, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+2, o+9, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+4, o+6, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+6, o+4, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+9, o+2, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+11, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+12, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+14, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+15, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+16, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+0, o+13, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+1, o+10, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+2, o+8, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+3, o+7, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+4, o+5, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+5, o+4, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+7, o+3, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+8, o+2, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+10, o+1, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+13, o+0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+1, o+13, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+2, o+10, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+3, o+8, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+4, o+7, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+5, o+5, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+7, o+4, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+8, o+3, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+10, o+2, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+13, o+1, Corner::RightBack);
      };

      apply_bc(QQxx);
      apply_bc(QQyx);
      apply_bc(phi);
      break;
    }
    //reservoir with chimney exit
    case 31:
    {
      //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeConcaveCorner(LX-d-1, 0, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.ApplyNeumannFreeWall(d+1, 0, LX-2*d-2, Wall::Front);
        field.ApplyNeumannFreeWall(d, 1, o, Wall::Left);
        field.ApplyNeumannFreeWall(LX-d-1, 1, o, Wall::Right);
        field.ApplyNeumannFreeConvexCorner(LX-d-1, o+1, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(d, o+1, Corner::RightBack);
        field.ApplyNeumannFreeWall(LX-d, o+1, d, Wall::Front);
        field.ApplyNeumannFreeWall(0, o+1, d, Wall::Front);
        field.CopyDerivativeFreeWall(0, LY-1, LX, Wall::Back);
      };

      apply_bc(QQxx);
      apply_bc(QQyx);
      apply_bc(phi);
      break;
    }

    //reservoir with chimney exit, radius corners=4
    case 32:
    {
      //automatically generated code
      auto apply_bc = [this](LBField& field) {

        field.ApplyNoSlipFreeWall(d+0, o+1-4, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-4, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeConvexCorner(d+0, o+1-3, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-1, o+1-2, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-2, o+1-1, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-3, o+1+0, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(d-1, o+1-3, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-2, o+1-2, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-3, o+1-1, Corner::LeftFront);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+1-4, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+4, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+3, o+1+0, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+2, o+1-1, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+1, o+1-2, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+0, o+1-3, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+3, o+1-1, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+2, o+1-2, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+1, o+1-3, Corner::RightFront);
        field.ApplyNoSlipFreeWall(d+1, 0, LX-2*d-2, Wall::Front);
        field.ApplyNoSlipFreeWall(d, 1, o-4, Wall::Left);
        field.ApplyNoSlipFreeWall(LX-d-1, 1, o-4, Wall::Right);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1, 0, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.ApplyNoSlipFreeWall(LX-d+4, o+1, d-4, Wall::Front);
        field.ApplyNoSlipFreeWall(0, o+1, d-4, Wall::Front);
        field.ApplyPressureOutletFreeWall(0, LY-1, LX, Wall::Back, n, ux, uy);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    //reservoir with chimney exit, radius corners=7
    case 33:
    {
      //automatically generated code
      auto apply_bc = [this](LBField& field) {

        field.ApplyNoSlipFreeWall(d+0, o+1-7, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+1-6, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-1, o+1-4, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-4, o+1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-6, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-7, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeConvexCorner(d+0, o+1-5, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-1, o+1-3, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-2, o+1-2, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-3, o+1-1, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-5, o+1+0, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(d-1, o+1-5, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-2, o+1-3, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-3, o+1-2, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-5, o+1-1, Corner::LeftFront);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o+1-4, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+1-6, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+1-7, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+7, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+6, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+4, o+1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+5, o+1+0, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+3, o+1-1, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+2, o+1-2, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+1, o+1-3, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+0, o+1-5, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+5, o+1-1, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+3, o+1-2, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+2, o+1-3, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+1, o+1-5, Corner::RightFront);
        field.ApplyNoSlipFreeWall(d+1, 0, LX-2*d-2, Wall::Front);
        field.ApplyNoSlipFreeWall(d, 1, o-7, Wall::Left);
        field.ApplyNoSlipFreeWall(LX-d-1, 1, o-7, Wall::Right);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1, 0, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.ApplyNoSlipFreeWall(LX-d+7, o+1, d-7, Wall::Front);
        field.ApplyNoSlipFreeWall(0, o+1, d-7, Wall::Front);
        field.ApplyPressureOutletFreeWall(0, LY-1, LX, Wall::Back, n, ux, uy);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    //reservoir with chimney exit, radius corners=10
    case 34:
    {
      //automatically generated code
      auto apply_bc = [this](LBField& field) {

        field.ApplyNoSlipFreeWall(d+0, o+1-10, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+1-9, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+1-8, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-1, o+1-6, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-6, o+1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-8, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-9, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-10, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeConvexCorner(d+0, o+1-7, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-1, o+1-5, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-2, o+1-4, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-3, o+1-3, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-4, o+1-2, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-5, o+1-1, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-7, o+1+0, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(d-1, o+1-7, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-2, o+1-5, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-3, o+1-4, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-4, o+1-3, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-5, o+1-2, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-7, o+1-1, Corner::LeftFront);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o+1-6, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+1-8, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+1-9, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+1-10, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+10, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+9, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+8, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+6, o+1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+7, o+1+0, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+5, o+1-1, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+4, o+1-2, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+3, o+1-3, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+2, o+1-4, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+1, o+1-5, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+0, o+1-7, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+7, o+1-1, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+5, o+1-2, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+4, o+1-3, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+3, o+1-4, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+2, o+1-5, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+1, o+1-7, Corner::RightFront);
        field.ApplyNoSlipFreeWall(d+1, 0, LX-2*d-2, Wall::Front);
        field.ApplyNoSlipFreeWall(d, 1, o-10, Wall::Left);
        field.ApplyNoSlipFreeWall(LX-d-1, 1, o-10, Wall::Right);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1, 0, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.ApplyNoSlipFreeWall(LX-d+10, o+1, d-10, Wall::Front);
        field.ApplyNoSlipFreeWall(0, o+1, d-10, Wall::Front);
        field.ApplyPressureOutletFreeWall(0, LY-1, LX, Wall::Back, n, ux, uy);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    //reservoir with chimney, radius corners=13
    case 35:
    {
      //automatically generated code
      auto apply_bc = [this](LBField& field) {

        field.ApplyNoSlipFreeWall(d+0, o+1-13, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+1-12, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+1-11, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-1, o+1-9, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-1, o+1-8, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-8, o+1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-9, o+1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-11, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-12, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-13, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeConvexCorner(d+0, o+1-10, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-1, o+1-7, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-2, o+1-6, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-3, o+1-5, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-4, o+1-4, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-5, o+1-3, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-6, o+1-2, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-7, o+1-1, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-10, o+1+0, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(d-1, o+1-10, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-2, o+1-7, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-3, o+1-6, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-4, o+1-5, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-5, o+1-4, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-6, o+1-3, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-7, o+1-2, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-10, o+1-1, Corner::LeftFront);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o+1-8, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o+1-9, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+1-11, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+1-12, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+1-13, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+13, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+12, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+11, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+9, o+1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+8, o+1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+10, o+1+0, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+7, o+1-1, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+6, o+1-2, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+5, o+1-3, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+4, o+1-4, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+3, o+1-5, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+2, o+1-6, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+1, o+1-7, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+0, o+1-10, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+10, o+1-1, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+7, o+1-2, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+6, o+1-3, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+5, o+1-4, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+4, o+1-5, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+3, o+1-6, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+2, o+1-7, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+1, o+1-10, Corner::RightFront);
        field.ApplyNoSlipFreeWall(d+1, 0, LX-2*d-2, Wall::Front);
        field.ApplyNoSlipFreeWall(d, 1, o-13, Wall::Left);
        field.ApplyNoSlipFreeWall(LX-d-1, 1, o-13, Wall::Right);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1, 0, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.ApplyNoSlipFreeWall(LX-d+13, o+1, d-13, Wall::Front);
        field.ApplyNoSlipFreeWall(0, o+1, d-13, Wall::Front);
        field.ApplyPressureOutletFreeWall(0, LY-1, LX, Wall::Back, n, ux, uy);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
    //reservoir with chimney exit, radius corners=16
    case 36:
    {
      //automatically generated code
      auto apply_bc = [this](LBField& field) {

        field.ApplyNoSlipFreeWall(d+0, o+1-16, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+1-15, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d+0, o+1-14, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-1, o+1-12, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-1, o+1-11, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-2, o+1-9, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-4, o+1-6, 1, Wall::Left);
        field.ApplyNoSlipFreeWall(d-6, o+1-4, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-9, o+1-2, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-11, o+1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-12, o+1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-14, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-15, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(d-16, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeConvexCorner(d+0, o+1-13, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-1, o+1-10, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-2, o+1-8, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-3, o+1-7, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-4, o+1-5, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-5, o+1-4, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-7, o+1-3, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-8, o+1-2, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-10, o+1-1, Corner::RightBack);
        field.ApplyNoSlipFreeConvexCorner(d-13, o+1+0, Corner::RightBack);
        field.ApplyNoSlipFreeConcaveCorner(d-1, o+1-13, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-2, o+1-10, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-3, o+1-8, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-4, o+1-7, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-5, o+1-5, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-7, o+1-4, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-8, o+1-3, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-10, o+1-2, Corner::LeftFront);
        field.ApplyNoSlipFreeConcaveCorner(d-13, o+1-1, Corner::LeftFront);
        field.ApplyNoSlipFreeWall(LX-d-1+4, o+1-6, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+2, o+1-9, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o+1-11, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+1, o+1-12, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+1-14, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+1-15, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+0, o+1-16, 1, Wall::Right);
        field.ApplyNoSlipFreeWall(LX-d-1+16, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+15, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+14, o+1+0, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+12, o+1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+11, o+1-1, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+9, o+1-2, 1, Wall::Front);
        field.ApplyNoSlipFreeWall(LX-d-1+6, o+1-4, 1, Wall::Front);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+13, o+1+0, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+10, o+1-1, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+8, o+1-2, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+7, o+1-3, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+5, o+1-4, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+4, o+1-5, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+3, o+1-7, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+2, o+1-8, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+1, o+1-10, Corner::LeftBack);
        field.ApplyNoSlipFreeConvexCorner(LX-d-1+0, o+1-13, Corner::LeftBack);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+13, o+1-1, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+10, o+1-2, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+8, o+1-3, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+7, o+1-4, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+5, o+1-5, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+4, o+1-7, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+3, o+1-8, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+2, o+1-10, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1+1, o+1-13, Corner::RightFront);
        field.ApplyNoSlipFreeWall(d+1, 0, LX-2*d-2, Wall::Front);
        field.ApplyNoSlipFreeWall(d, 1, o-16, Wall::Left);
        field.ApplyNoSlipFreeWall(LX-d-1, 1, o-16, Wall::Right);
        field.ApplyNoSlipFreeConcaveCorner(LX-d-1, 0, Corner::RightFront);
        field.ApplyNoSlipFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.ApplyNoSlipFreeWall(LX-d+16, o+1, d-16, Wall::Front);
        field.ApplyNoSlipFreeWall(0, o+1, d-16, Wall::Front);
        field.ApplyPressureOutletFreeWall(0, LY-1, LX, Wall::Back, n, ux, uy);
      };

      apply_bc(ff);
      apply_bc(fn);

      break;
    }
  }
}

void LyotropicFreeBoundary::BoundaryConditionsFields2()
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

      auto apply_uphi_bc = [this](ScalarField& field) {
        //walls
        field.ApplyDirichletFreeWall(0, 1, LY-2, Wall::Left,0);
        field.ApplyDirichletFreeWall(LX-1, 1, LY-2, Wall::Right,0);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front,0);
        field.ApplyDirichletFreeWall(1, LY-1, LX-2, Wall::Back,0);
        // corners (don't matter, but this depends on implementation of derivative)
        field.ApplyDirichletFreeConcaveCorner(LX-1,LY-1, Corner::RightBack,0);
        field.ApplyDirichletFreeConcaveCorner(LX-1,0, Corner::RightFront,0);
        field.ApplyDirichletFreeConcaveCorner(0,LY-1, Corner::LeftBack,0);
        field.ApplyDirichletFreeConcaveCorner(0,0, Corner::LeftFront,0);
        // walls
        field.ApplyDirichletFreeWall(SqX_left, SqY_front+1, SqWallHeight, Wall::Right,0);
        field.ApplyDirichletFreeWall(SqX_right, SqY_front+1 , SqWallHeight, Wall::Left,0);
        field.ApplyDirichletFreeWall(SqX_left+1, SqY_back , SqWallLength, Wall::Front,0);
        field.ApplyDirichletFreeWall(SqX_left+1,  SqY_front, SqWallLength, Wall::Back,0);
        // corners
        field.ApplyDirichletFreeConvexCorner(SqX_right , SqY_back ,Corner::RightBack,0);
        field.ApplyDirichletFreeConvexCorner(SqX_right ,SqY_front ,Corner::RightFront,0);
        field.ApplyDirichletFreeConvexCorner(SqX_left ,SqY_back,Corner::LeftBack,0);
        field.ApplyDirichletFreeConvexCorner(SqX_left ,SqY_front ,Corner::LeftFront,0);

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
      apply_uphi_bc(ux_phi);
      apply_uphi_bc(uy_phi);
      apply_bc(phi);
      apply_bc(MU);
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

      auto apply_uxphi_bc = [this](ScalarField& field) {
        //walls
        field.ApplyDirichletFreeWall(0, 1, LY-2, Wall::Left,0);
        field.ApplyDirichletFreeWall(LX-1, 1, LY-2, Wall::Right,0);
        field.ApplyNeumannFreeWall(1, 0, LX-2, Wall::Front);
        field.ApplyNeumannFreeWall(1, LY-1, LX-2, Wall::Back);
        // corners (don't matter, but this depends on implementation of derivative)
        field.ApplyDirichletFreeConcaveCorner(LX-1,LY-1, Corner::RightBack,0);
        field.ApplyDirichletFreeConcaveCorner(LX-1,0, Corner::RightFront,0);
        field.ApplyDirichletFreeConcaveCorner(0,LY-1, Corner::LeftBack,0);
        field.ApplyDirichletFreeConcaveCorner(0,0, Corner::LeftFront,0);
        // walls
        field.ApplyDirichletFreeWall(SqX_left, SqY_front+1, SqWallHeight, Wall::Right,0);
        field.ApplyDirichletFreeWall(SqX_right, SqY_front+1 , SqWallHeight, Wall::Left,0);
        field.ApplyNeumannFreeWall(SqX_left+1, SqY_back , SqWallLength, Wall::Front);
        field.ApplyNeumannFreeWall(SqX_left+1,  SqY_front, SqWallLength, Wall::Back);
        // corners
        field.ApplyDirichletFreeConvexCorner(SqX_right , SqY_back ,Corner::RightBack,0);
        field.ApplyDirichletFreeConvexCorner(SqX_right ,SqY_front ,Corner::RightFront,0);
        field.ApplyDirichletFreeConvexCorner(SqX_left ,SqY_back,Corner::LeftBack,0);
        field.ApplyDirichletFreeConvexCorner(SqX_left ,SqY_front ,Corner::LeftFront,0);

      };

      auto apply_uyphi_bc = [this](ScalarField& field) {
        //walls
        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front,0);
        field.ApplyDirichletFreeWall(1, LY-1, LX-2, Wall::Back,0);

        field.ApplyDirichletFreeConcaveCorner(LX-1,LY-1, Corner::RightBack,0);
        field.ApplyDirichletFreeConcaveCorner(LX-1,0, Corner::RightFront,0);
        field.ApplyDirichletFreeConcaveCorner(0,LY-1, Corner::LeftBack,0);
        field.ApplyDirichletFreeConcaveCorner(0,0, Corner::LeftFront,0);
        // walls
        field.ApplyNeumannFreeWall(SqX_left, SqY_front+1, SqWallHeight, Wall::Right);
        field.ApplyNeumannFreeWall(SqX_right, SqY_front+1 , SqWallHeight, Wall::Left);
        field.ApplyDirichletFreeWall(SqX_left+1, SqY_back , SqWallLength, Wall::Front,0);
        field.ApplyDirichletFreeWall(SqX_left+1,  SqY_front, SqWallLength, Wall::Back,0);
        // corners
        field.ApplyDirichletFreeConvexCorner(SqX_right , SqY_back ,Corner::RightBack,0);
        field.ApplyDirichletFreeConvexCorner(SqX_right ,SqY_front ,Corner::RightFront,0);
        field.ApplyDirichletFreeConvexCorner(SqX_left ,SqY_back,Corner::LeftBack,0);
        field.ApplyDirichletFreeConvexCorner(SqX_left ,SqY_front ,Corner::LeftFront,0);


      };

      auto apply_ux_bc = [this](ScalarField& field) {
        // walls
        field.CopyDerivativeFreeWall(0, 1, LY-2, Wall::Left);
        field.CopyDerivativeFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyNeumannFreeWall(1, 0, LX-2, Wall::Front);
        field.ApplyNeumannFreeWall(1, LY-1, LX-2, Wall::Back);

        field.CopyDerivativeFreeConcaveCorner(LX-1,LY-1, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-1,0, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(0,LY-1, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(0,0, Corner::LeftFront);
        // walls
        field.CopyDerivativeFreeWall(SqX_left, SqY_front+1, SqWallHeight, Wall::Right);
        field.CopyDerivativeFreeWall(SqX_right, SqY_front+1 , SqWallHeight, Wall::Left);
        field.ApplyNeumannFreeWall(SqX_left+1, SqY_front , SqWallLength, Wall::Back);
        field.ApplyNeumannFreeWall(SqX_left+1,  SqY_back, SqWallLength, Wall::Front);
        // corners
        field.CopyDerivativeFreeConvexCorner(SqX_right , SqY_back ,Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(SqX_right ,SqY_front ,Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(SqX_left ,SqY_back ,Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(SqX_left ,SqY_front ,Corner::LeftFront);


      };

      auto apply_uy_bc = [this](ScalarField& field) {
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

      apply_ux_bc(ux);
      apply_uy_bc(uy);
      apply_uxphi_bc(ux_phi);
      apply_uyphi_bc(uy_phi);
      apply_bc(phi);
      apply_bc(MU);
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

      auto apply_uphi_bc = [this](ScalarField& field) {
        //walls
        field.ApplyDirichletFreeWall(0, 1, LY-2, Wall::Left,0);
        field.ApplyDirichletFreeWall(LX-1, 1, LY-2, Wall::Right,0);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front,0);
        field.ApplyDirichletFreeWall(1, LY-1, LX-2, Wall::Back,0);
        // corners (don't matter, but this depends on implementation of derivative)
        field.ApplyDirichletFreeConcaveCorner(LX-1,LY-1, Corner::RightBack,0);
        field.ApplyDirichletFreeConcaveCorner(LX-1,0, Corner::RightFront,0);
        field.ApplyDirichletFreeConcaveCorner(0,LY-1, Corner::LeftBack,0);
        field.ApplyDirichletFreeConcaveCorner(0,0, Corner::LeftFront,0);
        // walls
        field.ApplyDirichletFreeWall(SqX_left, SqY_front+1, SqWallHeight, Wall::Right,0);
        field.ApplyDirichletFreeWall(SqX_right, SqY_front+1 , SqWallHeight, Wall::Left,0);
        field.ApplyDirichletFreeWall(SqX_left+1, SqY_back , SqWallLength, Wall::Front,0);
        field.ApplyDirichletFreeWall(SqX_left+1,  SqY_front, SqWallLength, Wall::Back,0);
        // corners
        field.ApplyDirichletFreeConvexCorner(SqX_right , SqY_back ,Corner::RightBack,0);
        field.ApplyDirichletFreeConvexCorner(SqX_right ,SqY_front ,Corner::RightFront,0);
        field.ApplyDirichletFreeConvexCorner(SqX_left ,SqY_back,Corner::LeftBack,0);
        field.ApplyDirichletFreeConvexCorner(SqX_left ,SqY_front ,Corner::LeftFront,0);

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
      apply_uphi_bc(ux_phi);
      apply_uphi_bc(uy_phi);
      apply_bc(phi);
      apply_bc(MU);
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
      };

      auto apply_uxphi_bc = [this](ScalarField& field) {
        //walls
        field.ApplyDirichletFreeWall(0, 1, LY-2, Wall::Left,0);
        field.ApplyDirichletFreeWall(LX-1, 1, LY-2, Wall::Right,0);
        field.ApplyNeumannFreeWall(1, 0, LX-2, Wall::Front);
        field.ApplyNeumannFreeWall(1, LY-1, LX-2, Wall::Back);
      };

      auto apply_uyphi_bc = [this](ScalarField& field) {
        //walls
        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front,0);
        field.ApplyDirichletFreeWall(1, LY-1, LX-2, Wall::Back,0);
      };

      auto apply_ux_bc = [this](ScalarField& field) {
        // walls
        field.CopyDerivativeFreeWall(0, 1, LY-2, Wall::Left);
        field.CopyDerivativeFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyNeumannFreeWall(1, 0, LX-2, Wall::Front);
        field.ApplyNeumannFreeWall(1, LY-1, LX-2, Wall::Back);
      };

      auto apply_uy_bc = [this](ScalarField& field) {
        // walls
        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.CopyDerivativeFreeWall(1, 0, LX-2, Wall::Front);
        field.CopyDerivativeFreeWall(1, LY-1, LX-2, Wall::Back);
      };

      apply_ux_bc(ux);
      apply_uy_bc(uy);
      apply_uxphi_bc(ux_phi);
      apply_uyphi_bc(uy_phi);
      apply_bc(phi);
      apply_bc(MU);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }
    //Free-slip box
    case 4:
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

      auto apply_uphi_bc = [this](ScalarField& field) {
        //walls
        field.ApplyDirichletFreeWall(0, 1, LY-2, Wall::Left,0);
        field.ApplyDirichletFreeWall(LX-1, 1, LY-2, Wall::Right,0);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front,0);
        field.ApplyDirichletFreeWall(1, LY-1, LX-2, Wall::Back,0);
        // corners (don't matter, but this depends on implementation of derivative)
        field.ApplyDirichletFreeConcaveCorner(LX-1,LY-1, Corner::RightBack,0);
        field.ApplyDirichletFreeConcaveCorner(LX-1,0, Corner::RightFront,0);
        field.ApplyDirichletFreeConcaveCorner(0,LY-1, Corner::LeftBack,0);
        field.ApplyDirichletFreeConcaveCorner(0,0, Corner::LeftFront,0);
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
      apply_uphi_bc(ux_phi);
      apply_uphi_bc(uy_phi);
      apply_bc(phi);
      apply_bc(MU);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }
    // von neumann channel
    case 2:
    case 5:
    {
      auto apply_bc = [this](ScalarField& field) {
        // pbc on the left and right walls
        field.ApplyNeumannFreeWall(0, 0, LX, Wall::Front);
        field.ApplyNeumannFreeWall(0, LY-1, LX, Wall::Back);
        // no corners
      };

      auto apply_uphi_bc = [this](ScalarField& field) {
        // pbc on the left and right walls
        field.ApplyDirichletFreeWall(0, 0, LX, Wall::Front, 0);
        field.ApplyDirichletFreeWall(0, LY-1, LX, Wall::Back, 0);
        // no corners
      };

      auto apply_u_bc = [this](ScalarField& field) {
        // pbc on the left and right walls
        field.CopyDerivativeFreeWall(0, 0, LX, Wall::Front);
        field.CopyDerivativeFreeWall(0, LY-1, LX, Wall::Back);
        // no corners
      };
      apply_u_bc(ux);
      apply_u_bc(uy);
      apply_uphi_bc(ux_phi);
      apply_uphi_bc(uy_phi);
      apply_bc(phi);
      apply_bc(MU);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }
    // von neumann channel in y-direction
    case 201:
    case 501:
    {
      auto apply_bc = [this](ScalarField& field) {
        // pbc on the left and right walls
        field.ApplyNeumannFreeWall(0, 0, LY, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 0, LY, Wall::Right);
        // no corners
      };

      auto apply_uphi_bc = [this](ScalarField& field) {
        // pbc on the left and right walls
        field.ApplyDirichletFreeWall(0, 0, LY, Wall::Left,0);
        field.ApplyDirichletFreeWall(LX-1, 0, LY, Wall::Right,0);
        // no corners
      };

      auto apply_u_bc = [this](ScalarField& field) {
        // pbc on the left and right walls
        field.CopyDerivativeFreeWall(0, 0, LY, Wall::Left);
        field.CopyDerivativeFreeWall(LX-1, 0, LY, Wall::Right);
        //field.ApplyDirichletFreeWall(0, 0, LY, Wall::Left,0);
        //field.ApplyDirichletFreeWall(0, LX-1, LY, Wall::Right,0);
        // no corners
      };
      apply_u_bc(ux);
      apply_u_bc(uy);
      apply_uphi_bc(ux_phi);
      apply_uphi_bc(uy_phi);
      apply_bc(phi);
      apply_bc(MU);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }
    // not implemented at the moment
    case 3:
    case 6:
    {
      break;
    }
    // channel that tightens and widens again with von Neumann
    case 7:
    case 8:
    {
      auto apply_bc = [this](ScalarField& field) {
        // walls
        field.ApplyNeumannFreeWall(0, 1,   o-1, Wall::Left);
        field.ApplyNeumannFreeWall(d, o+1, l-2, Wall::Left);
        field.ApplyNeumannFreeWall(0, o2,  r-1, Wall::Left);

        field.ApplyNeumannFreeWall(LX-1,   1,   o-1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1, o+1, l-2, Wall::Right);
        field.ApplyNeumannFreeWall(LX-1,   o2,  r-1, Wall::Right);

        field.ApplyNeumannFreeWall(1,    0,    LX-2, Wall::Front);
        field.ApplyNeumannFreeWall(1,    o2-1, d-1,  Wall::Front);
        field.ApplyNeumannFreeWall(LX-d, o2-1, d-1,  Wall::Front);

        field.ApplyNeumannFreeWall(1,    LY-1, LX-2, Wall::Back);
        field.ApplyNeumannFreeWall(1,    o,    d-1,  Wall::Back);
        field.ApplyNeumannFreeWall(LX-d, o,    d-1,  Wall::Back);

        // concave corners (don't matter, but this depends on implementation of derivative)
        field.ApplyNeumannFreeConcaveCorner(0, 0,    Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, o2-1, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, o,    Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(0, LY-1, Corner::LeftBack);

        field.ApplyNeumannFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, o,    Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, 0,    Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-1, o2-1, Corner::RightFront);

        //convex corners
        field.ApplyNeumannFreeConvexCorner(LX-d-1, o,    Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(d,      o,    Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1, o2-1, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(d,      o2-1, Corner::RightBack);
      };

      auto apply_uphi_bc = [this](ScalarField& field) {
        // walls
        field.ApplyDirichletFreeWall(0, 1,   o-1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d, o+1, l-2, Wall::Left, 0);
        field.ApplyDirichletFreeWall(0, o2,  r-1, Wall::Left, 0);

        field.ApplyDirichletFreeWall(LX-1,   1,   o-1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1, o+1, l-2, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-1,   o2,  r-1, Wall::Right, 0);

        field.ApplyDirichletFreeWall(1,    0,    LX-2, Wall::Front, 0);
        field.ApplyDirichletFreeWall(1,    o2-1, d-1,  Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d, o2-1, d-1,  Wall::Front, 0);

        field.ApplyDirichletFreeWall(1,    LY-1, LX-2, Wall::Back, 0);
        field.ApplyDirichletFreeWall(1,    o,    d-1,  Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d, o,    d-1,  Wall::Back, 0);

        // concave corners (don't matter, but this depends on implementation of derivative)
        field.ApplyDirichletFreeConcaveCorner(0, 0,    Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(0, o2-1, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(0, o,    Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(0, LY-1, Corner::LeftBack, 0);

        field.ApplyDirichletFreeConcaveCorner(LX-1, LY-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, o,    Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, 0,    Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, o2-1, Corner::RightFront, 0);

        //convex corners
        field.ApplyDirichletFreeConvexCorner(LX-d-1, o,    Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(d,      o,    Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1, o2-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(d,      o2-1, Corner::RightBack, 0);
      };

      auto apply_u_bc = [this](ScalarField& field) {
        // walls
        field.CopyDerivativeFreeWall(0, 1,   o-1, Wall::Left);
        field.CopyDerivativeFreeWall(d, o+1, l-2, Wall::Left);
        field.CopyDerivativeFreeWall(0, o2,  r-1, Wall::Left);

        field.CopyDerivativeFreeWall(LX-1,   1,   o-1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1, o+1, l-2, Wall::Right);
        field.CopyDerivativeFreeWall(LX-1,   o2,  r-1, Wall::Right);

        field.CopyDerivativeFreeWall(1,    0,    LX-2, Wall::Front);
        field.CopyDerivativeFreeWall(1,    o2-1, d-1,  Wall::Front);
        field.CopyDerivativeFreeWall(LX-d, o2-1, d-1,  Wall::Front);

        field.CopyDerivativeFreeWall(1,    LY-1, LX-2, Wall::Back);
        field.CopyDerivativeFreeWall(1,    o,    d-1,  Wall::Back);
        field.CopyDerivativeFreeWall(LX-d, o,    d-1,  Wall::Back);

        // concave corners (don't matter, but this depends on implementation of derivative)
        field.CopyDerivativeFreeConcaveCorner(0, 0,    Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(0, o2-1, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(0, o,    Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(0, LY-1, Corner::LeftBack);

        field.CopyDerivativeFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-1, o,    Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-1, 0,    Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-1, o2-1, Corner::RightFront);

        //convex corners
        field.CopyDerivativeFreeConvexCorner(LX-d-1, o,    Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(d,      o,    Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1, o2-1, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(d,      o2-1, Corner::RightBack);
      };
      apply_u_bc(ux);
      apply_u_bc(uy);
      apply_uphi_bc(ux_phi);
      apply_uphi_bc(uy_phi);
      apply_bc(phi);
      apply_bc(MU);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }

    // channel that widens with von Neumann.
    // In principle as BC==7/BC==8 but without the lower reservoir.
    case 9:
    case 10:
    {
      auto apply_bc = [this](ScalarField& field) {
        // walls
        field.ApplyNeumannFreeWall(d, 1,   l-1, Wall::Left);
        field.ApplyNeumannFreeWall(0, l+1, r-1, Wall::Left);

        field.ApplyNeumannFreeWall(LX-d-1, 1,   l-1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-1,   l+1, r-1, Wall::Right);

        field.ApplyNeumannFreeWall(d+1,  0, LX-2-2*d, Wall::Front);
        field.ApplyNeumannFreeWall(1,    l, d-1,      Wall::Front);
        field.ApplyNeumannFreeWall(LX-d, l, d-1,      Wall::Front);

        field.ApplyNeumannFreeWall(1,    LY-1, LX-2, Wall::Back);

        // concave corners (don't matter, but this depends on implementation of derivative)
        field.ApplyNeumannFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, l, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, LY-1, Corner::LeftBack);

        field.ApplyNeumannFreeConcaveCorner(LX-1,   LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1, 0,    Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-1,   l,    Corner::RightFront);

        //convex corners
        field.ApplyNeumannFreeConvexCorner(LX-d-1, l, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(d,      l, Corner::RightBack);
      };

      auto apply_uphi_bc = [this](ScalarField& field) { //velocity normal to boundary zero, parallel Neumann
        // walls
        field.ApplyDirichletFreeWall(d, 1,   l-1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(0, l+1, r-1, Wall::Left, 0);

        field.ApplyDirichletFreeWall(LX-d-1, 1,   l-1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-1,   l+1, r-1, Wall::Right, 0);

        field.ApplyDirichletFreeWall(d+1,  0, LX-2-2*d, Wall::Front, 0);
        field.ApplyDirichletFreeWall(1,    l, d-1,      Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d, l, d-1,      Wall::Front, 0);

        field.ApplyDirichletFreeWall(1,    LY-1, LX-2, Wall::Back, 0);

        // concave corners (don't matter, but this depends on implementation of derivative)
        field.ApplyDirichletFreeConcaveCorner(d, 0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(0, l, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(0, LY-1, Corner::LeftBack, 0);

        field.ApplyDirichletFreeConcaveCorner(LX-1,   LY-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1, 0,    Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1,   l,    Corner::RightFront, 0);

        //convex corners (one has to make a choice here, we choose dirichlet)
        field.ApplyDirichletFreeConvexCorner(LX-d-1, l, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(d,      l, Corner::RightBack, 0);
      };

      auto apply_u_bc = [this](ScalarField& field) { //velocity normal to boundary zero, parallel Neumann
        // walls
        field.CopyDerivativeFreeWall(d, 1,   l-1, Wall::Left);
        field.CopyDerivativeFreeWall(0, l+1, r-1, Wall::Left);

        field.CopyDerivativeFreeWall(LX-d-1, 1,   l-1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-1,   l+1, r-1, Wall::Right);

        field.CopyDerivativeFreeWall(d+1,  0, LX-2-2*d, Wall::Front);
        field.CopyDerivativeFreeWall(1,    l, d-1,      Wall::Front);
        field.CopyDerivativeFreeWall(LX-d, l, d-1,      Wall::Front);

        field.CopyDerivativeFreeWall(1,    LY-1, LX-2, Wall::Back);

        // concave corners (don't matter, but this depends on implementation of derivative)
        field.CopyDerivativeFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(0, l, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(0, LY-1, Corner::LeftBack);

        field.CopyDerivativeFreeConcaveCorner(LX-1,   LY-1, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1, 0,    Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-1,   l,    Corner::RightFront);

        //convex corners (one has to make a choice here, we choose dirichlet)
        field.CopyDerivativeFreeConvexCorner(LX-d-1, l, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(d,      l, Corner::RightBack);
      };

      apply_u_bc(ux);
      apply_u_bc(uy);
      apply_uphi_bc(ux_phi);
      apply_uphi_bc(uy_phi);
      apply_bc(phi);
      apply_bc(MU);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }
    // channel that tightens and widens again with von Neumann
    // Similar to C==7/BC==8 , but periodic at the reservoirs
    case 11:
    case 12:
    {
      auto apply_bc = [this](ScalarField& field) {
        // walls
        field.ApplyNeumannFreeWall(d, o+1, l-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-d-1, o+1, l-2, Wall::Right);

        field.ApplyNeumannFreeWall(0,    0,    LX, Wall::Front);
        field.ApplyNeumannFreeWall(0,    o2-1, d ,  Wall::Front);
        field.ApplyNeumannFreeWall(LX-d, o2-1, d ,  Wall::Front);

        field.ApplyNeumannFreeWall(0,    LY-1, LX, Wall::Back);
        field.ApplyNeumannFreeWall(0,    o,    d ,  Wall::Back);
        field.ApplyNeumannFreeWall(LX-d, o,    d ,  Wall::Back);

        //convex corners
        field.ApplyNeumannFreeConvexCorner(LX-d-1, o,    Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(d,      o,    Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1, o2-1, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(d,      o2-1, Corner::RightBack);
      };

      auto apply_uphi_bc = [this](ScalarField& field) {
        // walls
        field.ApplyDirichletFreeWall(d, o+1, l-2, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-d-1, o+1, l-2, Wall::Right, 0);

        field.ApplyDirichletFreeWall(0,    0,    LX, Wall::Front, 0);
        field.ApplyDirichletFreeWall(0,    o2-1, d ,  Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d, o2-1, d ,  Wall::Front, 0);

        field.ApplyDirichletFreeWall(0,    LY-1, LX, Wall::Back, 0);
        field.ApplyDirichletFreeWall(0,    o,    d ,  Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d, o,    d ,  Wall::Back, 0);

        //convex corners
        field.ApplyDirichletFreeConvexCorner(LX-d-1, o,    Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(d,      o,    Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1, o2-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(d,      o2-1, Corner::RightBack, 0);
      };

      auto apply_u_bc = [this](ScalarField& field) {
        // walls
        field.CopyDerivativeFreeWall(d, o+1, l-2, Wall::Left);
        field.CopyDerivativeFreeWall(LX-d-1, o+1, l-2, Wall::Right);

        field.CopyDerivativeFreeWall(0,    0,    LX, Wall::Front);
        field.CopyDerivativeFreeWall(0,    o2-1, d,  Wall::Front);
        field.CopyDerivativeFreeWall(LX-d, o2-1, d,  Wall::Front);

        field.CopyDerivativeFreeWall(0,    LY-1, LX, Wall::Back);
        field.CopyDerivativeFreeWall(0,    o,    d,  Wall::Back);
        field.CopyDerivativeFreeWall(LX-d, o,    d,  Wall::Back);

        //convex corners
        field.CopyDerivativeFreeConvexCorner(LX-d-1, o,    Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(d,      o,    Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1, o2-1, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(d,      o2-1, Corner::RightBack);
      };

      apply_u_bc(ux);
      apply_u_bc(uy);
      apply_uphi_bc(ux_phi);
      apply_uphi_bc(uy_phi);
      apply_bc(phi);
      apply_bc(MU);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }
    // channel that tightens and widens again with von Neumann
    // and with rounded corners
    case 13:
    case 14:
    {
      auto apply_bc = [this](ScalarField& field) {
        // walls
        field.ApplyNeumannFreeWall(0, 1,   o-1, Wall::Left);
        field.ApplyNeumannFreeWall(d, o+6, l-12, Wall::Left);
        field.ApplyNeumannFreeWall(0, o2,  r-1, Wall::Left);

        field.ApplyNeumannFreeWall(LX-1,   1,   o-1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1, o+6, l-12, Wall::Right);
        field.ApplyNeumannFreeWall(LX-1,   o2,  r-1, Wall::Right);

        field.ApplyNeumannFreeWall(1,      0,    LX-2, Wall::Front);
        field.ApplyNeumannFreeWall(1,      o2-1, d-6,  Wall::Front);
        field.ApplyNeumannFreeWall(LX-d+5, o2-1, d-6,  Wall::Front);

        field.ApplyNeumannFreeWall(1,      LY-1, LX-2, Wall::Back);
        field.ApplyNeumannFreeWall(1,      o,    d-6,  Wall::Back);
        field.ApplyNeumannFreeWall(LX-d+5, o,    d-6,  Wall::Back);

        // concave corners
        field.ApplyNeumannFreeConcaveCorner(0, 0,    Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, o2-1, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, o,    Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(0, LY-1, Corner::LeftBack);

        field.ApplyNeumannFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, o,    Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, 0,    Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-1, o2-1, Corner::RightFront);

        //rounded convex corners
        //Lower left corner
        field.ApplyNeumannFreeConvexCorner(d-5,      o,     Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-4,      o+1,   Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-2,      o+2,   Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-1,      o+4,   Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d  ,      o+5,   Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(d-5,      o+1,   Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-4,      o+2,   Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-2,      o+4,   Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-1,      o+5,   Corner::LeftBack);
        field.ApplyNeumannFreeWall(d-3, o+2, 1,  Wall::Back);
        field.ApplyNeumannFreeWall(d-2, o+3, 1,  Wall::Left);

        //Lower right corner
        field.ApplyNeumannFreeConvexCorner(LX-d+4, o,     Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d+3, o+1,   Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d+1, o+2,   Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d  , o+4,   Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1, o+5,   Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d+4, o+1,   Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d+3, o+2,   Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d+1, o+4,   Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d  , o+5,   Corner::RightBack);
        field.ApplyNeumannFreeWall(LX-d+2, o+2, 1,  Wall::Back);
        field.ApplyNeumannFreeWall(LX-d+1, o+3, 1,  Wall::Right);

        //Upper left corner
        field.ApplyNeumannFreeConvexCorner(d-5, o2-1, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-4, o2-2, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-2, o2-3, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-1, o2-5, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d  , o2-6, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(d-5, o2-2, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-4, o2-3, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-2, o2-5, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-1, o2-6, Corner::LeftFront);
        field.ApplyNeumannFreeWall(d-3, o2-3, 1,  Wall::Front);
        field.ApplyNeumannFreeWall(d-2, o2-4, 1,  Wall::Left);

        //Upper right corner
        field.ApplyNeumannFreeConvexCorner(LX-d+4, o2-1, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d+3, o2-2, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d+1, o2-3, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d  , o2-5, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1, o2-6, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d+4, o2-2, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d+3, o2-3, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d+1, o2-5, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d  , o2-6, Corner::RightFront);
        field.ApplyNeumannFreeWall(LX-d+2, o2-3, 1,  Wall::Front);
        field.ApplyNeumannFreeWall(LX-d+1, o2-4, 1,  Wall::Right);
      };

      auto apply_uphi_bc = [this](ScalarField& field) {
        // walls
        field.ApplyDirichletFreeWall(0, 1,   o-1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d, o+6, l-12, Wall::Left, 0);
        field.ApplyDirichletFreeWall(0, o2,  r-1, Wall::Left, 0);

        field.ApplyDirichletFreeWall(LX-1,   1,   o-1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1, o+6, l-12, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-1,   o2,  r-1, Wall::Right, 0);

        field.ApplyDirichletFreeWall(1,      0,    LX-2, Wall::Front, 0);
        field.ApplyDirichletFreeWall(1,      o2-1, d-6,  Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d+5, o2-1, d-6,  Wall::Front, 0);

        field.ApplyDirichletFreeWall(1,      LY-1, LX-2, Wall::Back, 0);
        field.ApplyDirichletFreeWall(1,      o,    d-6,  Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d+5, o,    d-6,  Wall::Back, 0);

        // concave corners
        field.ApplyDirichletFreeConcaveCorner(0, 0,    Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(0, o2-1, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(0, o,    Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(0, LY-1, Corner::LeftBack, 0);

        field.ApplyDirichletFreeConcaveCorner(LX-1, LY-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, o,    Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, 0,    Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, o2-1, Corner::RightFront, 0);

        //rounded convex corners
        //Lower left corner
        field.ApplyDirichletFreeConvexCorner(d-5,      o,     Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-4,      o+1,   Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-2,      o+2,   Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-1,      o+4,   Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d  ,      o+5,   Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-5,      o+1,   Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-4,      o+2,   Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-2,      o+4,   Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-1,      o+5,   Corner::LeftBack, 0);
        field.ApplyDirichletFreeWall(d-3, o+2, 1,  Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-2, o+3, 1,  Wall::Left, 0);

        //Lower right corner
        field.ApplyDirichletFreeConvexCorner(LX-d+4, o,     Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d+3, o+1,   Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d+1, o+2,   Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d  , o+4,   Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1, o+5,   Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d+4, o+1,   Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d+3, o+2,   Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d+1, o+4,   Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d  , o+5,   Corner::RightBack, 0);
        field.ApplyDirichletFreeWall(LX-d+2, o+2, 1,  Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d+1, o+3, 1,  Wall::Right, 0);

        //Upper left corner
        field.ApplyDirichletFreeConvexCorner(d-5, o2-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-4, o2-2, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-2, o2-3, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-1, o2-5, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d  , o2-6, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-5, o2-2, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-4, o2-3, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-2, o2-5, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-1, o2-6, Corner::LeftFront, 0);
        field.ApplyDirichletFreeWall(d-3, o2-3, 1,  Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-2, o2-4, 1,  Wall::Left, 0);

        //Upper right corner
        field.ApplyDirichletFreeConvexCorner(LX-d+4, o2-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d+3, o2-2, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d+1, o2-3, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d  , o2-5, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1, o2-6, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d+4, o2-2, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d+3, o2-3, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d+1, o2-5, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d  , o2-6, Corner::RightFront, 0);
        field.ApplyDirichletFreeWall(LX-d+2, o2-3, 1,  Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d+1, o2-4, 1,  Wall::Right, 0);
      };

      auto apply_u_bc = [this](ScalarField& field) {
        // walls
        field.CopyDerivativeFreeWall(0, 1,   o-1, Wall::Left);
        field.CopyDerivativeFreeWall(d, o+6, l-12, Wall::Left);
        field.CopyDerivativeFreeWall(0, o2,  r-1, Wall::Left);

        field.CopyDerivativeFreeWall(LX-1,   1,   o-1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1, o+6, l-12, Wall::Right);
        field.CopyDerivativeFreeWall(LX-1,   o2,  r-1, Wall::Right);

        field.CopyDerivativeFreeWall(1,      0,    LX-2, Wall::Front);
        field.CopyDerivativeFreeWall(1,      o2-1, d-6,  Wall::Front);
        field.CopyDerivativeFreeWall(LX-d+5, o2-1, d-6,  Wall::Front);

        field.CopyDerivativeFreeWall(1,      LY-1, LX-2, Wall::Back);
        field.CopyDerivativeFreeWall(1,      o,    d-6,  Wall::Back);
        field.CopyDerivativeFreeWall(LX-d+5, o,    d-6,  Wall::Back);

        // concave corners
        field.CopyDerivativeFreeConcaveCorner(0, 0,    Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(0, o2-1, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(0, o,    Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(0, LY-1, Corner::LeftBack);

        field.CopyDerivativeFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-1, o,    Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-1, 0,    Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-1, o2-1, Corner::RightFront);

        //rounded convex corners
        //Lower left corner
        field.CopyDerivativeFreeConvexCorner(d-5,      o,     Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-4,      o+1,   Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-2,      o+2,   Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-1,      o+4,   Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d  ,      o+5,   Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(d-5,      o+1,   Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-4,      o+2,   Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-2,      o+4,   Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-1,      o+5,   Corner::LeftBack);
        field.CopyDerivativeFreeWall(d-3, o+2, 1,  Wall::Back);
        field.CopyDerivativeFreeWall(d-2, o+3, 1,  Wall::Left);

        //Lower right corner
        field.CopyDerivativeFreeConvexCorner(LX-d+4, o,     Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d+3, o+1,   Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d+1, o+2,   Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d  , o+4,   Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1, o+5,   Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d+4, o+1,   Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d+3, o+2,   Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d+1, o+4,   Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d  , o+5,   Corner::RightBack);
        field.CopyDerivativeFreeWall(LX-d+2, o+2, 1,  Wall::Back);
        field.CopyDerivativeFreeWall(LX-d+1, o+3, 1,  Wall::Right);

        //Upper left corner
        field.CopyDerivativeFreeConvexCorner(d-5, o2-1, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-4, o2-2, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-2, o2-3, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-1, o2-5, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d  , o2-6, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(d-5, o2-2, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-4, o2-3, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-2, o2-5, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-1, o2-6, Corner::LeftFront);
        field.CopyDerivativeFreeWall(d-3, o2-3, 1,  Wall::Front);
        field.CopyDerivativeFreeWall(d-2, o2-4, 1,  Wall::Left);

        //Upper right corner
        field.CopyDerivativeFreeConvexCorner(LX-d+4, o2-1, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d+3, o2-2, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d+1, o2-3, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d  , o2-5, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1, o2-6, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d+4, o2-2, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d+3, o2-3, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d+1, o2-5, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d  , o2-6, Corner::RightFront);
        field.CopyDerivativeFreeWall(LX-d+2, o2-3, 1,  Wall::Front);
        field.CopyDerivativeFreeWall(LX-d+1, o2-4, 1,  Wall::Right);
      };

      apply_u_bc(ux);
      apply_u_bc(uy);
      apply_uphi_bc(ux_phi);
      apply_uphi_bc(uy_phi);
      apply_bc(phi);
      apply_bc(MU);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }
    //other type of rounded corners
    case 15:
    case 16:
    {
      //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(0, 1, o-1, Wall::Left);
        field.ApplyNeumannFreeWall(d, o+16, l-32, Wall::Left);
        field.ApplyNeumannFreeWall(0, o2, r-1, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, o-1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1, o+16, l-32, Wall::Right);
        field.ApplyNeumannFreeWall(LX-1, o2, r-1, Wall::Right);
        field.ApplyNeumannFreeWall(1, 0, LX-2, Wall::Front);
        field.ApplyNeumannFreeWall(1, o2-1, d-14, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d+15, o2-1, d-14, Wall::Front);
        field.ApplyNeumannFreeWall(1, LY-1, LX-2, Wall::Back);
        field.ApplyNeumannFreeWall(1, o, d-14, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d+15, o, d-14, Wall::Back);
        field.ApplyNeumannFreeWall(d-2, o+8, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o+10, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o+11, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+13, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+14, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+15, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-15, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-14, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-13, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-11, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-10, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-8, o+2, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(d-12, o+0, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-9, o+1, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-7, o+2, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-6, o+3, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-5, o+4, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-4, o+5, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-3, o+6, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-2, o+7, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-1, o+9, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d+0, o+12, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(d-12, o+1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-9, o+2, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-7, o+3, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-6, o+4, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-5, o+5, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-4, o+6, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-3, o+7, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-2, o+9, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-1, o+12, Corner::LeftBack);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+15, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+14, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+13, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o+11, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o+10, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+2, o+8, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+8, o+2, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+10, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+11, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+13, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+14, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+15, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+0, o+12, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+1, o+9, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+2, o+7, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+3, o+6, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+4, o+5, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+5, o+4, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+6, o+3, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+7, o+2, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+9, o+1, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+12, o+0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+1, o+12, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+2, o+9, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+3, o+7, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+4, o+6, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+5, o+5, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+6, o+4, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+7, o+3, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+9, o+2, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+12, o+1, Corner::RightBack);
        field.ApplyNeumannFreeWall(d+0, o2-1-15, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o2-1-14, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o2-1-13, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o2-1-11, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o2-1-10, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-2, o2-1-8, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-8, o2-1-2, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-10, o2-1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-11, o2-1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-13, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-14, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-15, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeConvexCorner(d+0, o2-1-12, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-1, o2-1-9, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-2, o2-1-7, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-3, o2-1-6, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-4, o2-1-5, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-5, o2-1-4, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-6, o2-1-3, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-7, o2-1-2, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-9, o2-1-1, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-12, o2-1+0, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(d-1, o2-1-12, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-2, o2-1-9, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-3, o2-1-7, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-4, o2-1-6, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-5, o2-1-5, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-6, o2-1-4, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-7, o2-1-3, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-9, o2-1-2, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-12, o2-1-1, Corner::LeftFront);
        field.ApplyNeumannFreeWall(LX-d-1+2, o2-1-8, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o2-1-10, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o2-1-11, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o2-1-13, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o2-1-14, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o2-1-15, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+15, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+14, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+13, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+11, o2-1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+10, o2-1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+8, o2-1-2, 1, Wall::Front);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+12, o2-1+0, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+9, o2-1-1, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+7, o2-1-2, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+6, o2-1-3, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+5, o2-1-4, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+4, o2-1-5, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+3, o2-1-6, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+2, o2-1-7, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+1, o2-1-9, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+0, o2-1-12, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+12, o2-1-1, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+9, o2-1-2, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+7, o2-1-3, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+6, o2-1-4, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+5, o2-1-5, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+4, o2-1-6, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+3, o2-1-7, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+2, o2-1-9, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+1, o2-1-12, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(0, 0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, o2-1, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, o, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(0, LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, o, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, 0, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-1, o2-1, Corner::RightFront);
      };

      auto apply_uphi_bc = [this](ScalarField& field) {

        field.ApplyDirichletFreeWall(0, 1, o-1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d, o+16, l-32, Wall::Left, 0);
        field.ApplyDirichletFreeWall(0, o2, r-1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-1, 1, o-1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1, o+16, l-32, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-1, o2, r-1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front, 0);
        field.ApplyDirichletFreeWall(1, o2-1, d-14, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d+15, o2-1, d-14, Wall::Front, 0);
        field.ApplyDirichletFreeWall(1, LY-1, LX-2, Wall::Back, 0);
        field.ApplyDirichletFreeWall(1, o, d-14, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d+15, o, d-14, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-2, o+8, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-1, o+10, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-1, o+11, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o+13, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o+14, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o+15, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-15, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-14, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-13, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-11, o+1, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-10, o+1, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-8, o+2, 1, Wall::Back, 0);
        field.ApplyDirichletFreeConvexCorner(d-12, o+0, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-9, o+1, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-7, o+2, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-6, o+3, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-5, o+4, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-4, o+5, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-3, o+6, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-2, o+7, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-1, o+9, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d+0, o+12, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-12, o+1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-9, o+2, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-7, o+3, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-6, o+4, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-5, o+5, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-4, o+6, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-3, o+7, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-2, o+9, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-1, o+12, Corner::LeftBack, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+15, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+14, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+13, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+1, o+11, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+1, o+10, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+2, o+8, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+8, o+2, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+10, o+1, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+11, o+1, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+13, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+14, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+15, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+0, o+12, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+1, o+9, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+2, o+7, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+3, o+6, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+4, o+5, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+5, o+4, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+6, o+3, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+7, o+2, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+9, o+1, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+12, o+0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+1, o+12, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+2, o+9, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+3, o+7, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+4, o+6, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+5, o+5, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+6, o+4, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+7, o+3, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+9, o+2, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+12, o+1, Corner::RightBack, 0);
        field.ApplyDirichletFreeWall(d+0, o2-1-15, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o2-1-14, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o2-1-13, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-1, o2-1-11, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-1, o2-1-10, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-2, o2-1-8, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-8, o2-1-2, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-10, o2-1-1, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-11, o2-1-1, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-13, o2-1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-14, o2-1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-15, o2-1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeConvexCorner(d+0, o2-1-12, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-1, o2-1-9, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-2, o2-1-7, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-3, o2-1-6, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-4, o2-1-5, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-5, o2-1-4, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-6, o2-1-3, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-7, o2-1-2, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-9, o2-1-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-12, o2-1+0, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-1, o2-1-12, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-2, o2-1-9, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-3, o2-1-7, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-4, o2-1-6, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-5, o2-1-5, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-6, o2-1-4, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-7, o2-1-3, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-9, o2-1-2, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-12, o2-1-1, Corner::LeftFront, 0);
        field.ApplyDirichletFreeWall(LX-d-1+2, o2-1-8, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+1, o2-1-10, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+1, o2-1-11, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o2-1-13, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o2-1-14, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o2-1-15, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+15, o2-1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+14, o2-1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+13, o2-1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+11, o2-1-1, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+10, o2-1-1, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+8, o2-1-2, 1, Wall::Front, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+12, o2-1+0, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+9, o2-1-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+7, o2-1-2, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+6, o2-1-3, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+5, o2-1-4, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+4, o2-1-5, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+3, o2-1-6, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+2, o2-1-7, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+1, o2-1-9, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+0, o2-1-12, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+12, o2-1-1, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+9, o2-1-2, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+7, o2-1-3, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+6, o2-1-4, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+5, o2-1-5, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+4, o2-1-6, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+3, o2-1-7, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+2, o2-1-9, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+1, o2-1-12, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(0, 0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(0, o2-1, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(0, o, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(0, LY-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, LY-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, o, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, 0, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, o2-1, Corner::RightFront, 0);
      };

      auto apply_u_bc = [this](ScalarField& field) {

        field.CopyDerivativeFreeWall(0, 1, o-1, Wall::Left);
        field.CopyDerivativeFreeWall(d, o+16, l-32, Wall::Left);
        field.CopyDerivativeFreeWall(0, o2, r-1, Wall::Left);
        field.CopyDerivativeFreeWall(LX-1, 1, o-1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1, o+16, l-32, Wall::Right);
        field.CopyDerivativeFreeWall(LX-1, o2, r-1, Wall::Right);
        field.CopyDerivativeFreeWall(1, 0, LX-2, Wall::Front);
        field.CopyDerivativeFreeWall(1, o2-1, d-14, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d+15, o2-1, d-14, Wall::Front);
        field.CopyDerivativeFreeWall(1, LY-1, LX-2, Wall::Back);
        field.CopyDerivativeFreeWall(1, o, d-14, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d+15, o, d-14, Wall::Back);
        field.CopyDerivativeFreeWall(d-2, o+8, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-1, o+10, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-1, o+11, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o+13, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o+14, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o+15, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-15, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-14, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-13, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-11, o+1, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-10, o+1, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-8, o+2, 1, Wall::Back);
        field.CopyDerivativeFreeConvexCorner(d-12, o+0, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-9, o+1, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-7, o+2, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-6, o+3, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-5, o+4, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-4, o+5, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-3, o+6, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-2, o+7, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-1, o+9, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d+0, o+12, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(d-12, o+1, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-9, o+2, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-7, o+3, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-6, o+4, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-5, o+5, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-4, o+6, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-3, o+7, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-2, o+9, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-1, o+12, Corner::LeftBack);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+15, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+14, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+13, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+1, o+11, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+1, o+10, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+2, o+8, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+8, o+2, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+10, o+1, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+11, o+1, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+13, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+14, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+15, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+0, o+12, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+1, o+9, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+2, o+7, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+3, o+6, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+4, o+5, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+5, o+4, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+6, o+3, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+7, o+2, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+9, o+1, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+12, o+0, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+1, o+12, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+2, o+9, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+3, o+7, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+4, o+6, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+5, o+5, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+6, o+4, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+7, o+3, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+9, o+2, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+12, o+1, Corner::RightBack);
        field.CopyDerivativeFreeWall(d+0, o2-1-15, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o2-1-14, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o2-1-13, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-1, o2-1-11, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-1, o2-1-10, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-2, o2-1-8, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-8, o2-1-2, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-10, o2-1-1, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-11, o2-1-1, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-13, o2-1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-14, o2-1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-15, o2-1+0, 1, Wall::Front);
        field.CopyDerivativeFreeConvexCorner(d+0, o2-1-12, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-1, o2-1-9, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-2, o2-1-7, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-3, o2-1-6, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-4, o2-1-5, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-5, o2-1-4, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-6, o2-1-3, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-7, o2-1-2, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-9, o2-1-1, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-12, o2-1+0, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(d-1, o2-1-12, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-2, o2-1-9, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-3, o2-1-7, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-4, o2-1-6, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-5, o2-1-5, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-6, o2-1-4, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-7, o2-1-3, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-9, o2-1-2, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-12, o2-1-1, Corner::LeftFront);
        field.CopyDerivativeFreeWall(LX-d-1+2, o2-1-8, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+1, o2-1-10, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+1, o2-1-11, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o2-1-13, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o2-1-14, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o2-1-15, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+15, o2-1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+14, o2-1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+13, o2-1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+11, o2-1-1, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+10, o2-1-1, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+8, o2-1-2, 1, Wall::Front);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+12, o2-1+0, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+9, o2-1-1, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+7, o2-1-2, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+6, o2-1-3, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+5, o2-1-4, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+4, o2-1-5, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+3, o2-1-6, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+2, o2-1-7, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+1, o2-1-9, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+0, o2-1-12, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+12, o2-1-1, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+9, o2-1-2, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+7, o2-1-3, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+6, o2-1-4, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+5, o2-1-5, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+4, o2-1-6, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+3, o2-1-7, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+2, o2-1-9, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+1, o2-1-12, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(0, 0, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(0, o2-1, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(0, o, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(0, LY-1, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-1, o, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-1, 0, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-1, o2-1, Corner::RightFront);
      };

      apply_u_bc(ux);
      apply_u_bc(uy);
      apply_uphi_bc(ux_phi);
      apply_uphi_bc(uy_phi);
      apply_bc(phi);
      apply_bc(MU);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }

    //other type of rounded corners with periodic reservoir
    case 17:
    case 18:
    {
     //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(d, o+16, l-32, Wall::Left);
        field.ApplyNeumannFreeWall(LX-d-1, o+16, l-32, Wall::Right);
        field.ApplyNeumannFreeWall(0, 0, LX, Wall::Front);
        field.ApplyNeumannFreeWall(0, o2-1, d-15, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d+15, o2-1, d-15, Wall::Front);
        field.ApplyNeumannFreeWall(0, LY-1, LX, Wall::Back);
        field.ApplyNeumannFreeWall(0, o, d-15, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d+15, o, d-15, Wall::Back);
        field.ApplyNeumannFreeWall(d-2, o+8, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o+10, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o+11, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+13, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+14, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+15, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-15, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-14, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-13, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-11, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-10, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-8, o+2, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(d-12, o+0, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-9, o+1, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-7, o+2, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-6, o+3, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-5, o+4, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-4, o+5, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-3, o+6, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-2, o+7, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-1, o+9, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d+0, o+12, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(d-12, o+1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-9, o+2, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-7, o+3, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-6, o+4, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-5, o+5, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-4, o+6, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-3, o+7, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-2, o+9, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-1, o+12, Corner::LeftBack);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+15, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+14, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+13, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o+11, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o+10, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+2, o+8, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+8, o+2, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+10, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+11, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+13, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+14, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+15, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+0, o+12, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+1, o+9, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+2, o+7, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+3, o+6, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+4, o+5, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+5, o+4, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+6, o+3, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+7, o+2, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+9, o+1, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+12, o+0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+1, o+12, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+2, o+9, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+3, o+7, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+4, o+6, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+5, o+5, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+6, o+4, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+7, o+3, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+9, o+2, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+12, o+1, Corner::RightBack);
        field.ApplyNeumannFreeWall(d+0, o2-1-15, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o2-1-14, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o2-1-13, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o2-1-11, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o2-1-10, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-2, o2-1-8, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-8, o2-1-2, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-10, o2-1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-11, o2-1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-13, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-14, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-15, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeConvexCorner(d+0, o2-1-12, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-1, o2-1-9, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-2, o2-1-7, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-3, o2-1-6, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-4, o2-1-5, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-5, o2-1-4, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-6, o2-1-3, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-7, o2-1-2, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-9, o2-1-1, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-12, o2-1+0, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(d-1, o2-1-12, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-2, o2-1-9, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-3, o2-1-7, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-4, o2-1-6, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-5, o2-1-5, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-6, o2-1-4, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-7, o2-1-3, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-9, o2-1-2, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-12, o2-1-1, Corner::LeftFront);
        field.ApplyNeumannFreeWall(LX-d-1+2, o2-1-8, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o2-1-10, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o2-1-11, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o2-1-13, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o2-1-14, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o2-1-15, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+15, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+14, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+13, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+11, o2-1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+10, o2-1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+8, o2-1-2, 1, Wall::Front);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+12, o2-1+0, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+9, o2-1-1, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+7, o2-1-2, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+6, o2-1-3, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+5, o2-1-4, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+4, o2-1-5, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+3, o2-1-6, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+2, o2-1-7, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+1, o2-1-9, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+0, o2-1-12, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+12, o2-1-1, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+9, o2-1-2, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+7, o2-1-3, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+6, o2-1-4, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+5, o2-1-5, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+4, o2-1-6, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+3, o2-1-7, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+2, o2-1-9, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+1, o2-1-12, Corner::RightFront);
      };

      auto apply_uphi_bc = [this](ScalarField& field) {

        field.ApplyDirichletFreeWall(d, o+16, l-32, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-d-1, o+16, l-32, Wall::Right, 0);
        field.ApplyDirichletFreeWall(0, 0, LX, Wall::Front, 0);
        field.ApplyDirichletFreeWall(0, o2-1, d-15, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d+15, o2-1, d-15, Wall::Front, 0);
        field.ApplyDirichletFreeWall(0, LY-1, LX, Wall::Back, 0);
        field.ApplyDirichletFreeWall(0, o, d-15, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d+15, o, d-15, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-2, o+8, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-1, o+10, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-1, o+11, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o+13, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o+14, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o+15, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-15, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-14, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-13, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-11, o+1, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-10, o+1, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-8, o+2, 1, Wall::Back, 0);
        field.ApplyDirichletFreeConvexCorner(d-12, o+0, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-9, o+1, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-7, o+2, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-6, o+3, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-5, o+4, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-4, o+5, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-3, o+6, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-2, o+7, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-1, o+9, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d+0, o+12, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-12, o+1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-9, o+2, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-7, o+3, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-6, o+4, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-5, o+5, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-4, o+6, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-3, o+7, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-2, o+9, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-1, o+12, Corner::LeftBack, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+15, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+14, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+13, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+1, o+11, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+1, o+10, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+2, o+8, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+8, o+2, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+10, o+1, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+11, o+1, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+13, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+14, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+15, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+0, o+12, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+1, o+9, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+2, o+7, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+3, o+6, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+4, o+5, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+5, o+4, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+6, o+3, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+7, o+2, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+9, o+1, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+12, o+0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+1, o+12, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+2, o+9, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+3, o+7, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+4, o+6, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+5, o+5, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+6, o+4, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+7, o+3, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+9, o+2, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+12, o+1, Corner::RightBack, 0);
        field.ApplyDirichletFreeWall(d+0, o2-1-15, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o2-1-14, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o2-1-13, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-1, o2-1-11, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-1, o2-1-10, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-2, o2-1-8, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-8, o2-1-2, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-10, o2-1-1, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-11, o2-1-1, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-13, o2-1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-14, o2-1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-15, o2-1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeConvexCorner(d+0, o2-1-12, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-1, o2-1-9, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-2, o2-1-7, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-3, o2-1-6, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-4, o2-1-5, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-5, o2-1-4, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-6, o2-1-3, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-7, o2-1-2, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-9, o2-1-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-12, o2-1+0, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-1, o2-1-12, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-2, o2-1-9, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-3, o2-1-7, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-4, o2-1-6, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-5, o2-1-5, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-6, o2-1-4, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-7, o2-1-3, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-9, o2-1-2, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-12, o2-1-1, Corner::LeftFront, 0);
        field.ApplyDirichletFreeWall(LX-d-1+2, o2-1-8, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+1, o2-1-10, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+1, o2-1-11, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o2-1-13, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o2-1-14, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o2-1-15, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+15, o2-1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+14, o2-1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+13, o2-1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+11, o2-1-1, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+10, o2-1-1, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+8, o2-1-2, 1, Wall::Front, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+12, o2-1+0, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+9, o2-1-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+7, o2-1-2, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+6, o2-1-3, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+5, o2-1-4, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+4, o2-1-5, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+3, o2-1-6, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+2, o2-1-7, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+1, o2-1-9, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+0, o2-1-12, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+12, o2-1-1, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+9, o2-1-2, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+7, o2-1-3, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+6, o2-1-4, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+5, o2-1-5, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+4, o2-1-6, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+3, o2-1-7, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+2, o2-1-9, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+1, o2-1-12, Corner::RightFront, 0);
      };

      auto apply_u_bc = [this](ScalarField& field) {

        field.CopyDerivativeFreeWall(d, o+16, l-32, Wall::Left);
        field.CopyDerivativeFreeWall(LX-d-1, o+16, l-32, Wall::Right);
        field.CopyDerivativeFreeWall(0, 0, LX, Wall::Front);
        field.CopyDerivativeFreeWall(0, o2-1, d-15, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d+15, o2-1, d-15, Wall::Front);
        field.CopyDerivativeFreeWall(0, LY-1, LX, Wall::Back);
        field.CopyDerivativeFreeWall(0, o, d-15, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d+15, o, d-15, Wall::Back);
        field.CopyDerivativeFreeWall(d-2, o+8, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-1, o+10, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-1, o+11, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o+13, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o+14, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o+15, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-15, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-14, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-13, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-11, o+1, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-10, o+1, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-8, o+2, 1, Wall::Back);
        field.CopyDerivativeFreeConvexCorner(d-12, o+0, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-9, o+1, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-7, o+2, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-6, o+3, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-5, o+4, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-4, o+5, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-3, o+6, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-2, o+7, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-1, o+9, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d+0, o+12, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(d-12, o+1, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-9, o+2, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-7, o+3, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-6, o+4, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-5, o+5, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-4, o+6, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-3, o+7, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-2, o+9, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-1, o+12, Corner::LeftBack);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+15, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+14, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+13, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+1, o+11, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+1, o+10, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+2, o+8, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+8, o+2, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+10, o+1, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+11, o+1, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+13, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+14, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+15, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+0, o+12, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+1, o+9, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+2, o+7, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+3, o+6, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+4, o+5, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+5, o+4, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+6, o+3, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+7, o+2, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+9, o+1, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+12, o+0, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+1, o+12, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+2, o+9, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+3, o+7, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+4, o+6, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+5, o+5, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+6, o+4, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+7, o+3, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+9, o+2, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+12, o+1, Corner::RightBack);
        field.CopyDerivativeFreeWall(d+0, o2-1-15, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o2-1-14, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o2-1-13, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-1, o2-1-11, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-1, o2-1-10, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-2, o2-1-8, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-8, o2-1-2, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-10, o2-1-1, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-11, o2-1-1, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-13, o2-1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-14, o2-1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-15, o2-1+0, 1, Wall::Front);
        field.CopyDerivativeFreeConvexCorner(d+0, o2-1-12, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-1, o2-1-9, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-2, o2-1-7, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-3, o2-1-6, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-4, o2-1-5, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-5, o2-1-4, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-6, o2-1-3, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-7, o2-1-2, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-9, o2-1-1, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-12, o2-1+0, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(d-1, o2-1-12, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-2, o2-1-9, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-3, o2-1-7, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-4, o2-1-6, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-5, o2-1-5, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-6, o2-1-4, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-7, o2-1-3, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-9, o2-1-2, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-12, o2-1-1, Corner::LeftFront);
        field.CopyDerivativeFreeWall(LX-d-1+2, o2-1-8, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+1, o2-1-10, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+1, o2-1-11, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o2-1-13, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o2-1-14, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o2-1-15, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+15, o2-1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+14, o2-1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+13, o2-1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+11, o2-1-1, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+10, o2-1-1, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+8, o2-1-2, 1, Wall::Front);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+12, o2-1+0, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+9, o2-1-1, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+7, o2-1-2, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+6, o2-1-3, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+5, o2-1-4, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+4, o2-1-5, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+3, o2-1-6, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+2, o2-1-7, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+1, o2-1-9, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+0, o2-1-12, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+12, o2-1-1, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+9, o2-1-2, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+7, o2-1-3, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+6, o2-1-4, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+5, o2-1-5, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+4, o2-1-6, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+3, o2-1-7, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+2, o2-1-9, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+1, o2-1-12, Corner::RightFront);
      };

      apply_u_bc(ux);
      apply_u_bc(uy);
      apply_uphi_bc(ux_phi);
      apply_uphi_bc(uy_phi);
      apply_bc(phi);
      apply_bc(MU);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }

    //other type of rounded corners with periodic reservoir, but different bc at top
    case 19:
    {
     //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(d, o+16, l-32, Wall::Left);
        field.ApplyNeumannFreeWall(LX-d-1, o+16, l-32, Wall::Right);
        field.ApplyNeumannFreeWall(0, 0, LX, Wall::Front);
        field.ApplyNeumannFreeWall(0, o2-1, d-15, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d+15, o2-1, d-15, Wall::Front);
        field.ApplyNeumannFreeWall(0, LY-1, LX, Wall::Back);
        field.ApplyNeumannFreeWall(0, o, d-15, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d+15, o, d-15, Wall::Back);
        field.ApplyNeumannFreeWall(d-2, o+8, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o+10, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o+11, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+13, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+14, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+15, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-15, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-14, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-13, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-11, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-10, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-8, o+2, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(d-12, o+0, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-9, o+1, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-7, o+2, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-6, o+3, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-5, o+4, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-4, o+5, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-3, o+6, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-2, o+7, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-1, o+9, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d+0, o+12, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(d-12, o+1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-9, o+2, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-7, o+3, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-6, o+4, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-5, o+5, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-4, o+6, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-3, o+7, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-2, o+9, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-1, o+12, Corner::LeftBack);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+15, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+14, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+13, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o+11, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o+10, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+2, o+8, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+8, o+2, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+10, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+11, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+13, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+14, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+15, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+0, o+12, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+1, o+9, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+2, o+7, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+3, o+6, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+4, o+5, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+5, o+4, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+6, o+3, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+7, o+2, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+9, o+1, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+12, o+0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+1, o+12, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+2, o+9, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+3, o+7, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+4, o+6, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+5, o+5, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+6, o+4, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+7, o+3, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+9, o+2, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+12, o+1, Corner::RightBack);
        field.ApplyNeumannFreeWall(d+0, o2-1-15, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o2-1-14, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o2-1-13, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o2-1-11, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o2-1-10, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-2, o2-1-8, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-8, o2-1-2, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-10, o2-1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-11, o2-1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-13, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-14, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-15, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeConvexCorner(d+0, o2-1-12, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-1, o2-1-9, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-2, o2-1-7, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-3, o2-1-6, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-4, o2-1-5, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-5, o2-1-4, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-6, o2-1-3, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-7, o2-1-2, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-9, o2-1-1, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-12, o2-1+0, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(d-1, o2-1-12, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-2, o2-1-9, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-3, o2-1-7, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-4, o2-1-6, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-5, o2-1-5, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-6, o2-1-4, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-7, o2-1-3, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-9, o2-1-2, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-12, o2-1-1, Corner::LeftFront);
        field.ApplyNeumannFreeWall(LX-d-1+2, o2-1-8, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o2-1-10, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o2-1-11, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o2-1-13, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o2-1-14, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o2-1-15, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+15, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+14, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+13, o2-1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+11, o2-1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+10, o2-1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+8, o2-1-2, 1, Wall::Front);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+12, o2-1+0, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+9, o2-1-1, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+7, o2-1-2, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+6, o2-1-3, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+5, o2-1-4, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+4, o2-1-5, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+3, o2-1-6, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+2, o2-1-7, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+1, o2-1-9, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+0, o2-1-12, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+12, o2-1-1, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+9, o2-1-2, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+7, o2-1-3, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+6, o2-1-4, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+5, o2-1-5, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+4, o2-1-6, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+3, o2-1-7, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+2, o2-1-9, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+1, o2-1-12, Corner::RightFront);
      };

      auto apply_uphi_bc = [this](ScalarField& field) {

        field.ApplyDirichletFreeWall(d, o+16, l-32, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-d-1, o+16, l-32, Wall::Right, 0);
        field.ApplyDirichletFreeWall(0, 0, LX, Wall::Front, 0);
        field.ApplyDirichletFreeWall(0, o2-1, d-15, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d+15, o2-1, d-15, Wall::Front, 0);
        field.ApplyNeumannFreeWall(0, LY-1, LX, Wall::Back);
        field.ApplyDirichletFreeWall(0, o, d-15, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d+15, o, d-15, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-2, o+8, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-1, o+10, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-1, o+11, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o+13, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o+14, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o+15, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-15, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-14, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-13, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-11, o+1, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-10, o+1, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-8, o+2, 1, Wall::Back, 0);
        field.ApplyDirichletFreeConvexCorner(d-12, o+0, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-9, o+1, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-7, o+2, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-6, o+3, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-5, o+4, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-4, o+5, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-3, o+6, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-2, o+7, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-1, o+9, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d+0, o+12, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-12, o+1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-9, o+2, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-7, o+3, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-6, o+4, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-5, o+5, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-4, o+6, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-3, o+7, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-2, o+9, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-1, o+12, Corner::LeftBack, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+15, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+14, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+13, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+1, o+11, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+1, o+10, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+2, o+8, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+8, o+2, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+10, o+1, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+11, o+1, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+13, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+14, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+15, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+0, o+12, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+1, o+9, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+2, o+7, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+3, o+6, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+4, o+5, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+5, o+4, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+6, o+3, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+7, o+2, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+9, o+1, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+12, o+0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+1, o+12, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+2, o+9, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+3, o+7, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+4, o+6, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+5, o+5, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+6, o+4, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+7, o+3, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+9, o+2, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+12, o+1, Corner::RightBack, 0);
        field.ApplyDirichletFreeWall(d+0, o2-1-15, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o2-1-14, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o2-1-13, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-1, o2-1-11, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-1, o2-1-10, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-2, o2-1-8, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-8, o2-1-2, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-10, o2-1-1, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-11, o2-1-1, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-13, o2-1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-14, o2-1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-15, o2-1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeConvexCorner(d+0, o2-1-12, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-1, o2-1-9, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-2, o2-1-7, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-3, o2-1-6, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-4, o2-1-5, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-5, o2-1-4, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-6, o2-1-3, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-7, o2-1-2, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-9, o2-1-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-12, o2-1+0, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-1, o2-1-12, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-2, o2-1-9, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-3, o2-1-7, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-4, o2-1-6, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-5, o2-1-5, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-6, o2-1-4, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-7, o2-1-3, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-9, o2-1-2, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-12, o2-1-1, Corner::LeftFront, 0);
        field.ApplyDirichletFreeWall(LX-d-1+2, o2-1-8, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+1, o2-1-10, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+1, o2-1-11, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o2-1-13, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o2-1-14, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o2-1-15, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+15, o2-1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+14, o2-1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+13, o2-1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+11, o2-1-1, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+10, o2-1-1, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+8, o2-1-2, 1, Wall::Front, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+12, o2-1+0, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+9, o2-1-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+7, o2-1-2, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+6, o2-1-3, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+5, o2-1-4, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+4, o2-1-5, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+3, o2-1-6, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+2, o2-1-7, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+1, o2-1-9, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+0, o2-1-12, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+12, o2-1-1, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+9, o2-1-2, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+7, o2-1-3, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+6, o2-1-4, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+5, o2-1-5, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+4, o2-1-6, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+3, o2-1-7, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+2, o2-1-9, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+1, o2-1-12, Corner::RightFront, 0);
      };

      auto apply_u_bc = [this](ScalarField& field) {

        field.CopyDerivativeFreeWall(d, o+16, l-32, Wall::Left);
        field.CopyDerivativeFreeWall(LX-d-1, o+16, l-32, Wall::Right);
        field.CopyDerivativeFreeWall(0, 0, LX, Wall::Front);
        field.CopyDerivativeFreeWall(0, o2-1, d-15, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d+15, o2-1, d-15, Wall::Front);
        field.ApplyNeumannFreeWall(0, LY-1, LX, Wall::Back);
        field.CopyDerivativeFreeWall(0, o, d-15, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d+15, o, d-15, Wall::Back);
        field.CopyDerivativeFreeWall(d-2, o+8, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-1, o+10, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-1, o+11, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o+13, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o+14, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o+15, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-15, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-14, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-13, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-11, o+1, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-10, o+1, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-8, o+2, 1, Wall::Back);
        field.CopyDerivativeFreeConvexCorner(d-12, o+0, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-9, o+1, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-7, o+2, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-6, o+3, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-5, o+4, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-4, o+5, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-3, o+6, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-2, o+7, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-1, o+9, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d+0, o+12, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(d-12, o+1, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-9, o+2, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-7, o+3, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-6, o+4, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-5, o+5, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-4, o+6, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-3, o+7, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-2, o+9, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-1, o+12, Corner::LeftBack);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+15, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+14, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+13, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+1, o+11, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+1, o+10, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+2, o+8, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+8, o+2, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+10, o+1, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+11, o+1, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+13, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+14, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+15, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+0, o+12, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+1, o+9, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+2, o+7, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+3, o+6, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+4, o+5, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+5, o+4, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+6, o+3, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+7, o+2, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+9, o+1, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+12, o+0, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+1, o+12, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+2, o+9, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+3, o+7, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+4, o+6, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+5, o+5, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+6, o+4, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+7, o+3, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+9, o+2, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+12, o+1, Corner::RightBack);
        field.CopyDerivativeFreeWall(d+0, o2-1-15, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o2-1-14, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o2-1-13, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-1, o2-1-11, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-1, o2-1-10, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-2, o2-1-8, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-8, o2-1-2, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-10, o2-1-1, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-11, o2-1-1, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-13, o2-1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-14, o2-1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-15, o2-1+0, 1, Wall::Front);
        field.CopyDerivativeFreeConvexCorner(d+0, o2-1-12, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-1, o2-1-9, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-2, o2-1-7, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-3, o2-1-6, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-4, o2-1-5, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-5, o2-1-4, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-6, o2-1-3, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-7, o2-1-2, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-9, o2-1-1, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-12, o2-1+0, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(d-1, o2-1-12, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-2, o2-1-9, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-3, o2-1-7, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-4, o2-1-6, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-5, o2-1-5, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-6, o2-1-4, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-7, o2-1-3, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-9, o2-1-2, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-12, o2-1-1, Corner::LeftFront);
        field.CopyDerivativeFreeWall(LX-d-1+2, o2-1-8, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+1, o2-1-10, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+1, o2-1-11, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o2-1-13, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o2-1-14, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o2-1-15, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+15, o2-1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+14, o2-1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+13, o2-1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+11, o2-1-1, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+10, o2-1-1, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+8, o2-1-2, 1, Wall::Front);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+12, o2-1+0, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+9, o2-1-1, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+7, o2-1-2, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+6, o2-1-3, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+5, o2-1-4, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+4, o2-1-5, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+3, o2-1-6, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+2, o2-1-7, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+1, o2-1-9, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+0, o2-1-12, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+12, o2-1-1, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+9, o2-1-2, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+7, o2-1-3, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+6, o2-1-4, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+5, o2-1-5, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+4, o2-1-6, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+3, o2-1-7, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+2, o2-1-9, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+1, o2-1-12, Corner::RightFront);
      };

      apply_u_bc(ux);
      apply_u_bc(uy);
      apply_uphi_bc(ux_phi);
      apply_uphi_bc(uy_phi);
      apply_bc(phi);
      apply_bc(MU);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }
    //reservoir with chimney
    case 20:
    {
      //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(d, o+1, LY-o-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-d-1, o+1, LY-o-2, Wall::Right);
        field.ApplyNeumannFreeWall(0, 0, LX, Wall::Front);
        field.ApplyNeumannFreeWall(d+1, LY-1, LX-2*d-2, Wall::Back);
        //field.CopyDerivativeFreeWall(d+1, LY-1, LX-2*d-2, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d, o, d, Wall::Back);
        field.ApplyNeumannFreeWall(0, o, d, Wall::Back);
        field.ApplyNeumannFreeConcaveCorner(d, LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1, o, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(d, o, Corner::RightFront);
      };

      auto apply_uphi_bc = [this](ScalarField& field) {

        field.ApplyDirichletFreeWall(d, o+1, LY-o-2, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-d-1, o+1, LY-o-2, Wall::Right, 0);
        field.ApplyDirichletFreeWall(0, 0, LX, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d+1, LY-1, LX-2*d-2, Wall::Back, 0);
        //field.CopyDerivativeFreeWall(d+1, LY-1, LX-2*d-2, Wall::Back);
        field.ApplyDirichletFreeWall(LX-d, o, d, Wall::Back, 0);
        field.ApplyDirichletFreeWall(0, o, d, Wall::Back, 0);
        field.ApplyDirichletFreeConcaveCorner(d, LY-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1, o, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(d, o, Corner::RightFront, 0);
      };

      auto apply_u_bc = [this](ScalarField& field) {

        field.CopyDerivativeFreeWall(d, o+1, LY-o-2, Wall::Left);
        field.CopyDerivativeFreeWall(LX-d-1, o+1, LY-o-2, Wall::Right);
        field.CopyDerivativeFreeWall(0, 0, LX, Wall::Front);
        field.CopyDerivativeFreeWall(d+1, LY-1, LX-2*d-2, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d, o, d, Wall::Back);
        field.CopyDerivativeFreeWall(0, o, d, Wall::Back);
        field.CopyDerivativeFreeConcaveCorner(d, LY-1, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1, o, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(d, o, Corner::RightFront);
      };

      apply_u_bc(ux);
      apply_u_bc(uy);
      apply_uphi_bc(ux_phi);
      apply_uphi_bc(uy_phi);
      apply_bc(phi);
      apply_bc(MU);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }
    //box with outlet at top
    case 21:
    {
      //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyNeumannFreeWall(1, 0, LX-2, Wall::Front);
        field.ApplyNeumannFreeWall(1, LY-1, LX-2, Wall::Back);
        field.ApplyNeumannFreeConcaveCorner(0, 0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, 0, Corner::RightFront);
      };

      auto apply_uphi_bc = [this](ScalarField& field) {

        field.ApplyDirichletFreeWall(0, 1, LY-2, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-1, 1, LY-2, Wall::Right, 0);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front, 0);
        field.ApplyNeumannFreeWall(1, LY-1, LX-2, Wall::Back);
        field.ApplyDirichletFreeConcaveCorner(0, 0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(0, LY-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, LY-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, 0, Corner::RightFront, 0);
      };

      auto apply_u_bc = [this](ScalarField& field) {

        field.CopyDerivativeFreeWall(0, 1, LY-2, Wall::Left);
        field.CopyDerivativeFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.CopyDerivativeFreeWall(1, 0, LX-2, Wall::Front);
        field.ApplyNeumannFreeWall(1, LY-1, LX-2, Wall::Back);
        field.CopyDerivativeFreeConcaveCorner(0, 0, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(0, LY-1, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-1, 0, Corner::RightFront);
      };
      apply_u_bc(ux);
      apply_u_bc(uy);
      apply_uphi_bc(ux_phi);
      apply_uphi_bc(uy_phi);
      apply_bc(phi);
      apply_bc(MU);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }
    //box with inlet at bottom and outlet at top
    case 22:
    {
      //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyNeumannFreeWall(1, 0, LX-2, Wall::Front);
        field.ApplyNeumannFreeWall(1, LY-1, LX-2, Wall::Back);
        field.ApplyNeumannFreeConcaveCorner(0, 0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, 0, Corner::RightFront);
      };

      auto apply_phi_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front,(1+noise*random_real())*conc);
        field.ApplyNeumannFreeWall(1, LY-1, LX-2, Wall::Back);
        field.ApplyNeumannFreeConcaveCorner(0, 0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, 0, Corner::RightFront);
      };

      auto apply_uxphi_bc = [this](ScalarField& field) {

        field.ApplyDirichletFreeWall(0, 1, LY-2, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-1, 1, LY-2, Wall::Right, 0);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front, 0);
        field.ApplyVxPhiOutletFreeWall(1, LY-1, LX-2, Wall::Back, uin, ff, phi);
        field.ApplyDirichletFreeConcaveCorner(0, 0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(0, LY-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, LY-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, 0, Corner::RightFront, 0);
      };

      auto apply_uyphi_bc = [this](ScalarField& field) {

        field.ApplyDirichletFreeWall(0, 1, LY-2, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-1, 1, LY-2, Wall::Right, 0);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front, (1+noise*random_real())*conc*uin);
        field.ApplyVyPhiOutletFreeWall(1, LY-1, LX-2, Wall::Back, uin, ff, phi);
        field.ApplyDirichletFreeConcaveCorner(0, 0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(0, LY-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, LY-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, 0, Corner::RightFront, 0);
      };

      auto apply_ux_bc = [this](ScalarField& field) {

        field.CopyDerivativeFreeWall(0, 1, LY-2, Wall::Left);
        field.CopyDerivativeFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front, 0);
        field.ApplyVxOutletFreeWall(1, LY-1, LX-2, Wall::Back, uin, ff);
        field.CopyDerivativeFreeConcaveCorner(0, 0, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(0, LY-1, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-1, 0, Corner::RightFront);
      };

      auto apply_uy_bc = [this](ScalarField& field) {

        field.CopyDerivativeFreeWall(0, 1, LY-2, Wall::Left);
        field.CopyDerivativeFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front, uin);
        field.ApplyVyOutletFreeWall(1, LY-1, LX-2, Wall::Back, uin, ff);
        field.CopyDerivativeFreeConcaveCorner(0, 0, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(0, LY-1, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-1, 0, Corner::RightFront);
      };

      apply_ux_bc(ux);
      apply_uy_bc(uy);
      apply_uxphi_bc(ux_phi);
      apply_uyphi_bc(uy_phi);
      apply_phi_bc(phi);
      apply_bc(MU);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }
    // box with pressure boundary on top
    case 23:
    {
      //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyNeumannFreeWall(1, 0, LX-2, Wall::Front);
        field.CopyDerivativeFreeWall(1, LY-1, LX-2, Wall::Back);
        field.ApplyNeumannFreeConcaveCorner(0, 0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, 0, Corner::RightFront);
      };

      auto apply_uphi_bc = [this](ScalarField& field) {

        field.ApplyDirichletFreeWall(0, 1, LY-2, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-1, 1, LY-2, Wall::Right, 0);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front, 0);
        field.CopyDerivativeFreeWall(1, LY-1, LX-2, Wall::Back);
        field.ApplyDirichletFreeConcaveCorner(0, 0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(0, LY-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, LY-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, 0, Corner::RightFront, 0);
      };

      auto apply_u_bc = [this](ScalarField& field) {

        field.CopyDerivativeFreeWall(0, 1, LY-2, Wall::Left);
        field.CopyDerivativeFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.CopyDerivativeFreeWall(1, 0, LX-2, Wall::Front);
        field.CopyDerivativeFreeWall(1, LY-1, LX-2, Wall::Back);
        field.CopyDerivativeFreeConcaveCorner(0, 0, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(0, LY-1, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-1, 0, Corner::RightFront);
      };

      apply_u_bc(ux);
      apply_u_bc(uy);
      apply_uphi_bc(ux_phi);
      apply_uphi_bc(uy_phi);
      apply_bc(phi);
      apply_bc(MU);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }
    // box with pressure boundary on top and diffuvise phi inlet at bottom
    case 24:
    {
      //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyNeumannFreeWall(1, 0, LX-2, Wall::Front);
        field.CopyDerivativeFreeWall(1, LY-1, LX-2, Wall::Back);
        field.ApplyNeumannFreeConcaveCorner(0, 0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, 0, Corner::RightFront);
      };

      auto apply_phi_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front, conc);
        field.CopyDerivativeFreeWall(1, LY-1, LX-2, Wall::Back);
        field.ApplyNeumannFreeConcaveCorner(0, 0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, 0, Corner::RightFront);
      };

      auto apply_uphi_bc = [this](ScalarField& field) {

        field.ApplyDirichletFreeWall(0, 1, LY-2, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-1, 1, LY-2, Wall::Right, 0);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front, 0);
        field.CopyDerivativeFreeWall(1, LY-1, LX-2, Wall::Back);
        field.ApplyDirichletFreeConcaveCorner(0, 0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(0, LY-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, LY-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, 0, Corner::RightFront, 0);
      };

      auto apply_u_bc = [this](ScalarField& field) {

        field.CopyDerivativeFreeWall(0, 1, LY-2, Wall::Left);
        field.CopyDerivativeFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.CopyDerivativeFreeWall(1, 0, LX-2, Wall::Front);
        field.CopyDerivativeFreeWall(1, LY-1, LX-2, Wall::Back);
        field.CopyDerivativeFreeConcaveCorner(0, 0, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(0, LY-1, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-1, 0, Corner::RightFront);
      };

      apply_u_bc(ux);
      apply_u_bc(uy);
      apply_uphi_bc(ux_phi);
      apply_uphi_bc(uy_phi);
      apply_phi_bc(phi);
      apply_bc(MU);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }
    //box with outlet at top and influx at bottom
    case 25:
    {
      //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyNeumannFreeWall(1, 0, LX-2, Wall::Front);
        field.CopyDerivativeFreeWall(1, LY-1, LX-2, Wall::Back);
        field.ApplyNeumannFreeConcaveCorner(0, 0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(0, LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1, 0, Corner::RightFront);
      };

      auto apply_phi_bc = [this](ScalarField& field) {
        // walls
        field.ApplyNeumannFreeWall(0, 1, LY-2, Wall::Left);
        field.ApplyNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front, conc);
        field.CopyDerivativeFreeWall(1, LY-1, LX-2, Wall::Back);
        // corners
        field.ApplyNeumannFreeConcaveCorner(LX-1,LY-1, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-1,0, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(0,LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(0,0, Corner::LeftFront);
      };

      auto apply_uxphi_bc = [this](ScalarField& field) {

        field.ApplyDirichletFreeWall(0, 1, LY-2, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-1, 1, LY-2, Wall::Right, 0);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front, 0);
        field.ApplyDirichletFreeWall(1, LY-1, LX-2, Wall::Back, 0);
        field.ApplyDirichletFreeConcaveCorner(0, 0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(0, LY-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, LY-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, 0, Corner::RightFront, 0);
      };

      auto apply_uyphi_bc = [this](ScalarField& field) {

        field.ApplyDirichletFreeWall(0, 1, LY-2, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-1, 1, LY-2, Wall::Right, 0);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front, conc*uin);
        field.ApplyDirichletFreeWall(1, LY-1, LX-2, Wall::Back, 0);
        field.ApplyDirichletFreeConcaveCorner(0, 0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(0, LY-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, LY-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, 0, Corner::RightFront, 0);
      };

      auto apply_sigmaXX_bc = [this](ScalarField& field) {

        field.ApplyDirichletFreeWall(0, 1, LY-2, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-1, 1, LY-2, Wall::Right, 0);
        field.ApplyStressNeumannInletFreeWall(1, 0, LX-2, Wall::Front, QQxx, QQyx, phi, AA, CC, KK, LL, zeta, xi, TensorComponent::XX);
        field.ApplyDirichletFreeWall(1, LY-1, LX-2, Wall::Back, 0);
        field.ApplyDirichletFreeConcaveCorner(0, 0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(0, LY-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, LY-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, 0, Corner::RightFront, 0);
      };
      auto apply_sigmaYY_bc = [this](ScalarField& field) {

        field.ApplyDirichletFreeWall(0, 1, LY-2, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-1, 1, LY-2, Wall::Right, 0);
        field.ApplyStressNeumannInletFreeWall(1, 0, LX-2, Wall::Front, QQxx, QQyx, phi, AA, CC, KK, LL, zeta, xi, TensorComponent::YY);
        field.ApplyDirichletFreeWall(1, LY-1, LX-2, Wall::Back, 0);
        field.ApplyDirichletFreeConcaveCorner(0, 0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(0, LY-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, LY-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, 0, Corner::RightFront, 0);
      };
      auto apply_sigmaYX_bc = [this](ScalarField& field) {

        field.ApplyDirichletFreeWall(0, 1, LY-2, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-1, 1, LY-2, Wall::Right, 0);
        field.ApplyStressNeumannInletFreeWall(1, 0, LX-2, Wall::Front, QQxx, QQyx, phi, AA, CC, KK, LL, zeta, xi, TensorComponent::YX);
        field.ApplyDirichletFreeWall(1, LY-1, LX-2, Wall::Back, 0);
        field.ApplyDirichletFreeConcaveCorner(0, 0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(0, LY-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, LY-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, 0, Corner::RightFront, 0);
      };
      auto apply_sigmaXY_bc = [this](ScalarField& field) {

        field.ApplyDirichletFreeWall(0, 1, LY-2, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-1, 1, LY-2, Wall::Right, 0);
        field.ApplyStressNeumannInletFreeWall(1, 0, LX-2, Wall::Front, QQxx, QQyx, phi, AA, CC, KK, LL, zeta, xi, TensorComponent::XY);
        field.ApplyDirichletFreeWall(1, LY-1, LX-2, Wall::Back, 0);
        field.ApplyDirichletFreeConcaveCorner(0, 0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(0, LY-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, LY-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, 0, Corner::RightFront, 0);
      };
      auto apply_ux_bc = [this](ScalarField& field) {

        field.CopyDerivativeFreeWall(0, 1, LY-2, Wall::Left);
        field.CopyDerivativeFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front,0);
        field.CopyDerivativeFreeWall(1, LY-1, LX-2, Wall::Back);
        field.CopyDerivativeFreeConcaveCorner(0, 0, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(0, LY-1, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-1, 0, Corner::RightFront);

        /*
        field.ApplyDirichletFreeWall(0, 1, LY-2, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-1, 1, LY-2, Wall::Right, 0);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front, 0);
        field.ApplyDirichletFreeWall(1, LY-1, LX-2, Wall::Back, 0);
        field.ApplyDirichletFreeConcaveCorner(0, 0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(0, LY-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, LY-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, 0, Corner::RightFront, 0);
        */
      };
      auto apply_uy_bc = [this](ScalarField& field) {

        field.CopyDerivativeFreeWall(0, 1, LY-2, Wall::Left);
        field.CopyDerivativeFreeWall(LX-1, 1, LY-2, Wall::Right);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front, uin);
        field.CopyDerivativeFreeWall(1, LY-1, LX-2, Wall::Back);
        field.CopyDerivativeFreeConcaveCorner(0, 0, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(0, LY-1, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(LX-1, LY-1, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-1, 0, Corner::RightFront);

        /*
        field.ApplyDirichletFreeWall(0, 1, LY-2, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-1, 1, LY-2, Wall::Right, 0);
        field.ApplyDirichletFreeWall(1, 0, LX-2, Wall::Front, uin);
        field.ApplyDirichletFreeWall(1, LY-1, LX-2, Wall::Back, 0);
        field.ApplyDirichletFreeConcaveCorner(0, 0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(0, LY-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, LY-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-1, 0, Corner::RightFront, 0);
        */
      };
      apply_ux_bc(ux);
      apply_uy_bc(uy);
      apply_uxphi_bc(ux_phi);
      apply_uyphi_bc(uy_phi);
      apply_phi_bc(phi);
      apply_bc(MU);
      apply_sigmaXX_bc(sigmaXX);
      apply_sigmaYY_bc(sigmaYY);
      apply_sigmaYX_bc(sigmaYX);
      apply_sigmaXY_bc(sigmaXY);
      break;
    }
    //reservoir with chimney, radius corners=4
    case 26:
    {
      //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(d, o+5, LY-o-6, Wall::Left);
        field.ApplyNeumannFreeWall(LX-d-1, o+5, LY-o-6, Wall::Right);
        field.ApplyNeumannFreeWall(0, 0, LX, Wall::Front);
        field.CopyDerivativeFreeWall(d+1, LY-1, LX-2*d-2, Wall::Back);
        field.ApplyNeumannFreeConcaveCorner(d, LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeWall(LX-d+4, o, d-4, Wall::Back);
        field.ApplyNeumannFreeWall(0, o, d-4, Wall::Back);
        field.ApplyNeumannFreeWall(d+0, o+4, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-4, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(d-3, o+0, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-2, o+1, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-1, o+2, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d+0, o+3, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(d-3, o+1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-2, o+2, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-1, o+3, Corner::LeftBack);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+4, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+4, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+0, o+3, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+1, o+2, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+2, o+1, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+3, o+0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+1, o+3, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+2, o+2, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+3, o+1, Corner::RightBack);
      };

      auto apply_uphi_bc = [this](ScalarField& field) {

        field.ApplyDirichletFreeWall(d, o+5, LY-o-6, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-d-1, o+5, LY-o-6, Wall::Right, 0);
        field.ApplyDirichletFreeWall(0, 0, LX, Wall::Front, 0);
        field.CopyDerivativeFreeWall(d+1, LY-1, LX-2*d-2, Wall::Back);
        field.ApplyDirichletFreeConcaveCorner(d, LY-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeWall(LX-d+4, o, d-4, Wall::Back, 0);
        field.ApplyDirichletFreeWall(0, o, d-4, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d+0, o+4, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-4, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeConvexCorner(d-3, o+0, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-2, o+1, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-1, o+2, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d+0, o+3, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-3, o+1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-2, o+2, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-1, o+3, Corner::LeftBack, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+4, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+4, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+0, o+3, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+1, o+2, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+2, o+1, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+3, o+0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+1, o+3, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+2, o+2, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+3, o+1, Corner::RightBack, 0);
      };

      auto apply_u_bc = [this](ScalarField& field) {

        field.CopyDerivativeFreeWall(d, o+5, LY-o-6, Wall::Left);
        field.CopyDerivativeFreeWall(LX-d-1, o+5, LY-o-6, Wall::Right);
        field.CopyDerivativeFreeWall(0, 0, LX, Wall::Front);
        field.CopyDerivativeFreeWall(d+1, LY-1, LX-2*d-2, Wall::Back);
        field.CopyDerivativeFreeConcaveCorner(d, LY-1, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack);
        field.CopyDerivativeFreeWall(LX-d+4, o, d-4, Wall::Back);
        field.CopyDerivativeFreeWall(0, o, d-4, Wall::Back);
        field.CopyDerivativeFreeWall(d+0, o+4, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-4, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeConvexCorner(d-3, o+0, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-2, o+1, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-1, o+2, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d+0, o+3, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(d-3, o+1, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-2, o+2, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-1, o+3, Corner::LeftBack);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+4, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+4, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+0, o+3, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+1, o+2, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+2, o+1, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+3, o+0, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+1, o+3, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+2, o+2, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+3, o+1, Corner::RightBack);
      };

      apply_u_bc(ux);
      apply_u_bc(uy);
      apply_uphi_bc(ux_phi);
      apply_uphi_bc(uy_phi);
      apply_bc(phi);
      apply_bc(MU);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }
    //reservoir with chimney, radius corners=7
    case 27:
    {
      //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(d, o+8, LY-o-9, Wall::Left);
        field.ApplyNeumannFreeWall(LX-d-1, o+8, LY-o-9, Wall::Right);
        field.ApplyNeumannFreeWall(0, 0, LX, Wall::Front);
        field.CopyDerivativeFreeWall(d+1, LY-1, LX-2*d-2, Wall::Back);
        field.ApplyNeumannFreeConcaveCorner(d, LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeWall(LX-d+7, o, d-7, Wall::Back);
        field.ApplyNeumannFreeWall(0, o, d-7, Wall::Back);
        field.ApplyNeumannFreeWall(d-1, o+4, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+6, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+7, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-7, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-6, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-4, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(d-5, o+0, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-3, o+1, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-2, o+2, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-1, o+3, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d+0, o+5, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(d-5, o+1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-3, o+2, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-2, o+3, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-1, o+5, Corner::LeftBack);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+7, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+6, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o+4, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+4, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+6, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+7, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+0, o+5, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+1, o+3, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+2, o+2, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+3, o+1, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+5, o+0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+1, o+5, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+2, o+3, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+3, o+2, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+5, o+1, Corner::RightBack);
      };

      auto apply_uphi_bc = [this](ScalarField& field) {

        field.ApplyDirichletFreeWall(d, o+8, LY-o-9, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-d-1, o+8, LY-o-9, Wall::Right, 0);
        field.ApplyDirichletFreeWall(0, 0, LX, Wall::Front, 0);
        field.CopyDerivativeFreeWall(d+1, LY-1, LX-2*d-2, Wall::Back);
        field.ApplyDirichletFreeConcaveCorner(d, LY-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeWall(LX-d+7, o, d-7, Wall::Back, 0);
        field.ApplyDirichletFreeWall(0, o, d-7, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-1, o+4, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o+6, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o+7, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-7, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-6, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-4, o+1, 1, Wall::Back, 0);
        field.ApplyDirichletFreeConvexCorner(d-5, o+0, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-3, o+1, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-2, o+2, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-1, o+3, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d+0, o+5, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-5, o+1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-3, o+2, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-2, o+3, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-1, o+5, Corner::LeftBack, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+7, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+6, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+1, o+4, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+4, o+1, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+6, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+7, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+0, o+5, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+1, o+3, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+2, o+2, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+3, o+1, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+5, o+0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+1, o+5, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+2, o+3, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+3, o+2, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+5, o+1, Corner::RightBack, 0);
      };

      auto apply_u_bc = [this](ScalarField& field) {

        field.CopyDerivativeFreeWall(d, o+8, LY-o-9, Wall::Left);
        field.CopyDerivativeFreeWall(LX-d-1, o+8, LY-o-9, Wall::Right);
        field.CopyDerivativeFreeWall(0, 0, LX, Wall::Front);
        field.CopyDerivativeFreeWall(d+1, LY-1, LX-2*d-2, Wall::Back);
        field.CopyDerivativeFreeConcaveCorner(d, LY-1, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack);
        field.CopyDerivativeFreeWall(LX-d+7, o, d-7, Wall::Back);
        field.CopyDerivativeFreeWall(0, o, d-7, Wall::Back);
        field.CopyDerivativeFreeWall(d-1, o+4, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o+6, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o+7, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-7, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-6, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-4, o+1, 1, Wall::Back);
        field.CopyDerivativeFreeConvexCorner(d-5, o+0, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-3, o+1, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-2, o+2, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-1, o+3, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d+0, o+5, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(d-5, o+1, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-3, o+2, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-2, o+3, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-1, o+5, Corner::LeftBack);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+7, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+6, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+1, o+4, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+4, o+1, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+6, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+7, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+0, o+5, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+1, o+3, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+2, o+2, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+3, o+1, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+5, o+0, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+1, o+5, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+2, o+3, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+3, o+2, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+5, o+1, Corner::RightBack);
      };

      apply_u_bc(ux);
      apply_u_bc(uy);
      apply_uphi_bc(ux_phi);
      apply_uphi_bc(uy_phi);
      apply_bc(phi);
      apply_bc(MU);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }
    //reservoir with chimney, radius corners=10
    case 28:
    {
      //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(d, o+11, LY-o-12, Wall::Left);
        field.ApplyNeumannFreeWall(LX-d-1, o+11, LY-o-12, Wall::Right);
        field.ApplyNeumannFreeWall(0, 0, LX, Wall::Front);
        field.CopyDerivativeFreeWall(d+1, LY-1, LX-2*d-2, Wall::Back);
        field.ApplyNeumannFreeConcaveCorner(d, LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeWall(LX-d+10, o, d-10, Wall::Back);
        field.ApplyNeumannFreeWall(0, o, d-10, Wall::Back);
        field.ApplyNeumannFreeWall(d-1, o+6, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+8, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+9, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+10, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-10, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-9, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-8, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-6, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(d-7, o+0, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-5, o+1, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-4, o+2, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-3, o+3, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-2, o+4, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-1, o+5, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d+0, o+7, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(d-7, o+1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-5, o+2, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-4, o+3, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-3, o+4, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-2, o+5, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-1, o+7, Corner::LeftBack);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+10, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+9, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+8, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o+6, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+6, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+8, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+9, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+10, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+0, o+7, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+1, o+5, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+2, o+4, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+3, o+3, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+4, o+2, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+5, o+1, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+7, o+0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+1, o+7, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+2, o+5, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+3, o+4, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+4, o+3, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+5, o+2, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+7, o+1, Corner::RightBack);
      };

      auto apply_uphi_bc = [this](ScalarField& field) {

        field.ApplyDirichletFreeWall(d, o+11, LY-o-12, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-d-1, o+11, LY-o-12, Wall::Right, 0);
        field.ApplyDirichletFreeWall(0, 0, LX, Wall::Front, 0);
        field.CopyDerivativeFreeWall(d+1, LY-1, LX-2*d-2, Wall::Back);
        field.ApplyDirichletFreeConcaveCorner(d, LY-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeWall(LX-d+10, o, d-10, Wall::Back, 0);
        field.ApplyDirichletFreeWall(0, o, d-10, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-1, o+6, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o+8, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o+9, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o+10, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-10, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-9, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-8, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-6, o+1, 1, Wall::Back, 0);
        field.ApplyDirichletFreeConvexCorner(d-7, o+0, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-5, o+1, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-4, o+2, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-3, o+3, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-2, o+4, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-1, o+5, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d+0, o+7, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-7, o+1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-5, o+2, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-4, o+3, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-3, o+4, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-2, o+5, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-1, o+7, Corner::LeftBack, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+10, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+9, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+8, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+1, o+6, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+6, o+1, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+8, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+9, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+10, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+0, o+7, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+1, o+5, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+2, o+4, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+3, o+3, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+4, o+2, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+5, o+1, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+7, o+0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+1, o+7, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+2, o+5, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+3, o+4, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+4, o+3, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+5, o+2, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+7, o+1, Corner::RightBack, 0);
      };

      auto apply_u_bc = [this](ScalarField& field) {

        field.CopyDerivativeFreeWall(d, o+11, LY-o-12, Wall::Left);
        field.CopyDerivativeFreeWall(LX-d-1, o+11, LY-o-12, Wall::Right);
        field.CopyDerivativeFreeWall(0, 0, LX, Wall::Front);
        field.CopyDerivativeFreeWall(d+1, LY-1, LX-2*d-2, Wall::Back);
        field.CopyDerivativeFreeConcaveCorner(d, LY-1, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack);
        field.CopyDerivativeFreeWall(LX-d+10, o, d-10, Wall::Back);
        field.CopyDerivativeFreeWall(0, o, d-10, Wall::Back);
        field.CopyDerivativeFreeWall(d-1, o+6, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o+8, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o+9, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o+10, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-10, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-9, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-8, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-6, o+1, 1, Wall::Back);
        field.CopyDerivativeFreeConvexCorner(d-7, o+0, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-5, o+1, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-4, o+2, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-3, o+3, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-2, o+4, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-1, o+5, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d+0, o+7, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(d-7, o+1, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-5, o+2, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-4, o+3, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-3, o+4, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-2, o+5, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-1, o+7, Corner::LeftBack);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+10, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+9, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+8, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+1, o+6, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+6, o+1, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+8, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+9, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+10, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+0, o+7, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+1, o+5, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+2, o+4, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+3, o+3, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+4, o+2, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+5, o+1, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+7, o+0, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+1, o+7, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+2, o+5, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+3, o+4, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+4, o+3, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+5, o+2, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+7, o+1, Corner::RightBack);
      };

      apply_u_bc(ux);
      apply_u_bc(uy);
      apply_uphi_bc(ux_phi);
      apply_uphi_bc(uy_phi);
      apply_bc(phi);
      apply_bc(MU);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }
    //reservoir with chimney, radius corners=13
    case 29:
    {
      //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(d, o+14, LY-o-15, Wall::Left);
        field.ApplyNeumannFreeWall(LX-d-1, o+14, LY-o-15, Wall::Right);
        field.ApplyNeumannFreeWall(0, 0, LX, Wall::Front);
        field.CopyDerivativeFreeWall(d+1, LY-1, LX-2*d-2, Wall::Back);
        field.ApplyNeumannFreeConcaveCorner(d, LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeWall(LX-d+13, o, d-13, Wall::Back);
        field.ApplyNeumannFreeWall(0, o, d-13, Wall::Back);
        field.ApplyNeumannFreeWall(d-1, o+8, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o+9, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+11, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+12, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+13, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-13, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-12, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-11, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-9, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-8, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(d-10, o+0, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-7, o+1, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-6, o+2, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-5, o+3, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-4, o+4, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-3, o+5, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-2, o+6, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-1, o+7, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d+0, o+10, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(d-10, o+1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-7, o+2, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-6, o+3, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-5, o+4, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-4, o+5, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-3, o+6, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-2, o+7, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-1, o+10, Corner::LeftBack);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+13, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+12, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+11, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o+9, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o+8, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+8, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+9, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+11, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+12, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+13, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+0, o+10, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+1, o+7, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+2, o+6, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+3, o+5, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+4, o+4, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+5, o+3, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+6, o+2, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+7, o+1, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+10, o+0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+1, o+10, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+2, o+7, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+3, o+6, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+4, o+5, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+5, o+4, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+6, o+3, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+7, o+2, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+10, o+1, Corner::RightBack);
      };

      auto apply_uphi_bc = [this](ScalarField& field) {

        field.ApplyDirichletFreeWall(d, o+14, LY-o-15, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-d-1, o+14, LY-o-15, Wall::Right, 0);
        field.ApplyDirichletFreeWall(0, 0, LX, Wall::Front, 0);
        field.CopyDerivativeFreeWall(d+1, LY-1, LX-2*d-2, Wall::Back);
        field.ApplyDirichletFreeConcaveCorner(d, LY-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeWall(LX-d+13, o, d-13, Wall::Back, 0);
        field.ApplyDirichletFreeWall(0, o, d-13, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-1, o+8, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-1, o+9, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o+11, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o+12, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o+13, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-13, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-12, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-11, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-9, o+1, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-8, o+1, 1, Wall::Back, 0);
        field.ApplyDirichletFreeConvexCorner(d-10, o+0, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-7, o+1, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-6, o+2, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-5, o+3, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-4, o+4, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-3, o+5, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-2, o+6, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-1, o+7, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d+0, o+10, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-10, o+1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-7, o+2, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-6, o+3, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-5, o+4, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-4, o+5, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-3, o+6, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-2, o+7, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-1, o+10, Corner::LeftBack, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+13, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+12, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+11, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+1, o+9, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+1, o+8, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+8, o+1, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+9, o+1, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+11, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+12, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+13, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+0, o+10, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+1, o+7, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+2, o+6, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+3, o+5, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+4, o+4, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+5, o+3, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+6, o+2, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+7, o+1, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+10, o+0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+1, o+10, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+2, o+7, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+3, o+6, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+4, o+5, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+5, o+4, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+6, o+3, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+7, o+2, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+10, o+1, Corner::RightBack, 0);
      };

      auto apply_u_bc = [this](ScalarField& field) {

        field.CopyDerivativeFreeWall(d, o+14, LY-o-15, Wall::Left);
        field.CopyDerivativeFreeWall(LX-d-1, o+14, LY-o-15, Wall::Right);
        field.CopyDerivativeFreeWall(0, 0, LX, Wall::Front);
        field.CopyDerivativeFreeWall(d+1, LY-1, LX-2*d-2, Wall::Back);
        field.CopyDerivativeFreeConcaveCorner(d, LY-1, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack);
        field.CopyDerivativeFreeWall(LX-d+13, o, d-13, Wall::Back);
        field.CopyDerivativeFreeWall(0, o, d-13, Wall::Back);
        field.CopyDerivativeFreeWall(d-1, o+8, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-1, o+9, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o+11, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o+12, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o+13, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-13, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-12, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-11, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-9, o+1, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-8, o+1, 1, Wall::Back);
        field.CopyDerivativeFreeConvexCorner(d-10, o+0, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-7, o+1, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-6, o+2, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-5, o+3, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-4, o+4, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-3, o+5, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-2, o+6, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-1, o+7, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d+0, o+10, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(d-10, o+1, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-7, o+2, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-6, o+3, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-5, o+4, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-4, o+5, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-3, o+6, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-2, o+7, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-1, o+10, Corner::LeftBack);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+13, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+12, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+11, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+1, o+9, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+1, o+8, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+8, o+1, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+9, o+1, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+11, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+12, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+13, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+0, o+10, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+1, o+7, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+2, o+6, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+3, o+5, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+4, o+4, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+5, o+3, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+6, o+2, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+7, o+1, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+10, o+0, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+1, o+10, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+2, o+7, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+3, o+6, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+4, o+5, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+5, o+4, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+6, o+3, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+7, o+2, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+10, o+1, Corner::RightBack);
      };

      apply_u_bc(ux);
      apply_u_bc(uy);
      apply_uphi_bc(ux_phi);
      apply_uphi_bc(uy_phi);
      apply_bc(phi);
      apply_bc(MU);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }
    //reservoir with chimney, radius corners=16
    case 30:
    {
      //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(d, o+17, LY-o-18, Wall::Left);
        field.ApplyNeumannFreeWall(LX-d-1, o+17, LY-o-18, Wall::Right);
        field.ApplyNeumannFreeWall(0, 0, LX, Wall::Front);
        field.CopyDerivativeFreeWall(d+1, LY-1, LX-2*d-2, Wall::Back);
        field.ApplyNeumannFreeConcaveCorner(d, LY-1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack);
        field.ApplyNeumannFreeWall(LX-d+16, o, d-16, Wall::Back);
        field.ApplyNeumannFreeWall(0, o, d-16, Wall::Back);
        field.ApplyNeumannFreeWall(d-4, o+6, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-2, o+9, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o+11, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o+12, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+14, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+15, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+16, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-16, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-15, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-14, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-12, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-11, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-9, o+2, 1, Wall::Back);
        field.ApplyNeumannFreeWall(d-6, o+4, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(d-13, o+0, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-10, o+1, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-8, o+2, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-7, o+3, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-5, o+4, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-4, o+5, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-3, o+7, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-2, o+8, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d-1, o+10, Corner::RightFront);
        field.ApplyNeumannFreeConvexCorner(d+0, o+13, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(d-13, o+1, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-10, o+2, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-8, o+3, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-7, o+4, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-5, o+5, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-4, o+7, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-3, o+8, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-2, o+10, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(d-1, o+13, Corner::LeftBack);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+16, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+15, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+14, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o+12, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o+11, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+2, o+9, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+4, o+6, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+6, o+4, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+9, o+2, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+11, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+12, o+1, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+14, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+15, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeWall(LX-d-1+16, o+0, 1, Wall::Back);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+0, o+13, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+1, o+10, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+2, o+8, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+3, o+7, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+4, o+5, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+5, o+4, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+7, o+3, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+8, o+2, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+10, o+1, Corner::LeftFront);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+13, o+0, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+1, o+13, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+2, o+10, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+3, o+8, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+4, o+7, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+5, o+5, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+7, o+4, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+8, o+3, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+10, o+2, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+13, o+1, Corner::RightBack);
      };

      auto apply_uphi_bc = [this](ScalarField& field) {

        field.ApplyDirichletFreeWall(d, o+17, LY-o-18, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-d-1, o+17, LY-o-18, Wall::Right, 0);
        field.ApplyDirichletFreeWall(0, 0, LX, Wall::Front, 0);
        field.CopyDerivativeFreeWall(d+1, LY-1, LX-2*d-2, Wall::Back);
        field.ApplyDirichletFreeConcaveCorner(d, LY-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeWall(LX-d+16, o, d-16, Wall::Back, 0);
        field.ApplyDirichletFreeWall(0, o, d-16, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-4, o+6, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-2, o+9, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-1, o+11, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-1, o+12, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o+14, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o+15, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o+16, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-16, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-15, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-14, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-12, o+1, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-11, o+1, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-9, o+2, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(d-6, o+4, 1, Wall::Back, 0);
        field.ApplyDirichletFreeConvexCorner(d-13, o+0, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-10, o+1, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-8, o+2, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-7, o+3, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-5, o+4, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-4, o+5, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-3, o+7, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-2, o+8, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d-1, o+10, Corner::RightFront, 0);
        field.ApplyDirichletFreeConvexCorner(d+0, o+13, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-13, o+1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-10, o+2, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-8, o+3, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-7, o+4, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-5, o+5, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-4, o+7, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-3, o+8, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-2, o+10, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-1, o+13, Corner::LeftBack, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+16, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+15, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+14, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+1, o+12, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+1, o+11, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+2, o+9, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+4, o+6, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+6, o+4, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+9, o+2, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+11, o+1, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+12, o+1, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+14, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+15, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeWall(LX-d-1+16, o+0, 1, Wall::Back, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+0, o+13, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+1, o+10, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+2, o+8, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+3, o+7, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+4, o+5, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+5, o+4, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+7, o+3, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+8, o+2, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+10, o+1, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+13, o+0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+1, o+13, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+2, o+10, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+3, o+8, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+4, o+7, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+5, o+5, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+7, o+4, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+8, o+3, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+10, o+2, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+13, o+1, Corner::RightBack, 0);
      };

      auto apply_u_bc = [this](ScalarField& field) {

        field.CopyDerivativeFreeWall(d, o+17, LY-o-18, Wall::Left);
        field.CopyDerivativeFreeWall(LX-d-1, o+17, LY-o-18, Wall::Right);
        field.CopyDerivativeFreeWall(0, 0, LX, Wall::Front);
        field.CopyDerivativeFreeWall(d+1, LY-1, LX-2*d-2, Wall::Back);
        field.CopyDerivativeFreeConcaveCorner(d, LY-1, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1, LY-1, Corner::RightBack);
        field.CopyDerivativeFreeWall(LX-d+16, o, d-16, Wall::Back);
        field.CopyDerivativeFreeWall(0, o, d-16, Wall::Back);
        field.CopyDerivativeFreeWall(d-4, o+6, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-2, o+9, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-1, o+11, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-1, o+12, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o+14, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o+15, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o+16, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-16, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-15, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-14, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-12, o+1, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-11, o+1, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-9, o+2, 1, Wall::Back);
        field.CopyDerivativeFreeWall(d-6, o+4, 1, Wall::Back);
        field.CopyDerivativeFreeConvexCorner(d-13, o+0, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-10, o+1, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-8, o+2, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-7, o+3, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-5, o+4, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-4, o+5, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-3, o+7, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-2, o+8, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d-1, o+10, Corner::RightFront);
        field.CopyDerivativeFreeConvexCorner(d+0, o+13, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(d-13, o+1, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-10, o+2, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-8, o+3, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-7, o+4, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-5, o+5, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-4, o+7, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-3, o+8, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-2, o+10, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(d-1, o+13, Corner::LeftBack);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+16, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+15, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+14, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+1, o+12, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+1, o+11, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+2, o+9, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+4, o+6, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+6, o+4, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+9, o+2, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+11, o+1, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+12, o+1, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+14, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+15, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeWall(LX-d-1+16, o+0, 1, Wall::Back);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+0, o+13, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+1, o+10, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+2, o+8, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+3, o+7, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+4, o+5, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+5, o+4, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+7, o+3, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+8, o+2, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+10, o+1, Corner::LeftFront);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+13, o+0, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+1, o+13, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+2, o+10, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+3, o+8, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+4, o+7, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+5, o+5, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+7, o+4, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+8, o+3, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+10, o+2, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+13, o+1, Corner::RightBack);
      };

      apply_u_bc(ux);
      apply_u_bc(uy);
      apply_uphi_bc(ux_phi);
      apply_uphi_bc(uy_phi);
      apply_bc(phi);
      apply_bc(MU);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }
    //reservoir with chimney exit
    case 31:
    {
      //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeConcaveCorner(LX-d-1, 0, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.ApplyNeumannFreeWall(d+1, 0, LX-2*d-2, Wall::Front);
        field.ApplyNeumannFreeWall(d, 1, o, Wall::Left);
        field.ApplyNeumannFreeWall(LX-d-1, 1, o, Wall::Right);
        field.ApplyNeumannFreeConvexCorner(LX-d-1, o+1, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(d, o+1, Corner::RightBack);
        field.ApplyNeumannFreeWall(LX-d, o+1, d, Wall::Front);
        field.ApplyNeumannFreeWall(0, o+1, d, Wall::Front);
        field.CopyDerivativeFreeWall(0, LY-1, LX, Wall::Back);
      };

      auto apply_uphi_bc = [this](ScalarField& field) {

        field.ApplyDirichletFreeConcaveCorner(LX-d-1, 0, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d, 0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeWall(d+1, 0, LX-2*d-2, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d, 1, o, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-d-1, 1, o, Wall::Right, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1, o+1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(d, o+1, Corner::RightBack, 0);
        field.ApplyDirichletFreeWall(LX-d, o+1, d, Wall::Front, 0);
        field.ApplyDirichletFreeWall(0, o+1, d, Wall::Front, 0);
        field.CopyDerivativeFreeWall(0, LY-1, LX, Wall::Back);
      };

      auto apply_u_bc = [this](ScalarField& field) {

        field.CopyDerivativeFreeConcaveCorner(LX-d-1, 0, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.CopyDerivativeFreeWall(d+1, 0, LX-2*d-2, Wall::Front);
        field.CopyDerivativeFreeWall(d, 1, o, Wall::Left);
        field.CopyDerivativeFreeWall(LX-d-1, 1, o, Wall::Right);
        field.CopyDerivativeFreeConvexCorner(LX-d-1, o+1, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(d, o+1, Corner::RightBack);
        field.CopyDerivativeFreeWall(LX-d, o+1, d, Wall::Front);
        field.CopyDerivativeFreeWall(0, o+1, d, Wall::Front);
        field.CopyDerivativeFreeWall(0, LY-1, LX, Wall::Back);
      };

      apply_u_bc(ux);
      apply_u_bc(uy);
      apply_uphi_bc(ux_phi);
      apply_uphi_bc(uy_phi);
      apply_bc(phi);
      apply_bc(MU);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }



    //reservoir with chimney exit, radius corners=4
    case 32:
    {
      //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(d+0, o+1-4, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-4, o+1+0, 1, Wall::Front);
        field.ApplyNeumannFreeConvexCorner(d+0, o+1-3, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-1, o+1-2, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-2, o+1-1, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-3, o+1+0, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(d-1, o+1-3, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-2, o+1-2, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-3, o+1-1, Corner::LeftFront);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+1-4, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+4, o+1+0, 1, Wall::Front);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+3, o+1+0, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+2, o+1-1, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+1, o+1-2, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+0, o+1-3, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+3, o+1-1, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+2, o+1-2, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+1, o+1-3, Corner::RightFront);
        field.ApplyNeumannFreeWall(d+1, 0, LX-2*d-2, Wall::Front);
        field.ApplyNeumannFreeWall(d, 1, o-4, Wall::Left);
        field.ApplyNeumannFreeWall(LX-d-1, 1, o-4, Wall::Right);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1, 0, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.ApplyNeumannFreeWall(LX-d+4, o+1, d-4, Wall::Front);
        field.ApplyNeumannFreeWall(0, o+1, d-4, Wall::Front);
        field.CopyDerivativeFreeWall(0, LY-1, LX, Wall::Back);
      };

      auto apply_uphi_bc = [this](ScalarField& field) {

        field.ApplyDirichletFreeWall(d+0, o+1-4, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-4, o+1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeConvexCorner(d+0, o+1-3, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-1, o+1-2, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-2, o+1-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-3, o+1+0, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-1, o+1-3, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-2, o+1-2, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-3, o+1-1, Corner::LeftFront, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+1-4, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+4, o+1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+3, o+1+0, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+2, o+1-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+1, o+1-2, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+0, o+1-3, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+3, o+1-1, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+2, o+1-2, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+1, o+1-3, Corner::RightFront, 0);
        field.ApplyDirichletFreeWall(d+1, 0, LX-2*d-2, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d, 1, o-4, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-d-1, 1, o-4, Wall::Right, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1, 0, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d, 0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeWall(LX-d+4, o+1, d-4, Wall::Front, 0);
        field.ApplyDirichletFreeWall(0, o+1, d-4, Wall::Front, 0);
        field.CopyDerivativeFreeWall(0, LY-1, LX, Wall::Back);
      };

      auto apply_u_bc = [this](ScalarField& field) {

        field.CopyDerivativeFreeWall(d+0, o+1-4, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-4, o+1+0, 1, Wall::Front);
        field.CopyDerivativeFreeConvexCorner(d+0, o+1-3, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-1, o+1-2, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-2, o+1-1, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-3, o+1+0, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(d-1, o+1-3, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-2, o+1-2, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-3, o+1-1, Corner::LeftFront);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+1-4, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+4, o+1+0, 1, Wall::Front);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+3, o+1+0, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+2, o+1-1, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+1, o+1-2, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+0, o+1-3, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+3, o+1-1, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+2, o+1-2, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+1, o+1-3, Corner::RightFront);
        field.CopyDerivativeFreeWall(d+1, 0, LX-2*d-2, Wall::Front);
        field.CopyDerivativeFreeWall(d, 1, o-4, Wall::Left);
        field.CopyDerivativeFreeWall(LX-d-1, 1, o-4, Wall::Right);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1, 0, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.CopyDerivativeFreeWall(LX-d+4, o+1, d-4, Wall::Front);
        field.CopyDerivativeFreeWall(0, o+1, d-4, Wall::Front);
        field.CopyDerivativeFreeWall(0, LY-1, LX, Wall::Back);
      };

      apply_u_bc(ux);
      apply_u_bc(uy);
      apply_uphi_bc(ux_phi);
      apply_uphi_bc(uy_phi);
      apply_bc(phi);
      apply_bc(MU);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }
    //reservoir with chimney exit, radius corners=7
    case 33:
    {
      //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(d+0, o+1-7, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+1-6, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o+1-4, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-4, o+1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-6, o+1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-7, o+1+0, 1, Wall::Front);
        field.ApplyNeumannFreeConvexCorner(d+0, o+1-5, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-1, o+1-3, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-2, o+1-2, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-3, o+1-1, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-5, o+1+0, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(d-1, o+1-5, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-2, o+1-3, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-3, o+1-2, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-5, o+1-1, Corner::LeftFront);
        field.ApplyNeumannFreeWall(LX-d-1+1, o+1-4, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+1-6, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+1-7, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+7, o+1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+6, o+1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+4, o+1-1, 1, Wall::Front);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+5, o+1+0, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+3, o+1-1, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+2, o+1-2, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+1, o+1-3, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+0, o+1-5, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+5, o+1-1, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+3, o+1-2, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+2, o+1-3, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+1, o+1-5, Corner::RightFront);
        field.ApplyNeumannFreeWall(d+1, 0, LX-2*d-2, Wall::Front);
        field.ApplyNeumannFreeWall(d, 1, o-7, Wall::Left);
        field.ApplyNeumannFreeWall(LX-d-1, 1, o-7, Wall::Right);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1, 0, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.ApplyNeumannFreeWall(LX-d+7, o+1, d-7, Wall::Front);
        field.ApplyNeumannFreeWall(0, o+1, d-7, Wall::Front);
        field.CopyDerivativeFreeWall(0, LY-1, LX, Wall::Back);
      };

      auto apply_uphi_bc = [this](ScalarField& field) {

        field.ApplyDirichletFreeWall(d+0, o+1-7, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o+1-6, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-1, o+1-4, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-4, o+1-1, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-6, o+1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-7, o+1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeConvexCorner(d+0, o+1-5, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-1, o+1-3, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-2, o+1-2, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-3, o+1-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-5, o+1+0, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-1, o+1-5, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-2, o+1-3, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-3, o+1-2, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-5, o+1-1, Corner::LeftFront, 0);
        field.ApplyDirichletFreeWall(LX-d-1+1, o+1-4, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+1-6, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+1-7, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+7, o+1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+6, o+1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+4, o+1-1, 1, Wall::Front, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+5, o+1+0, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+3, o+1-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+2, o+1-2, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+1, o+1-3, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+0, o+1-5, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+5, o+1-1, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+3, o+1-2, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+2, o+1-3, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+1, o+1-5, Corner::RightFront, 0);
        field.ApplyDirichletFreeWall(d+1, 0, LX-2*d-2, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d, 1, o-7, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-d-1, 1, o-7, Wall::Right, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1, 0, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d, 0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeWall(LX-d+7, o+1, d-7, Wall::Front, 0);
        field.ApplyDirichletFreeWall(0, o+1, d-7, Wall::Front, 0);
        field.CopyDerivativeFreeWall(0, LY-1, LX, Wall::Back);
      };

      auto apply_u_bc = [this](ScalarField& field) {

        field.CopyDerivativeFreeWall(d+0, o+1-7, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o+1-6, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-1, o+1-4, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-4, o+1-1, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-6, o+1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-7, o+1+0, 1, Wall::Front);
        field.CopyDerivativeFreeConvexCorner(d+0, o+1-5, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-1, o+1-3, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-2, o+1-2, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-3, o+1-1, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-5, o+1+0, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(d-1, o+1-5, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-2, o+1-3, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-3, o+1-2, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-5, o+1-1, Corner::LeftFront);
        field.CopyDerivativeFreeWall(LX-d-1+1, o+1-4, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+1-6, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+1-7, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+7, o+1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+6, o+1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+4, o+1-1, 1, Wall::Front);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+5, o+1+0, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+3, o+1-1, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+2, o+1-2, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+1, o+1-3, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+0, o+1-5, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+5, o+1-1, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+3, o+1-2, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+2, o+1-3, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+1, o+1-5, Corner::RightFront);
        field.CopyDerivativeFreeWall(d+1, 0, LX-2*d-2, Wall::Front);
        field.CopyDerivativeFreeWall(d, 1, o-7, Wall::Left);
        field.CopyDerivativeFreeWall(LX-d-1, 1, o-7, Wall::Right);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1, 0, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.CopyDerivativeFreeWall(LX-d+7, o+1, d-7, Wall::Front);
        field.CopyDerivativeFreeWall(0, o+1, d-7, Wall::Front);
        field.CopyDerivativeFreeWall(0, LY-1, LX, Wall::Back);
      };

      apply_u_bc(ux);
      apply_u_bc(uy);
      apply_uphi_bc(ux_phi);
      apply_uphi_bc(uy_phi);
      apply_bc(phi);
      apply_bc(MU);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }
    //reservoir with chimney exit, radius corners=10
    case 34:
    {
      //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(d+0, o+1-10, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+1-9, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+1-8, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o+1-6, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-6, o+1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-8, o+1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-9, o+1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-10, o+1+0, 1, Wall::Front);
        field.ApplyNeumannFreeConvexCorner(d+0, o+1-7, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-1, o+1-5, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-2, o+1-4, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-3, o+1-3, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-4, o+1-2, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-5, o+1-1, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-7, o+1+0, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(d-1, o+1-7, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-2, o+1-5, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-3, o+1-4, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-4, o+1-3, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-5, o+1-2, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-7, o+1-1, Corner::LeftFront);
        field.ApplyNeumannFreeWall(LX-d-1+1, o+1-6, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+1-8, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+1-9, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+1-10, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+10, o+1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+9, o+1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+8, o+1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+6, o+1-1, 1, Wall::Front);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+7, o+1+0, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+5, o+1-1, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+4, o+1-2, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+3, o+1-3, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+2, o+1-4, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+1, o+1-5, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+0, o+1-7, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+7, o+1-1, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+5, o+1-2, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+4, o+1-3, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+3, o+1-4, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+2, o+1-5, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+1, o+1-7, Corner::RightFront);
        field.ApplyNeumannFreeWall(d+1, 0, LX-2*d-2, Wall::Front);
        field.ApplyNeumannFreeWall(d, 1, o-10, Wall::Left);
        field.ApplyNeumannFreeWall(LX-d-1, 1, o-10, Wall::Right);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1, 0, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.ApplyNeumannFreeWall(LX-d+10, o+1, d-10, Wall::Front);
        field.ApplyNeumannFreeWall(0, o+1, d-10, Wall::Front);
        field.CopyDerivativeFreeWall(0, LY-1, LX, Wall::Back);
      };

      auto apply_uphi_bc = [this](ScalarField& field) {

        field.ApplyDirichletFreeWall(d+0, o+1-10, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o+1-9, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o+1-8, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-1, o+1-6, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-6, o+1-1, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-8, o+1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-9, o+1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-10, o+1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeConvexCorner(d+0, o+1-7, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-1, o+1-5, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-2, o+1-4, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-3, o+1-3, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-4, o+1-2, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-5, o+1-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-7, o+1+0, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-1, o+1-7, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-2, o+1-5, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-3, o+1-4, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-4, o+1-3, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-5, o+1-2, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-7, o+1-1, Corner::LeftFront, 0);
        field.ApplyDirichletFreeWall(LX-d-1+1, o+1-6, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+1-8, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+1-9, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+1-10, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+10, o+1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+9, o+1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+8, o+1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+6, o+1-1, 1, Wall::Front, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+7, o+1+0, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+5, o+1-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+4, o+1-2, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+3, o+1-3, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+2, o+1-4, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+1, o+1-5, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+0, o+1-7, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+7, o+1-1, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+5, o+1-2, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+4, o+1-3, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+3, o+1-4, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+2, o+1-5, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+1, o+1-7, Corner::RightFront, 0);
        field.ApplyDirichletFreeWall(d+1, 0, LX-2*d-2, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d, 1, o-10, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-d-1, 1, o-10, Wall::Right, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1, 0, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d, 0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeWall(LX-d+10, o+1, d-10, Wall::Front, 0);
        field.ApplyDirichletFreeWall(0, o+1, d-10, Wall::Front, 0);
        field.CopyDerivativeFreeWall(0, LY-1, LX, Wall::Back);
      };

      auto apply_u_bc = [this](ScalarField& field) {

        field.CopyDerivativeFreeWall(d+0, o+1-10, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o+1-9, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o+1-8, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-1, o+1-6, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-6, o+1-1, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-8, o+1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-9, o+1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-10, o+1+0, 1, Wall::Front);
        field.CopyDerivativeFreeConvexCorner(d+0, o+1-7, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-1, o+1-5, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-2, o+1-4, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-3, o+1-3, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-4, o+1-2, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-5, o+1-1, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-7, o+1+0, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(d-1, o+1-7, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-2, o+1-5, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-3, o+1-4, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-4, o+1-3, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-5, o+1-2, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-7, o+1-1, Corner::LeftFront);
        field.CopyDerivativeFreeWall(LX-d-1+1, o+1-6, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+1-8, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+1-9, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+1-10, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+10, o+1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+9, o+1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+8, o+1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+6, o+1-1, 1, Wall::Front);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+7, o+1+0, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+5, o+1-1, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+4, o+1-2, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+3, o+1-3, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+2, o+1-4, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+1, o+1-5, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+0, o+1-7, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+7, o+1-1, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+5, o+1-2, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+4, o+1-3, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+3, o+1-4, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+2, o+1-5, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+1, o+1-7, Corner::RightFront);
        field.CopyDerivativeFreeWall(d+1, 0, LX-2*d-2, Wall::Front);
        field.CopyDerivativeFreeWall(d, 1, o-10, Wall::Left);
        field.CopyDerivativeFreeWall(LX-d-1, 1, o-10, Wall::Right);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1, 0, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.CopyDerivativeFreeWall(LX-d+10, o+1, d-10, Wall::Front);
        field.CopyDerivativeFreeWall(0, o+1, d-10, Wall::Front);
        field.CopyDerivativeFreeWall(0, LY-1, LX, Wall::Back);
      };

      apply_u_bc(ux);
      apply_u_bc(uy);
      apply_uphi_bc(ux_phi);
      apply_uphi_bc(uy_phi);
      apply_bc(phi);
      apply_bc(MU);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }
    //reservoir with chimney exit, radius corners=13
    case 35:
    {
      //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(d+0, o+1-13, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+1-12, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+1-11, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o+1-9, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o+1-8, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-8, o+1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-9, o+1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-11, o+1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-12, o+1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-13, o+1+0, 1, Wall::Front);
        field.ApplyNeumannFreeConvexCorner(d+0, o+1-10, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-1, o+1-7, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-2, o+1-6, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-3, o+1-5, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-4, o+1-4, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-5, o+1-3, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-6, o+1-2, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-7, o+1-1, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-10, o+1+0, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(d-1, o+1-10, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-2, o+1-7, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-3, o+1-6, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-4, o+1-5, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-5, o+1-4, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-6, o+1-3, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-7, o+1-2, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-10, o+1-1, Corner::LeftFront);
        field.ApplyNeumannFreeWall(LX-d-1+1, o+1-8, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o+1-9, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+1-11, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+1-12, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+1-13, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+13, o+1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+12, o+1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+11, o+1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+9, o+1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+8, o+1-1, 1, Wall::Front);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+10, o+1+0, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+7, o+1-1, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+6, o+1-2, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+5, o+1-3, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+4, o+1-4, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+3, o+1-5, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+2, o+1-6, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+1, o+1-7, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+0, o+1-10, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+10, o+1-1, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+7, o+1-2, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+6, o+1-3, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+5, o+1-4, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+4, o+1-5, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+3, o+1-6, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+2, o+1-7, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+1, o+1-10, Corner::RightFront);
        field.ApplyNeumannFreeWall(d+1, 0, LX-2*d-2, Wall::Front);
        field.ApplyNeumannFreeWall(d, 1, o-13, Wall::Left);
        field.ApplyNeumannFreeWall(LX-d-1, 1, o-13, Wall::Right);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1, 0, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.ApplyNeumannFreeWall(LX-d+13, o+1, d-13, Wall::Front);
        field.ApplyNeumannFreeWall(0, o+1, d-13, Wall::Front);
        field.CopyDerivativeFreeWall(0, LY-1, LX, Wall::Back);
      };

      auto apply_uphi_bc = [this](ScalarField& field) {

        field.ApplyDirichletFreeWall(d+0, o+1-13, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o+1-12, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o+1-11, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-1, o+1-9, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-1, o+1-8, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-8, o+1-1, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-9, o+1-1, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-11, o+1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-12, o+1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-13, o+1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeConvexCorner(d+0, o+1-10, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-1, o+1-7, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-2, o+1-6, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-3, o+1-5, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-4, o+1-4, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-5, o+1-3, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-6, o+1-2, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-7, o+1-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-10, o+1+0, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-1, o+1-10, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-2, o+1-7, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-3, o+1-6, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-4, o+1-5, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-5, o+1-4, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-6, o+1-3, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-7, o+1-2, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-10, o+1-1, Corner::LeftFront, 0);
        field.ApplyDirichletFreeWall(LX-d-1+1, o+1-8, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+1, o+1-9, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+1-11, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+1-12, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+1-13, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+13, o+1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+12, o+1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+11, o+1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+9, o+1-1, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+8, o+1-1, 1, Wall::Front, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+10, o+1+0, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+7, o+1-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+6, o+1-2, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+5, o+1-3, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+4, o+1-4, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+3, o+1-5, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+2, o+1-6, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+1, o+1-7, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+0, o+1-10, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+10, o+1-1, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+7, o+1-2, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+6, o+1-3, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+5, o+1-4, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+4, o+1-5, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+3, o+1-6, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+2, o+1-7, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+1, o+1-10, Corner::RightFront, 0);
        field.ApplyDirichletFreeWall(d+1, 0, LX-2*d-2, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d, 1, o-13, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-d-1, 1, o-13, Wall::Right, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1, 0, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d, 0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeWall(LX-d+13, o+1, d-13, Wall::Front, 0);
        field.ApplyDirichletFreeWall(0, o+1, d-13, Wall::Front, 0);
        field.CopyDerivativeFreeWall(0, LY-1, LX, Wall::Back);
      };

      auto apply_u_bc = [this](ScalarField& field) {

        field.CopyDerivativeFreeWall(d+0, o+1-13, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o+1-12, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o+1-11, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-1, o+1-9, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-1, o+1-8, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-8, o+1-1, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-9, o+1-1, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-11, o+1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-12, o+1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-13, o+1+0, 1, Wall::Front);
        field.CopyDerivativeFreeConvexCorner(d+0, o+1-10, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-1, o+1-7, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-2, o+1-6, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-3, o+1-5, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-4, o+1-4, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-5, o+1-3, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-6, o+1-2, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-7, o+1-1, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-10, o+1+0, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(d-1, o+1-10, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-2, o+1-7, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-3, o+1-6, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-4, o+1-5, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-5, o+1-4, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-6, o+1-3, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-7, o+1-2, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-10, o+1-1, Corner::LeftFront);
        field.CopyDerivativeFreeWall(LX-d-1+1, o+1-8, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+1, o+1-9, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+1-11, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+1-12, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+1-13, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+13, o+1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+12, o+1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+11, o+1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+9, o+1-1, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+8, o+1-1, 1, Wall::Front);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+10, o+1+0, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+7, o+1-1, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+6, o+1-2, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+5, o+1-3, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+4, o+1-4, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+3, o+1-5, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+2, o+1-6, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+1, o+1-7, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+0, o+1-10, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+10, o+1-1, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+7, o+1-2, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+6, o+1-3, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+5, o+1-4, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+4, o+1-5, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+3, o+1-6, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+2, o+1-7, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+1, o+1-10, Corner::RightFront);
        field.CopyDerivativeFreeWall(d+1, 0, LX-2*d-2, Wall::Front);
        field.CopyDerivativeFreeWall(d, 1, o-13, Wall::Left);
        field.CopyDerivativeFreeWall(LX-d-1, 1, o-13, Wall::Right);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1, 0, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.CopyDerivativeFreeWall(LX-d+13, o+1, d-13, Wall::Front);
        field.CopyDerivativeFreeWall(0, o+1, d-13, Wall::Front);
        field.CopyDerivativeFreeWall(0, LY-1, LX, Wall::Back);
      };

      apply_u_bc(ux);
      apply_u_bc(uy);
      apply_uphi_bc(ux_phi);
      apply_uphi_bc(uy_phi);
      apply_bc(phi);
      apply_bc(MU);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }
    //reservoir with chimney exit, radius corners=16
    case 36:
    {
      //automatically generated code
      auto apply_bc = [this](ScalarField& field) {

        field.ApplyNeumannFreeWall(d+0, o+1-16, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+1-15, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d+0, o+1-14, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o+1-12, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-1, o+1-11, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-2, o+1-9, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-4, o+1-6, 1, Wall::Left);
        field.ApplyNeumannFreeWall(d-6, o+1-4, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-9, o+1-2, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-11, o+1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-12, o+1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-14, o+1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-15, o+1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(d-16, o+1+0, 1, Wall::Front);
        field.ApplyNeumannFreeConvexCorner(d+0, o+1-13, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-1, o+1-10, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-2, o+1-8, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-3, o+1-7, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-4, o+1-5, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-5, o+1-4, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-7, o+1-3, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-8, o+1-2, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-10, o+1-1, Corner::RightBack);
        field.ApplyNeumannFreeConvexCorner(d-13, o+1+0, Corner::RightBack);
        field.ApplyNeumannFreeConcaveCorner(d-1, o+1-13, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-2, o+1-10, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-3, o+1-8, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-4, o+1-7, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-5, o+1-5, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-7, o+1-4, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-8, o+1-3, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-10, o+1-2, Corner::LeftFront);
        field.ApplyNeumannFreeConcaveCorner(d-13, o+1-1, Corner::LeftFront);
        field.ApplyNeumannFreeWall(LX-d-1+4, o+1-6, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+2, o+1-9, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o+1-11, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+1, o+1-12, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+1-14, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+1-15, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+0, o+1-16, 1, Wall::Right);
        field.ApplyNeumannFreeWall(LX-d-1+16, o+1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+15, o+1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+14, o+1+0, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+12, o+1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+11, o+1-1, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+9, o+1-2, 1, Wall::Front);
        field.ApplyNeumannFreeWall(LX-d-1+6, o+1-4, 1, Wall::Front);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+13, o+1+0, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+10, o+1-1, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+8, o+1-2, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+7, o+1-3, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+5, o+1-4, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+4, o+1-5, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+3, o+1-7, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+2, o+1-8, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+1, o+1-10, Corner::LeftBack);
        field.ApplyNeumannFreeConvexCorner(LX-d-1+0, o+1-13, Corner::LeftBack);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+13, o+1-1, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+10, o+1-2, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+8, o+1-3, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+7, o+1-4, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+5, o+1-5, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+4, o+1-7, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+3, o+1-8, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+2, o+1-10, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1+1, o+1-13, Corner::RightFront);
        field.ApplyNeumannFreeWall(d+1, 0, LX-2*d-2, Wall::Front);
        field.ApplyNeumannFreeWall(d, 1, o-16, Wall::Left);
        field.ApplyNeumannFreeWall(LX-d-1, 1, o-16, Wall::Right);
        field.ApplyNeumannFreeConcaveCorner(LX-d-1, 0, Corner::RightFront);
        field.ApplyNeumannFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.ApplyNeumannFreeWall(LX-d+16, o+1, d-16, Wall::Front);
        field.ApplyNeumannFreeWall(0, o+1, d-16, Wall::Front);
        field.CopyDerivativeFreeWall(0, LY-1, LX, Wall::Back);
      };

      auto apply_uphi_bc = [this](ScalarField& field) {

        field.ApplyDirichletFreeWall(d+0, o+1-16, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o+1-15, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d+0, o+1-14, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-1, o+1-12, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-1, o+1-11, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-2, o+1-9, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-4, o+1-6, 1, Wall::Left, 0);
        field.ApplyDirichletFreeWall(d-6, o+1-4, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-9, o+1-2, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-11, o+1-1, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-12, o+1-1, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-14, o+1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-15, o+1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d-16, o+1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeConvexCorner(d+0, o+1-13, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-1, o+1-10, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-2, o+1-8, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-3, o+1-7, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-4, o+1-5, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-5, o+1-4, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-7, o+1-3, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-8, o+1-2, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-10, o+1-1, Corner::RightBack, 0);
        field.ApplyDirichletFreeConvexCorner(d-13, o+1+0, Corner::RightBack, 0);
        field.ApplyDirichletFreeConcaveCorner(d-1, o+1-13, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-2, o+1-10, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-3, o+1-8, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-4, o+1-7, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-5, o+1-5, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-7, o+1-4, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-8, o+1-3, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-10, o+1-2, Corner::LeftFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d-13, o+1-1, Corner::LeftFront, 0);
        field.ApplyDirichletFreeWall(LX-d-1+4, o+1-6, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+2, o+1-9, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+1, o+1-11, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+1, o+1-12, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+1-14, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+1-15, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+0, o+1-16, 1, Wall::Right, 0);
        field.ApplyDirichletFreeWall(LX-d-1+16, o+1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+15, o+1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+14, o+1+0, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+12, o+1-1, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+11, o+1-1, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+9, o+1-2, 1, Wall::Front, 0);
        field.ApplyDirichletFreeWall(LX-d-1+6, o+1-4, 1, Wall::Front, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+13, o+1+0, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+10, o+1-1, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+8, o+1-2, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+7, o+1-3, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+5, o+1-4, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+4, o+1-5, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+3, o+1-7, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+2, o+1-8, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+1, o+1-10, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConvexCorner(LX-d-1+0, o+1-13, Corner::LeftBack, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+13, o+1-1, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+10, o+1-2, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+8, o+1-3, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+7, o+1-4, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+5, o+1-5, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+4, o+1-7, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+3, o+1-8, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+2, o+1-10, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1+1, o+1-13, Corner::RightFront, 0);
        field.ApplyDirichletFreeWall(d+1, 0, LX-2*d-2, Wall::Front, 0);
        field.ApplyDirichletFreeWall(d, 1, o-16, Wall::Left, 0);
        field.ApplyDirichletFreeWall(LX-d-1, 1, o-16, Wall::Right, 0);
        field.ApplyDirichletFreeConcaveCorner(LX-d-1, 0, Corner::RightFront, 0);
        field.ApplyDirichletFreeConcaveCorner(d, 0, Corner::LeftFront, 0);
        field.ApplyDirichletFreeWall(LX-d+16, o+1, d-16, Wall::Front, 0);
        field.ApplyDirichletFreeWall(0, o+1, d-16, Wall::Front, 0);
        field.CopyDerivativeFreeWall(0, LY-1, LX, Wall::Back);
      };

      auto apply_u_bc = [this](ScalarField& field) {

        field.CopyDerivativeFreeWall(d+0, o+1-16, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o+1-15, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d+0, o+1-14, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-1, o+1-12, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-1, o+1-11, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-2, o+1-9, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-4, o+1-6, 1, Wall::Left);
        field.CopyDerivativeFreeWall(d-6, o+1-4, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-9, o+1-2, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-11, o+1-1, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-12, o+1-1, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-14, o+1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-15, o+1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(d-16, o+1+0, 1, Wall::Front);
        field.CopyDerivativeFreeConvexCorner(d+0, o+1-13, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-1, o+1-10, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-2, o+1-8, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-3, o+1-7, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-4, o+1-5, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-5, o+1-4, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-7, o+1-3, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-8, o+1-2, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-10, o+1-1, Corner::RightBack);
        field.CopyDerivativeFreeConvexCorner(d-13, o+1+0, Corner::RightBack);
        field.CopyDerivativeFreeConcaveCorner(d-1, o+1-13, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-2, o+1-10, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-3, o+1-8, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-4, o+1-7, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-5, o+1-5, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-7, o+1-4, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-8, o+1-3, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-10, o+1-2, Corner::LeftFront);
        field.CopyDerivativeFreeConcaveCorner(d-13, o+1-1, Corner::LeftFront);
        field.CopyDerivativeFreeWall(LX-d-1+4, o+1-6, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+2, o+1-9, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+1, o+1-11, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+1, o+1-12, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+1-14, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+1-15, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+0, o+1-16, 1, Wall::Right);
        field.CopyDerivativeFreeWall(LX-d-1+16, o+1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+15, o+1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+14, o+1+0, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+12, o+1-1, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+11, o+1-1, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+9, o+1-2, 1, Wall::Front);
        field.CopyDerivativeFreeWall(LX-d-1+6, o+1-4, 1, Wall::Front);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+13, o+1+0, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+10, o+1-1, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+8, o+1-2, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+7, o+1-3, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+5, o+1-4, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+4, o+1-5, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+3, o+1-7, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+2, o+1-8, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+1, o+1-10, Corner::LeftBack);
        field.CopyDerivativeFreeConvexCorner(LX-d-1+0, o+1-13, Corner::LeftBack);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+13, o+1-1, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+10, o+1-2, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+8, o+1-3, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+7, o+1-4, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+5, o+1-5, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+4, o+1-7, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+3, o+1-8, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+2, o+1-10, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1+1, o+1-13, Corner::RightFront);
        field.CopyDerivativeFreeWall(d+1, 0, LX-2*d-2, Wall::Front);
        field.CopyDerivativeFreeWall(d, 1, o-16, Wall::Left);
        field.CopyDerivativeFreeWall(LX-d-1, 1, o-16, Wall::Right);
        field.CopyDerivativeFreeConcaveCorner(LX-d-1, 0, Corner::RightFront);
        field.CopyDerivativeFreeConcaveCorner(d, 0, Corner::LeftFront);
        field.CopyDerivativeFreeWall(LX-d+16, o+1, d-16, Wall::Front);
        field.CopyDerivativeFreeWall(0, o+1, d-16, Wall::Front);
        field.CopyDerivativeFreeWall(0, LY-1, LX, Wall::Back);
      };

      apply_u_bc(ux);
      apply_u_bc(uy);
      apply_uphi_bc(ux_phi);
      apply_uphi_bc(uy_phi);
      apply_bc(phi);
      apply_bc(MU);
      apply_bc(sigmaXX);
      apply_bc(sigmaYY);
      apply_bc(sigmaYX);
      apply_bc(sigmaXY);
      break;
    }
  }
}

void LyotropicFreeBoundary::RuntimeChecks()
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
      case 4:
      {
        test[0]=phi.CheckNeumannFreeWall(0, 1, LY-2, Wall::Left);
        test[1]=phi.CheckNeumannFreeWall(LX-1, 1, LY-2, Wall::Right);
        test[2]=phi.CheckNeumannFreeWall(1, 0, LX-2, Wall::Front);
        test[3]=phi.CheckNeumannFreeWall(1, LY-1, LX-2, Wall::Back);

        std::cout << "\nCheck of bc for field (incl corners):\nLeft: "
                  << test[0] << ", Right: " << test[1] << ", Front: "
                  << test[2] << ", Back: " << test[3];
      break;
      }
      case 2:
      case 5:
      {
        test[0]=uy.CheckNeumannFreeWall(0, 0, LX, Wall::Front);
        test[1]=uy.CheckNeumannFreeWall(0, LY-1, LX, Wall::Back);
        test[2]=phi.CheckNeumannFreeWall(0, 0, LX, Wall::Front);
        test[3]=phi.CheckNeumannFreeWall(0, LY-1, LX, Wall::Back);

        std::cout << "\nCheck of bc for field:\nuy front: "
                  << test[0] << ", uy back: " << test[1] << ", phi front: "
                  << test[2] << ", phi back: " << test[3];
      break;
      }
      case 201:
      case 501:
      {
        test[0]=uy.CheckNeumannFreeWall(0, 0, LY, Wall::Left);
        test[1]=uy.CheckNeumannFreeWall(LX-1, 0, LY, Wall::Right);
        test[2]=phi.CheckNeumannFreeWall(0, 0, LY, Wall::Left);
        test[3]=phi.CheckNeumannFreeWall(LX-1, 0, LY, Wall::Right);

        std::cout << "\nCheck of bc for field:\nuy front: "
                  << test[0] << ", uy back: " << test[1] << ", phi front: "
                  << test[2] << ", phi back: " << test[3];
      break;
      }
      case 3:
      case 6:
      {
        test[0]=QQxx.CheckNeumannFreeWall(40, 10, 21, Wall::Left);
        test[1]=QQxx.CheckNeumannFreeWall(10, 10, 21, Wall::Right);
        test[2]=QQxx.CheckNeumannFreeWall(10, 30, 31, Wall::Front);
        test[3]=QQxx.CheckNeumannFreeWall(10, 10, 31, Wall::Back);

        std::cout << "\nCheck of bc for field (incl corners):\nLeft: "
                  << test[0] << ", Right: " << test[1] << ", Front: "
                  << test[2] << ", Back: " << test[3];
      break;
      }
      case 7:
      case 8:
      {
        test[0]=QQxx.CheckNeumannFreeWall(40, 10, 21, Wall::Left);
        test[1]=QQxx.CheckNeumannFreeWall(10, 10, 21, Wall::Right);
        test[2]=QQxx.CheckNeumannFreeWall(10, 30, 31, Wall::Front);
        test[3]=QQxx.CheckNeumannFreeWall(10, 10, 31, Wall::Back);


        test[0]= phi.CheckNeumannFreeWall(0, 1,    o-1, Wall::Left)
                +phi.CheckNeumannFreeWall(d, o,    l,   Wall::Left)
                +phi.CheckNeumannFreeWall(0, o2+1, r-1, Wall::Left);

        test[1]= phi.CheckNeumannFreeWall(LX-1,   1,    o-1, Wall::Right)
                +phi.CheckNeumannFreeWall(LX-d-1, o,    l,   Wall::Right)
                +phi.CheckNeumannFreeWall(LX-1,   o2+1, r-1, Wall::Right);

        test[2]= phi.CheckNeumannFreeWall(1,      0,    LX-2, Wall::Front)
                +phi.CheckNeumannFreeWall(1,      o2-1, d,    Wall::Front)
                +phi.CheckNeumannFreeWall(LX-d-1, o2-1, d,    Wall::Front);

        test[3]= phi.CheckNeumannFreeWall(1,      LY-1, LX-2, Wall::Back)
                +phi.CheckNeumannFreeWall(1,      o,    d,    Wall::Back)
                +phi.CheckNeumannFreeWall(LX-d-1, o,    d,    Wall::Back);


        std::cout << "\nCheck of bc for field (incl corners):\nLeft: "
                  << test[0] << ", Right: " << test[1] << ", Front: "
                  << test[2] << ", Back: " << test[3];
      break;
      }
      case 9:
      case 10:
      {
        test[0]= phi.CheckNeumannFreeWall(d, 1,   l, Wall::Left)
                +phi.CheckNeumannFreeWall(0, l+1, r-1, Wall::Left);

        test[1]= phi.CheckNeumannFreeWall(LX-1-d,   1, l, Wall::Right)
                +phi.CheckNeumannFreeWall(LX-1,   l+1, r-1, Wall::Right);

        test[2]= phi.CheckNeumannFreeWall(d+1,    0, LX-2-2*d, Wall::Front)
                +phi.CheckNeumannFreeWall(1,      l, d,      Wall::Front)
                +phi.CheckNeumannFreeWall(LX-d-1, l, d,      Wall::Front);

        test[3]= phi.CheckNeumannFreeWall(1,      LY-1, LX-2, Wall::Back);


        std::cout << "\nCheck of bc for field (incl corners):\nLeft: "
                  << test[0] << ", Right: " << test[1] << ", Front: "
                  << test[2] << ", Back: " << test[3];
      break;
      }
      case 11:
      case 12:
      {
        test[0]=QQxx.CheckNeumannFreeWall(40, 10, 21, Wall::Left);
        test[1]=QQxx.CheckNeumannFreeWall(10, 10, 21, Wall::Right);
        test[2]=QQxx.CheckNeumannFreeWall(10, 30, 31, Wall::Front);
        test[3]=QQxx.CheckNeumannFreeWall(10, 10, 31, Wall::Back);


        test[0]= phi.CheckNeumannFreeWall(d, o,    l,   Wall::Left);
        test[1]= phi.CheckNeumannFreeWall(LX-d-1, o,    l,   Wall::Right);

        test[2]= phi.CheckNeumannFreeWall(0,      0,    LX, Wall::Front)
                +phi.CheckNeumannFreeWall(0,      o2-1, d+1,    Wall::Front)
                +phi.CheckNeumannFreeWall(LX-d-1, o2-1, d+1,    Wall::Front);

        test[3]= phi.CheckNeumannFreeWall(0,      LY-1, LX, Wall::Back)
                +phi.CheckNeumannFreeWall(0,      o,    d+1,    Wall::Back)
                +phi.CheckNeumannFreeWall(LX-d-1, o,    d+1,    Wall::Back);


        std::cout << "\nCheck of bc for field (incl corners):\nLeft: "
                  << test[0] << ", Right: " << test[1] << ", Front: "
                  << test[2] << ", Back: " << test[3];
      break;
      }
    }
  }
#endif //DEBUG

  // check that the sum of f is constant
  {
    double fcheck = 0;
    for(unsigned k=0; k<DomainSize; ++k)
      if(!outside[k]) fcheck = accumulate(begin(ff[k]), end(ff[k]), fcheck);
    std::cout << "\nftot=" << ftot << ", difference: " << ftot-fcheck;
    if(abs(ftot-fcheck)>1 and (BC!=19 and BC!=20 and BC!=21 and BC!=22 and BC!=23 and BC!=24 and BC!=25 and BC!=26 and BC!=27 and BC!=28 and BC!=29 and BC!=30 and BC!=31 and BC!=32 and BC!=33 and BC!=34 and BC!=35 and BC!=36))
      throw error_msg("f is not conserved (", ftot, "/", fcheck, ")");
  }

  // check that phi is conserved
  {
    double pcheck = 0;
    for(unsigned k=0; k<DomainSize; ++k)
      if(!outside[k])pcheck += phi[k];
    std::cout << "\nptot=" << ptot << ", difference: " << ptot-pcheck;
#ifdef DEBUG
     std::cout << " ,integral over flux: "<<countphi_flux << ", integral over laplacian: "<<countphi_laplace;
#endif //DEBUG
    //if(abs(ptot-pcheck)>1)
    //  throw error_msg("phi is not conserved (", ptot, "/", pcheck, ")");
  }

  cout << "\n";
}
