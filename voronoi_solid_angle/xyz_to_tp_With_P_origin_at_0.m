function [ t, p ] = xyz_to_tp_With_P_origin_at_0 ( x, y, z )

%*****************************************************************************80
%
%% XYZ_TO_TP converts (X,Y,Z) to (Theta,Phi) coordinates on the unit sphere.
%
%  Modified:
%
%    09 September 2010
%
%  Author:
%
%    Dmitri Laikov
%
%  Reference:
%
%    Vyacheslav Lebedev, Dmitri Laikov,
%    A quadrature formula for the sphere of the 131st
%    algebraic order of accuracy,
%    Russian Academy of Sciences Doklady Mathematics,
%    Volume 59, Number 3, 1999, pages 477-481.
%
%  Parameters:
%
%    Input, real X, Y, Z, the Cartesian coordinates of a point
%    on the unit sphere.
%
%    Output, real T, P, the Theta and Phi coordinates of
%    the point.
%
  p = acos ( z );

  fact = sqrt ( x * x + y * y );

  if ( 0.0 < fact )
    t = acos ( x / fact );
  else
    t = acos ( x );
  end

  if ( y < 0.0 )
    t = - t;
  end

%  Convert to degrees. Azimuth going from 0 to 360. elevation from -90 to
%  90. ////TM adjustment 
  t = t * 180.0 / pi;
  t = t - 90;
  if t < 0
      t = t + 360;
  end
  
  p = (p * 180.0 / pi) - 90; 

  return
end
