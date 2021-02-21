function cvt_2d_sampling ( g_num, it_num, s_num )

%*****************************************************************************80
%
%% CVT_2D_SAMPLING carries out the Lloyd algorithm in a 2D unit box.
%
%  Discussion:
%
%    This program is a variation of the CVT_2D_LLOYD method.
%
%    Instead of using an exact technique to determine the Voronoi
%    regions, it uses sampling.
%
%    MATLAB highly inconvenienced me (let me put it nicely) by removing
%    the DSEARCH function, making all my CVT codes fail until I found
%    the replacement function.  It would have killed them to leave the
%    old thing available, so my codes would still run?
%    JVB, 06 June 2016.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    06 June 2016
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer G_NUM, the number of generators.
%    A value of 50 is reasonable.
%
%    Input, integer IT_NUM, the number of CVT iterations.
%    A value of 20 or 50 might be reasonable.
%
%    Input, integer S_NUM, the number of sample points to use
%    when estimating the Voronoi regions.
%    A value of 1,000 is too low.  A value of 1,000,000 is somewhat high.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'CVT_2D_SAMPLING\n' );
  fprintf ( 1, '  MATLAB version\n' );
  fprintf ( 1, '  Use sampling to approximate Lloyd''s algorithm\n' );
  fprintf ( 1, '  in the 2D unit square.\n' );

  if ( nargin < 1 )
    g_num = input ( '  Enter number of generators: ' );
  elseif ( ischar ( g_num ) )
    g_num = str2num ( g_num );
  end

  if ( nargin < 2 ) 
    it_num = input ( '  Enter number of iterations: ' );
  elseif ( ischar ( it_num ) )
    it_num = str2num ( it_num );
  end

  if ( nargin < 3 ) 
    s_num = input ( '  Enter number of sample points: ' );
  elseif ( ischar ( s_num ) )
    s_num = str2num ( s_num );
  end

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Number of generators is %d\n', g_num );
  fprintf ( 1, '  Number of iterations is %d\n', it_num );
  fprintf ( 1, '  Number of samples is %d\n', s_num );
%
%  Initialize the generators.
%
  g = rand ( g_num, 2 );
%
%  Carry out the iteration.
%
  step = 1 : it_num;
  e = nan ( it_num, 1 );
  gm = nan ( it_num, 1 );

  for it = 1 : it_num
%
%  Compute the Delaunay triangle information T for the current nodes.
%
    t = delaunay ( g(:,1), g(:,2) );
%
%  Display the Voronoi cells.
%
    figure ( 1 );
    voronoi ( g(:,1), g(:,2), t );
    title_string = sprintf ( 'Voronoi, step %d', it );
    title ( title_string );
    axis equal
    axis ( [  0.0, 1.0, 0.0, 1.0 ] )
    axis square
    drawnow
%
%  Display the Delaunay triangulation.
%
    figure ( 2 );
    trimesh ( t, g(:,1), g(:,2), zeros(g_num,1), 'EdgeColor', 'r' )
    hold on
    scatter ( g(:,1), g(:,2), 'k.' )
    title_string = sprintf ( 'Delaunay, step %d', it );
    title ( title_string );
    axis ( [  0.0, 1.0, 0.0, 1.0 ] )
    axis square
    view ( 2 )
    hold off
%
%  Generate sample points.
%  New option for fixed grid sampling.
%
    if ( false )
      s = rand ( s_num, 2 );
    else
      s2 = floor ( sqrt ( s_num ) );
      s2 = max ( s2, 2 );
      s_num = s2 * s2;
      [ sx, sy ] = meshgrid ( linspace ( 0.0, 1.0, s2 ) );
      sx = reshape ( sx, s2 * s2, 1 );
      sy = reshape ( sy, s2 * s2, 1 );
      s = [ sx, sy ];
    end
%
%  For each sample point, find K, the index of the nearest generator.
%  We do this efficiently by using the Delaunay information with
%  Matlab's DSEARCH command, rather than a brute force nearest neighbor
%  computation.  Also, DSEARCH has been removed, so we need DSEARCHN.
%  
    k = dsearchn ( g, t, s );

    m(:,1) = accumarray ( k, ones(s_num,1) );

    g_new(:,1) = accumarray ( k, s(:,1) ) ./ m(:,1);
    g_new(:,2) = accumarray ( k, s(:,2) ) ./ m(:,1);
%
%  Compute the energy.
%
    e(it,1) = sum ( ( s(:,1) - g(k(:,1),1) ).^2 ...
                  + ( s(:,2) - g(k(:,1),2) ).^2 ) / s_num;
%
%  Display the energy.
%
    figure ( 3 )
    plot ( step, log ( e ), 'm-*' )
    title ( 'Log (Energy)' )
    xlabel ( 'Step' )
    ylabel ( 'Energy' )
    grid
%
%  Compute the generator motion.
%
    gm(it,1) = sum ( ( g_new(:,1) - g(:,1) ).^2 ...
                   + ( g_new(:,2) - g(:,2) ).^2 ) / g_num;
%
%  Display the generator motion.
%
    figure ( 4 )
    plot ( step, log ( gm ), 'm-*' )
    title ( 'Log (Average generator motion)' )
    xlabel ( 'Step' )
    ylabel ( 'Energy' )
    grid
%
%  Update the generators.
%
    g = g_new;

  end
%
%  Save final graphics.
%
  figure ( 1 )
  filename = 'voronoi.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );

  figure ( 2 )
  filename = 'delaunay.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );

  figure ( 3 )
  filename = 'energy.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );

  figure ( 4 )
  filename = 'motion.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );
%
%  Terminate.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'CVT_2D_SAMPLING\n' );
  fprintf ( 1, '  Normal end of execution.\n' );

  return
end
