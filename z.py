'''
  Finds location of maximum probability in a 2D space defined by three seperate PDFs.
  
  | z.py: teaser at https://tech.zalando.com/jobs/65946/?gh_jid=65946
  
  | author: Dominik Leier
  
  | Requirements: python libs: numpy, scipy, matplotlib, mpl_toolkits.basemap
  
  | Usage: python z.py spree.dat 100
  |   For the program to run an ascii file containing coordinates (column1=latitude, column2=longitude)
  |   of the course of the river Spree must be provided!
  |   The number is the number of grid cells in one dimension, i.e. 100 creates a 100x100 grid
  
  | Output: a png file (zmap.png) and on-screen log
  
  | Instructions for the quest are as follows:

  | The candidate is likely to be close to the river Spree.
  | The probability at any point is given by a Gaussian function of its shortest distance to the river.
  | The function peaks at zero and has 95 per cent of its total integral within +/-2730m.
  
  | A probability distribution centered around the Brandenburg Gate also informs us
  | of the candidate's location. The distribution's radial profile is log-normal with
  | a mean of 4700m and a mode of 3877m in every direction.
  | A satellite offers further information:
  | with 95 per cent probability she is located within 2400 m distance of the satellite's path
  | (assuming a normal probability distribution)
  
'''
__version__ = "0.1.0"


import numpy as np
import matplotlib.pyplot as plt
import sys
from pylab import *
from scipy.special import erf
from scipy.special import erfinv
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LogNorm


# global variables
# points on great circle
la_sat=np.array([52.590117,52.437385])
lo_sat=np.array([13.39915,13.553989])
# coords brandenburg gate
la_bgate=52.516288
lo_bgate=13.377689
# spree dat from file in top-level script at EOF



def main():
  """Finds location of maximum probability in a 2D space defined by three seperate PDFs."""
  
  # plotting
  fig=plt.figure(11,figsize=(6,6))
  plt.clf()
  plt.subplots_adjust(left=0.01, bottom=0.01, right=0.99, top=0.99, wspace=0., hspace=0.)
  ax = plt.subplot(1,1,1)
  xticklabels = getp(gca(), 'xticklabels')
  yticklabels = getp(gca(), 'yticklabels')
  
  # create map and set size
  boxsize=6000
  lat, lon=52.52437, 13.41053 # centre berlin map
  m = Basemap(projection='ortho',lon_0=lon,lat_0=lat,resolution='l',area_thresh='l',llcrnrx=-boxsize,llcrnry=-boxsize,urcrnrx=+boxsize,urcrnry=+boxsize)
  
  # read and plot street map / water areas from shapefile
  m.readshapefile('berlin_germany/berlin_germany_osm_roads_gen1', 'Streets',color='#DDDDDD')
  m.readshapefile('berlin_germany/berlin_germany_osm_waterareas', 'Water',color='#add8e6')
  
  # compute and plot points on great circle
  x,y=m.gcpoints(lo_sat[0],la_sat[0], lo_sat[1],la_sat[1], 100)
  m.plot(x,y,'r-',zorder=7)
  
  # plot spree list of coords
  x,y=m(lo_spree,la_spree)
  m.plot(x,y,'b.',markersize=2)
  
  # plot brandenburg gate
  gx,gy=m(lo_bgate,la_bgate) # coord transform
  m.plot(gx,gy,'bx')
  
  # create search grid
  bestY,bestX=m(13.4497839044, 52.5142566528 ) # grid centre
  minX=bestX-8000
  maxX=bestX+8000
  minY=bestY-8000
  maxY=bestY+8000
  
  rangeX=np.linspace(minX,maxX,nbins)
  rangeY=np.linspace(minY,maxY,nbins)
  npoints=None
  
  meshX, meshY = meshgrid(rangeX, rangeY)
  mzz=np.zeros(shape=(len(rangeY),len(rangeX)))
  maxmzz=0
  
  # compute probabilities for grid and location with max. prob.
  for ii in range(0, len(rangeY)):
    s = str(int(float(ii+1)/nbins*100.)) + '% computed'
    print s,
    sys.stdout.flush()
    backspace(len(s))
    for jj in range(0, len(rangeX)):
      mY,mX=m(rangeY[ii],rangeX[jj],inverse=True)
      prob_bgate,prob_spree,prob_sat,x_bgate,x_spree,x_sat=z_prob(m,mX,mY,npoints)
      mzz[ii][jj]=prob_bgate*prob_spree*prob_sat
      if mzz[ii][jj]>maxmzz:
        maxmzz=mzz[ii][jj]
        maxmX=mX
        maxmY=mY
        best_x_bgate=x_bgate
        best_x_spree=x_spree
        best_x_sat=x_sat
  
  # plot location with max. prob.
  loc_x,loc_y= m(maxmY,maxmX)
  m.plot(loc_x,loc_y,'r.',zorder=8)
  
  # output
  print "\nDone!"
  print "+++ search log +++"
  print "max prob.:", maxmzz, " found at: (lat,lon)=(",maxmX,",",maxmY,")"
  print "with dist. to spree: ", best_x_spree
  print "with dist. to bgate: ", best_x_bgate
  print "with dist. to sat: ", best_x_sat
  print "resolution in m: (dx,dy)=",(maxX-minX)/nbins,(maxY-minY)/nbins #
  
  # plotting: color map and contours
  logcrange=arange(-5,-2.1,0.1)
  crange=10**logcrange
  cs = m.contourf(meshY,meshX,mzz,levels=crange,cmap=cm.Blues,zorder=5, alpha=0.3)
  
  # plotting legend
  prop = matplotlib.font_manager.FontProperties(size=10)
  b, = plot([10000,20000],'bx', linewidth=1)
  c, = plot([10000,20000],'b.', linewidth=1)
  d, = plot([10000,20000],'r-', linewidth=1)
  e, = plot([10000,20000],'r.', linewidth=1)
  legend([b,c,d,e,], ['Brandenburg Gate','Spree (file)','projected satellite path','coordinates with max. probability'],loc='upper left', ncol=1, shadow=False, fancybox=True, numpoints=1, prop=prop,labelspacing=-0.0)
  
  # plotting scale
  shift=6700
  sx=6368777.0109+shift
  sy=6370097.04832+shift
  m.plot([sx,sx+1000],[sy,sy],'k-',linewidth=5)
  s1lo,s1la=m(sx+1000,sy,inverse=True)
  s2lo,s2la=m(sx,sy,inverse=True)
  text(sx,sy-400,r' $1$ $km$')
  
  #plotting colorbar
  cax = axes([0.05,0.05,0.9,0.04]) # cbar location and dimension
  cbar = plt.colorbar(cs,boundaries=logcrange,cax=cax,orientation="horizontal",norm=LogNorm(vmin=1E-6, vmax=7E-3),ticks=[5e-5,1e-4,5e-4,1e-3,5e-3])
  cbar.ax.tick_params(labelsize=9)
  cbar.set_label(r'$\Pi_i P_i$                                                                                           ',labelpad=-25)

  savefig('zmap.png')
  
  return


def log_norm(x,mu,sigma):
  return exp(-(log(x)-mu)**2./(2*sigma**2))/(x*sigma*sqrt(2*pi))

def gauss(x,mu,sigma):
  return exp(-0.5*(x-mu)**2./sigma**2.)/sigma/sqrt(2*pi)

def shortest_distance_spree(m,la, lo):
  """Shortest distance between point and list of points
    m - basemap;
    calculates the shortest distance between coords la,lo
    and the connecting lines of two consecutive coordinates in the list
    of coordinates la_spree,lo_spree
  """
  min_dist=1E12
  for ii in range(0,len(la_spree)-1):
    dist_tmp= distance_line_point_la_lo(m,la,lo,la_spree[ii],lo_spree[ii],la_spree[ii+1],lo_spree[ii+1])
    segment,xl,yl=projection_in_segment(m,la,lo,la_spree[ii],lo_spree[ii],la_spree[ii+1],lo_spree[ii+1])
    if (segment== True and dist_tmp<min_dist):
      min_dist=dist_tmp
  dist_tmp=dist_closest_support_point(m,la,lo,la_spree,lo_spree)
  if dist_tmp<min_dist:
    min_dist=dist_tmp
  return min_dist

def dist_closest_support_point(m,la,lo,la_list,lo_list):
  """used within shortest_distance_spree
    m - basemap;
    calculates the shortest distance between point la,lo
    and points in the list of coordinates la_spree,lo_spree
  """
  min_dist=1E12
  for i in range(0,len(lo_spree)):
    dist_tmp= distance_on_earth(la_list[i], lo_list[i], la, lo)
    if dist_tmp<min_dist:
      min_dist=dist_tmp
  return min_dist

def projection_in_segment(m,la0,lo0,la1,lo1,la2,lo2):
  """used within shortest_distance_spree
    m - basemap;
    returns boolean value and coordinates of orthogonal projection point of la0,lo0
    onto the line going through la1,lo1 and la2,lo2
    boolean value True: orthogonal projection lies between the two points 
    boolean value False: orthogonal projection does not lie between the two points
  """
  x0,y0=m(lo0,la0)
  x1,y1=m(lo1,la1)
  x2,y2=m(lo2,la2)
  r=((x1-x2)*(x0-x2)+(y1-y2)*(y0-y2))/(((y2-y1)**2.+(x2-x1)**2.))
  xl=r*(x1-x2)+x2
  yl=r*(y1-y2)+y2
  dist_1=sqrt((x1-xl)**2.+(y1-yl)**2.)
  dist_2=sqrt((x2-xl)**2.+(y2-yl)**2.)
  dist_12=sqrt((x2-x1)**2.+(y2-y1)**2.)
  if (dist_1 <= dist_12+0.01 and dist_2 <= dist_12+0.01):
    return True, xl,yl
  else:
    return False, xl,yl

def shortest_distance_sat(m,la, lo, npoints=None):
  """Shortest distance between point and great circle
    m - basemap;
    npoints is optional:
    if npoints is not specified shortest_distance_sat1 is called.
    if npoints is specified shortest_distance_sat2 is called.
  """
  if npoints==None:
    min_dist=shortest_distance_sat1(m,la, lo)
  if npoints!=None:
    min_dist=shortest_distance_sat2(m,la, lo,npoints)
  return min_dist

def shortest_distance_sat1(m,la, lo):
  """Shortest distance between point and line
    m - basemap;
    calculates the shortest distance between coords la,lo
    and a line through two points defining the satellite path
  """
  min_dist=distance_line_point_la_lo(m,la,lo,la_sat[0],lo_sat[0],la_sat[1],lo_sat[1])
  return min_dist

def shortest_distance_sat2(m,la, lo, npoints):
  """Substitute function for shortest_distance_sat1
    m - basemap;
    calculates the shortest distance between coords la,lo and a great circle
    with number of equidistant points npoints.
  """
  x,y=m.gcpoints(lo_sat[0],la_sat[0], lo_sat[1],la_sat[1], npoints)
  lo_list,la_list=m(x,y,inverse=True)
  min_dist=1E12
  for i in range(0,len(la_list)):
    dist=distance_on_earth(la_list[i], lo_list[i], la, lo)
    if dist<min_dist:
      min_dist=dist
  return min_dist


def z_prob(m,la,lo,npoints=None):
  """Computes distances to point (la,lo) and probabilities according to the three PDFs
    m - basemap;
    npoints - optional, if
  """
  #part1 distance from river and prob
  sigma=2.730/sqrt(2)/erfinv(2*0.95-1.)
  x_spree=shortest_distance_spree(m,la, lo)
  prob_spree=gauss(x_spree,0,sigma)
  #part2 distance from brandenburg gate and prob
  mean=4.700
  mode=3.877
  x_bgate=distance_on_earth(la_bgate, lo_bgate, la, lo)
  mu=(2.*log(mean)+log(mode))/3.
  sigma=sqrt((log(mean)-log(mode))*2./3.)
  if x_bgate!=0:
    prob_bgate=log_norm(x_bgate,mu,sigma)
  else:
    prob_bgate=0
  #part3 distance from great circle and prob
  sigma=2.400/sqrt(2)/erfinv(2*0.95-1.)
  x_sat=shortest_distance_sat(m,la,lo,npoints)
  prob_sat=gauss(x_sat,0,sigma)
  return prob_bgate,prob_spree,prob_sat,x_bgate,x_spree,x_sat


def distance_line_point_la_lo(m,la0,lo0,la1,lo1,la2,lo2):
  """Computes distance between point (la0,lo0)
    and line through (la1,lo1) and (la2,lo2) given in longitude and latitude in km
    m - basemap;
  """
  x0,y0=m(lo0,la0)
  x1,y1=m(lo1,la1)
  x2,y2=m(lo2,la2)
  return distance_line_point(x0,y0,x1,y1,x2,y2)*0.001


def distance_line_point(x0,y0,x1,y1,x2,y2):
  """Computes distance between point (x0,y0) and line through (x1,y1) and (x2,y2)"""
  return (abs((y2-y1)*x0-(x2-x1)*y0+x2*y1-y2*x1))/(sqrt((y2-y1)**2.+(x2-x1)**2.))

def distance_on_earth(lat1, long1, lat2, long2):
  """ distance between two points given in long and lat on a sphere with radius 6371"""
  radius=6371.
  degrees_to_radians = math.pi/180.0
  
  # phi = 90 - latitude
  phi1 = (90.0 - lat1)*degrees_to_radians
  phi2 = (90.0 - lat2)*degrees_to_radians
  
  # theta = longitude
  theta1 = long1*degrees_to_radians
  theta2 = long2*degrees_to_radians
  
  cos = (math.sin(phi1)*math.sin(phi2)*math.cos(theta1 - theta2) + math.cos(phi1)*math.cos(phi2))
  arc = math.acos( cos )
  
  return arc*radius



############## for calculating uncertainties ###########
def bgate_error():
  return distance_on_earth(la_bgate,lo_bgate,52.516272, 13.377722) #latter coords taken from GeoHack

def sat_error():
  boxsize=6000
  lat,lon=52.52437,13.41053 # centre Berlin map
  m = Basemap(projection='ortho',lon_0=lon,lat_0=lat,resolution='l',area_thresh='l',llcrnrx=-boxsize,llcrnry=-boxsize,urcrnrx=+boxsize,urcrnry=+boxsize)
  
  x0,y0=m(lo_sat[0],la_sat[0])
  x1,y1=m(lo_sat[1],la_sat[1])
  dog=sqrt((x1-x0)**2.+(y1-y0)**2.)
  print x0+0.5*(x1-x0),y0+0.5*(y1-y0)
  print distance_on_earth(la_sat[0],lo_sat[0],la_sat[1],lo_sat[1])
  
  x,y=m.gcpoints(lo_sat[0],la_sat[0], lo_sat[1],la_sat[1], 10000)
  
  log,lag=m(x[5000],y[5000],inverse=True)
  lop,lap=m(x0+0.5*(x1-x0),y0+0.5*(y1-y0),inverse=True)
  print lag,log,lap,lop
  return distance_on_earth(lag,log,lap,lop)

def backspace(n):
  print '\r' * n,
  return





if __name__ == "__main__":

  if len(sys.argv) < 3 or len(sys.argv) > 3:
    print "type: python z.py <file> <resolution>"
    print "e.g. python z.py spree.dat 100"
    print "continue with default: spree.dat 100"
    z_data=np.loadtxt('spree.dat',delimiter=",")
    nbins=100
  
  else:
    file=sys.argv[1]
    nbins=int(sys.argv[2])
    z_data=np.loadtxt(file,delimiter=",")

  la_spree=z_data[:,0]
  lo_spree=z_data[:,1]

  main()

else:

  file="spree.dat"
  nbins=100
  z_data=np.loadtxt(file,delimiter=",")
  la_spree=z_data[:,0]
  lo_spree=z_data[:,1]

