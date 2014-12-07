# Program of ECE 569 Project B
# Authors :- Debasmit Das, Tongyang Liu
# Input :- DH Parameters and Locus Parameters
# Output :-Drawing a Circle and Lines

from openravepy import * # This is imported for OpenRave Methods/Functions
import time # This is imported for invoking delay in the Program
import numpy as np # This is Matrix and Array Operations
from numpy import *
from sympy import * # This is imported for symbolic variables
from sympy import Matrix 
from numpy import linalg # this is useful for Linear ALgebra Applications

p=raw_input("Enter Path for DH File\n")
DH=numpy.genfromtxt(p) # Matrix of D-H Parameters and ranges read from a file
th=DH[:,0] # Separating the Joint Angles
thr=th*(np.pi/180) # Converting to Radians 
al=DH[:,1] # Separating the Link Twist
alr=al*(np.pi/180) # Converting to Radians 
a=DH[:,2] # Separating the Link Length
d=DH[:,3] # Separating the Joint Offset
r1=DH[:,4]*(np.pi/180) # Converting lower limit of range to radians
r2=DH[:,5]*(np.pi/180) # Converting upper limit of range to radians

# Here each joint angle which are the variables are stored 
q=[]
for i in range(1,7):
  q.append(Symbol('q'+repr(i))) # A Matrix of Symbolic Joint Variables is created

A=[[0.0],[0.0],[0.0],[0.0],[0.0],[0.0]];#Create a list of zero matrix This list will later fill up with symbolic matrices

# The Link transformation matrix is computed

for i in range(0,6):
  A[i]=Matrix([[cos(q[i]),-cos(alr[i])*sin(q[i]),sin(alr[i])*sin(q[i]),a[i]*cos(q[i])],[sin(q[i]),cos(alr[i])*cos(q[i]),-sin(alr[i])*cos(q[i]),a[i]*sin(q[i])],[0,sin(alr[i]),cos(alr[i]),d[i]],[0,0,0,1]])
  print "\n" 

# Now the Geometric Inverse Kinmeatics Subroutine is written that is to be used for position control

def GeomIK(T,ARM,ELB,WRI,OHM):
  px=T[0,3];py=T[1,3];pz=T[2,3];nx=T[0,0]; ny=T[1,0] ; nz=T[2,0] ; sx=T[0,1] ; sy=T[1,1] ; sz=T[2,1] ; ax=T[0,2] ; ay=T[1,2] ; az=T[2,2];
  Px=px - d[5]*ax ; Py=py - d[5]*ay ; Pz=pz - d[5]*az ;
  rg = (Px**2 + Py**2 - d[1]**2)**0.5 ; R1 = (Px**2 + Py**2)**0.5
  Qg=np.array([0.0,0.0,0.0,0.0,0.0,0.0])
  Qg[0]=atan2(-ARM*Py*rg - Px*d[1] , -ARM*Px*rg + Py*d[1]) 
  R2=(Px**2 + Py**2 + Pz**2 - d[1]**2)**0.5
  sina = -Pz/R2 ; cosa = -ARM*rg/R2 ; cosb=(a[1]**2 + R2**2 - d[3]**2 - a[2]**2)/(2*a[1]*R2) ; sinb = (1-cosb**2)**0.5
  Qg[1]=atan2(sina*cosb + ARM*ELB*cosa*sinb , cosa*cosb - ARM*ELB*sina*sinb) 
  cosphi = (a[1]**2 + d[3]**2 + a[2]**2 - R2**2)/(2*a[1]*(d[3]**2 + a[2]**2)**0.5)
  sinphi = ARM*ELB*(1 - cosphi**2)**0.5
  sinbeta = d[3]/((d[3]**2 + a[2]**2)**0.5)
  cosbeta = abs(a[2])/((d[3]**2 + a[2]**2))**0.5
  Qg[2] = atan2(sinphi*cosbeta - cosphi*sinbeta, cosphi*cosbeta + sinphi*sinbeta)
  M=WRI*OHM;
  Qg[3] = atan2(M*(cos(Qg[0])*ay - sin(Qg[0])*ax) , M*(cos(Qg[0])*cos(Qg[1]+Qg[2])*ax + sin(Qg[0])*cos(Qg[1]+Qg[2])*ay - sin(Qg[1]+Qg[2])*az))
  Qg[4] = atan2((cos(Qg[0])*cos(Qg[1]+Qg[2])*cos(Qg[3]) - sin(Qg[0])*sin(Qg[3]))*ax + (sin(Qg[0])*cos(Qg[1]+Qg[2])*cos(Qg[3]) + cos(Qg[0])*sin(Qg[3]))*ay - cos(Qg      [3])*sin(Qg[1]+Qg[2])*az , cos(Qg[0])*sin(Qg[1]+Qg[2])*ax + sin(Qg[0])*sin(Qg[1]+Qg[2])*ay + cos(Qg[1]+Qg[2])*az)
  Qg[5] = atan2((-sin(Qg[0])*cos(Qg[3]) - cos(Qg[0])*cos(Qg[1]+Qg[2])*sin(Qg[3]))*nx + (cos(Qg[0])*cos(Qg[3]) - sin(Qg[0])*cos(Qg[1]+Qg[2])*sin(Qg[3]))*ny + sin(Qg[1]+Qg[2])*sin(Qg[3])*nz, (-sin(Qg[0])*cos(Qg[3]) - cos(Qg[0])*cos(Qg[1]+Qg[2])*sin(Qg[3]))*sx + (cos(Qg[0])*cos(Qg[3]) - sin(Qg[0])*cos(Qg[1]+Qg[2])*sin(Qg[3]))*sy + sin(Qg[1]+Qg[2])*sin(Qg[3])*sz)
  if(Qg[1]<=np.pi and Qg[1]>=np.pi/2):
    Qg[1]=Qg[1]-2*np.pi
  if(Qg[2]<=-np.pi/2 and Qg[2]>=-np.pi):
    Qg[2]=Qg[2]+2*np.pi
  return Qg
  
#Now Lets try finding the Jacobian
# First Let us try the cross-product method

P0=(A[0]*A[1]*A[2]*A[3]*A[4]*A[5])[0:3,3]
Z0=Matrix([[0],[0],[1]])
J0=Z0.cross(P0).transpose().col_join(Z0)

P1=(A[0]*A[1]*A[2]*A[3]*A[4]*A[5]-A[0])[0:3,3]
Z1=A[0][0:3,2]
J1=Z1.cross(P1).transpose().col_join(Z1)

P2=(A[0]*A[1]*A[2]*A[3]*A[4]*A[5]-A[0]*A[1])[0:3,3]
Z2=(A[0]*A[1])[0:3,2]
J2=Z2.cross(P2).transpose().col_join(Z2)

P3=(A[0]*A[1]*A[2]*A[3]*A[4]*A[5]-A[0]*A[1]*A[2])[0:3,3]
Z3=(A[0]*A[1]*A[2])[0:3,2]
J3=Z3.cross(P3).transpose().col_join(Z3)

P4=(A[0]*A[1]*A[2]*A[3]*A[4]*A[5]-A[0]*A[1]*A[2]*A[3])[0:3,3]
Z4=(A[0]*A[1]*A[2]*A[3])[0:3,2]
J4=Z4.cross(P4).transpose().col_join(Z4)

P5=(A[0]*A[1]*A[2]*A[3]*A[4]*A[5]-A[0]*A[1]*A[2]*A[3]*A[4])[0:3,3]
Z5=(A[0]*A[1]*A[2]*A[3]*A[4])[0:3,2]
J5=Z5.cross(P5).transpose().col_join(Z5)

Jv=J0.row_join(J1).row_join(J2).row_join(J3).row_join(J4).row_join(J5)

#Now, Let us try the Differential Translation and Rotation Approach
U0=A[0]*A[1]*A[2]*A[3]*A[4]*A[5]
J0=Matrix([[U0[0,3]*U0[1,0]-U0[1,3]*U0[0,0]],[U0[0,3]*U0[1,1]-U0[1,3]*U0[0,1]],[U0[0,3]*U0[1,2]-U0[1,3]*U0[0,2]],[U0[2,0]],[U0[2,1]],[U0[2,2]]])

U1=A[1]*A[2]*A[3]*A[4]*A[5]
J1=Matrix([[U1[0,3]*U1[1,0]-U1[1,3]*U1[0,0]],[U1[0,3]*U1[1,1]-U0[1,3]*U1[0,1]],[U1[0,3]*U1[1,2]-U1[1,3]*U1[0,2]],[U1[2,0]],[U1[2,1]],[U1[2,2]]])

U2=A[2]*A[3]*A[4]*A[5]
J2=Matrix([[U2[0,3]*U2[1,0]-U2[1,3]*U2[0,0]],[U2[0,3]*U2[1,1]-U2[1,3]*U2[0,1]],[U2[0,3]*U2[1,2]-U2[1,3]*U2[0,2]],[U2[2,0]],[U2[2,1]],[U2[2,2]]])

U3=A[3]*A[4]*A[5]
J3=Matrix([[U3[0,3]*U3[1,0]-U3[1,3]*U3[0,0]],[U3[0,3]*U3[1,1]-U3[1,3]*U3[0,1]],[U3[0,3]*U3[1,2]-U3[1,3]*U3[0,2]],[U3[2,0]],[U3[2,1]],[U3[2,2]]])

U4=A[4]*A[5]
J4=Matrix([[U4[0,3]*U4[1,0]-U4[1,3]*U4[0,0]],[U4[0,3]*U4[1,1]-U4[1,3]*U4[0,1]],[U4[0,3]*U4[1,2]-U4[1,3]*U4[0,2]],[U4[2,0]],[U4[2,1]],[U4[2,2]]])

U5=A[5]
J5=Matrix([[U5[0,3]*U5[1,0]-U5[1,3]*U5[0,0]],[U5[0,3]*U5[1,1]-U5[1,3]*U5[0,1]],[U5[0,3]*U5[1,2]-U5[1,3]*U5[0,2]],[U5[2,0]],[U5[2,1]],[U5[2,2]]])

Jd=J0.row_join(J1).row_join(J2).row_join(J3).row_join(J4).row_join(J5)

# The vector cross product sub-routine

def JacoV(QJ):
 return np.array(Jv.subs(q[0],QJ[0]).subs(q[1],QJ[1]).subs(q[2],QJ[2]).subs(q[3],QJ[3]).subs(q[4],QJ[4]).subs(q[5],QJ[5]))

# Defining Camera Transformation
Tstart=np.array([[ 0.9983323 ,  0.03651881,  0.04471023, -0.12004546],[-0.02196539, -0.47593777,  0.87920462, -4.18794775],[ 0.0533868 , -0.87872044, -0.47434189,  2.22462606],[ 0.        ,  0.        ,  0.        ,  1.        ]])

Tcirc=np.array([[ 0.20604119,  0.33831506, -0.9181993 ,  3.84545088],[ 0.97824604, -0.09434114,  0.18475504, -0.79177779],[-0.02411856, -0.93629198, -0.35039353,  1.40385485],[ 0.        ,  0.        ,  0.        ,  1.        ]])

Tline=np.array([[ 0.02187936,  0.99890052,  0.04146134,  0.02088733],[ 0.98965908, -0.01575924, -0.14257124,  0.30318534],[-0.14176108,  0.04415196, -0.98891577,  2.93201566],[ 0.        ,  0.        ,  0.        ,  1.        ]])

p=raw_input('Enter Path for XML File\n') # Change the path here
env=Environment();
time.sleep(1);
env.SetViewer('qtcoin');
time.sleep(1)
env.Load(p); # Load the simple robot scene of PUMA 560
time.sleep(1)
env.GetViewer().SetCamera(Tstart)
robot = env.GetRobots()[0] # Get the first robot
Tz=np.array([[1,0,0,0],[0,1,0,0],[0,0,1,-1.37],[0,0,0,1]])
robot.SetTransform(Tz) 

time.sleep(2)

print 'METHOD : INVERSE KINEMATICS/ POSITION CONTROL'
print ' '

#Drawing circle using Inverse Kinematics Method :
env.GetViewer().SetCamera(Tcirc)
xcir=float(raw_input("Enter distance along x axis where you want to draw circle in m(0.508 is a good choice)\n"))
rcir=float(raw_input("Enter Radius of Circle in m (0.254 is a good choice)\n"))
v=float(raw_input("Enter velocity of Motion in m/s(0.1 is a good choice)\n"))
env.GetViewer().SetTitle('Drawing Circle Using Inverse Kinematics/Position Control Method')
time.sleep(2)



sam=0.005
handles=[]
ycir=0.0;
while (ycir<=rcir): # First Quadrant
  T=np.array([[1,0,0,xcir],[0,1,0,ycir],[0,0,1,(rcir**2 - ycir**2)**0.5],[0,0,0,1]])
  Q=GeomIK(T,-1,1,1,1)
  robot.SetDOFValues(Q)
  handles.append(env.plot3(robot.GetLinks()[6].GetTransform()[0:3,3],5))
  ycir=ycir+sam;
  time.sleep(sam/v);

ycir=rcir;
while (ycir>=0): # Second Quadrant
  T=np.array([[1,0,0,xcir],[0,1,0,ycir],[0,0,1,-(rcir**2 - ycir**2)**0.5],[0,0,0,1]])
  Q=GeomIK(T,-1,1,1,1)
  robot.SetDOFValues(Q)
  handles.append(env.plot3(robot.GetLinks()[6].GetTransform()[0:3,3],5))
  ycir=ycir-sam;
  time.sleep(sam/v);

ycir=0.0;
while (ycir>=-rcir): # Third Quadrant
  T=np.array([[1,0,0,xcir],[0,1,0,ycir],[0,0,1,-(rcir**2 - ycir**2)**0.5],[0,0,0,1]])
  Q=GeomIK(T,-1,1,1,1)
  robot.SetDOFValues(Q)
  handles.append(env.plot3(robot.GetLinks()[6].GetTransform()[0:3,3],5))
  ycir=ycir-sam;
  time.sleep(sam/v);

ycir=-rcir;
while (ycir<=0): # Fourth Quadrant
  T=np.array([[1,0,0,xcir],[0,1,0,ycir],[0,0,1,(rcir**2 - ycir**2)**0.5],[0,0,0,1]])
  Q=GeomIK(T,-1,1,1,1)
  robot.SetDOFValues(Q)
  handles.append(env.plot3(robot.GetLinks()[6].GetTransform()[0:3,3],5))
  ycir=ycir+sam;
  time.sleep(sam/v);

time.sleep(2)

# Drawing Line using inverse Kinematics Method 
env.GetViewer().SetCamera(Tline)
zline=float(raw_input("Enter distance along z axis where you want to draw line in m (0.434 is a good choice)\n"))
aline=float(raw_input("Enter slope of the line (-1 is a good choice)\n"))
bline=float(raw_input("Enter y-intercept of the line in m (0.50 is a good choice)\n"))
v=float(raw_input("Enter velocity of Motion in m/s(0.1 is a good choice)\n"))
env.GetViewer().SetTitle('Drawing Line Using Inverse Kinematics/Position Control Method')
time.sleep(2)

i=0;sam=0.005  # sam for sample
xline=0;
handles=[]
while (xline<0.4):
  T=np.array([[1,0,0,xline],[0,1,0,aline*xline+bline],[0,0,1,zline],[0,0,0,1]]) 
  Q=GeomIK(T,-1,-1,1,1)
  robot.SetDOFValues(Q)
  handles.append(env.plot3(robot.GetLinks()[5].GetTransform()[0:3,3],5))
  time.sleep(sam/v)
  xline=xline+ sam*cos(atan(aline))


time.sleep(2)
print 'METHOD : INVERSE JACOBIAN/RATE CONTROL'
print ' '

# Drawing Circle using Inverse Jacobian method
env.GetViewer().SetCamera(Tcirc)
xcir=float(raw_input("Enter distance along x axis where you want to draw circle in m(0.508 is a good choice)\n"))
rcir=float(raw_input("Enter Radius of Circle in m (0.254 is a good choice)\n"))
v=float(raw_input("Enter velocity of Motion in m/s(0.1 is a good choice)\n"))
env.GetViewer().SetTitle('Drawing Circle Using Inverse Jacobian/Rate Control Method')
time.sleep(2)

w=v/rcir; 
tH=0;
Tstart=np.array([[1,0,0,xcir],[0,1,0,rcir*sin(tH)],[0,0,1,rcir*cos(tH)],[0,0,0,1]]) 
Qstart=GeomIK(Tstart,-1,1,1,1)
robot.SetDOFValues(Qstart)
handles=[];
sam=0.01
Qset=Qstart
while(tH<=2*np.pi):
  vtask=np.array([[0],[v*cos(tH)],[-v*sin(tH)],[0],[0],[0]])
  #T=np.array([[1,0,0,xcir],[0,1,0,rcir*sin(tH)],[0,0,1,rcir*cos(tH)],[0,0,0,1]])
  #Q=GeomIK(T,-1,-np.sign(cos(tH)),1,1)
  Qdot=np.dot(linalg.inv(JacoV(Qset)),vtask)
  Qdotm1=np.vstack(Qdot[:, 0]).astype(np.float)
  Qdotm2=Qdotm1.transpose()[0]
  Qset=Qset + (Qdotm2*sam/w);
  robot.SetDOFValues(Qset)
  handles.append(env.plot3(robot.GetLinks()[6].GetTransform()[0:3,3],5))
  time.sleep(sam/w)
  tH=tH+sam;

time.sleep(2)

#Drawing Line using Inverse Jacobian Method
env.GetViewer().SetCamera(Tline)
env.GetViewer().SetCamera(Tline)
zline=float(raw_input("Enter distance along z axis where you want to draw line in m (0.434 is a good choice)\n"))
aline=float(raw_input("Enter slope of the line (-1 is a good choice)\n"))
bline=float(raw_input("Enter y-intercept of the line in m (0.50 is a good choice)\n"))
v=float(raw_input("Enter velocity of Motion in m/s(0.1 is a good choice)\n"))
env.GetViewer().SetTitle('Drawing Line Using Inverse Jacobian/Rate Control Method')
time.sleep(2)

i=0;sam=0.005 # sam for sample
xstart=0.0;
Tstart=np.array([[1,0,0,xstart],[0,1,0,aline*xstart + bline],[0,0,1,zline],[0,0,0,1]])
Qstart=GeomIK(Tstart,-1,-1,1,1)
robot.SetDOFValues(Qstart);xline=0.0;
handles=[];Qset=Qstart
while (xline<0.4):
  xdot=v*cos(atan(aline))
  ydot=v*sin(atan(aline))
  vtask=np.array([[xdot],[ydot],[0],[0],[0],[0]])
  #T=np.array([[1,0,0,xline],[0,1,0,aline*xline+bline],[0,0,1,zline],[0,0,0,1]]) 
  #Q=GeomIK(T,-1,-1,1,1)
  Qdot=np.dot(linalg.inv(JacoV(Qset)),vtask) 
  Qdotm1=np.vstack(Qdot[:, 0]).astype(np.float)
  Qdotm2=Qdotm1.transpose()[0]
  Qset=Qset + (Qdotm2*sam);
  robot.SetDOFValues(Qset)
  handles.append(env.plot3(robot.GetLinks()[5].GetTransform()[0:3,3],5))
  time.sleep(sam)
  xline=xline+v*sam*cos(atan(aline))
  
  


  







