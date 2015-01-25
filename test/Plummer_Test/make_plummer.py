import numpy as np
import os
import sys
import random

# Code adpted from code obtained from Denis Erkal.

###########################################################################
#                                 Options                                 #
##################################r#########################################

PlummerMass = float(sys.argv[1]);
PlummerScale= float(sys.argv[2]);
PlummerNBody = int(sys.argv[3]);
PlummerCOMPos = np.array(sys.argv[4]);
PlummerCOMVel = np.array(sys.argv[5]);
out_name = sys.argv[6];

def GetUnitVector():
    CosTheta = random.uniform(-1.0, 1.0);
    theta = np.arccos(CosTheta);
    phi = random.uniform(0.0, 2.0 * np.pi);
    x = np.sin(theta) * np.cos(phi);
    y = np.sin(theta) * np.sin(phi);
    z = CosTheta;
    return np.array([x,y,z]);

###########################################################################
#                           Get x,y,z positions                           #
###########################################################################

#G = 1 (Henon Units)
def profile(scaledR, mass):
    val = mass * (1 + scaledR**2)**(-2.5);
    return float(val);


def differentialEnclosedMass(r, scale, mass):
    scaledR = r/scale;
    val = 4.0 * np.pi * r**(2.0) * profile(scaledR, mass);
    return float(val);


#############################################################################
#
# The maximum occurs where the derivative of the enclosed mass integrand is 0   # The root of interest occuers at scaledR = sqrt(2/3)  
#
#############################################################################

def PosAcceptRejectMC(scale, mass, N):
    radii = [];
    MAX_RHO =  differentialEnclosedMass(np.sqrt(2.0/3.0)*scale, scale, mass);
    wasAccepted = False;
    for i in range(0,N):
        while (wasAccepted == False):
            unifVal = random.uniform(0.0,1.0);
            rScale = random.uniform(0.0,5.0) #truncate at 5 scale radii
            relProbabilityDistr = differentialEnclosedMass(rScale * scale, scale,mass);
            acceptParam = relProbabilityDistr / MAX_RHO;
            print acceptParam
            if (acceptParam > unifVal):
                radii.append(rScale * scale);
                wasAccepted = True; 
        wasAccepted = False
    return radii;

radii = PosAcceptRejectMC(PlummerScale, PlummerMass, PlummerNBody);
posx = [];
posy = [];
posz = [];
for r in radii:
    vector = (r *np.array( GetUnitVector()));
    posx.append(vector[0]);
    posy.append(vector[1]);
    posz.append(vector[2]);

pos = [posx, posy, posz];
pos = np.transpose(pos);
    


###########################################################################
#                             Get Velocities                              #
###########################################################################

def MassEnclosed(rScale, scale, mass):
    r = rScale * scale;
    val = mass * r**(3.0) * (r**(2.0) + scale**(2.0))**(-3.0/2.0);
    return val;

# Temporary for debug: assign circular orbits.
# G = 1

def escapeVelMagnitude(rScale, scale, mass):
    G = 6.67 * 10.0**(-20.0); # length in km
    # length to km
    scale = scale * 3.086 * 10**16.0;
    r = rScale * scale;
    mass = mass * 1.989 * 10.**30.; # solar masses
    #massEnclosed = MassEnclosed(rScale, scale, mass)
    val = np.sqrt(2 * G * mass / np.sqrt(r**2 + scale**2));
    return val;

def selectFromProfile(): # von Neumann
    x = random.uniform(0.0,1.0);
    y = random.uniform(0.0, 0.1);
    
    while(y > (x**2 *(1 - x)**(3.5))):
        x = random.uniform(0.0,1.0);
        y = random.uniform(0.0, 0.1);
    return x;


#radii = PosAccepRejectMC(PlummerScale, PlummerMass, PlummerNBody);

velx = [];
vely = [];
velz = [];
for r in radii:
    rScale = r / PlummerScale;
    velProfile = selectFromProfile();
    velMag = velProfile * escapeVelMagnitude(r, PlummerScale, PlummerMass);
    v = velMag * GetUnitVector();
    velx.append(v[0]);
    vely.append(v[1]);
    velz.append(v[2]);

#print velx
#print vely
#print velz

vel = [velx, vely, velz];
vel = np.transpose(vel);

print vel;
###########################################################################
#                             CoM Correction                              #
###########################################################################

for p in pos:
    #print p
    #print PlummerCOMPos
    p = np.add(p, PlummerCOMPos);

for v in vel:
    v = np.add(v,PlummerCOMVel);

#print pos
#print len(pos)
#print sys.getsizeof(pos);
#print v
#print len(vel)
#print sys.getsizeof(vel)

###########################################################################
#                  Adding Orbital Position and Velocity                   #
###########################################################################

###########################################################################
#                           Output to text file                           #
###########################################################################

f = open(out_name,'w');
for i in range(0, len(pos)):
    line = str(i) + " 1 " + str(PlummerMass/float(PlummerNBody)) + " " +str(pos[i][0]) + " " + str(pos[i][1]) + " " + str(pos[i][2]) \
           + " " + str(vel[i][0]) + " " + str(vel[i][1]) + " " + str(vel[i][2]) + "\n";
    f.write(line);

f.close();

###########################################################################
#                            Output HDF5 file                             #
###########################################################################


#print "Outputting Binary file -", out_name
#print

#try:
#    with open( out_name ) as f:
#        pass;
#    print "Output file exists - deleting it.";
#    os.remove( out_name );
#except IOError:
#    pass;

#file = open(out_name,'wb');

#buffer_size = np.array([256],np.uint32);
#buffer_size.tofile(file,"",format="%i");

#Write Header

#nparts = PlummerNBody;
#npart_array = np.array([0,nparts,0,0,0,0],np.uint32);
#npart_array.tofile(file,"",format="%i");

#24 bytes
#mpart = PlummerMass / PlummerNBody;
#mass_array = np.array([0,mpart,0,0,0,0],np.float64); # masses for each type
#mass_array.tofile(file,"",format="%d");

#72 bytes

#time = np.array([0],np.float64);
#time.tofile(file,"",format="%d");

#80 bytes

#redshift = np.array([0],np.float64)
#redshift.tofile(file,"",format="%d")

#flagsfr = np.array([0],np.uint32)
#flagsfr.tofile(file,"",format="%i")

#flagFB = np.array([0],np.uint32)
#flagFB.tofile(file,"",format="%i")

#npart_array.tofile(file,"",format="%i")

#flagcooling = np.array([0],np.uint32)
#flagcooling.tofile(file,"",format="%i")

#numfiles = np.array([0],np.uint32)
#numfiles.tofile(file,"",format="%i")

#boxsize = np.array([0],np.float64)
#boxsize.tofile(file,"",format="%d")

#Omega0 = np.array([0],np.float64)
#Omega0.tofile(file,"",format="%d")

#OmegaL0 = np.array([0],np.float64)
#OmegaL0.tofile(file,"",format="%d")

#H0 = np.array([0],np.float64)
#H0.tofile(file,"",format="%d")

#unused_buffer = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],np.uint32)
#unused_buffer.tofile(file,"",format="%d")

#buffer_size.tofile(file,"",format="%i")

#Write Positions

#value = nparts*3*4
#buffer_size = np.array([value],np.uint32)

#buffer_size.tofile(file,"",format="%i")
#coords = np.array(pos, np.float32)
#coords.tofile(file,"",format="%f")
#buffer_size.tofile(file,"",format="%i")

#Write Velocities

#buffer_size.tofile(file,"",format="%i")
#vels   = np.transpose(np.array(vel, np.float32 ))
#vels.tofile(file,"",format="%f")
#buffer_size.tofile(file,"",format="%i")

#Write IDs

#value = nparts*4
#buffer_size = np.array([value],np.uint32)

#buffer_size.tofile(file,"",format="%i")
#ids    = np.arange( 1, nparts + 1, dtype=np.uint32 )
#ids.tofile(file,"",format="%i")
#buffer_size.tofile(file,"",format="%i")

#file.close()

#HDF5

#of = snapHDF5.openfile( out_name )

#npart  = np.array([0, 0, nparts, 0, 0, 0], np.uint32)
#marray = np.array([0, 0, mpart, 0, 0, 0], np.float64)
#header = snapHDF5.snapshot_header(npart=npart, nall=npart, massarr=marray)

#coords = np.array( [x, y, z], np.float32 ).T
#vels   = np.array( [vx, vy, vz], np.float32 ).T
#ids    = np.arange( 1, nparts + 1, dtype=np.uint32 )

#snapHDF5.writeheader(of, header)
#snapHDF5.write_block(of, "POS ", 2, coords)
#snapHDF5.write_block(of, "VEL ", 2, vels  )
#snapHDF5.write_block(of, "ID  ", 2, ids   )

#snapHDF5.closefile(of)
print "IC finished"

