# Drex #

4-17-15 Christopher Thissen, Yale University
christopher.thissen@yale.edu

This code simulates the development of LPO textures in a deforming
olivine aggregate. This is a modified version of the popular fortran code
D-Rex, originally released by Ed Kaminski, Neil Ribe, and others. 

The code has been modified from the original fortran version. This
version of the code only tracks the development of crystal sizes and
orientations. All path-line integration and seismic anistropy
calculations present in the original fortran version have been removed.
All calculations relating to enstatite have also been removed. The
code is suitable for simulating olivine fabrics, and can easily be
modified for calculating fabric development for a given velocity field.

Other modifications include:
1. The initially random LPO is now generated using a random stream, and
follows the method in Morawiec, A. (2003). Orientations and rotations.
Springer-Verlag.
2. The code warns the user when the rotation matrix used to update a
crystal orientation is not a proper rotation matrix. 
3. An indexing error in the original fortran code related to the
calculation of dislocation densities has be fixed. 

The theoretical development of the model is described in the following
publications: 
Kaminski, E., & Ribe, N. M. (2001). A kinematic model for 
recrystallization and texture development in olivine polycrystals.
Earth and Planetary Science Letters, 189(3), 253-267.

Kaminski, E., Ribe, N. M., & Browaeys, J. T. (2004). D-Rex, a program 
for calculation of seismic anisotropy due to crystal lattice preferred 
orientation in the convective upper mantle. Geophysical Journal
International, 158(2), 744-752.

### USAGE ###

#### INPUT: ####
The code does not require input arguments and can be run by
calling the program the commandline or pushing the green "Run" button.
The user specifies two types of parameters in the following sections. The
first type of parameters relates to the olivine aggregate, such as the
number of olivine crystals and the stress exponent. These parameters are
defined in the next section as part of the "Grain" structure. Default
values and alternative options are also discussed.

The second type of input parameter relates to the imposed deformation.
The user may choose from among the pre-defined deformation gradient
tensors by setting the "deformationSymmetry" variable to one of the
following strings: axisymmetricCompression, axisymmetricExtension,
orthorhombic, simpleShear, triclinic, or noDeformation. Custom
deformation gradient tensors can also be added in the section labelled
"Define Deformation Gradient Tensor". The number of integration steps
must also be specified.

#### OUTPUT: ####
The code outputs several text files and figures. Five text files are
created. The first is an info file that lists important parameters and
final finite strain results for the aggregtate. The second text file
lists the final orientations (ZXZ, Bunge convention, in radians),
and the volume fraction of each crystal. Three other text files are
also output that contain only euler angles. The unweighted file includes
a single measurement for all crystals that have volume fraction greater
than zero. This is akin to making a single measurement of each grain. The
volume weighted file includes repeated measurements for the same grain scaled by the final volume fraction. 
This is akin to making measurements on a predefined grid. The inverse
volume weighted file includes repreate measurements for each grain

The code also outputs a figure containing pole figures for the
[100],[010], and [001] axes. The default coordinate frame for the pole
figures is defined as follows:
<pre>                 
            %%%     %%%
       %%%              %%%
 
   %%%                      %%%
 
  %%%                         %%%
 %              +Z             +X % 
  %%%                         %%%
 
  %%%                        %%%
 
     %%%                  %%%
 
           %%%  +Y  %%%
</pre>               
The pole figures are saved using the export_fig toolbox. A version of this
toolbox is included.
