<B>Tahoe</B>

Tahoe is a research-oriented platform for the development of numerical methods and material models. The goal of the work surrounding Tahoe is the simulation of stresses and deformations for situations that cannot be treated by standard continuum simulation techniques. These situations include material fracture or failure, interfacial adhesion and debonding, shear banding, length-scale dependent elasticity and plasticity, and deformation in small-scale structures. Aside from a collection of standard finite elements, Tahoe includes meshfree simulation capability and particle methods. Tahoe includes a number of "cohesive" approaches for modeling fracture. These include both surface and bulk constitutive models that incorporate cohesive behavior. Tahoe is capable of performing static and transient dynamic coupled-physics analysis in two and three dimensions. Many capabilities support parallel execution.

sourceforge page: <a href="http://tahoe.sourceforge.net/">http://tahoe.sourceforge.net/</a>

<B>Building Tahoe</B>
There are a straightforwrd method to install Tahoe in your unix machine. Simply go to <a href="http://sourceforge.net/projects/tahoe/">http://sourceforge.net/projects/tahoe/</a> and download the files. Then execute the following command line:

<code>$~: tahoe-manager init buil</code>

for more details take a look at README file in this repository.

<B>Pre Built Tahoe</B>

This is a pre built version of Tahoe specifically configurated for developing Dielectric Elastomer models.

<B>Linux Ubuntu 12.04</B>

It seems for linux (Ubuntu 12.04 and 14.04 is tested) the best choose for compiler is GNU-GCC4.6. Before building the Tahoe make sure you have GCC 4.6 installed in your machine. Defualt GCC compiler on Ubuntu 12.04 is 4.6. However, ubuntu 14.04 uses GCC 4.8. However, it's possible to install GCC 4.6 using Synaptic Package Manager or directly from ubuntu repository. Following command:

<code>$~: sudo apt-get install gcc-4.6 gfortran-4.6 cpp-4.6</code>

Now some changes should be made in macro files. Simply go to Tahoe/macros and edit the file GNU-GCC4.6.macros and change all gcc, g++ and gfortran to gcc-4.6, g++-4.6 and gfortran-4.6. This should work.
